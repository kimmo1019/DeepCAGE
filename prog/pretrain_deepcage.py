#multi-gpu-model
import argparse
import random,os,sys
import numpy as np
from scipy import stats
from pyfasta import Fasta
import time
from multiprocessing import Pool 
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_auc_score
from sklearn import preprocessing
import pandas as pd
import keras.backend as K
from keras.models import Model, Sequential
from keras.models import load_model
from keras.layers import Input,InputLayer,Multiply,ZeroPadding2D
from keras.layers import Conv2D, MaxPooling2D
from keras.layers import Dense,Activation,Dropout,Flatten,Concatenate
from keras.layers import BatchNormalization
from keras import optimizers,utils
from keras.constraints import max_norm
from keras import regularizers
from keras.callbacks import ModelCheckpoint,Callback,EarlyStopping,History
from keras.utils import multi_gpu_model
from keras.optimizers import Adam, SGD
from keras.models import model_from_json
import tensorflow as tf
from sklearn.metrics import average_precision_score

######################## Usage #######################
usage='''
Usage: python 4.Regression.py [FOLD_ID] [RATIO]
-- a program for regression (cross cell types)
[GPU_ID] : GPU ID (from 0 to 7)
[FOLD_ID]: which fold (0-4)
'''
# if len(sys.argv)!=3:
#     print usage
#     sys.exit(1)
######################## Global Settings #######################
os.environ["CUDA_VISIBLE_DEVICES"] = sys.argv[1]
SEQ_LEN = 1000
#NUM_TF = 769
NUM_TF = 678

fold_idx = 4
ratio = 1
suffix = 'fold%s_ratio%s'%(str(fold_idx),str(ratio))
###################### Model Implementation #####################
#gene expression quantitile normalization
def quantile_norm(matrix):
    rank_mean = matrix.stack().groupby(matrix.rank(method='first').stack().astype(int)).mean()
    return matrix.rank(method='min').stack().astype(int).map(rank_mean).unstack()

#one-hot coding
def one_hot_seq(row_name):
    chrom = row_name.split(':')[0]
    start = int(row_name.split(':')[-1].split('-')[0])
    end = int(row_name.split(':')[-1].split('-')[-1])
    assert end - start == 200
    midpoint = (start+end)/2
    seq = genome[chrom][midpoint-SEQ_LEN/2:midpoint+SEQ_LEN/2]
    d = {'a':0, 'A':0, 'c':1, 'C':1, 'g':2, 'G':2, 't':3, 'T':3, 'N':4, 'n':4}
    mat = np.zeros((5,SEQ_LEN,1))  
    for i in range(len(seq)):
        mat[d[seq[i]],i,0] = 1
    mat = mat[:4,:,:]
    return mat

def conv_block(x, stage, branch, nb_filter, dropout_rate=None, weight_decay=1e-4):
    eps = 1.1e-5
    conv_name_base = 'conv' + str(stage) + '_' + str(branch)
    relu_name_base = 'relu' + str(stage) + '_' + str(branch)

    # 1x1 Convolution (Bottleneck layer)
    inter_channel = nb_filter
    x = BatchNormalization(epsilon=eps, name=conv_name_base+'_x1_bn')(x)
    x = Activation('relu', name=relu_name_base+'_x1')(x)
    x = Conv2D(inter_channel, (1, 1), name=conv_name_base+'_x1',use_bias=False)(x)
    if dropout_rate:
        x = Dropout(dropout_rate)(x)

    # 3x3 Convolution
    x = BatchNormalization(epsilon=eps, name=conv_name_base+'_x2_bn')(x)
    x = Activation('relu', name=relu_name_base+'_x2')(x)
    x = ZeroPadding2D((0, 1), name=conv_name_base+'_x2_zeropadding')(x)
    x = Conv2D(nb_filter, (1, 3), name=conv_name_base+'_x2',use_bias=False)(x)
    if dropout_rate:
        x = Dropout(dropout_rate)(x)

    return x

def dense_block(x, stage, nb_layers, nb_filter, growth_rate, dropout_rate=None, weight_decay=1e-4, grow_nb_filters=False):
    eps = 1.1e-5
    concat_feat = x

    for i in range(nb_layers):
        branch = i+1
        x = conv_block(concat_feat, stage, branch, growth_rate, dropout_rate, weight_decay)
        concat_feat = Concatenate(axis = -1,name='concat_'+str(stage)+'_'+str(branch))([concat_feat, x])
        if grow_nb_filters:
            nb_filter += growth_rate

    return concat_feat, nb_filter

def precision(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision
def recall(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall
def f1_score(y_true, y_pred):
    prec = precision(y_true, y_pred)
    recal = recall(y_true, y_pred)
    return 2.0*prec*recal/(prec+recal+K.epsilon())
def average_precision(y_true, y_pred):
    return tf.py_func(average_precision_score, (y_true, y_pred), tf.double)    

class Mycallback(Callback)            :
    def __init__(self,model):
        self.single_model = model
    def on_epoch_end(self,epoch,logs={}):
        lr = self.model.optimizer.lr
        decay = self.model.optimizer.decay
        iterations = self.model.optimizer.iterations
        lr_with_decay = lr / (1. + decay * K.cast(iterations, K.dtype(decay)))
        print(K.eval(lr_with_decay))
        self.single_model.save_weights('%s/checkpoint/regression_multi_gpu_model_weights_at_epoch_%d.h5'%(DPATH,epoch))

class roc_callback(Callback):
    def __init__(self,validation_data,patience):
        self.x_val = validation_data[0]
        self.y_val = validation_data[1]
        self.patience = patience
        self.best_weight = None
        self._time = time.time()

    def on_train_begin(self, logs={}):
        self.wait = 0
        self.stopped_epoch = 0
        self.best = -np.Inf
        return

    def on_train_end(self, logs={}):
        self.model.set_weights(self.best_weight)
        if self.stopped_epoch > 0 :
            print('Epoch %05d: early stopping' % (self.stopped_epoch + 1))
        return

    def on_epoch_begin(self, epoch, logs={}):
        return

    def on_epoch_end(self, epoch, logs={}):
        y_pred_val = self.model.predict(self.x_val)
        roc_val = roc_auc_score(self.y_val, y_pred_val)
        precision,recall,_, = metrics.precision_recall_curve(self.y_val,y_pred_val)
        pr_val = -np.trapz(precision,recall)
        print '\nroc-auc_val: %s' % str(round(roc_val,4))
        print 'pr-auc_val: %s' % str(round(pr_val,4))
        print 'epoch training time: %d s' % (time.time()-self._time)
        if pr_val > self.best:
            self.best = pr_val
            self.wait = 0
            self.best_weight = self.model.get_weights()
        else:
            self.wait+=1
            if self.wait >= self.patience:
                self.stopped_epoch = epoch
                self.model.stop_training = True
        return

    def on_batch_begin(self, batch, logs={}):
        return

    def on_batch_end(self, batch, logs={}):
        return



def model_construct(output_dim):
    #paramters define
    conv1_params = {'nb_filter':160,'filter_size':(4,15),'pool_size':(1,4)}
    conv2_params = {'nb_filter':160,'filter_size':(1,12),'pool_size':(1,4)}
    conv3_params = {'nb_filter':160,'filter_size':(1,12),'pool_size':(1,4)}
    seq_input = Input(shape=(4,SEQ_LEN,1),name='seq_input')
    #gexp_input = Input(shape=(NUM_TF,),name='tf_gexp_input')
    #motif_input = Input(shape=(NUM_TF,),name='tf_motifscore_input')
    x = Conv2D(filters=conv1_params['nb_filter'], kernel_size=conv1_params['filter_size'],strides=(1, 1), activation = 'relu',padding='valid',name='conv1')(seq_input)
    x = BatchNormalization()(x)
    x = MaxPooling2D(pool_size=conv1_params['pool_size'],name='pooling1')(x)
    x = Dropout(0.2)(x)
    
    #Denseblock1
    x, nb_filter = dense_block(x, stage=1, nb_layers=5, nb_filter=64, growth_rate=32, dropout_rate=0.0, weight_decay=1e-4)

    x = Conv2D(filters=conv2_params['nb_filter'], kernel_size=conv2_params['filter_size'],strides=(1, 1), activation = 'relu',padding='valid',name='conv2')(x)
    x = BatchNormalization()(x)
    x = MaxPooling2D(pool_size=conv2_params['pool_size'],name='pooling2')(x)
    x = Dropout(0.2)(x)
    
    #Denseblock2
    x, nb_filter = dense_block(x, stage=2, nb_layers=5, nb_filter=64, growth_rate=32, dropout_rate=0.0, weight_decay=1e-4)

    x = Conv2D(filters=conv3_params['nb_filter'], kernel_size=conv3_params['filter_size'],strides=(1, 1), activation = 'relu',padding='valid',name='conv3')(x)
    x = BatchNormalization()(x)
    x = MaxPooling2D(pool_size=conv3_params['pool_size'],name='pooling3')(x)
    x = Dropout(0.2)(x)
    
    #Denseblock3
    x, nb_filter = dense_block(x, stage=3, nb_layers=5, nb_filter=64, growth_rate=32, dropout_rate=0.0, weight_decay=1e-4)
    
    x = Flatten()(x)
    #att_vec = Dense(402, activation='softmax', name='attention_vec')(gexp_input)
    #x1 = Multiply()([gexp_input,motif_input])
    #x1 = gexp_input
    #x = Concatenate()([x,x1])
    x = Dense(256)(x)
    x = Activation('relu')(x)
    output = Dense(output_dim,name='output')(x)
    x = Activation('relu')(x)
    #model  = Model(inputs=[seq_input,gexp_input,motif_input],outputs=output)    
    model  = Model(inputs=seq_input,outputs=output)    
    #model.summary()
    return model



################################## Model Training##################################
def model_training(model,batch_size,epochs):
    #callbacks = [roc_callback(validation_data=(X_valid,Y_valid),patience=6)]
    callbacks = [EarlyStopping(monitor='val_loss', patience=5, verbose=1),ModelCheckpoint('../data/encode/checkpoint/pretrain_full.h5', monitor='val_loss',save_best_only=True,save_weights_only=False)]
    model.compile(optimizer='adam',loss='mean_squared_error',metrics=['mse'])
    model.fit(X_seq_train,Y_train,batch_size=batch_size,epochs=epochs,callbacks=callbacks,validation_data=(X_seq_valid,Y_valid))
    return model

def model_evaluation(model,batch_size):
    Y_pred = model.predict(X_seq_valid,batch_size=batch_size)
    pcc=[]
    #f_out = open('%s/log/pretrain.log','w')
    for i in range(Y_valid.shape[-1]):
        pearsonr = stats.pearsonr(Y_pred[:,i],Y_valid[:,i])[0]
        pcc.append(pearsonr)
        print pearsonr
    print(np.mean(pcc),np.median(pcc))

       
        
if  __name__ == "__main__" :
    DPATH='/home/liuqiao/software/DeepCAGE/data/encode'
    genome = Fasta('/home/liuqiao/openness_pre/genome.fa')
    readscount_file='%s/readscount_normalized_filtered.csv'%DPATH
    readscount = pd.read_csv(readscount_file,header=0,index_col=[0],sep='\t')
    readscount = np.log(readscount+1)
    exps = list(readscount.columns)
    exps_train = exps[:int(len(exps)*0.8)]
    exps_test = exps[int(len(exps)*0.8):]
    readscount_train = readscount[exps_train]
    readscount_test = readscount[exps_test]
    nb_exps_train = len(exps_train)
    if True:
    #if False:
        #X_seq = np.load('../data/encode/X_seq.npy')
        #X_seq = X_seq[...,np.newaxis]
        X_seq = np.load('../data/encode/X_seq_full.npy')
    else:
        pool = Pool()
        #X_seq=np.stack(pool.map(one_hot_seq,list(readscount_train.index)))
        X_seq=np.stack(pool.map(one_hot_seq,list(readscount.index)))
        pool.close()
        pool.join()
        np.save('../data/encode/X_seq_full.npy',X_seq)
    Y = readscount.values
    random.seed(0)
    idx = np.arange(X_seq.shape[0])
    random.shuffle(idx)
    X_seq_train = X_seq[idx[:int(0.95*len(idx))]]
    Y_train = Y[idx[:int(0.95*len(idx))]]
    X_seq_valid = X_seq[idx[int(0.95*len(idx)):]]
    Y_valid = Y[idx[int(0.95*len(idx)):]]

    model = model_construct(Y.shape[-1])
    optimizer = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0000)
    #model = model_training(model,batch_size=128,epochs=100)
    model = load_model('../data/encode/checkpoint/pretrain_full.h5')
    #model.save('%s/CNN/models/final/classification_fold%d_seed%d.h5'%(DPATH,fold_idx,int(sys.argv[4])))
    model_evaluation(model,batch_size=1024)
    
    
    
    
    
    

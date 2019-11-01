#multi-gpu-model
import argparse
import random,os,sys
import numpy as np
from scipy import stats
from pyfasta import Fasta
import time
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
os.environ["CUDA_VISIBLE_DEVICES"] = '0,1,2,3,4,5,6,7'
SEQ_LEN = 1000
#NUM_TF = 769
NUM_TF = 678

fold_idx = int(sys.argv[1])
ratio = float(sys.argv[2])
suffix = 'fold%s_ratio%s'%(str(fold_idx),str(ratio))
###################### Model Implementation #####################
#gene expression quantitile normalization
def quantile_norm(matrix):
    rank_mean = matrix.stack().groupby(matrix.rank(method='first').stack().astype(int)).mean()
    return matrix.rank(method='min').stack().astype(int).map(rank_mean).unstack()

#one-hot coding
def seq_to_mat(seq):
    d = {'a':0, 'A':0, 'g':1, 'G':1, 'c':2, 'C':2, 't':3, 'T':3, 'N':4, 'n':4}
    mat = np.zeros((5,SEQ_LEN))  
    for i in range(len(seq)):
        mat[d[seq[i]],i] = 1
    mat = mat[:4,:]
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



def model_construct():
    #paramters define
    conv1_params = {'nb_filter':160,'filter_size':(4,15),'pool_size':(1,4)}
    conv2_params = {'nb_filter':160,'filter_size':(1,12),'pool_size':(1,4)}
    conv3_params = {'nb_filter':160,'filter_size':(1,12),'pool_size':(1,4)}
    seq_input = Input(shape=(4,SEQ_LEN,1),name='seq_input')
    gexp_input = Input(shape=(NUM_TF,),name='tf_gexp_input')
    motif_input = Input(shape=(NUM_TF,),name='tf_motifscore_input')
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
    x1 = Multiply()([gexp_input,motif_input])
    #x1 = gexp_input
    x = Concatenate()([x,x1])
    x = Dense(256)(x)
    x = BatchNormalization()(x)
    x = Activation('relu')(x)
    x = Dropout(0.5)(x)
    output = Dense(1, kernel_regularizer=regularizers.l2(0.01),name='output')(x)
    model  = Model(inputs=[seq_input,gexp_input,motif_input],outputs=output)    
    #model.summary()
    return model

def generate_batch_data(training_queue,batch_size,padding=400,bootstrap=False):
    X_batch_mat = np.zeros((batch_size,4,SEQ_LEN,1))
    X_batch_gexp = np.zeros((batch_size,NUM_TF))
    X_batch_motif = np.zeros((batch_size,NUM_TF))
    Y_batch = np.zeros(batch_size)
    if bootstrap:
        while True:
            present_queue = random.sample(training_queue,batch_size)
            for i in range(len(present_queue)):
                region_idx,cellid = present_queue[i]
                info = readscount.index[region_idx]
                chrom,start,end = info.split(':')[0],int(info.split(':')[1].split('-')[0]),int(info.split(':')[1].split('-')[1])
                sequence = genome[chrom][start:end]
                motif_idx = chrom+':'+str(start)+'-'+str(end)
                X_batch_mat[i,:,:,0] = seq_to_mat(sequence)
                X_batch_gexp[i,:] = tf_gexp_norm.loc[cellid].values.reshape((1,-1))
                X_batch_motif[i,:] = motif.loc[motif_idx].values.reshape((1,-1)) if motif_idx in motif.index else np.zeros((1,motif.shape[1]))
                Y_batch[i] = readscount.iloc[region_idx].loc[str(cellid)]
            yield [X_batch_mat,X_batch_gexp,X_batch_motif], Y_batch            
    else:
        while True:
            for step in range(len(training_queue)/batch_size):
                present_queue = training_queue[step*batch_size:(step+1)*batch_size]
                for i in range(len(present_queue)):
                    region_idx,cellid = present_queue[i]
                    info = readscount.index[region_idx]
                    chrom,start,end = info.split(':')[0],int(info.split(':')[1].split('-')[0]),int(info.split(':')[1].split('-')[1])
                    sequence = genome[chrom][start:end]
                    motif_idx = chrom+':'+str(start)+'-'+str(end)
                    X_batch_mat[i,:,:,0] = seq_to_mat(sequence)
                    X_batch_gexp[i,:] = tf_gexp_norm.loc[cellid].values.reshape((1,-1))
                    X_batch_motif[i,:] = motif.loc[motif_idx].values.reshape((1,-1))
                    Y_batch[i] = readscount.iloc[region_idx].loc[str(cellid)]
                yield [X_batch_mat,X_batch_gexp,X_batch_motif], Y_batch


################################## Model Training##################################
def model_training(parallel_model,model,batch_size,epochs):
    data_queue=[]
    for idx in region_idx:
        for cellid in train_cell_idx:
            data_queue.append((idx,cellid))
    random.shuffle(data_queue)
    train_queue = data_queue[:int(0.95*len(data_queue))]
    valid_queue = data_queue[int(0.95*len(data_queue)):]
    print '%d training and %d validation examples'%(len(train_queue),len(valid_queue))
    steps_per_epoch = len(train_queue)/batch_size
    validation_steps = len(valid_queue)/batch_size
    callbacks = [Mycallback(model),
                 EarlyStopping(monitor='val_loss', patience=5, verbose=1,mode='auto'),
                 ModelCheckpoint('../data/checkpoint/regression_best_weights_multigpu.h5',monitor='val_loss',save_best_only=True,save_weights_only=True,mode='auto')]
    training_generator = generate_batch_data(train_queue,batch_size=batch_size,padding=400,bootstrap=False)
    validation_generator = generate_batch_data(valid_queue,batch_size=batch_size,padding=400,bootstrap=False)
    parallel_model.fit_generator(generator=training_generator,steps_per_epoch=steps_per_epoch,epochs=epochs,verbose=1,callbacks=callbacks,validation_data=validation_generator,validation_steps=validation_steps,max_queue_size=100,workers=28,use_multiprocessing=True)
    return parallel_model, model

################################## Model Evaluation#################################
def model_evaluation(model,region_idx,padding=400):
    corrs=[]
    for cellid in test_cell_idx:
        print 'Predicting cell %d'%cellid
        X_test_mat = np.zeros((len(region_idx),4,SEQ_LEN,1))
        X_test_gexp = np.zeros((len(region_idx),NUM_TF))
        X_test_motif = np.zeros((len(region_idx),NUM_TF))
        Y_test = np.zeros(len(region_idx))
        for i in range(len(region_idx)):
            info = readscount.index[region_idx[i]]
            chrom,start,end = info.split(':')[0],int(info.split(':')[1].split('-')[0]),int(info.split(':')[1].split('-')[1])
            sequence = genome[chrom][start:end]
            motif_idx = chrom+':'+str(start)+'-'+str(end)
            X_test_mat[i,:,:,0] = seq_to_mat(sequence)
            X_test_gexp[i,:] = tf_gexp_norm.loc[cellid].values.reshape((1,-1))
            X_test_motif[i,:] = motif.loc[motif_idx].values.reshape((1,-1))
            Y_test[i] = readscount.iloc[region_idx[i]].loc[str(cellid)]
        Y_pred = model.predict([X_test_mat,X_test_gexp,X_test_motif],batch_size=1024)
        pearsonr = stats.pearsonr(Y_pred[:,0],Y_test)[0]
        corrs.append(pearsonr)
        f_out = open('%s/regression_%s.log'%(DPATH,suffix),'a+')
        f_out.write('%d\t%.5f\n'%(cellid,pearsonr))
        f_out.close()
    print('Median pearson r: %.4f'%np.median(corrs))
    print('Mean pearson r: %.4f'%np.mean(corrs))


        
        
if  __name__ == "__main__" :
    DPATH='/home/liuqiao/openness_pre/expression/data'
    genome = Fasta('/home/liuqiao/openness_pre/genome.fa')
    readscount_file='%s/processed_RNA_DNase/union.peaks.readscount.norm'%DPATH
    tf_gexp_file = '%s/processed_RNA_DNase/tf_gexp_matrix_filtered.txt'%DPATH
    motifscore_file = '%s/motif_db/HOCOMOCO/motif_score_mat_pvalue0.001_filtered.txt' %DPATH
    tf_gexp = pd.read_csv(tf_gexp_file ,sep='\t',header=0,index_col=[0])
    tf_gexp_log = np.log(tf_gexp+1)
    tf_gexp_norm = pd.DataFrame.transpose(quantile_norm(pd.DataFrame.transpose(tf_gexp_log)))
    motif = pd.read_csv(motifscore_file,sep='\t',header=0,index_col=[0])
    readscount = pd.read_csv(readscount_file,header=0,index_col=[0],sep='\t')
    readscount = np.log(readscount+1)
    if ratio<0 or ratio>1:
        print 'Input ratio between 0 and 1'
        sys.exit()    
    cellid_info = open('%s/processed_RNA_DNase/five_fold_cells.txt'%DPATH).readlines()[fold_idx]
    train_cell_idx = [int(cellid) for cellid in cellid_info.split('\t')[2].split(' ')]
    test_cell_idx = [int(cellid) for cellid in cellid_info.strip().split('\t')[4].split(' ')]
    random.seed(1234)
    region_idx = random.sample(list(range(readscount.shape[0])),int(readscount.shape[0]*ratio))
    model = model_construct()
    optimizer = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0000)
    parallel_model = multi_gpu_model(model,gpus = 8)
    parallel_model.compile(optimizer=optimizer,loss='mean_squared_error',metrics=['mse'])
    print 'Start training...'
    parallel_model, model = model_training(parallel_model,model,batch_size=128*8,epochs=20)
    #model.save('%s/CNN/models/final/classification_fold%d_seed%d.h5'%(DPATH,fold_idx,int(sys.argv[4])))
    model_evaluation(parallel_model,region_idx,padding=400)
    
    
    
    
    
    

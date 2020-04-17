#multi-gpu-model, use 1672 tf gexpr from fantom database, we don't use motif score now
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
Usage: python fine_tune_deepcage_fantom.py [GPU_ID]
-- a program for regression (cross cell types)
[GPU_ID] : GPU ID (from 0 to 7)
'''

######################## Global Settings #######################
os.environ["CUDA_VISIBLE_DEVICES"] = sys.argv[1]
SEQ_LEN = 1000
#NUM_TF = 711
NUM_TF = 1672

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

def extent_region(df,padding = 400):
    region_info = df.index
    region_extend_info=[]
    for each in region_info:
        chrom,start,end = each.split(':')[0],int(each.split(':')[1].split('-')[0])-padding,int(each.split(':')[1].split('-')[1])+padding
        region_extend_info.append(chrom+':'+str(start)+'-'+str(end))
    df.index = region_extend_info
    return df
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

class Mycallback(Callback)            :
    def __init__(self,model):
        self.model = model
    def on_epoch_end(self,epoch,logs={}):
        #lr = self.model.optimizer.lr
        #decay = self.model.optimizer.decay
        #iterations = self.model.optimizer.iterations
        #lr_with_decay = lr / (1. + decay * K.cast(iterations, K.dtype(decay)))
        #print(K.eval(lr_with_decay))
        #self.model.save('../data/encode/checkpoint/model_at_epoch_%d.h5'%epoch)
        model.save('%s/checkpoint/model_at_epoch_%d_fantom.h5'%(DPATH,epoch))


def model_construct():
    #paramters define
    conv1_params = {'nb_filter':160,'filter_size':(4,15),'pool_size':(1,4)}
    conv2_params = {'nb_filter':160,'filter_size':(1,12),'pool_size':(1,4)}
    conv3_params = {'nb_filter':160,'filter_size':(1,12),'pool_size':(1,4)}
    seq_input = Input(shape=(4,SEQ_LEN,1),name='seq_input')
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

    gexp_input = Input(shape=(NUM_TF,),name='tf_gexp_input')#add1
    motif_input = Input(shape=(NUM_TF,),name='tf_motifscore_input')#add2
    #x1 = Multiply()([gexp_input,motif_input])#add3
    #x = Concatenate()([x,x1])#add4
    x = Concatenate()([x,gexp_input])

    x = Dense(256)(x)
    x = Activation('relu')(x)
    #output = Dense(output_dim,name='output')(x)
    output = Dense(1,name='output')(x)#modify
    x = Activation('relu')(x)
    model  = Model(inputs=[seq_input,gexp_input,motif_input],outputs=output)
    #model  = Model(inputs=seq_input,outputs=output)    
    #model.summary()
    return model



def generate_batch_data(training_queue,batch_size,use_paral=False):
    if not use_paral:
        X_batch_mat = np.zeros((batch_size,4,SEQ_LEN,1))
        X_batch_gexp = np.zeros((batch_size,NUM_TF))
        X_batch_motif = np.zeros((batch_size,NUM_TF))
        Y_batch = np.zeros(batch_size)
    def meta2data(info):
        cell_idx,region_idx = present_queue[i]
        x_onehot = X_seq[region_idx]
        x_gexp = tf_gexp_mat[cell_idx,:]
        x_motifscore = tf_motifscore_mat[region_idx,:]
        y = readscount_mat[region_idx,cell_idx]
        return [x_onehot,x_gexp,x_motifscore,y]

    while True:
        for step in range(len(training_queue)/batch_size):
            present_queue = training_queue[step*batch_size:(step+1)*batch_size]
            if use_paral:
                pool = Pool()
                batch_data =pool.map(meta2data,present_queue)
                X_batch_mat = np.stack([item[0] for item in batch_data])
                X_batch_gexp = np.stack([item[1] for item in batch_data])
                X_batch_motif = np.stack([item[2] for item in batch_data])
                Y_batch = np.array([item[3] for item in batch_data])
                yield [X_batch_mat,X_batch_gexp,X_batch_motif], Y_batch
            else:
                for i in range(len(present_queue)):
                    cell_idx,region_idx = present_queue[i]
                    #X_batch_mat[i,:,:,:] = one_hot_seq(region_info[region_idx])
                    X_batch_mat[i,:,:,:] = X_seq[region_idx]
                    X_batch_gexp[i,:] = tf_gexp_mat[cell_idx,:]
                    #X_batch_motif[i,:] = tf_motifscore_mat[region_idx,:]
                    Y_batch[i] = readscount_mat[region_idx,cell_idx]
                yield [X_batch_mat,X_batch_gexp,X_batch_motif], Y_batch

################################## Model Training##################################
def model_training(data_queue,batch_size,epochs):
    print('Training size %d'%len(data_queue))
    train_queue = data_queue[:int(0.95*len(data_queue))]
    valid_queue = data_queue[int(0.95*len(data_queue)):]
    steps_per_epoch = len(train_queue)/batch_size
    validation_steps = len(valid_queue)/batch_size
    print steps_per_epoch,validation_steps
    callbacks = [Mycallback(model),
        EarlyStopping(monitor='val_loss', patience=2, verbose=1,mode='auto'),
        ModelCheckpoint('%s/checkpoint/best_fine_tune_weights_fantom.h5'%(DPATH),monitor='val_loss',save_best_only=True,save_weights_only=True,mode='auto')]
    training_generator = generate_batch_data(train_queue,batch_size=batch_size,use_paral=False)
    validation_generator = generate_batch_data(valid_queue,batch_size=batch_size,use_paral=False)
    print 'Start training...'
    parallel_model.fit_generator(generator=training_generator,steps_per_epoch=steps_per_epoch,epochs=epochs,verbose=1,callbacks=callbacks,validation_data=validation_generator,validation_steps=validation_steps,max_queue_size=50,workers=14,use_multiprocessing=True)



def model_evaluation(batch_size):
    corrs=[]
    nb_regions = X_seq.shape[0]
    for cell_idx in cell_idx_test:
        print 'Predicting cell %d'%cell_idx
        X_test_mat = np.zeros((nb_regions,4,SEQ_LEN,1))
        X_test_gexp = np.zeros((nb_regions,NUM_TF))
        X_test_motif = np.zeros((nb_regions,NUM_TF))
        Y_test = np.zeros(nb_regions)
        for i in range(nb_regions):
            X_test_mat[i,:,:,:] = X_seq[i]
            X_test_gexp[i,:] = tf_gexp_mat[cell_idx,:]
            X_test_motif[i,:] = tf_motifscore_mat[i,:]
            Y_test[i] = readscount_mat[i,cell_idx]
        #Y_pred = model.predict([X_test_mat,X_test_gexp,X_test_motif],batch_size=batch_size)
        Y_pred = model.predict([X_test_mat,X_test_gexp],batch_size=batch_size)
        pearsonr = stats.pearsonr(Y_pred[:,0],Y_test)[0]
        corrs.append(pearsonr)
        f_out = open('%s/log/regression_fantom.log'%(DPATH),'a+')
        f_out.write('%d\t%.5f\n'%(cell_idx,pearsonr))
        f_out.close()
    print('Median pearson r: %.4f'%np.median(corrs))
    print('Mean pearson r: %.4f'%np.mean(corrs))  

        
        
if  __name__ == "__main__" :
    DPATH='/home/liuqiao/software/DeepCAGE/data/encode'
    genome = Fasta('/home/liuqiao/openness_pre/genome.fa')
    readscount_file = '%s/readscount_normalized_filtered.csv'%DPATH
    #tf_gexp_file = '%s/preprocessed/tf_gexp.csv'%DPATH
    tf_gexp_file = '%s/preprocessed/fantom_tf_gexp.csv'%DPATH
    tf_motifscore_file = '%s/preprocessed/tf_motif_score.csv'%DPATH
    readscount = pd.read_csv(readscount_file,header=0,index_col=[0],sep='\t')
    readscount = np.log(readscount+1)
    readscount = extent_region(readscount)
    tf_gexp = pd.read_csv(tf_gexp_file ,sep='\t',header=0,index_col=[0])
    tf_gexp_log = np.log(tf_gexp+1)
    tf_gexp_log = pd.DataFrame.transpose(quantile_norm(pd.DataFrame.transpose(tf_gexp_log)))
    tf_motifscore = pd.read_csv(tf_motifscore_file,sep='\t',header=0,index_col=[0])
    readscount_mat = readscount.values
    tf_gexp_mat = tf_gexp_log.values
    tf_motifscore_mat = tf_motifscore.values
    tf_motifscore_mat = np.zeros((tf_motifscore_mat.shape[0],NUM_TF))
    region_info = readscount.index

    random.seed(0)
    cell_idx_all = np.arange(tf_gexp_log.shape[0])
    random.shuffle(cell_idx_all)
    cell_idx_train = cell_idx_all[:int(tf_gexp_log.shape[0]*0.8)]
    cell_idx_test = cell_idx_all[int(tf_gexp_log.shape[0]*0.8):]
    train_data_queue,test_data_queue = [],[]

    for i in range(tf_motifscore.shape[0]):#i denotes region idx
        for j in cell_idx_train:
            train_data_queue.append((j,i))
        for j in cell_idx_test:
            test_data_queue.append((j,i))
    random.shuffle(train_data_queue)
    random.shuffle(test_data_queue)

    if True:
    #if False:
        X_seq = np.load('%s/X_seq_full.npy'%DPATH)
    else:
        pool = Pool()
        #X_seq=np.stack(pool.map(one_hot_seq,list(readscount_train.index)))
        X_seq=np.stack(pool.map(one_hot_seq,list(readscount.index)))
        pool.close()
        pool.join()
        np.save('%s/X_seq_full.npy'%DPATH,X_seq)
    #pretrain_model = load_model('../data/encode/checkpoint/pretrain.h5')
    model = model_construct()

    model = load_model('%s/checkpoint/single_model_fantom.h5'%(DPATH))
    model_evaluation(batch_size=1024)
    sys.exit()
    #for i in range(132):
    #    model.layers[i].trainable = False    
    
    nb_gpus = 4
    parallel_model = multi_gpu_model(model,gpus = nb_gpus)

    #if os.path.isfile('%s/checkpoint/best_fine_tune_weights.h5'%DPATH):
    if False:
        parallel_model.load_weights('%s/checkpoint/best_fine_tune_weights.h5'%DPATH)
    else:
        pretrain_model = load_model('%s/checkpoint/pretrain_full.h5'%DPATH)
        for i in range(132):
            model.layers[i].set_weights(pretrain_model.layers[i].get_weights())
            #model.layers[i].trainable = False

    optimizer = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0000)
    ##parallel_model.compile(optimizer=optimizer,loss='mean_squared_error',metrics=['mse'])
    ##model_training(train_data_queue,batch_size=nb_gpus*128,epochs=100)
    #parallel_model.load_weights('%s/checkpoint/best_fine_tune_weights.h5'%DPATH)
    ##parallel_model.load_weights('%s/checkpoint/best_fine_tune_weights_fantom.h5'%(DPATH))
    ##single_model = parallel_model.layers[-2]
    ##single_model.save('%s/checkpoint/single_model_fantom.h5'%(DPATH))
    model = load_model('%s/checkpoint/single_model_fantom.h5'%(DPATH))
    model_evaluation(batch_size=1024)
    #model.save('%s/CNN/models/final/classification_fold%d_seed%d.h5'%(DPATH,fold_idx,int(sys.argv[4])))
    #model_evaluation(parallel_model,region_idx,padding=400)
    
    
    
    
    
    

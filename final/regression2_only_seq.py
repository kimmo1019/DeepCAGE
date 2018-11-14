#regression 2 difference from 1 is line 301 (gene quantile norm)
#remove attention module
import random,os,sys
import numpy as np
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
from scipy.stats import pearsonr
######################## Global Settings #######################
os.environ["CUDA_VISIBLE_DEVICES"] = sys.argv[1]
SEQ_LEN = 1000
NUM_TF = 402


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

class EarlyStoppingByLossVal(Callback):
    def __init__(self, monitor='loss', value=0.01, verbose=0):
        super(Callback, self).__init__()
        self.monitor = monitor
        self.value = value
        self.verbose = verbose
        self.wait = 0
        self.stopped_epoch = 0
    def on_epoch_end(self, epoch, logs={}):
        current = logs.get(self.monitor)
        if current is None:
            print("Early stopping requires %s available!" % self.monitor)
            exit()

        if current < self.value:
            if self.verbose > 0:
                print("Epoch %05d: early stopping THR" % epoch)
            self.model.stop_training = True
            
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
        self._time = time.time()
        return

    def on_epoch_end(self, epoch, logs={}):
        y_pred_val = self.model.predict(self.x_val)
        corr = pearsonr(self.y_val, y_pred_val[:,0])[0]
        print '\npearson coefficient: %s' % str(round(corr,4))
        print 'epoch training time: %d s' % (time.time()-self._time)
        if corr > self.best:
            self.best = corr
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


# Dense-Attention layer visulization
def get_att_vec(model,inputs,layer_name):
    output = [layer.output for layer in model.layers if layer.name==layer_name]
    assert len(output)==1
    func = K.function(model.input,output)
    return func(inputs)[0]

def DenseAttNet():
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
    #x1 = Multiply()([att_vec,gexp_input])
    #x1 = Multiply()([att_vec,motif_input])
    #x = Concatenate()([x,x1,motif_input])
    x = Dense(512)(x)
    x = BatchNormalization()(x)
    x = Activation('relu')(x)
    x = Dropout(0.5)(x)
    output = Dense(1,kernel_regularizer=regularizers.l2(0.01), name='output')(x)
    model  = Model(inputs=[seq_input,gexp_input,motif_input],outputs=output)    
    #model.summary()
    return model

def generate_batch_data(readscount,tf_gexp_norm,motif,train_cell_idx,training_queue,batch_size=128,padding=400,bootstrap=True):
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
                chrom,start,end = info.split(':')[0],int(info.split(':')[1].split('-')[0])-padding,int(info.split(':')[1].split('-')[1])+padding
                sequence = genome[chrom][start:end]
                motif_idx = chrom+':'+str(start)+'-'+str(end)
                X_batch_mat[i,:,:,0] = seq_to_mat(sequence)
                X_batch_gexp[i,:] = tf_gexp_norm.loc[cellid].values.reshape((1,-1))
                X_batch_motif[i,:] = motif.loc[motif_idx].values.reshape((1,-1))
                Y_batch[i] = readscount.iloc[region_idx].loc[str(cellid)]
            yield [X_batch_mat,X_batch_gexp,X_batch_motif], Y_batch            
    else:
        while True:
            for step in range(len(training_queue)/batch_size):
                present_queue = training_queue[step*batch_size:(step+1)*batch_size]
                for i in range(len(present_queue)):
                    region_idx,cellid = present_queue[i]
                    info = readscount.index[region_idx]
                    chrom,start,end = info.split(':')[0],int(info.split(':')[1].split('-')[0])-padding,int(info.split(':')[1].split('-')[1])+padding
                    sequence = genome[chrom][start:end]
                    motif_idx = chrom+':'+str(start)+'-'+str(end)
                    X_batch_mat[i,:,:,0] = seq_to_mat(sequence)
                    X_batch_gexp[i,:] = tf_gexp_norm.loc[cellid].values.reshape((1,-1))
                    X_batch_motif[i,:] = motif.loc[motif_idx].values.reshape((1,-1))
                    Y_batch[i] = readscount.iloc[region_idx].loc[str(cellid)]
                yield [X_batch_mat,X_batch_gexp,X_batch_motif], Y_batch

            
def generate_validation_data(readscount,tf_gexp_norm,motif, test_cell_idx, region_idx, padding=400):
    region_sample_idx = random.sample(region_idx,min(20000,len(region_idx)))
    X_valid_mat = np.zeros((len(region_sample_idx)*len(test_cell_idx),4,SEQ_LEN,1))
    X_valid_gexp = np.zeros((len(region_sample_idx)*len(test_cell_idx),NUM_TF))
    X_valid_motif = np.zeros((len(region_sample_idx)*len(test_cell_idx),NUM_TF))
    Y_valid = np.zeros(len(region_sample_idx)*len(test_cell_idx))
    idx = 0
    for cellid in test_cell_idx:
        for i in range(len(region_sample_idx)):
            info = readscount.index[region_sample_idx[i]]
            chrom,start,end = info.split(':')[0],int(info.split(':')[1].split('-')[0])-padding,int(info.split(':')[1].split('-')[1])+padding
            sequence = genome[chrom][start:end]
            motif_idx = chrom+':'+str(start)+'-'+str(end)
            X_valid_mat[idx,:,:,0] = seq_to_mat(sequence)
            X_valid_gexp[idx,:] = tf_gexp_norm.loc[cellid].values.reshape((1,-1))
            X_valid_motif[idx,:] = motif.loc[motif_idx].values.reshape((1,-1))
            Y_valid[idx] = readscount.iloc[region_sample_idx[i]].loc[str(cellid)]
            idx+=1
    return [X_valid_mat,X_valid_gexp,X_valid_motif], Y_valid

################################## Model Training##################################
def model_training(model,region_idx,batch_size=128,epochs=100,steps_per_epoch = 10000):
    training_queue=[]
    for idx in region_idx:
        for cellid in train_cell_idx:
            training_queue.append((idx,cellid))
    random.shuffle(training_queue)
    steps_per_epoch = len(training_queue)/batch_size
    X_valid,Y_valid = generate_validation_data(readscount,tf_gexp_norm,motif,test_cell_idx,region_idx,padding=400)
    callbacks = [roc_callback(validation_data=(X_valid,Y_valid),patience=6)]
    #callbacks = [EarlyStopping(monitor='val_loss', patience=3, verbose=1)]
    training_generator = generate_batch_data(readscount,tf_gexp_norm,motif,train_cell_idx, training_queue,batch_size=batch_size,padding=400,bootstrap=False)
    model.compile(optimizer='adam',loss='mean_squared_error',metrics=['mae'])
    model.fit_generator(generator=training_generator,steps_per_epoch = steps_per_epoch,epochs=epochs,verbose=1,callbacks=callbacks,validation_data=(X_valid,Y_valid),max_queue_size=8,use_multiprocessing=True)

################################## Model Evaluation#################################
def model_evaluation(model,region_idx,padding=400):
    for cellid in test_cell_idx:
        X_test_mat = np.zeros((len(region_idx),4,SEQ_LEN,1))
        X_test_gexp = np.zeros((len(region_idx),NUM_TF))
        X_test_motif = np.zeros((len(region_idx),NUM_TF))
        Y_test = np.zeros(len(region_idx))
        for i in range(len(region_idx)):
            info = readscount.index[region_idx[i]]
            chrom,start,end = info.split(':')[0],int(info.split(':')[1].split('-')[0])-padding,int(info.split(':')[1].split('-')[1])+padding
            sequence = genome[chrom][start:end]
            motif_idx = chrom+':'+str(start)+'-'+str(end)
            X_test_mat[i,:,:,0] = seq_to_mat(sequence)
            X_test_gexp[i,:] = tf_gexp_norm.loc[cellid].values.reshape((1,-1))
            X_test_motif[i,:] = motif.loc[motif_idx].values.reshape((1,-1))
            Y_test[i] = readscount.iloc[region_idx[i]].loc[str(cellid)]
        Y_pred = model.predict([X_test_mat,X_test_gexp,X_test_motif])
        corr = pearsonr(Y_test, Y_pred[:,0])[0]
        f_out = open('%s/CNN/outcome/final/regression2_fold%d_seed%d_%s_only_seq.log'%(DPATH,fold_idx,int(sys.argv[4]),str(ratio)),'a+')
        f_out.write('%d\t%.5f\n'%(cellid,corr))
        f_out.close()


        
        
if  __name__ == "__main__" :
    DPATH='/home/liuqiao/openness_pre/expression/data'
    genome = Fasta('/home/liuqiao/openness_pre/genome.fa')
    readscount_file='%s/processed_RNA_DNase/newdata/union.peaks.train.readscount'%DPATH
    tf_gexp_file = '%s/processed_RNA_DNase/tf_gexp.txt'%DPATH
    motifscore_file = '%s/processed_RNA_DNase/newdata/union.peaks.train.motifscores.txt' %DPATH
    tf_gexp = pd.read_csv(tf_gexp_file ,sep='\t',header=0,index_col=[0])
    tf_gexp_log = np.log(tf_gexp+1)
    tf_gexp_norm = pd.DataFrame.transpose(quantile_norm(pd.DataFrame.transpose(tf_gexp_log)))
    motif = pd.read_csv(motifscore_file,sep='\t',header=0,index_col=[0])
    readscount = pd.read_csv(readscount_file,header=0,index_col=[0],sep='\t')
    readscount = np.log(readscount+1)
    fold_idx = int(sys.argv[2])
    ratio = float(sys.argv[3])
    if ratio<0 or ratio>1:
        print 'Input ratio between 0 and 1'
        sys.exit()
    cellid_info = open('%s/processed_RNA_DNase/five_fold_cells.txt'%DPATH).readlines()[fold_idx]
    train_cell_idx = [int(cellid) for cellid in cellid_info.split('\t')[2].split(' ')]
    test_cell_idx = [int(cellid) for cellid in cellid_info.strip().split('\t')[4].split(' ')]
    random.seed(int(sys.argv[4]))
    region_idx = random.sample(list(range(readscount.shape[0])),int(readscount.shape[0]*ratio))
    model = DenseAttNet()
    model_training(model,region_idx,batch_size=128,epochs=100,steps_per_epoch = 5000)
    #t = time.strftime("%m-%d-%H", time.localtime()) 
    model.save('%s/CNN/models/final/regression2_fold%d_seed%d_%s_only_seq.h5'%(DPATH,fold_idx,int(sys.argv[4]),str(ratio)))
    model_evaluation(model,region_idx,padding=400)
    
           
        
        
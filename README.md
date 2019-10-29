# DeepCAGE
 A deep densely connected convolutional network for predicting chromatin accessibility with attention on gene expression
 
 ![model](https://github.com/kimmo1019/DeepCAGE/blob/master/model.png)
 
 DeepCAGE contains a deep densely connected convolutional network and a joint module for incorporating TF gene expression and motif score.
 
 # Requirements
- Keras==2.1.4
- TensorFlow==1.13.1
- hickle >= 2.1.0

# Installation
DeepCAGE can be downloaded by
```shell
git clone https://github.com/kimmo1019/DeepCAGE
```
Installation has been tested in a Linux/MacOS platform.

# Instructions
We provide detailed step-by-step instructions for running DeepCAGE model including data preprocessing, model training, and model test.

**Step 1**: Download raw DNase-seq and RNA-seq data


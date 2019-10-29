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

We provided `1.Download_raw_data.sh` for download RNA-seq data (.tsv) and DNase-seq data (.narrowPeak and .bam) from the ENCODE project
We pre-defined cell type ID from 1-60. After downloading the meta data from ENCODE website (`20180113.txt.gz`), one can run the following script:

```shell
bash 1.Download_raw_data.bash  -c <CELL_ID> -r -p -b
-c  CELLID: pre-defined cell ID (from 1 to 60)
-r  download RNA-seq data (.tsv)
-p  download chromatin accessible peaks from DNase-seq data (.narrowPeak)
-b  download chromatin accessible readscount from DNase-seq data (.bam)
```
one can also run ```bash 1.Download_raw_data.bash  -h``` to show the script instructions. Note that `.bam` files downloading may take time. After downloading the raw data, the raw data folder will be organized by `cell-assay-experiment-file` order. See an example of the folder tree:

```
data/
    |-- raw_data/
    |   |-- 1/
    |   |   |-- dseq/
    |   |   |   |-- ENCSR000EIE/
    |   |   |   |   |-- ENCFF953HEA.bed.gz
    |   |   |   |   |-- ENCFF983PML.bam
    |   |   |   |-- ENCSR000ELW/
    |   |   |   |   |...
    |   |   |-- rseq/
    |   |   |   |-- ENCSR000BXY/
    |   |   |   |   |-- ENCFF110IED.tsv
    |   |   |   |   |-- ENCFF219FVQ.tsv
    |   |   |   |-- ENCSR000BYH/
    |   |   |   |   |...
```











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
one can also run ```bash 1.Download_raw_data.bash  -h``` to show the script instructions. Note that `.bam` files downloading may take time. After downloading the raw data, the raw data folder will be organized by `cell-assay-experiment-file` order. Note that each experiment may contain multiple replicates. See an example of the folder tree:

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

**Step 2**: Merge multiple replicates of DNase-seq and RNA-seq data

We merge multiple replicate of RNA-seq data by taking the average expression of each gene across replicates in a cell type. As for DNase-seq data, we only keep bins that appear in more than half of the replicates with respect to a cell type. One can run the following scripts to merge relicates of both DNase-seq and RNA-seq data. Note that the referece genome (`hg19`) will be automatically downloaded.

```shell
python 2.Merge_multi_rep_data  -c <CELL_ID> 
-c  CELLID: pre-defined cell ID (from 1 to 60)
```
The merged data (`e.g. 1.TPM.tsv and 1.peak.bins.bed`) will be located in `data/processed_RNA_DNase` folder.

**Step 3**: Loci filtering and candidate regulatory regions selection

Please refer to `Supplementary Figure 1` for candidate regulatory regions selection strategy. Directly run `bash 3.0.Generate_peak_bin.sh` to generate candidate regulatory regions set (`union.peaks.bed` and `union.peaks.pad1k.bed`)


















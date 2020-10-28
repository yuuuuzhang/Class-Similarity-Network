# Class Similarity Network
We propose a novel network named Class Similarity Network, and we show its application in coding RNA (mRNA) and long non-coding RNA (lncRNA) classification.

The overview of the Class Similarity Network is illustrated in Figure 1. 
![alt text](https://github.com/yuuuuzhang/Class-Similarity-Network/blob/main/fig1.png)
Figure 1. An overview of the Class Similarity Network. Class Similarity Network contains three modules: Class Similarity Measurement module, Fully Connected module, and Decision module.


## EXPLANATION
This repository contains four folders, code, data, input and output.

### Code folder:
contains the python codes.  
```
utils_csn.py -- functions will be used in csn_RNA.py  
csn_RNA.py -- main file for mRNA and lncRNA classification  
csn_built.py -- how the Class Similarity Network built
```
### rawdata folder:
This folder contains the test dataset II in [1].

### input_files folder:
The model and files to be used for mRNA and lncRNA classification.

### output_files folder:
This folder contains the predicted results of input file.

## USAGE:
Tested on the Linux system. 

### Python Dependencies:
Based on python2 and tensorflow 1.14 
  
### Python modules:  
```
numpy  
pandas  
csv  
tensorflow
keras
Bio
math
argparse
```
### Installation and run:
```
git clone https://github.com/yuuuuzhang/Class-Similarity-Network.git
cd Class-Similarity-Network
```

Input: require .fa or .fasta file

output: the results will re recorded in a .csv or .txt file as specified. 
```
python code/csv_RNA.py -input inputfile -output output_file_name
```
for example:
```
python code/csv_RNA.py -input input/test2_lncrna.fa -output output/test2_lncrna.csv

```

More details can be found from [1]

## REFERANCE
[1] 

## CONTACT
If you have any inqueries, please contact YU007@e.ntu.edu.sg


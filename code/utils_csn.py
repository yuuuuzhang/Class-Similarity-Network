from __future__ import division

import numpy as np
import pandas as pd
import math
import csv
from Bio import SeqIO
import tensorflow as tf
from keras.models import load_model

def countnum(seq,nuacid):
    return len([1 for i in range(len(seq)) if seq.startswith(nuacid,i)])

def construct_kmer():
    ntarr = ("A","C","G","T")
    
    kmerArray = []
    
    
    for n in range(4):
        kmerArray.append(ntarr[n])
    
    for n in range(4):
        str1 = ntarr[n]
        for m in range(4):
            str2 = str1 + ntarr[m]
            kmerArray.append(str2)
    #############################################
    for n in range(4):
        str1 = ntarr[n]
        for m in range(4):
            str2 = str1 + ntarr[m]
            for x in range(4):
                str3 = str2 + ntarr[x]
                kmerArray.append(str3)
    #############################################
#change this part for 3mer or 4mer
    for n in range(4):
        str1 = ntarr[n]
        for m in range(4):
            str2 = str1 + ntarr[m]
            for x in range(4):
                str3 = str2 + ntarr[x]
                for y in range(4):
                    str4 = str3 + ntarr[y]
                    kmerArray.append(str4)
    ############################################
    for n in range(4):
        str1 = ntarr[n]
        for m in range(4):
            str2 = str1 + ntarr[m]
            for x in range(4):
                str3 = str2 + ntarr[x]
                for y in range(4):
                    str4 = str3 + ntarr[y]
                    for z in range(4):
                        str5 = str4 + ntarr[z]
                        kmerArray.append(str5)
    ####################### 6-mer ##############
    for n in range(4):
        str1 = ntarr[n]
        for m in range(4):
            str2 = str1 + ntarr[m]
            for x in range(4):
                str3 = str2 + ntarr[x]
                for y in range(4):
                    str4 = str3 + ntarr[y]
                    for z in range(4):
                        str5 = str4 + ntarr[z]
                        for t in range(4):
                            str6 = str5 + ntarr[t]
                            kmerArray.append(str6)

    return kmerArray

def kmer_encode1(seq,kmerarray):
    result = np.zeros((len(seq),len(kmerarray)))
    for i in range(len(seq)):
        for j in range(len(kmerarray)):
            result[i,j] = seq[i].count(kmerarray[j])/(len(seq[i]))
    return result

def kmer_encode3(seq,kmerarray):
    result = np.zeros((len(seq),len(kmerarray)))
    for i in range(len(seq)):
        for j in range(len(kmerarray)):
            result[i,j] = seq[i].count(kmerarray[j])/(len(seq[i])-2)
    return result

def hexamer_sin(seq,nc_m,c_m,kmerarray,x):
    l = len(seq)-x+1
    log_r = np.zeros((l))
    for i in range(l):
        tempseq = seq[i:i+x]
        idx = kmerarray.index(tempseq)
        Fc = c_m[int(idx)]
        Fnc = nc_m[int(idx)]
        if Fc==0 and Fnc==0:
            log_r[i]=0
        elif Fc==0 and Fnc!=0:
            log_r[i]=-1
        elif Fnc==0 and Fc!=0:
            log_r[i]=1
        else:
            log_r[i] = math.log(Fc/Fnc)
    miu = sum(log_r)/l
    
    return miu

def hexamer_score(seq,nc_m,c_m,kmerarray,x):
    miu = np.zeros((len(seq)))
    for i in range(len(seq)):
        if len(seq[i])>6:
            miu[i] = hexamer_sin(seq[i],nc_m,c_m,kmerarray,x)
    miu0 = np.expand_dims(miu, axis=1)
    return miu0

def orf_single(seq):
    startss = []
    stopss = []
    starts_ = []
    stop_s = []
    start = 0
    sslen = 0
    s_len1 = 0
    s_len2 = 0
    newseq = seq
    newseq_6 = []
    max_l = len(seq)
    l = len(seq)
    for i in range(len(seq)):
        if (seq[i:i+3]=="ATG"):
            start = 1 # has start codon
            for j in range(int((len(seq)-(i+3))/3)):
                if (seq[i+3+3*j:i+3+3*j+3]=="TAA") or (seq[i+3+3*j:i+3+3*j+3]=="TAG") or (seq[i+3+3*j:i+3+3*j+3]=="TGA"):
                    startss.append(i)
                    stopss.append(i+3+3*j+3)
                    break
            if len(startss)==0 :
                starts_.append(i)

    if start == 0:
        for k in range(len(seq)):
            if (seq[k:k+3]=="TAA") or (seq[k:k+3]=="TAG") or (seq[k:k+3]=="TGA"):
                stop_s.append(k+3)

    if len(startss)!=0:
        startss = np.array(startss)
        stopss = np.array(stopss)
        coding_len = stopss-startss
        max_len_position = np.argmax(coding_len)
        sslen = coding_len[max_len_position]
        newseq = seq[(startss[max_len_position]):(stopss[max_len_position])]
        max_l = sslen
        if (startss[max_len_position]-3)>=0 and (startss[max_len_position]+5)<l:
            newseq_6 = seq[(startss[max_len_position]-3):       (startss[max_len_position])]+seq[(startss[max_len_position]+3):(startss[max_len_position]+6)]

    elif len(starts_)!=0:
        starts_ = np.array(starts_)
        s_len1 = len(seq)-starts_[0]
        newseq = seq[(starts_[0]):len(seq)]
        max_l = s_len1
        if (starts_[0]-3)>=0 and (starts_[0]+5)<l:
            newseq_6 = seq[(starts_[0]-3):(starts_[0])]+seq[(starts_[0]+3):(starts_[0]+6)]

    elif len(stop_s)!=0:
        stop_s = np.array(stop_s)
        s_len1 = stop_s[-1]
        newseq = seq[0:(stop_s[-1])]
        max_l = s_len1
    
    orf_feature = (sslen/len(seq),s_len1/len(seq))
    
    return orf_feature,max_l,newseq,newseq_6

def orf_feature(seq):
    orf = []
    max_l = []
    newseq = []
    newseq_nu6 = []
    for i in range(len(seq)):
        orfsin,max_lsin,newseqsin,newseq_nu6sin = orf_single(seq[i])
        orf.append(orfsin)
        max_l.append(max_lsin)
        newseq.append(newseqsin)
        newseq_nu6.append(newseq_nu6sin)
    orf = np.array(orf)
    max_l = np.array(max_l)
    return orf,max_l,newseq,newseq_nu6

def nucleic_biasin(seq,arrlnc,arrpc):
    if len(seq)==6:
        bia = np.zeros((6))
        for i in range(6):
            if seq[i]=='A':
                bia[i]=math.log(arrpc[i,0]/arrlnc[i,0])
            elif seq[i]=='C':     
                bia[i]=math.log(arrpc[i,1]/arrlnc[i,1])
            elif seq[i]=='G':      
                bia[i]=math.log(arrpc[i,2]/arrlnc[i,2])
            elif seq[i]=='T':        
                bia[i]=math.log(arrpc[i,3]/arrlnc[i,3])
            else:
                print ('contain letter other than ACGT')
        bia0 = sum(bia)/6
    else:
        bia0=0
    return bia0
      
def nucleic_bia(seq,arrlnc,arrpc):
    bia=np.zeros((len(seq)))
    for i in range(len(seq)):
        bia[i] = nucleic_biasin(seq[i],arrlnc,arrpc)
    bia0 = np.expand_dims(bia, axis=1)
    return bia0

def encode(RNA_seq):
    # all
    
    maxl=107976
    lnc_arr=np.array([[0.26794686, 0.2211678 , 0.27436642, 0.23651892],
                      [0.33744557, 0.19448476, 0.25231662, 0.21575304],
                      [0.28714972, 0.26750028, 0.29887239, 0.14647762],
                      [0.32438316, 0.21430166, 0.28553087, 0.1757843 ],
                      [0.30328235, 0.23908675, 0.22323323, 0.23439768],
                      [0.24159875, 0.21809758, 0.29178296, 0.24852071]])
    pc_arr= np.array([[0.43128673 ,0.12353085, 0.35409525, 0.09108717],
                      [0.29952253, 0.34491308, 0.21244491, 0.14311949],
                      [0.20277302, 0.42718536, 0.28121939, 0.08882223],
                      [0.23105411, 0.15548482, 0.46709721, 0.14636386],
                      [0.28281097, 0.37524486, 0.17911361, 0.16283056],
                      [0.17296156, 0.24978575, 0.35489104, 0.22236166]])
    
    kmerArray = construct_kmer()
    
    RNA_orf,RNA_maxlen0,RNA_seqnew,nuc6 = orf_feature(RNA_seq)
    
    RNA_maxlen1 = RNA_maxlen0/maxl
    RNA_maxlen = np.expand_dims(RNA_maxlen1, axis=1)
    
    RNA_kmer1 = kmer_encode1(RNA_seq,kmerArray[0:4])
    RNA_kmer3 = kmer_encode3(RNA_seqnew,kmerArray[20:84])
    
    nucbia = nucleic_bia(nuc6,lnc_arr,pc_arr)
    
    nc_m = pd.read_csv('input/mer16_lncrnamean.csv',header=None,delimiter = ',')
    c_m = pd.read_csv('input/mer16_mrnamean.csv',header=None,delimiter = ',')
    pos = np.array(c_m)
    neg = np.array(nc_m)
    hexamer1 = hexamer_score(RNA_seq,pos[0:4],neg[0:4],kmerArray[0:4],1)
    hexamer2 = hexamer_score(RNA_seq,pos[4:20],neg[4:20],kmerArray[4:20],2)
    hexamer3 = hexamer_score(RNA_seq,pos[20:84],neg[20:84],kmerArray[20:84],3)
    hexamer4 = hexamer_score(RNA_seq,pos[84:340],neg[84:340],kmerArray[84:340],4)
    hexamer5 = hexamer_score(RNA_seq,pos[340:1364],neg[340:1364],kmerArray[340:1364],5)
    hexamer6 = hexamer_score(RNA_seq,pos[1364:],neg[1364:],kmerArray[1364:5460],6)

    nc_m1 = pd.read_csv('input/inmer16_lncrnamean.csv',header=None,delimiter = ',')
    c_m1 = pd.read_csv('input/inmer16_mrnamean.csv',header=None,delimiter = ',')
    pos1 = np.array(c_m1)
    neg1 = np.array(nc_m1)
    inhexamer1 = hexamer_score(RNA_seqnew,pos1[0:4],neg1[0:4],kmerArray[0:4],1)
    inhexamer2 = hexamer_score(RNA_seqnew,pos1[4:20],neg1[4:20],kmerArray[4:20],2)
    inhexamer3 = hexamer_score(RNA_seqnew,pos1[20:84],neg1[20:84],kmerArray[20:84],3)
    inhexamer4 = hexamer_score(RNA_seqnew,pos1[84:340],neg1[84:340],kmerArray[84:340],4)
    inhexamer5 = hexamer_score(RNA_seqnew,pos1[340:1364],neg1[340:1364],kmerArray[340:1364],5)
    inhexamer6 = hexamer_score(RNA_seqnew,pos1[1364:],neg1[1364:],kmerArray[1364:5460],6)
    
    
    RNA_data = np.concatenate((RNA_orf,RNA_maxlen,nucbia,RNA_kmer1,RNA_kmer3,hexamer1,hexamer2,hexamer3,hexamer4,hexamer5,hexamer6,inhexamer1,inhexamer2,inhexamer3,inhexamer4,inhexamer5,inhexamer6),axis=1) #
    return RNA_data


def test_model(input_name,output_name):
    records = list(SeqIO.parse(str(input_name), "fasta"))
    seq = []
    RNA_id = []
    for j in range(len(records)):
        seq.append(records[j].seq)
        RNA_id.append(records[j].id)
    print ('start encoding features')
    fea = encode(seq)
    print ('feature encoding finish')
    xnm=pd.read_csv('input/xnm.csv',header=None,delimiter = ',')
    xpm=pd.read_csv('input/xpm.csv',header=None,delimiter = ',')
    xnm = np.array(xnm)
    xpm = np.array(xpm)
    xpm = np.reshape(xpm,[1,len(xpm),1])
    xnm = np.reshape(xnm,[1,len(xnm),1])
    model_name=['input/m1.h5','input/m2.h5','input/m3.h5']

    X_test = fea.reshape(len(fea),84,1)
    Xp_t = np.repeat(xpm, len(X_test), axis=0)
    Xn_t = np.repeat(xnm, len(X_test), axis=0)
    
    bs=256
    probas = np.zeros((len(X_test),2))
    print ('start prediction')
    for filename in model_name:
        model = load_model(filename)
        pre_rst = model.predict([X_test,Xp_t,Xn_t]) # predict predict_classes
        probas = probas + pre_rst
        print ("finish model")
        ypredict = probas.argmax(axis=-1)
        rnaclass=[]
        for i in range(len(ypredict)):
            if ypredict[i]==0:
                rnaclass.append('mRNA')
            else:
                rnaclass.append('lncRNA')
    print ('finish prediction, wring to output...')
    with open(output_name, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerows(zip(RNA_id,rnaclass))
    
    return

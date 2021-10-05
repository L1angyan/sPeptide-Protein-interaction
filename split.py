# -*- coding:utf-8 -*-
#Find Protein Protein interactions from the 3rd sequencing data of Yeast 2 hybird
#Author:Liangyan Huazhong Agriculture University College of Plant Science and Technology

import pandas as pd
import sys

def read_fa(fa):
    fa_obj = open(fa,"r")
    sequence = {}
    while 1:
        line = fa_obj.readline().strip()
        if line == "":break
        if ">" in line:name = line[1:]
        sequence[name] = fa_obj.readline().strip()
    #将fa文件的序列读入一个字典，序列名为索引
    fa_obj.close()
    return(sequence)

def filter_txt(txt):
    df_in = pd.read_table(txt,sep="\t",header=None)
    df_in = df_in.iloc[:,[0,1,6,7]]
    df_out = pd.DataFrame()
    for i in df_in.groupby(0):
        df_read = i[1]
        l = list(df_read.iloc[:,1])
        if "AD" in l and "BD" in l and "aatL" in l and len(l) ==3:
            df_out = df_out.append(df_read.sort_values(1))
    return(df_out)

def split(df_read,seq_dic):
    df_read["mid_pos"] = (df_read.iloc[:,2]+df_read.iloc[0,3])/2
    read = df_read.iloc[0,0]
    ad , bd = df_read.iloc[0,4], df_read.iloc[1,4]
    attl_start, attl_end = df_read.iloc[2,2], df_read.iloc[2,3]
    if ad < bd:
        ad_seq = seq_dic[read][:attl_start]
        bd_seq = seq_dic[read][attl_end:]
    else:
        bd_seq = seq_dic[read][:attl_start]
        ad_seq = seq_dic[read][attl_end:]
    return(ad_seq,bd_seq)


fa = sys.argv[1]
txt = sys.argv[2]
head1 = sys.argv[3]+"_ad.fa"
head2 = sys.argv[3]+"_bd.fa"

seq_dic = read_fa(fa)
df = filter_txt(txt)
ad = open(head1,"w")
bd = open(head2,"w")
for i in df.groupby(0):
    read = i[0]
    df_read = i[1]
    seq = split(df_read,seq_dic)
    ad.write(">"+read+"\n"+seq[0]+"\n")
    bd.write(">"+read+"\n"+seq[1]+"\n")

ad.close()
bd.close()

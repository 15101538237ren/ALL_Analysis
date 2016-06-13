# -*- coding:utf-8 -*-
import csv
import random
import math
import os
import itertools
import datetime
from collections import namedtuple
import Queue
def read_csv(file, columns, type_name="Row"):
	try:
		row_type = namedtuple(type_name, columns)
	except ValueError:
		row_type = tuple
	rows = iter(csv.reader(open(file, "rb")))
	header = rows.next()
	mapping = [header.index(x) for x in columns]
	for row in rows:
		row = row_type(*[row[i] for i in mapping])
		yield row

#extract the corresponding mRNA according to its methylation gene name

def extract_mRNA_and_methylation_expression(in_m_RNA_file_path,in_methylation_file_path,out_dir_path,mr_Header=None,me_Header=None,mr_Sample=None,me_Sample=None):
    rna_lines = csv.reader(open(in_m_RNA_file_path, "rb"))
    mrna_dataset = list(rna_lines)
    rna_sample_id=mrna_dataset[0][0]
    rna_gene_names=[]
    rna_sample_names=[]
    methy_gene_names=[]
    methy_sample_names=[]
    if mr_Header:
        rna_gene_names=mrna_dataset[0][1:len(mrna_dataset[0])]
        del mrna_dataset[0]
    if mr_Sample:
        for row in read_csv(in_m_RNA_file_path,rna_sample_id.split()):
            rna_sample_names.append(row)
    methy_lines = csv.reader(open(in_methylation_file_path, "rb"))
    methy_dataset = list(methy_lines)
    methy_sample_id=methy_dataset[0][0]
    if me_Header:
        methy_gene_names=methy_dataset[0][1:len(methy_dataset[0])-1]
        del methy_dataset[0]
    if me_Sample:
        for row in read_csv(in_methylation_file_path,methy_sample_id.split()):
            methy_sample_names.append(row)
    mRNAqueue=Queue.Queue(maxsize = len(rna_gene_names)+10)
    methyqueue=Queue.Queue(maxsize=len(methy_gene_names)+10)

    for mrna_gene in rna_gene_names:
        mRNAqueue.put(mrna_gene)
    for methy_gene in methy_gene_names:
        methyqueue.put(methy_gene)
    mRNAtop=mRNAqueue.get()
    methytop=methyqueue.get()

    rna_index=0
    methy_index=0
    #匹配的methy,mRNA pair
    methy_mRNA_pair=[]
    while (not mRNAqueue.empty() and not methyqueue.empty()):
        #甲基化的基因名:ADAD2 > RNA基因名:A1CF,此时应将rna队列出队,更新index再次比较
        if methytop>mRNAtop:
            mRNAtop=mRNAqueue.get()
            rna_index=rna_index+1
            print str(mRNAqueue.qsize())
        #已经走过了直接对methy_queue出队即可
        elif methytop<mRNAtop:
            methytop=methyqueue.get()
            methy_index=methy_index+1
        elif methytop==mRNAtop:
            methy_mRNA_pair.append((methy_index,rna_index))
            methytop=methyqueue.get()
            methy_index=methy_index+1
            mRNAtop=mRNAqueue.get()
            rna_index=rna_index+1
    #print len(methy_mRNA_pair)
    hash_of_common_gene={}
    for methy_index,rna_index in methy_mRNA_pair:
        gene_name=rna_gene_names[rna_index]
        methy_expr=[]
        for i in range(len(methy_dataset)):
            methy_expr.append(methy_dataset[i][methy_index+1])
        mrna_expr=[]
        for i in range(len(mrna_dataset)):
            mrna_expr.append(mrna_dataset[i][rna_index+1])
        hash_of_common_gene[gene_name]=[methy_expr,mrna_expr]
    for (gene_name,expr_list) in hash_of_common_gene.items():
        methy_expr_list=expr_list[0]
        mrna_expr_list=expr_list[1]
        out_file_path=out_dir_path+os.sep+gene_name+".csv"
        out_file=open(out_file_path,"w")
        header_to_write="sample_name,methy_expr,mrna_expr\n"
        out_file.write(header_to_write)
        for i in range(len(methy_expr_list)):
            str_to_wrt=methy_sample_names[i][0]+","+str(methy_expr_list[i])+","+str(mrna_expr_list[i])+"\n"
            out_file.write(str_to_wrt)
        out_file.close()
if __name__=='__main__':
    starttime = datetime.datetime.now()

    in_m_RNA_file_path="mRNA_statistics2.csv"

    in_methylation_file_path="methy_gene_list_ALL_and_Normal_sorted.csv"

    out_dir_path="mrna_methy_compare"
    #抽取出甲基化文件和mRNA表达数据中共有的基因,行分别为:methylation,mRNA的表达水平,列为基因名,或相反,输出到out_dir_path
    extract_mRNA_and_methylation_expression(in_m_RNA_file_path,in_methylation_file_path,out_dir_path,True,True,True,True)
    endtime = datetime.datetime.now()
    print "running in "+str((endtime - starttime).seconds)+" seconds\n"
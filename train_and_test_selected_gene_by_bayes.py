# -*- coding:utf-8 -*-
#  Example of Naive Bayes implemented from Scratch in Python
import csv
import random
import math
import os
import itertools
import datetime
from Naive_Bayes import *
if __name__=='__main__':
    starttime = datetime.datetime.now()
    trainingtimes=20
    n_fold=10
    splitRatio = 1.0-1.0/n_fold
    input_dir="data"+os.sep+"methylation_data"
    #读取基因列表的base文件夹
    target_dir="target_gene_statistics"
    out_dir="train_test_result"
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    #每有一个methylation文件
    for filename in os.listdir(input_dir):
        if filename.startswith("."):
            continue
        filepath=os.path.join(input_dir,filename)
        if os.path.isfile(filepath):
            #甲基化原始数据文件
            dataset,gene_names = loadCsv(filepath,True,True)

            #基因列表文件
            out_gene_target_list_file_name=target_dir+os.sep+filename[0:len(filename)-4]+"gene_target"+os.sep+"overall_genes.csv"
            out_gene_list=open(out_gene_target_list_file_name,'r')

            out_file=open(out_dir+os.sep+filename,'w')
            line=out_gene_list.readline()
            while line:
                gene_split=line.split(",")
                gene_split=gene_split[0:len(gene_split)-1]
                #放组合基因的index的
                combination=[]
                #存放组合基因的甲基化数据的
                dataset_with_combination=[]
                for gene_name in gene_split:
                    gene_index=gene_names.index(gene_name)
                    combination.append(gene_index)
                #获取单行的列数
                column_len=len(dataset[0])
                #迭代行
                for j in range(0,len(dataset)):
                    #追加空list
                    dataset_with_combination.append([])
                    #迭代列
                    for k in range(0,column_len):
                        if k in combination:
                            dataset_with_combination[j].append(dataset[j][k])
                    dataset_with_combination[j].append(dataset[j][len(dataset[0])-1])


                evaluate_pridiction_sum=0.0
                for train_time in range(1,21):
                    evaluate_pridiction=0.0
                    i=1
                    while i<=n_fold:
                        trainingSet, testSet = splitDataset(dataset_with_combination, splitRatio)
                        # prepare model
                        summaries = summarizeByClass(trainingSet)
                        # test model
                        predictions = getPredictions(summaries, testSet)
                        accuracy = getAccuracy(testSet, predictions)
                        i=i+1
                        evaluate_pridiction=evaluate_pridiction+(1.0/n_fold)*(accuracy/100.0)
                    evaluate_pridiction_sum=evaluate_pridiction_sum+(1.0/trainingtimes)*evaluate_pridiction
                genes=(",").join(gene_split)
                print "file %s, genes:%s ,rate:%f" % (filename,genes,evaluate_pridiction_sum)
                line_to_write=genes+","+str(evaluate_pridiction_sum)+"\n"
                out_file.write(line_to_write)
                line=out_gene_list.readline()
            out_file.close()
            out_gene_list.close()
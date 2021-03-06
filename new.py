# -*- coding:utf-8 -*-
#  Example of Naive Bayes implemented from Scratch in Python
import csv,os
from numpy import *
IPBP=["B0303C","B1494C","B1919C","R2326C","R2236C","R2692C","R2812C","R2791C"]
IPBM=["B0585C","B0588C","B1589C","R1767C","R2533C","R2811C"]
IMBP=["B0202C","B0782C","B1452C","B1477C","R2710C","R2682C"]
IMBM=["B1493C","B1614C","B1657C","B1077C","B0201C","B06150C","R2026C","R2031C","R2372C","R2163C"]
NORMAL=["B9N1","B9N2","B9N3","B9N4","B9N5","R9N6","R9N7","R9N8","R9N9","R9N10"]

def process_methylation_file(file_path,out_file_path,out_2):
    ABNORMAL=[]
    ABNORMAL.extend(IPBP)
    ABNORMAL.extend(IPBM)
    ABNORMAL.extend(IMBP)
    ABNORMAL.extend(IMBM)
    f=file(file_path)
    out_file=open(out_file_path,'w')
    out_2_file=open(out_2,'w')
    reader = csv.reader(f)
    hash_table={}
    index=0
    sample_name_list=[]
    #read csv and store it as a big hashset[sample_name][gene_name][peak_item_list_id]
    for line in reader:
        if index==0:
            index=index+1
            continue
        sample_name=line[0]
        if sample_name not in hash_table.keys():
            hash_table[sample_name]={}
            sample_name_list.append(sample_name)
        gene_name=line[8]
        peak_M_value=line[7]
        if gene_name not in hash_table[sample_name].keys():
            hash_table[sample_name][gene_name]=[]
        hash_table[sample_name][gene_name].append(float(peak_M_value))

    sample_name_list.sort()
    #sort the hashtable and merge the gene_list
    sorted(hash_table.iteritems(), key = lambda asd:asd[0])
    old_set=set()
    for key in hash_table.keys():
        new_set=set(hash_table[key].keys())
        old_set=new_set|old_set
    gene_list=list(old_set)
    gene_list.sort()

    #crate a table [row:samples,column:genes],calc every cell's methylation level by avg(methylation_list)
    gene_mety_hash={}
    for gene in gene_list:
        gene_mety_hash[gene]=repeat([0.0],40).tolist()
    hash_table_keys_sort=hash_table.keys()
    hash_table_keys_sort.sort()
    for index,sample in enumerate(hash_table_keys_sort):
        for key in hash_table[sample].keys():
            methylation_list=hash_table[sample][key]
            count=0
            sum=0.0
            for i in range(0,len(methylation_list)):
                sum=sum+hash_table[sample][key][i]
                count=count+1
            avg=sum/count
            gene_mety_hash[key][index]=avg
    #sorted(gene_mety_hash.iteritems(), key = lambda asd:asd[0])
    #write_to outfile the header of table [row:samples,column:genes] header:gene1,gene2,.....

    # for index,gene in enumerate(gene_mety_hash.keys()):
    #     if index==0:
    #         out_file.write("Sample_Name")
    #     if index!=len(gene_mety_hash.keys()):
    #         out_file.write(","+gene)
    # out_file.write(",class")
    #     #print gene+",",
    # #print "\n",
    # out_file.write("\n")

    #write_to outfile the content of table [row:samples,column:genes] content:sample_name,value1,value2...

    len_keys=len(gene_mety_hash.keys())
    for i in range(0,len(hash_table.keys())):
        #out_file.write(sample_name_list[i])

        # out_2_file.write(sample_name_list[i]+",")
        # print sample_name_list[i]+",",

        for index,key in enumerate(gene_mety_hash.keys()):
            #if i ==0:
                #print key+",",
                out_file.write(str(gene_mety_hash[key][i])+",")

        class_normal = "Normal"
        if sample_name_list[i] in ABNORMAL:
            class_normal="Cancer"
        out_file.write(class_normal)#+str(class_normal))
        out_file.write("\n")
    #print "\n",
    # out_2_file.write("\n")
    f.close()
    out_file.close()
    out_2_file.close()
#将ALL和正常组甲基化的样本按照类别将样本进行划分,并将其输出到out_dir_类别.csv中
def all_vs_normal_split_by_type(input_file_path,out_dir):
    input_file=open(input_file_path, "r")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    ipbp_out_file=open(out_dir+os.sep+"i+b+_VS_normal.csv","w")
    ipbm_out_file=open(out_dir+os.sep+"i+b-_VS_normal.csv","w")
    imbp_out_file=open(out_dir+os.sep+"i-b+_VS_normal.csv","w")
    imbm_out_file=open(out_dir+os.sep+"i-b-_VS_normal.csv","w")
    list_file=[ipbp_out_file,ipbm_out_file,imbp_out_file,imbm_out_file]
    list_type=[IPBP,IPBM,IMBP,IMBM]
    line=input_file.readline()
    if line:
        for file in list_file:
            file.write(line)
    line=input_file.readline()
    while line:
        contents=line.split(",")
        type_name=contents[0].strip('\"').strip('\'')
        for index,typelist in enumerate(list_type):
            if type_name in typelist:
                list_file[index].write(line)
        if type_name in NORMAL:
            for file in list_file:
                file.write(line)
        line=input_file.readline()

    input_file.close()

    for file in list_file:
        file.close()
def reorder_samples_by_all_types(input_file_path,out_put_file_path):

    output_file=open(out_put_file_path,'w')

    list_type=[IPBP,IPBM,IMBP,IMBM,NORMAL]
    list_type_name=["IPBP","IPBM","IMBP","IMBM","NORMAL"]
    header=True
    for index,type in enumerate(list_type):
        input_file=open(input_file_path, "r")
        line=input_file.readline()
        if header:
            output_file.write(line)
            header=False
        line=input_file.readline()
        while line:
            contents=line.split(",")
            type_name=contents[0].strip('\"').strip('\'')

            if type_name in type:
                contents[len(contents)-1]=list_type_name[index]
                line_to_write=(",").join(contents)+"\n"
                output_file.write(line_to_write)
            line=input_file.readline()
        input_file.close()
if __name__ == '__main__':
    input_file="data"+os.sep+"methylation_data"+os.sep+"methy_gene_list_ALL_and_Normal.csv"
    out_dir="split_all_and_normal_by_type"
    #all_vs_normal_split_by_type(input_file,out_dir)
    reorder_samples_by_all_types(input_file,out_dir+os.sep+"reordered_all_and_normal.csv")
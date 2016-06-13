import csv
import random
import math
import os
import itertools

def statistic_gene_for_combination(input_dir,output_file_dir,hash_gene_names):
    if not os.path.exists(output_file_dir):
        os.makedirs(output_file_dir)
    gene_names=hash_gene_names
    global_static_hash={}
    for i in range(2,15):
        global_static_hash[i]={}
        for gene_name in gene_names:
             global_static_hash[i][gene_name]=0
    for filename in os.listdir(input_dir):
        if os.path.isfile(os.path.join(input_dir,filename)):
            input_file=open(input_dir+os.sep+filename, "r")
            output_file=open(output_file_dir+os.sep+filename,'w')
            line_cnt=1
            line=input_file.readline()
            statistic_hash={}
            while line:
                line=line[0:len(line)-1]
                chars=line.split(",")
                length=len(chars)-1
                if length not in statistic_hash.keys():
                    statistic_hash[length]={}
                    for gene_name in gene_names:
                        statistic_hash[length][gene_name]=0
                for index,val in enumerate(chars):
                    if index!=0:
                        statistic_hash[length][val]=statistic_hash[length][val]+1
                        global_static_hash[length][val]=global_static_hash[length][val]+1
                line=input_file.readline()
                line_cnt=line_cnt+1
            input_file.close()
        sorted(statistic_hash.items(), key=lambda d:d[0])

        # for (key,val) in statistic_hash.items():
        #     sorted(statistic_hash[key].items(), key=lambda d:d[0])
        gene_name_of_hash=[]
        for (key,val) in statistic_hash.items():
            gene_name_of_hash=statistic_hash[key].keys()
            break
        line=(',').join(gene_name_of_hash)
        line="len,"+line+"\n"
        output_file.write(line)

        for (key,val) in statistic_hash.items():
            val_list=[]
            for (gene_name,count) in statistic_hash[key].items():
                val_list.append(str(count))
            line=(',').join(val_list)
            line=str(key)+","+line+"\n"
            output_file.write(line)
        output_file.close()
    sorted(global_static_hash.items(), key=lambda d:d[0])
    output_file=open(output_file_dir+os.sep+"overall.csv",'w')
    gene_name_of_hash=global_static_hash[4].keys()
    line=(',').join(gene_name_of_hash)
    line="len,"+line+"\n"
    output_file.write(line)

    for (key,val) in global_static_hash.items():
        val_list=[]
        for (gene_name,count) in global_static_hash[key].items():
            number=count/float(len(global_static_hash))

            val_list.append(str(round(number,3)))
        line=(',').join(val_list)
        line=str(key)+","+line+"\n"
        output_file.write(line)
    output_file.close()
def from_over_all_file_find_top_genes(input_overall_file_path,out_file_path):

    #print "now handling %s" % input_overall_file_path
    input_file=open(input_overall_file_path, "r")
    output_file=open(out_file_path,'w')
    line=input_file.readline()
    line=line[0:len(line)-1]

    gene_names=line.split(",")
    gene_name=gene_names[1:len(gene_names)]

    line=input_file.readline()
    line_cnt=1
    while line:
        #print line_cnt
        line=line[0:len(line)-1]
        line_data=line.split(",")
        length=0
        val_hash={}
        sum_calc=0.0
        for index,data in enumerate(line_data):
            if index==0:
                length=int(data)
            else:
                val=float(data)
                sum_calc=sum_calc+val
                val_hash[index-1]=val
        val_hash_sorted=sorted(val_hash.iteritems(), key=lambda d:d[1],reverse = True)
        cnt_now=0
        gene_to_append=[]
        sum_taken=0.0
        for (key,val) in val_hash_sorted:
            if cnt_now<length:
                sum_taken=sum_taken+val
                gene_to_append.append(gene_name[key])
            else:
                break
            cnt_now=cnt_now+1

        if sum_calc==0:
            continue
        percent=sum_taken/sum_calc
        gene_to_append.append(str(round(percent,3)))
        line_to_write=(",").join(gene_to_append)
        line_to_write=line_to_write+"\n"
        output_file.write(line_to_write)

        line=input_file.readline()
        line_cnt=line_cnt+1
    input_file.close()
    output_file.close()
import datetime
if __name__=='__main__':
    all_genes_file_path="data"+os.sep+"methylation_data"+os.sep+"methy_gene_list_ALL_and_Normal.csv"
    starttime = datetime.datetime.now()

    input_dir_raw="target_gene"
    output_dir_raw="target_gene_statistics"
    if not os.path.exists(output_dir_raw):
        os.makedirs(output_dir_raw)
    for filename in os.listdir(input_dir_raw):
        if os.path.isdir(os.path.join(input_dir_raw,filename)):
            file_with_gene_names=filename[0:len(filename)-11]
            all_genes_file_path="data"+os.sep+"methylation_data"+os.sep+file_with_gene_names+".csv"
            all_genes_file=open(all_genes_file_path,'r')
            line=all_genes_file.readline().split(',')
            gene_names=line[1:len(line)-1]
            gene_names.sort()
            hash_gene_names={}
            for index,name in enumerate(gene_names):
                hash_gene_names[name]=index
            input_dir=input_dir_raw+os.sep+filename
            output_dir=output_dir_raw+os.sep+filename
            if filename!="methy_gene_list_I+B+_and_I+B-gene_target":
                statistic_gene_for_combination(input_dir,output_dir,hash_gene_names)
                from_over_all_file_find_top_genes(output_dir+os.sep+"overall.csv",output_dir+os.sep+"overall_genes.csv")
                print "hand over %s" % filename
    endtime = datetime.datetime.now()
    print "running in "+str((endtime - starttime).seconds)+" seconds\n"
from pre_process import *
import os
file_path="kangcheng_methylation.csv"
out_file_name="data"+os.sep+"split_methylation_peak"+os.sep+"out.csv"
out_2="sample_name.txt"
if __name__ == '__main__':
    process_methylation_file(file_path,out_file_name,out_2)

'''
Download and unpack data.
Set up .sex file
'''
import numpy as np
import subprocess

tracklet="atlas_T08-T3753359"
file_dl_loc="/Users/jrobinson/Downloads"

# file names
diff_files="{}/{}-diff.zip".format(file_dl_loc,tracklet)
red_files="{}/{}-reduced.zip".format(file_dl_loc,tracklet)
sex_settings="./extra_files/sextractor_files/default.sex"
# sex_save_file="{}/{}.cat".format(tracklet,tracklet)
sex_file1="{}/{}_comet".format(tracklet,tracklet)
sex_file2="{}/{}_stars".format(tracklet,tracklet)

# mkdir
mkdir_cmd="mkdir {}".format(tracklet)
print(mkdir_cmd)
subprocess.Popen(mkdir_cmd,shell=True)

# unpack data
for files in [diff_files,red_files]:

    mkdir_cmd="mkdir {}/{}".format(tracklet,files.split("/")[-1].split(".")[0])
    print(mkdir_cmd)
    subprocess.Popen(mkdir_cmd,shell=True)

    unzip_cmd="unzip {} -d {}/{}".format(files, tracklet,files.split("/")[-1].split(".")[0])
    print(unzip_cmd)
    subprocess.Popen(unzip_cmd,shell=True)

    # # remove zips
    # rm_cmd="rm {}".format(files)
    # print(rm_cmd)

# set up .sex files for comet stack and star stack
for sf in [sex_file1,sex_file2]:

    fin=open(sex_settings,"r")
    fout=open("{}.sex".format(sf),"w")

    keyword1="CATALOG_NAME"
    replacement1="CATALOG_NAME\t{}.cat\n".format(sf)

    for line in fin:
        if (keyword1 in line):
            line = line.replace(line, replacement1)
        print(line)
        fout.write(line)

    fout.close()
    fin.close()

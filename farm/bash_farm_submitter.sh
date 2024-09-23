#!/bin/bash


#this code is used to create the Table.dat file needed for a farm
#out_dir is where you want the new ntp to be created
#code is the exact location of the opticalizer.py file
#in_dir is the location of the cal files to be opticlaized
#file_farm is the location of your farm with the table.dat file to be filled

#if your cal files are not of this format  AmBe_mc_cal_[0-9]+.root please edit line 16 "files"

out_dir='/scratch/trevorh/rat6_ambe/ntp/'
code="/project/6004969/trevorh/trevor_branch/optical/opticalizer.py"
in_dir='/scratch/trevorh/rat6_ambe/cal/'

files=($(ls "$in_dir" | grep -E "AmBe_mc_cal_[0-9]+.root"))

num_jobs=${#files[@]}

file_farm='/scratch/trevorh/farmer/optical/table.dat'
count=1

for file in "${files[@]}"; do
    count=$((count + 1))
    outfile=$(echo "$file" | sed 's/.root/_ntp.root/')
    
    echo "python $code $out_dir$outfile $in_dir$file" >> "$file_farm"
    
   
done

#!/bin/bash

# delete previous outputs
rm *fasta.txt
rm meme_output_smaller.csv
rm meme_outputs_smaller.zip

# create output csv file
echo "Input file,K,Motif,E Value,Time,Objective Function" >>meme_output_smaller.csv

# loop over
for target in $(ls *fasta); do
    c_=0
    # meme
    for i in {8..15}; do
        for objf_n in cd de; do
            ~/bin/meme ${target} -nostatus -objfun ${objf_n} -dna -w ${i} -text >meme_out_s_k_${i}_${c_}__${target}.txt
            ((c_++))

        done

        # many command line combinations are not eligible for output, delete these 0 size files
        find -size 0 -delete

        # this should work on 5 * 2 * 3 * 3 = 90 files
        for file_target in $(ls meme_out_s*); do
            # source file name from output filename
            file_name=$(echo "$file_target" | grep -o -E "[a-z0-9]+\.fasta")

            # motif name(this is a combined name ,i.e. all combinations of neucleotides
            # have a name, N and X are wildcard) from output file
            motif=$(grep -o -G -m 1 "MOTIF [ATGCRYKMSWBDHVNX]*" $file_target | sed -r "s/MOTIF[ ]*//g")

            # width of motif from motif metadata
            w=$(grep -o -E "width =[ ]*[0-9]+" $file_target | sed -r "s/width[ ]*=[ ]*//g")

            # time is on another line
            time_spent=$(grep "Time" $file_target | sed -r -e "s/Time[ ]*//g" -e "s/[ ]*secs.//g" | tail -n 1)

            # there is usually no score in meme output, rather an E-Value
            evalue=$(grep -o -E "E\-value = [\-]*[0-9e\.]*[\+\-]?[0-9]*" $file_target | sed -r "s/E\-value[ ]*=[ ]*//g")

            # following 1 value are collected from the command line arguments line
            obj_fun=$(grep -o -E "\-objfun[ ]*[a-z]+" $file_target | sed -r "s/\-objfun[ ]*//g")

            echo "$file_name,$w,$motif,$evalue,$time_spent,$obj_fun" >>meme_output_smaller.csv
        done

        # zip all outputs and add to archive
        find . -name "meme_out_s*" -print | zip -q meme_outputs_smaller -@

        # delete meme outputs, clear folder a little bit for next step
        find -name "meme_out_s*" -delete

        #clear screen
        clear
    done
done

echo "done job"

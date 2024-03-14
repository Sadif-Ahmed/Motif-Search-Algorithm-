#!/bin/bash

# delete previous outputs
rm *fasta.txt
rm streme_output.csv
rm streme_outputs.zip

for target in $(ls *fasta); do
	for i in {8..15}; do
		for objf_n in cd de; do
			~/bin/streme --p ${target} --verbosity 1 --objfun ${objf_n} --w $i --nmotifs 1 --alph dna_alphabet.txt --text >streme_out_k_${i}_${objf_n}__${target}.txt
		done
	done
done

# many command line combinations are not eligible for output, delete these 0 size files
find -size 0 -delete

echo "Input file,K,Motif,Score,Time,Objective Function" >>streme_output.csv

for target in $(ls streme_out_*); do
	# source file name from output filename
	filename=$(echo "$target" | grep -o -E "[a-z0-9]+\.fasta")
	# motif name(this is a combined name ,i.e. all combinations of neucleotides
	# have a name, N and X are wildcard) from output file
	motif=$(grep -o -G "1\-[ATGCRYKMSWBDHVNX]*" $target | sed -n -r "s/1\-//p")
	# width of motif from the same line
	w=$(grep -o -E "w= [0-9]+" $target | sed -r -n "s/w=[ ]*//p")
	# time is in another line
	time_spent=$(grep "FINALTIME" $target | sed -n -r -e "s/FINALTIME[ ]*://p" -e "s/[ ]*seconds//p" | tail -n 1)
	# score is in same line of motif
	score=$(grep -o -E "S= [\-]*[0-9e\.]*[\+\-]?[0-9]*" $target | sed -n -r "s/S=[ ]*//p")
	# object function from output file name
	obj_fun=$(echo "$target" | grep -o -E "_[cde]+_" | sed "s/_//g")

	echo "$filename,$w,$motif,$score,$time_spent,$obj_fun" >>streme_output.csv
done

# zip all outputs
find . -name "streme_out_*" -print | zip -q streme_outputs -@

# delete streme outputs, clear folder
find -name "streme_out_*" -delete

echo "done job"

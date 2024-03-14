for file_target in $(ls meme_out_*); do
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

    # following 4 values are collected from the command line arguments line
    obj_fun=$(grep -o -E "\-objfun[ ]*[a-z]+" $file_target | sed -r "s/\-objfun[ ]*//g")

    prior=$(grep -o -E "\-prior[ ]*[a-z]+" $file_target | sed -r "s/\-prior[ ]*//g")

    stat_test=$(grep -o -E "\-test[ ]*[a-z]+" $file_target | sed -r "s/\-test[ ]*//g")

    mod_n=$(grep -o -E "\-mod[ ]*[a-z]*" $file_target | sed -r "s/\-mod[ ]*//g")

    echo "$file_name,$w,$motif,$evalue,$time_spent,$obj_fun,$prior,$stat_test,$mod_n" >>meme_output.csv
done

#!/bin/bash

for target in /*.fasta
do
        for i in {8..15}
        do
                ~/bin/streme --p hm03.fasta --objfun de --w $i --text > streme_out_de_k$i_$target.txt
        done

        for i in {8..15}
        do
                ~/bin/streme --p hm03.fasta --objfun cd --w $i --text > streme_out_cd_k$i_$target.txt
        done

        for 
done
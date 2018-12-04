seqkit stat $1|sed -E 's/ +/\t/g'|sed '1d'|awk '{print $4}'|sed 's/,//g'|awk '{sum += $1};END {print sum}' >>stat/stat.txt
seqkit stat step2/*|sed -E 's/ +/\t/g'|sed '1d'|awk '{print $4}'|sed 's/,//g'|awk '{sum += $1};END {print sum}' >>stat/stat.txt
seqkit stat step3/*|sed -E 's/ +/\t/g'|sed '1d'|awk '{print $4}'|sed 's/,//g'|awk '{sum += $1};END {print sum}' >>stat/stat.txt
seqkit stat $2|sed -E 's/ +/\t/g'|sed '1d'|awk '{print $4}'|sed 's/,//g'|awk '{sum += $1};END {print sum}' >>stat/stat.txt
#head -10 step6/map_Log.final.out |tail -1|cut -f 2 -d "|"|sed -E 's/\s+//g'>>stat/stat.txt
#wc -l step6/Tag_Reads_gene.tsv|awk '{print $1}' >>stat/stat.txt
paste script/stat_anno.txt stat/stat.txt >stat/stat_anno.tsv


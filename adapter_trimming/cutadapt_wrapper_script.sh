i=$1
cutadapt -a TCGCTGTCTCTTATACACATCT $i/${i}_R1.fq.gz -o ${i}/${i}_R1_trimmed.fq
cutadapt -a TCGCTGTCTCTTATACACATCT $i/${i}_R2.fq.gz -o ${i}/${i}_R2_trimmed.fq

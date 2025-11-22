
#----------------------------------------------------------------------------------------------------
# ChIP-seq


fastp -i pk15_input_1.fq.gz -o pk15_input_1.clean.fq.gz  -I pk15_input_2.fq.gz -O pk15_input_2.clean.fq.gz 
fastp -i pk15_ip_1.fq.gz -o pk15_ip_1.clean.fq.gz  -I pk15_ip_2.fq.gz -O pk15_ip_2.clean.fq.gz


# make index
/public/home/chenlijuan/bioinfo_software/bowtie2-master/bowtie2-build  --threads 20  T2T_PK15.fa   T2T_PK15.fa
samtools  faidx  T2T_PK15.fa

#02 mapping INPUT
/public/home/chenlijuan/bioinfo_software/bowtie2-master/bowtie2   -p 5   --very-sensitive --no-mixed --no-discordant -k 10   -x  T2T_PK15.fa   -1  pk15_input_1.clean.fq.gz  -2  pk15_input_2.clean.fq.gz  | samtools sort -O bam -@ 10 -o - > pk15.input.bam

#03 maping IP
/public/home/chenlijuan/bioinfo_software/bowtie2-master/bowtie2   -p 5  --very-sensitive --no-mixed --no-discordant -k 10    -x  T2T_PK15.fa   -1  pk15_ip_1.clean.fq.gz  -2  pk15_ip_2.clean.fq.gz  | samtools sort -O bam -@ 10 -o - >  pk15.CENP.bam


macs3    callpeak   -t  pk15.CENP.bam   -c  pk15.input.bam  -g   2.5e9  -f  BAM  -n PK15  -B -q 0.01
samtools index pk15.CENP.bam pk15.CENP.bam.bai
samtools index pk15.input.bam pk15.input.bam.bai
bamCoverage -b pk15.CENP.bam -o pk15.CENP.bw
bamCoverage -b pk15.input.bam -o pk15.input.bw


#----------------------------------------------------------------------------------------------------
# TRASH


/public/home/chenlijuan/bioinfo_software/TRASH/TRASH_run.sh T2T_PK15.fa \
--o /public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/TT_HiC/centromere/TRASH



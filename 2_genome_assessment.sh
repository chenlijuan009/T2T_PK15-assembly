
#----------------------------------------------------------------------------------------------------
# Calulate QC
source ~/anaconda3/bin/activate 
conda activate merqury

asm=/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/T2T_PK15.fa
hifi_reads=/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/HiFi/jellyfish_each/HiFi.fastq.gz
meryl k=21 count output hifi_reads.meryl $hifi_reads
merqury.sh  hifi_reads.meryl  $asm  T2T_PK15



#----------------------------------------------------------------------------------------------------
#Mummer & GenomeSyn 

ref=/public/agis/yiguoqiang_group/chenlijuan/database/Sus_scrofa.fa
asm=/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/Final_asm/T2T_PK15.fa


nucmer  -g 1000 -c 90 -l 40 -t 20  -p PK15 $ref $asm
delta-filter -r -q  -l 1000 PK15.delta >PK15.delta.filter
show-coords -TrHcl PK15.delta.filter >PK15.delta.filter.coords 
mummerplot -color  -postscript -R $ref_file -p PK15.delta.filter PK15.delta.filter

GenomeSyn -t 3 -g1 $ref -g2 $asm -cf1 11.coords



#----------------------------------------------------------------------------------------------------
#BUSCO

busco -i /public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/T2T_PK15.fa \
-l /public/agis/yiguoqiang_group/chenlijuan/database/busco_downloads/mammalia_odb10 -m genome -o busco_output --offline



#----------------------------------------------------------------------------------------------------
#mapping 

hifi_reads=/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/HiFi/jellyfish_each/HiFi.fastq.gz
ont_reads=/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/ONT/merge_cell/readfq/ONT_fastqpass_file60k.fq.gz
asm=/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/T2T_PK15.fa
minimap2 -ax map-hifi -t 20 $asm $hifi_reads|samtools sort -o hifi.map.sort.bam -
samtools index hifi.map.sort.bam
minimap2 -ax map-ont -t 20 $asm $ont_reads|samtools sort -o ont.map.sort.bam -
samtools index ont.map.sort.bam

samtools flagstat hifi.map.sort.bam >hifi_flagstat_final.txt
samtools flagstat ont.map.sort.bam >ont_flagstat_final.txt

samtools coverage hifi.map.sort.bam >hifi_cov.txt
samtools coverage ont.map.sort.bam >ont_cov.txt


samtools bedcov hifi.map.sort.bam >hifi_bedcov.txt
samtools bedcov ont.map.sort.bam >ont_bedcov.txt


#----------------------------------------------------------------------------------------------------
#gap region visulization

samtools view -@ 5 -b hifi.map.sort.bam Chr3:49002919-51002919 -o Chr3_hifi_gap3.bam
samtools index Chr3_hifi_gap3.bam
bamCoverage -b Chr3_hifi_gap3.bam -o Chr3_hifi_gap3.bw




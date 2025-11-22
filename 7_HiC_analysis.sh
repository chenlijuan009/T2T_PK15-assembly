
#HiC-pro
bowtie2-build -f T2T_PK15.fa T2T_PK15.fa
python /public/home/chenlijuan/bioinfo_software/HiC_dir/HiC-Pro-3.1.0/bin/utils/digest_genome.py T2T_PK15.fa -r mboi -o asm.mboi.bed
samtools faidx T2T_PK15.fa
awk '{print $1 "\t" $2}' T2T_PK15.fa.fai > chrom.sizes

HiCreads=/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/HIC/Reads

/public/home/chenlijuan/bioinfo_software/HiC_dir/HiC-Pro-3.1.0/bin/HiC-Pro -i $HiCreads  \
-o HIC_out -c config-hicpro.txt 

/public/home/chenlijuan/bioinfo_software/HiC_dir/HiC-Pro-3.1.0/bin/HiC-Pro \
-i /public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/HIC/HIC_out/hic_results/data/ \
-s build_contact_maps -o Hic_map -c config-hicpro.txt 







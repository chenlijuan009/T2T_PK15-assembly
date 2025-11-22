



#----------------------------------------------------------------------------------------------------
#Long reads variation detection

mkdir workdir

#***01 mapping
minimap2  --secondary=no -t 20   -R  '@RG\tID:BEG\tSM:${sample}'   -ax  map-hifi reference.fa  hifi.fa.gz  | samtools view --threads 20    -bS | samtools sort --threads 8 -m 50G -o ${sample}.sort.bam 
samtools index  -@ 10 ${sample}.sort.bam


ref=/public/agis/yiguoqiang_group/chenlijuan/database/Sus_scrofa.chr18_addX.fa
asm=/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/T2T_PK15.fa

minimap2  --secondary=no -t 20 -Y --MD  -R  '@RG\tID:BEG\tSM:sample'   -ax  map-hifi $asm  sample.fastq.gz  | samtools view --threads 20   -bS | samtools sort --threads 20  -o sample.asm.sort.bam
samtools index  -@ 10 sample.sort.bam



#***02 calling
#cute SV
cuteSV --max_cluster_bias_INS  1000  --diff_ratio_merging_INS  0.9   --max_cluster_bias_DEL  1000  --diff_ratio_merging_DEL  0.5  --genotype  -t 20  \
sample.sort.bam $asm     sample.vcf  ./workdir/

#pbsv
pbsv  discover  sample.sort.bam  sample.svsig.gz 
pbsv call  --ccs   --num-threads 20   $asm   sample.svsig.gz  sample.pbsv.vcf

#sniffles
sniffles   --threads  20   --input sample.sort.bam   --vcf  sample.vcf


#***03 merge all sample SV dataset
ls >allsample_vcf_list
SURVIVOR merge allsample_vcf_list 1000 2 1 0 0 50 allsample.asm.merge.vcf
cat allsample.asm.merge.vcf |wc -l


#----------------------------------------------------------------------------------------------------
#short reads variation detection
gatk=/public/home/chenlijuan/bioinfo_software/SNP/gatk-4.2.6.1/gatk

bwa index Sus_scrofa.chr18_addX.fa
bwa index T2T_PK15.fa
bwa mem -t 10 T2T_PK15.fa ./clean/merge_f1.clean.fastq.gz ./clean/merge_r2.clean.fastq.gz |samtools sort -@ 40 -o merge_sort.bam -

samtools index -@ 10 ./bam/${sample}.sort.bam 
samtools rmdup -s  ./bam/${sample}.sort.bam  ./bam/${sample}.sort.rmdup.bam
samtools index -@  10 ./bam/${sample}.sort.rmdup.bam

$gatk   HaplotypeCaller  --native-pair-hmm-threads 20   -R  T2T_PK15.fa   -I  ./bam/${sample}.sort.rmdup.bam     -ERC GVCF  -O  ${sample}.g.vcf.gz




for i in {1..18} X; do mkdir -p ./my_database/$SLURM_JOBID/tmp; \
$gatk GenomicsDBImport --reader-threads 5  --sample-name-map ./DBI.list -R Raw_np1_np2.fa --genomicsdb-workspace-path ./my_database/${SLURM_JOB_ID}/Chr${i}_gdb --tmp-dir "./my_database/$SLURM_JOBID/tmp" --intervals Chr${i};done;








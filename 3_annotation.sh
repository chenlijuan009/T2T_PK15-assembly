#----------------------------------------------------------------------------------------------------
# repeat annotation


#RepeatModeler
BuildDatabase -name Pig_asm -engine ncbi ../T2T_PK15.fa
RepeatModeler -database Pig_asm  -pa 30 

#RepeatMasker
RepeatMasker -pa 40 -nolow  -no_is -gff -norna -e ncbi -lib Pig_asm-families.fa -dir RepeatModeler_output1 ../T2T_PK15.fa

#DeepTE
mkdir output_dir
python DeepTE.py -d working_dir -o output_dir -i Pig_asm-families_RMRBSeq.fa -sp M  -m_dir working_dir/Metazoans_model 

#----------------------------------------------------------------------------------------------------
# non-coding annotation

#snRNA, miRNA

cmsearch --cut_ga --rfam --nohmmonly  --cpu 8 --tblout rfam_out.tab --clanin /public/agis/yiguoqiang_group/chenlijuan/database/Rfam.clanin\
--rfam /public/agis/yiguoqiang_group/chenlijuan/database/Rfam15.0.cm $asm >T2T_PK15.fa.cmscan

#tRNA
tRNAscan-SE --thread 4 -E -I -m tRNA.stats -o PK15.tRNA.tsv -f PK15.tRNA.structure.txt $asm

#rRNA

barrnap --kingdom euk  -o rrna.fa < genome.softmasked.fasta >rna.gff

#----------------------------------------------------------------------------------------------------
# gene annotation

#Iso-seq
asm=/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/genome.softmasked.fasta
minimap2 -ax splice -f 1000 --sam-hit-only --secondary=no --eqx -K 100M -t 30 -uf --cap-sw-mem=3g $asm sample_subreads.fastq.gz >sample.sam
sort -k3,3 -k4,4n ${i}.sam  >sample.sorted.sam 
~/mambaforge/envs/Cupcake/bin/collapse_isoforms_by_sam.py  --input sample_subreads.fastq --fq -s sample.sorted.sam -c 0.9 -i 0.9 --dun-merge-5-shorter  --max_3_diff 500 -o sample
collapse_isoforms_by_sam.py   --input sample_subreads.fastq --fq -s sample.sorted.sam -c 0.9 -i 0.9 --dun-merge-5-shorter  --max_3_diff 500 -o sample



#RNA-seq
hisat2-build -p 10 genome.softmasked.fasta genome.softmasked.fasta
hisat2 --dta --max-intronlen 20000 --min-intronlen 20 -x genome.softmasked.fasta -p 30 -1 ./clean/sample_1.clean.fq.gz -2 ./clean/sample_2.clean.fq.gz | samtools sort -@ 40 > ./bam/sample.sort.bam
stringtie -p 10 -o ./stringtie/sample.gtf ./bam/sample.sort.bam
gffread ./stringtie/sample.gtf -o- >./stringtie/sample.gff3


#PASA
/public/home/chenlijuan/mambaforge/envs/Gene_Anno/bin/gtf_genome_to_cdna_fasta.pl merge_iso_ngs.gtf.combined.gtf $asm >transcripts.fasta
/public/home/chenlijuan/mambaforge/envs/Gene_Anno/bin/gtf_to_alignment_gff3.pl merge_iso_ngs.gtf.combined.gtf >merge_iso_ngs.gtf.combined.gff
diamond makedb --in uniprot_sprot.fasta --db uniprot_sprot.fasta
diamond blastp --evalue 1e-5 --outfmt 6 -d uniprot_sprot.fasta 
TransDecoder.LongOrfs -t transcripts.fasta
diamond blastp --evalue 1e-5 --outfmt 6 -d uniprot_sprot.fasta -q transcripts.fasta.transdecoder_dir/longest_orfs.pep --max-target-seqs 1 > blastp.result
TransDecoder.Predict -t transcripts.fasta --retain_blastp_hits blastp.result 
/public/home/chenlijuan/mambaforge/envs/Gene_Anno/bin/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 \
merge_iso_ngs.gtf.combined.gff transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3
cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3


makeblastdb -in UniVec -dbtype nucl -input_type fasta -parse_seqids -out UniVec
seqclean transcripts.fasta -c 10 -n 10000 -v /public/agis/yiguoqiang_group/chenlijuan/database/UniVec/


DATABASE=/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/TT_HiC/merge_ISO_RNA/db.sqlite



perl /public/home/chenlijuan/mambaforge/envs/pasa/opt/pasa-2.4.1/Launch_PASA_pipeline.pl -c align.conf -C -R \
-g $asm -t transcripts.fasta.clean -T -u transcripts.fasta --ALIGNER gmap --CPU 20

/public/home/chenlijuan/mambaforge/envs/pasa/opt/pasa-2.4.1/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta transcripts.fasta --pasa_transcripts_gff3  db.sqlite.valid_blat_alignments.gff3



#Augustus
augustus --strand=both --genemodel=partial  --gff3=on --alternatives-from-sampling=true \
--AUGUSTUS_CONFIG_PATH=/public/home/chenlijuan/mambaforge/envs/Augustus/config --outfile=agugstus_PK15.gff3 --species=pig $asm

/public/home/chenlijuan/mambaforge/envs/evidencemodeler/bin/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl  agugstus_complete.gff3 >agugstus_complete_EVM.gff3


#GeMoMA

ourdir=/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/TT_HiC/Gene_Annotation/step5_homo_preidct 
Sus_fa=/public/agis/yiguoqiang_group/chenlijuan/database/Sus_scrofa.fa 
Sus_gff=/public/agis/yiguoqiang_group/chenlijuan/database/Sus_scrofa.Sscrofa11.1.112.gff3
hm_fa=/public/agis/yiguoqiang_group/chenlijuan/database/pig_variety_database/Homo_sapiens.GRCh38.dna.toplevel.fa
hm_gff=/public/agis/yiguoqiang_group/chenlijuan/database/pig_variety_database/Homo_sapiens.GRCh38.113.gff3
mouse_fa=/public/agis/yiguoqiang_group/chenlijuan/database/pig_variety_database/Mus_musculus.GRCm39.dna.toplevel.fa
mouse_gff=/public/agis/yiguoqiang_group/chenlijuan/database/pig_variety_database/Mus_musculus.GRCm39.113.gff3
catle_fa=/public/agis/yiguoqiang_group/chenlijuan/database/pig_variety_database/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa
catle_gff=/public/agis/yiguoqiang_group/chenlijuan/database/pig_variety_database/Bos_taurus.ARS-UCD1.3.113.gff3
dog_fa=/public/agis/yiguoqiang_group/chenlijuan/database/pig_variety_database/Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa
dog_gff=/public/agis/yiguoqiang_group/chenlijuan/database/pig_variety_database/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gff3
sheep_fa=/public/agis/yiguoqiang_group/chenlijuan/database/pig_variety_database/Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.dna.toplevel.fa
sheep_gff=/public/agis/yiguoqiang_group/chenlijuan/database/pig_variety_database/Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.113.gff3


GeMoMa -Xmx100g GeMoMaPipeline threads=5 \
outdir=$ourdir GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=$asm \
s=own i=1 a=$Sus_gff g=$Sus_fa s=own i=2  a=$hm_gff g=$hm_fa s=own i=3  a=$mouse_gff g=$mouse_fa \
s=own i=4  a=$catle_gff g=$catle_fa \
s=own i=5  a=$dog_gff g=$dog_fa \
s=own i=6  a=$sheep_gff g=$sheep_fa

GeMoMa_gff_to_gff3.pl final_annotation.gff  > GeMoMa.gff


#RUN EVM
/public/home/chenlijuan/mambaforge/envs/evidencemodeler/bin/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl  agugstus_complete.gff3 >agugstus_complete_EVM.gff3
/public/home/chenlijuan/mambaforge/envs/evidencemodeler/opt/evidencemodeler-2.1.0/EvmUtils/misc/GeMoMa_gff_to_gff3.pl final_annotation.gff >protein_alignments.gff3
/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/TT_HiC/step4_denovo_repdict/agugstus_complete_EVM.gff3 
/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/Gen_asm/TT_HiC/step5_homo_preidct/protein_alignments.gff3



/public/home/chenlijuan/mambaforge/envs/evidencemodeler/opt/evidencemodeler-2.1.0/EvmUtils/partition_EVM_inputs.pl --genome $asm \
--gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcript_alignments.gff3 \
--segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out


/public/home/chenlijuan/mambaforge/envs/evidencemodeler/opt/evidencemodeler-2.1.0/EvmUtils/write_EVM_commands.pl --genome $asm --weights weights.txt \
 --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcript_alignments.gff3 \
 --output_file_name evm.out --partitions partitions_list.out > commands.list

recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output_file_name evm.out --genome ./*.fasta











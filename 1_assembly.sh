
# hifiasm 
/public/home/chenlijuan/bioinfo_software/hifiasm/hifiasm -o HiFi_PK15 -t 20  HiFi.fastq.gz


#NextDenovo 
/public/home/chenlijuan/bioinfo_software/assmebly_dir/NextDenovo/nextDenovo run.cfg


#***run.cfg configure file
[General]
job_type = local
job_prefix = nextDenovo
task = all # 'all', 'correct', 'assemble'
rewrite = yes # yes/no
deltmp = yes
rerun = 3
parallel_jobs = 30
input_type = raw
read_type = ont  # clr, ont, hifi
input_fofn = ./input.fofn
workdir = ./01_work

[correct_option]
read_cutoff = 6k
genome_size = 2500M
pa_correction = 50
sort_options = -m 50g -t 20
minimap2_options_raw =  -t 20
correction_options = -p 30

[assemble_option]
minimap2_options_cns = -t 10
nextgraph_options = -a 1















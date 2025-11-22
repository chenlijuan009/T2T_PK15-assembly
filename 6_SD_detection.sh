

source ~/mambaforge/bin/activate 
conda activate BISER


asm=/public/agis/yiguoqiang_group/chenlijuan/repetamasker/genome.softmasked.fasta
biser --max-error 20 --max-edit-error 10  -o pig.sd -t 1 $asm --gc-heap 32G --kmer-size 31





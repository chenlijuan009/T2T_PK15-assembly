# T2T_PK15-assembly
A gap-free telomere-to-telomere pig genome (~2.6 Gb), including all autosomes and chromosome X.  A combination of PacBio HiFi and ultra-long ONT was emplyed to generate T2T-PK15 genome.


Descriptions
1_assembly.sh script for the pig T2T genome assembly 

2_genome_assessment.sh script for annotation
calculate QV, collinearity analysis, BUSCO completeness, and mapping rate 

3_annotation.sh script for annotation 
Repetitive sequence annotation, non-coding RNA annotation, gene structure prediction, and integration

4_variation_detection.sh script for variation detection
Long reads variation detection and  short reads variation detection

5_centromere.sh script for centoremere 
ChIP-seq: ChIP-seq reads were aligned to the reference genome, followed by peak calling to identify significant enrichment regions.
centromeric repeat unit: identifying centromeric repeat units using TRASH software

6_SD_detection.sh scripts for segments duplications identification

7_HiC_analysis.sh script for HiC analysis

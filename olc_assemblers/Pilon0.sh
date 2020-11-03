# Pilon is a software tool which can be used to:
# Automatically improve draft assemblies and find variation among strains, including large event detection
# RUN IN A SCREEN - step 1 can be longer than 3 hours so run in medium, long or himem partition
# srun --partition himem --time 1-12:00:00 --mem-per-cpu 20G --cpus-per-task 24 --pty bash (long was full)
srun --partition short --time 0-03:00:00 --mem-per-cpu 20G --cpus-per-task 24 --pty bash
# First install pilon into conda env (I used olc_assemblers env)
# conda install -c bioconda pilon

# activate conda env
# conda activate olc_assemblers

# Default usage
# Usage: pilon --genome genome.fasta [--frags frags.bam] [--jumps jumps.bam] [--unpaired unpaired.bam]

# Step 0: Index your assembled genome (do for each assembly)
# bwa index assembly/flye/F.oxysporum_fsp_lactucae/race_1/flye_raw/racon_10/medaka/assembly_r4_renamed_FR.fasta index_prefix
bwa index assembly/miniasm/F.oxysporum_fsp_lactucae/race_1/racon_10/medaka/race_1_miniasm_racon_round_2_renamed.fasta

# Step 1: Map your reads (PacBio or Oxford Nanopore) back to your (consensus) assembly using BWA. This will create a mapping file in SAM format
# did in pilon dir
# v1 bwa mem –t 2 –x ont2d -p assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/racon_10/medaka/race_1_smartdenovo_racon_round_10_renamed.fasta assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz > bwa_mapping_SDen.sam
# v2 bwa mem -t 4 -p race_1_smartdenovo_racon_round_10_renamed.fasta FAL_trim.fastq.gz >bwa_mapping_SDen.sam (worked)
# bwa mem -t 4 -p assembly/flye/F.oxysporum_fsp_lactucae/race_1/flye_raw/racon_10/medaka/assembly_r4_renamed_FR.fasta assembly/flye/F.oxysporum_fsp_lactucae/race_1/flye_raw/racon_10/medaka/FAL_trim.fastq.gz >bwa_mapping_flyeR.sam
bwa mem -t 4 -p assembly/miniasm/F.oxysporum_fsp_lactucae/race_1/racon_10/medaka/race_1_miniasm_racon_round_2_renamed.fasta assembly/miniasm/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz >bwa_mapping_mini.sam

# Step 2: Convert your SAM file to the binary BAM format using samtools
samtools view -Sb bwa_mapping_SDen.sam > bwa_mapping_SDen.bam
# samtools view -Sb bwa_mapping_flyeR.sam > bwa_mapping_flyeR.bam
# samtools view -Sb bwa_mapping_mini.sam > bwa_mapping_mini.bam

# Step 3:Sort and index the BAM file for faster access using samtools
samtools sort –o bwa_mapping.SDEN.sorted.bam bwa_mapping_SDen.bam
samtools index bwa_mapping.SDEN.sorted.bam

# Step 4: Run Pilon using the genome assembly and the sorted and indexed bam file
pilon --genome race_1_smartdenovo_racon_round_10_renamed.fasta --bam bwa_mapping.SDEN.sorted.bam --threads 2 # <-<- Doesn't work: causes memory errors- "Exception in thread "main" java.lang.OutOfMemoryError: Java heap space""

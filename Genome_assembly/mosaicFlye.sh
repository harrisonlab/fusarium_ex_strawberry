# mosaicFLye
# have clone to git_repos then run in a node in a screen

# screen -r
# srun --partition long --time 0-06:00:00 --mem 40G --cpus-per-task 24 --pty bash
# run in this dir /projects/fusarium_EX_Lactucae
# Use "ctrl+d to close python2"
# Run for flye, miniasm and SMARTdenovo outputs

# mosaic --reads <file> -o <dir> (--genome-size <int> | --flye-dir <dir> | --contigs <contigs>) [--threads <int>]
# removed --threads

# Examples of how to run mosaicFLye correctly - example 2 is best otherwise it reruns flye
# python2 /path/to/mosaic --reads path/to/reads/xyz.fastq.gz -o /path/to/output 60000000
# python2 /path/to/mosaic --reads path/to/reads/xyz.fastq.gz -o /path/to/output --flye-dir path/to/flye --threads 8
# python2 /path/to/mosaic --reads path/to/reads/xyz.fastq.gz -o /path/to/output --genome-size 60000000 --flye-dir path/to/flye --contigs path/to/contigs --threads 8

python2 /home/akinya/git_repos/assembly_fusarium_ex/mosaic/ --reads assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz -o assembly/mosaicFLye/flye --flye-dir /home/akinya/miniconda3/envs/olc_assemblers/bin/ --threads 8

# Examples of how to run mosaicFLye correctly
# python2 /path/to/mosaic --reads path/to/reads/xyz.fastq.gz -o /path/to/output 60000000
# python2 /path/to/mosaic --reads path/to/reads/xyz.fastq.gz -o /path/to/output --flye-dir path/to/flye --threads 8
# python2 /path/to/mosaic --reads path/to/reads/xyz.fastq.gz -o /path/to/output --genome-size 60000000 --flye-dir path/to/flye --contigs path/to/contigs --threads 8

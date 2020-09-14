# mosaicFLye
# have clone to git_repos then run in a node in a screen

# screen -r
# srun --partition long --time 0-06:00:00 --mem 40G --cpus-per-task 24 --pty bash
# run in this dir /projects/fusarium_EX_Lactucae
# Use "ctrl+d to close python2"
# Run for flye, miniasm and SMARTdenovo outputs

# mosaic --reads <file> -o <dir> (--genome-size <int> | --flye-dir <dir> | --contigs <contigs>) [--threads <int>]

python2 /home/akinya/git_repos/assembly_fusarium_ex/mosaic/mosaic --reads assembly/flye/F.oxysporum_fsp_lactucae/race_1/FAL_trim.fastq.gz -o /assembly/mosaicFLye/flye (--genome-size 60m | --flye-dir miniconda3/envs/olc_assemblers/bin/ | --contigs assembly.fasta) [--threads 8]

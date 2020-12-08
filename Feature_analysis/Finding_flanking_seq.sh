# Use to find region upstrean and downstream of gene of interesting i.e. SIX genes
# copy six gene hit sequence from blast analysis to "Six_4_seq.fasta" & "six_6_seq.fasta"

# 1) make a BLAST database for your genome
makeblastdb -in reference.fa -dbtype nucl # makeblastdb -in pilon_10_renamed.fasta -dbtype nucl

# 2) make a FASTA index of the genome
samtools faidx reference.fa # samtools faidx pilon_10_renamed.fasta

# 3) get the chromosome sizes
cut -f 1,2 reference.fa.fai > chrom.sizes # cut -f 1,2 pilon_10_renamed.fasta.fai > chrom.sizes

# 4) perform BLAST (note use your desired settings) in conda env
blastn -db reference -query query.fasta -outfmt 6 -out query.fasta-against-reference.blast

# 5) convert BLAST 6 format to BED
cut -f 2,9,10 query.fasta-against-reference.blast | awk '{if($2 > $3) {t = $2; $2 = $3; $3 = t; print;} else if($2 < $3) print; }' OFS='\t' |awk '{a=$2-1;print $1,a,$3;}' OFS='\t'|bedtools sort > query.fasta-against-reference.blast.bed

# 6) add on 500 bases upstream and downstream
bedtools slop -i query.fasta-against-reference.blast.bed -g chrom.sizes -b 500 > query.fasta-against-reference.blast.500.bed

# 7) get the accompanying bases for BLAST coordinates but 500 bases upstream and downstream
bedtools getfasta -fi reference.fa -bed query.fasta-against-reference.blast.500.bed -fo query.fasta-against-reference.blast.500.fasta

#!/home/akinya/miniconda3/envs/RMask/bin/perl
#!/usr/bin/perl
use strict;
use warnings;

# TE_to_gff.pl
# This script will accept the contigs & coordinates from a transposonPSI output file and relate it to gene annotations in a combined .gff3 file
# and will output a merged feature track for the hit gene and cand TE, using the features supplied.

# build a file of genes that match contig and overlap coordinates in TE file from final_genes_appended_renamed.gff3

# read each line of gff file
#       split the gff feature line
#       if the first overlaps and fourth & fifth column of files lie in range
#               then store this as the "transposon_gene"
#       if the third column contains 'gene'
#               then split the 9th column to get the transcript ID
#               if the transcript ID matches an entry in the contig & coord then
#                       print the previous gene as the parent feature (in gff3 format)
#                       print array of line
#                       ...and column 2 with the specifed program...
#                       ...and column 9 containing 'transcript_id "gene_name+transcript_no"; gene_id "gene_name"; parent "previous_gene_name";'
#       repeat for next gff line
#               add 10th column for TE ID and hit no from bestPerLocus file
#                       print array of line

#[prediction source] - AUGUSTUS/CodingQuarry_v2.0
my $usage = "TE_to_gff.pl <*.fa.TPSI.allHits.chains.bestPerLocus.gff3> <final_genes_appended.gff3>  [prediction source] [transcript_id_to_search_for_(ie.ID/Name)]

my $infile = shift or die $usage;
my $gene_models = shift or die $usage;
my $prediction_source = shift or die $usage;
my $feature = shift or die $usage;
my $model_switch = shift or die $usage;
my %contig_hash
my $TE_ID
my $OldCol1
my $OldCol2
my $OldCol3
my $OldCol4
my $OldCol5
my $OldCol6
my $OldCol7
my $OldCol8
my $OldCol9
my $gene_num;
my $Prev_num;
my $PrevCol1;
my $PrevCol2;
my $PrevCol3;
my $PrevCol4;
my $PrevCol5;
my $PrevCol6;
my $PrevCol7;
my $PrevCol8;
my $PrevCol9;

# This is a switch to turn on printing subfeatures with parents
# that are in the list of desired transcripts.

my $print_subfeatures = 0;

print "##gff-version 3\n";

open (INFILE, $infile) or die "can't open: $infile\n$usage\n\n";
open (GENE_MODELS, $gene_models) or die "can't open: $gene_models\n$usage\n\n";

while (my $contig_name = <INFILE>) {
        chomp $contig_name;
        $contig_hash{$contig_name} = "True";
}

while (my $gff_line = <GENE_MODELS>) {
        chomp $gff_line;
        if ($model_switch eq "Augustus") {
                my @split_line = split ('\t', $gff_line);
                if ($split_line[0] =~ m/gff-version 3/) { next;
                } elsif ($split_line[2] eq 'gene') {
                                $PrevCol1 = $split_line[0];
                                $PrevCol2 = $prediction_source;
                                $PrevCol3 = "gene";
                                $PrevCol4 = $split_line[3];
                                $PrevCol5 = $split_line[4];
                                $PrevCol6 = $split_line[5];
                                $PrevCol7 = $split_line[6];
                                $PrevCol8 = $split_line[7];
                                $PrevCol9 = $split_line[8];
                                $print_subfeatures = 0;
                } elsif ($split_line[2] eq 'transcript') {
                        $split_line[8] =~ m/$feature=[\S]+/;
                        $transcript_ID = $&;
                        $transcript_ID =~ s/^.*?=//;
                        $transcript_ID =~ s/;.*//;
                        if ($gene_hash{$transcript_ID} && $split_line[2] eq 'transcript') {
                                $transcript_ID =~ m/\.t[\d]\b/;
                                my $gene_num = $`;
                                my $col1 = $split_line[0];
                                my $col2 = $prediction_source;
                                my $col3 = $split_line[2];
                                my $col4 =      $split_line[3];
                                my $col5 =      $split_line[4];
                                my $col6 =      $split_line[5];
                                my $col7 =      $split_line[6];
                                my $col8 =      $split_line[7];
                                my $col9 =      $split_line[8];
                                print join("\t", $PrevCol1, $PrevCol2, $PrevCol3, $PrevCol4, $PrevCol5, $PrevCol6, $PrevCol7, $PrevCol8, $PrevCol9), "\n";
                                print join("\t", $col1, $col2, $col3, $col4, $col5, $col6, $col7, $col8, $col9), "\n";
                                $print_subfeatures = 1;

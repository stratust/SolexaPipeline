#!/usr/bin/perl
#
#              INGLÊS/ENGLISH
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  http://www.gnu.org/copyleft/gpl.html
#
#
#             PORTUGUÊS/PORTUGUESE
#  Este programa é distribuído na expectativa de ser útil aos seus
#  usuários, porém NÃO TEM NENHUMA GARANTIA, EXPLÍCITAS OU IMPLÍCITAS,
#  COMERCIAIS OU DE ATENDIMENTO A UMA DETERMINADA FINALIDADE.  Consulte
#  a Licença Pública Geral GNU para maiores detalhes.
#  http://www.gnu.org/copyleft/gpl.html
#
#  Copyright (C) 2008  Fundação Hemocentro de Ribeirão Preto
#
#  Laboratório de Bioinformática
#  BiT -  Bioinformatic Team
#  Fundação Hemocentro de Ribeirão Preto
#  Rua Tenente Catão Roxo, 2501
#  Ribeirão Preto - São Paulo
#  Brasil
#  CEP 14051-140
#  Fone: 55 16 2101-9300 Ramal 9365
#
#  Thiago Yukio Kikuchi Oliveira 
#  stratus@lgmb.fmrp.usp.br
#  http://lgmb.fmrp.usp.br
#
# $Id$

=head1 NAME 

=head1 SYNOPSIS

=head1 ABSTRACT

=head1 DESCRIPTION

=head1 AUTHOR

Thiago Yukio Kikuchi Oliveira E<lt>stratus@lgmb.fmrp.usp.brE<gt>

Copyright (c) 2006 Regional Blood Center of Ribeirao Preto

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

package Main;
use strict;
use warnings;
use Getopt::Long;
use Moose;
use FindBin qw($Bin);                # where was script installed?
use lib "$Bin/lib";      # use that dir for libs, too
use ProcessFastq;
use File::Basename;


#GetOptions variables
my ( $barcode, $file1, $file2, $verbose, $reference_genome_path, $np );

Usage("Too few arguments") if $#ARGV<1;
GetOptions( 
    "h|?|help" => sub { Usage(); },
    "barcode=s" => \$barcode,
    "f1=s" => \$file1,
    "f2=s" => \$file2,
    "refpath=s" => \$reference_genome_path,
    "np=s" => \$np,
    "verbose" => \$verbose,

    
)or Usage();

# Presetting;
$barcode = 'CGCGCCT' unless $barcode;
$np = 2 unless $np;


# STEP1
# Runing Filter
print "--------------------------------\n";
print "STEP1\n";
print "--------------------------------\n";
print "Generating Quality informations of the sequence...\n";
my $sequence_info_output = "sequence_info"; # Where the sequence info will be send

# creating directory if not exist
mkdir $sequence_info_output unless (-e $sequence_info_output);
system( "fastx_quality_stats -i $file1 -o $sequence_info_output/" . basename($file1) . ".stats" );
system( "fastx_quality_stats -i $file2 -o $sequence_info_output/" . basename($file2) . ".stats" );

print "Generating Quality Graph...\n";
system( "fastq_quality_boxplot_graph.sh -i $sequence_info_output/".basename($file1).".stats -o $sequence_info_output/".basename($file1).".png -t '".basename($file1)."'" );
system( "fastq_quality_boxplot_graph.sh -i $sequence_info_output/".basename($file2).".stats -o $sequence_info_output/".basename($file2).".png -t '".basename($file2)."'" );

print "Filtering and Spliting Barcode Sequences...\n\n";
# Creating object ProcessFastq;
my $barcode_filter_path = "STEP1-barcode_filtered"; # Where the filtered files will be send
my $self = ProcessFastq->new(file1 => $file1, file2=> $file2, barcode => $barcode, output_path => $barcode_filter_path);
$self->filter_fastq();

print "Generating Quality informations of the Filtered sequences...\n";
# creating directory if not exist
mkdir "$barcode_filter_path/$sequence_info_output" unless (-e "$barcode_filter_path/$sequence_info_output");
system( "fastx_quality_stats -i $barcode_filter_path/".basename($file1).".barcode_matched -o $barcode_filter_path/$sequence_info_output/" . basename($file1) . ".barcode_matched.stats" );
system( "fastx_quality_stats -i $barcode_filter_path/".basename($file2).".barcode_matched -o $barcode_filter_path/$sequence_info_output/" . basename($file2) . ".barcode_matched.stats" );

print "Generating Quality Graph...\n";
system( "fastq_quality_boxplot_graph.sh -i $barcode_filter_path/$sequence_info_output/".basename($file1).".barcode_matched.stats -o $barcode_filter_path/$sequence_info_output/".basename($file1).".png -t '".basename($file1).".barcode_matched'" );
system( "fastq_quality_boxplot_graph.sh -i $barcode_filter_path/$sequence_info_output/".basename($file2).".barcode_matched.stats -o $barcode_filter_path/$sequence_info_output/".basename($file2).".png -t '".basename($file2).".barcode_matched'" );



#STEP2 Alignment Paired-end
print "--------------------------------\n";
print "STEP2\n";
print "--------------------------------\n";

print "Align Paired-end Sequences...\n\n";
my $paired_alignment_output = "STEP2-paired_alignment"; # Where the paired alignments will be

# creating directory if not exist
mkdir $paired_alignment_output unless (-e $paired_alignment_output);

system("time bowtie -S --best --fr -p $np -X 1000 $reference_genome_path -1 $barcode_filter_path/".basename($file1).".barcode_matched -2 $barcode_filter_path/".basename($file2).".barcode_matched --un $paired_alignment_output/unaligned.fq --al $paired_alignment_output/aligned.fq $paired_alignment_output/valid_paired_alignments.sam");



#STEP4 Align Unmatched sequences as single-ends
print "--------------------------------\n";
print "STEP3\n";
print "--------------------------------\n\n";

my $single_alignment_output = "STEP3-single_alignment"; # Where the single alignments will be

# creating directory if not exist
mkdir $single_alignment_output unless (-e $single_alignment_output);

print "Align Pair One Sequences as single-ends...\n";
print "--------------------------------------------\n\n";
system("time bowtie -S --best -p $np $reference_genome_path $paired_alignment_output/unaligned_1.fq  --un $single_alignment_output/single_unalignment1.fq  $single_alignment_output/single_alignment1.sam");

print "Align Pair Two  Sequences as single-ends...\n";
print "--------------------------------------------\n\n";
system("time bowtie -S --best -p $np $reference_genome_path $paired_alignment_output/unaligned_2.fq  --un $single_alignment_output/single_unalignment2.fq  $single_alignment_output/single_alignment2.sam");



#STEP4 Order single alignment
print "--------------------------------\n";
print "STEP4\n";
print "--------------------------------\n";

my $ordered_single_alignment_output = "STEP4-ordered_single_alignment"; # Where the single alignments will be
mkdir $ordered_single_alignment_output unless (-e $ordered_single_alignment_output);

print "Creating BAM Alignment File One...\n";
system("samtools view -bS -o $ordered_single_alignment_output/single_alignment1.bam $single_alignment_output/single_alignment1.sam");
#print "Sorting BAM Alignment File One...\n\n";
#system("samtools sort -n $ordered_single_alignment_output/single_alignment1.bam $ordered_single_alignment_output/single_alignment1.sorted");

print "Creating BAM Alignment File Two...\n";
system("samtools view -bS -o $ordered_single_alignment_output/single_alignment2.bam $single_alignment_output/single_alignment2.sam");
#print "Sorting BAM Alignment File Two...\n\n";
#system("samtools sort -n $ordered_single_alignment_output/single_alignment2.bam $ordered_single_alignment_output/single_alignment2.sorted");

#STEP5 Merging single alignment
print "--------------------------------\n";
print "STEP5\n";
print "--------------------------------\n";

my $merged_single_alignment_output = "STEP5-merged_single_alignment"; # Where the merged single alignments will be
mkdir $merged_single_alignment_output unless (-e $merged_single_alignment_output);

print "Merging single_alignment1.bam and single_alignment2.bam ...\n";
system("samtools merge -n  $merged_single_alignment_output/merged.bam  $ordered_single_alignment_output/single_alignment1.bam $ordered_single_alignment_output/single_alignment2.bam");

print "Sorting Mate Pairs ...\n";
system("samtools sort -n $merged_single_alignment_output/merged.bam  $merged_single_alignment_output/merged.sorted");

print "Fixing Mate Pairs ...\n";
system("samtools fixmate $merged_single_alignment_output/merged.sorted.bam  $merged_single_alignment_output/merged.sorted.fixed.bam");


#STEP6 Translocated pairs
print "--------------------------------\n";
print "STEP6\n";
print "--------------------------------\n";

my $translocated_pairs_output = "STEP6-translocated_pairs"; # Where the tranlocated pairs will be
mkdir $translocated_pairs_output unless (-e $translocated_pairs_output);

print "Generating Translocated Pairs SAM list ...\n";
system("samtools view -h $merged_single_alignment_output/merged.sorted.fixed.bam  | processBAM.pl $translocated_pairs_output");

print "Convert Translocated Pairs to BAM list ...\n";
system("samtools view -bS -o $translocated_pairs_output/translocated_pairs.bam  $translocated_pairs_output/translocated_pairs.sam");

print "Sorting Translocated Pairs to BAM list ...\n";
system("samtools sort $translocated_pairs_output/translocated_pairs.bam  $translocated_pairs_output/translocated_pairs.sorted");






###############
# Subroutines #
###############
sub Usage {
	my($msg) = @_;
print STDERR "\nERR: $msg\n\n" if $msg;
print STDERR qq[$0  ].q[$Revision$].qq[\n];
print STDERR<<EOU;
Thiago Yukio Kikuchi Oliveira (stratus\@lgmb.fmrp.usp.br) 
(c)2008 Regional Blood Center of Ribeirao Preto

Usage 

    sample -refpath [reference_path] <-np [2]> -f1 [file1] -f2 [file2] -barcode ACTG

     Options:
       -help            brief help message
       -refpath         input reference genome path for the alignment. Required
       -np              the number of processors to use on alignment. Default: 2
       -f1              input solexa pair one. Required
       -f2              input solexa pair two. Required
       -barcode         the barcode string. Default: CGCGCCT


EOU
exit(1);
}

1;

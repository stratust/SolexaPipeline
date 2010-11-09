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


use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

Usage("Too few arguments") if $#ARGV<0;
GetOptions("h|?|help"=> sub { Usage(); } ) or Usage();
	   


use Bio::DB::Sam;

 my $sam = Bio::DB::Sam->new(
#                             -fasta => shift,
                             -bam  => shift,
                             -autoindex => 1,
-expand_flags  => 1,

                         );

                        my @targets = $sam->seq_ids;
foreach (@targets){
    print $_."\n";
 my @alignments = $sam->get_features_by_location(
                                                     -type => 'paired',
                                                 -seqid => $_,
                                                 -start  => 1,
                                                 -end    => 9999999800
 
                     );


print Dumper(@alignments);

=cut
for my $a (@alignments) {
    print $a->seq_id;
    my $start  = $a->start;
    my $end    = $a->end;
    my $strand = $a->strand;
    my $ref_dna= $a->dna;
    my $query_start  = $a->query->start;
    my $query_end    = $a->query->end;
    my $query_strand = $a->query->strand;
    my $query_dna    = $a->query->dna;
   
    my $cigar     = $a->cigar_str;
    my @scores    = $a->qscore;     # per-base quality scores
    my $match_qual= $a->qual;       # quality of the match

    my $paired = $a->get_tag_values('PAIRED');
 }
=cut
}

#print Dumper(@chromosomes);

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

EOU
exit(1);
}



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
use Bio::SeqIO;
use Data::Dumper;
use 5.10.0;
Usage("Too few arguments") if $#ARGV<0;
GetOptions("h|?|help"=> sub { Usage(); } ) or Usage();
	   
my $filename = shift;

my $id_file = shift;

my %hash;

open( my $id_in, "<", $id_file ) or die "Cannot open $id_file";

say "Building Hash...";

while ( <$id_in> ) {
    chomp;
	$hash{ $_ }++;
}

close( $id_in );

# grabs the FASTQ parser, specifies the Illumina variant
my $in = Bio::SeqIO->new(
	-format => 'fastq-illumina',
	-file   => $filename
);

my $out = Bio::SeqIO->new(
	-format => 'fastq-illumina',
	-file   => ">$filename"."_justbarcode"
);


# $seq is a Bio::Seq::Quality object
my $i=1;
while ( my $seq = $in->next_seq ) {
	#$in->write_seq($seq);  # convert Illumina 1.3 to Sanger format
    my $id =  $seq->display_id;
   
    $id =~ s/_read_2//g;
    $id =~ s/\/\d$//g;

	if ( $hash{ $id } == 1 ) {
		$out->write_seq( $seq );
		say $seq->display_id ."\t". $i;
	}
    $i++;

}




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



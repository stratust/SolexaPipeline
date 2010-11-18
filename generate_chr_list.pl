#!/data4/stratus/local/bin/perl
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
use 5.10.0;
use Moose;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Data::Dumper;

# My window size

my $w_size = 10000;
my $max_chr_size = 197195432;

my %part;

for (my $i = $w_size; $i <= $max_chr_size; $i = $i+ $w_size){
    $part{$i}=0;
}

 $part{$max_chr_size}=0;

my @chr_names =(qw/
    chr1
    chr2
    chr3
    chr4
    chr5
    chr6
    chr7
    chr8
    chr9
    chr10
    chr11
    chr12
    chr13
    chr14
    chr15
    chr16
    chr17
    chr18
    chr19
    chrM
    chrX
    chrY
    /);
 

    my %chr_final;
    foreach my $chr (@chr_names){
         
         my %aux = %part;
         open( my $in, '<', 'STEP8-Defining_regions/'.$chr.'.txt.sorted' ) or die "cannot open file";

         my @reads = <$in>;
         close($in);

        foreach my $read (@reads) {
            chomp $read;
            my $pos = $1 if ($read =~ /\t(\d+)$/);
            my $k_before = 0;
            foreach my $key ( sort {$a <=> $b} keys %part ) {
                if ( $pos > $k_before && $pos <= $key ) {
                    $aux{$key}++;
                    last;
                }
                $k_before = $key;

            }
        }
         
        $chr_final{$chr} = \%aux;    
    }

    my $header = "Position";
    foreach my $chr(@chr_names){
        $header .= "\t".$chr;
    }
    say $header;
    foreach my $key ( sort {$a <=> $b}keys %part ) {
        my $row = "$key";
         foreach my $chr(@chr_names){
                 $row .= "\t".$chr_final{$chr}->{$key};
        }
        say $row;
    }        


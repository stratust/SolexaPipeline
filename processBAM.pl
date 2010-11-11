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

my $dir_name = shift;

my $print_next = 0;
my $no_print   = 0;

open(my $sam,'>',$dir_name.'/translocated_pairs.sam');
open(my $out,'>',$dir_name.'/info.txt');

my $total_pairs;
my $total_nomapped_pairs;
my $total_translocated_pairs;
my $total_invalid_pairs;


while ( my $line = <STDIN> ) {

    if ($line =~ /^@/){
        print $sam $line;
        next;
    }

    $total_pairs++;

    my @f = split( "\t", $line );

    if ( $f[2] =~ /\*/ || $f[6] =~ /\*/ ){
        $total_nomapped_pairs++;
        next;
    }

    if ($print_next) {
        $print_next = 0;
        print $sam $line;
        $total_translocated_pairs++;
    }
    elsif ($no_print) {
        $no_print = 0;
        $total_invalid_pairs++;
        next;
    }
    else {

        if (
            (

                # For primer as mate1
                # igH Primer R
                (
                       $f[2] =~ /chr12/
                    && $f[3] >= 114664840
                    && ( $f[3] + 36 ) <= 114664886
                )
                ||

                # igH Primer F
                (
                       $f[2] =~ /chr12/
                    && $f[3] >= 114665161
                    && ( $f[3] + 36 ) <= 114665206
                )
            )
            ||

            # For primer as mate2
            # igH Primer R
            (
                (
                       $f[6] =~ /chr12/
                    && $f[7] >= 114664840
                    && ( $f[7] + 36 ) <= 114664886
                )
                ||

                # igH Primer F
                (
                       $f[6] =~ /chr12/
                    && $f[7] >= 114665161
                    && ( $f[7] + 36 ) <= 114665206
                )
            )
          )
        {
            $print_next = 1;
            $total_translocated_pairs++;
            print $sam $line;
        }
        else {
            $total_invalid_pairs++;
            $no_print = 1;
        }

    }
}
close($sam);

print $out "  Filter Information\n";
print $out "--------------------------------------------------------------\n";
print $out "Total of pairs:                  ".($total_pairs/2)."\n";
print $out "Total of nonmapped pairs:        ".($total_nomapped_pairs/2)."\n";
print $out "Total of invalid pairs:          ".($total_invalid_pairs/2)."\n";
print $out "Total of translocated pairs:     ".($total_translocated_pairs/2)."\n"; 
print $out "--------------------------------------------------------------\n";
close($out);

print  "  Filter Information\n";
print  "--------------------------------------------------------------\n";
print  "Total of pairs:                  ".($total_pairs/2)."\n";
print  "Total of nonmapped pairs:        ".($total_nomapped_pairs/2)."\n";
print  "Total of invalid pairs:          ".($total_invalid_pairs/2)."\n";
print  "Total of translocated pairs:     ".($total_translocated_pairs/2)."\n"; 
print  "--------------------------------------------------------------\n";
close();

1;

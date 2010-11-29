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

my $dir_name = shift;
my $ref_path_name =shift;
my $have_insert =shift;

my $print_next = 0;
my $invalid   = 0;

open(my $sam,'>',$dir_name.'/translocated_pairs.sam');

open(my $sam_invalid,'>',$dir_name.'/invalid_pairs.sam');

open(my $out,'>',$dir_name.'/info.txt');

my $total_pairs;
my $total_nomapped_pairs;
my $total_translocated_pairs;
my $total_invalid_pairs;


while ( my $line = <STDIN> ) {

    if ( $line =~ /^@/ ) {
        print $sam $line;
        next;
    }

    $total_pairs++;

    my @f = split( "\t", $line );

    if ( $f[2] =~ /\*/ || $f[6] =~ /\*/ ) {
        $total_nomapped_pairs++;
        next;
    }

    if ($print_next) {
        $print_next = 0;
        print $sam $line;
        $total_translocated_pairs++;
    }
    elsif ($invalid) {
        $invalid = 0;
        $total_invalid_pairs++;
        print $sam_invalid $line;
        next;
    }
    else {
        my $read_size;
        my $read = length $f[9];
        if ( $read < 40 ) {
            $read_size = 36;
        }
        else {
            $read_size = 54;
        }

        my ( $primerF_start, $primerR_end, $insert_size, $chr );

        if ( $ref_path_name !~ m/_igH/i ) {

            # Setting myc primers positions
            $chr           = "chr15";
            $primerF_start = 61818182;
            $primerR_end = 61818392;
            $insert_size   = 114;
            $insert_size   = 0  unless $have_insert;
            $primerR_end = $primerR_end + $insert_size;


        }
        else {
            # R primer always after the insert doesn't matter the name of the primer (R or F)
            # Setting igh primers positions
            $chr           = "chr12";
            $primerF_start = 114664845;
            $primerR_end = 114665029;
            $insert_size   = 137;
            $insert_size   = 0  unless $have_insert;
            $primerR_end = $primerR_end + $insert_size;
        }

        if (
            (

                # For primer as mate1
                # igH Primer R
                (
                       $f[2] =~ /$chr/ 
                    && $f[3] >= ( $primerR_end - $read_size - 10 )
                    && ( $f[3] + $read_size ) <=
                    ( $primerR_end + 5 )
                )
                ||

                # igH Primer F
                (
                       $f[2] =~ /$chr/
                    && $f[3] >= ( $primerF_start - 5 )
                    && ( $f[3] + $read_size ) <=
                    ( $primerF_start + $read_size + 10 )
                )
            )
            ||

            # For primer as mate2
            # igH Primer R
            (
                (
                       ($f[6] =~ /$chr/ || ($f[2] =~ /$chr/ && $f[6] eq '='))
                    && $f[7] >=
                    ( $primerR_end - $read_size - 10 )

                    && ( $f[7] + $read_size ) <=
                    ( $primerR_end + 5 )
                )
                ||

                # igH Primer F
                (
                       ($f[6] =~ /$chr/ || ($f[2] =~ /$chr/ && $f[6] eq '='))
                    && $f[7] >= ( $primerF_start - 5 )
                    && ( $f[7] + $read_size ) <=
                    ( $primerF_start + $read_size + 10 )
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
            $invalid = 1;
            print $sam_invalid $line;
        }

    }
}
close($sam);

print $out "  Filter Information\n";
print $out "--------------------------------------------------------------\n";
print $out "Total of pairs:                  " . ( $total_pairs / 2 ) . "\n";
print $out "Total of nonmapped pairs:        "
  . ( $total_nomapped_pairs / 2 ) . "\n";
print $out "Total of invalid pairs:          "
  . ( $total_invalid_pairs / 2 ) . "\n";
print $out "Total of translocated pairs:     "
  . ( $total_translocated_pairs / 2 ) . "\n";
print $out "--------------------------------------------------------------\n";
close($out);

print "  Filter Information\n";
print "--------------------------------------------------------------\n";
print "Total of pairs:                  " . ( $total_pairs / 2 ) . "\n";
print "Total of nonmapped pairs:        "
  . ( $total_nomapped_pairs / 2 ) . "\n";
print "Total of invalid pairs:          " . ( $total_invalid_pairs / 2 ) . "\n";
print "Total of translocated pairs:     "
  . ( $total_translocated_pairs / 2 ) . "\n";
print "--------------------------------------------------------------\n";
close();

1;

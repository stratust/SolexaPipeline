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
#  Copyright (C) 2005  Fundação Hemocentro de Ribeirão Preto
#
#  Laboratório de Bioinformática
#  BiT -  Bioinformatic Team
#  Fundação Hemocentro de Ribeirão Preto
#  Rua Tenente Catão Roxo, 2501
#  Ribeirão Preto - São Paulo
#  Brasil
#  CEP 14051-140
#  Fone: 55 16 39639300 Ramal 9603
#
#  Thiago Yukio Kikuchi Oliveira
#  stratus@lgmb.fmrp.usp.br
#  http://lgmb.fmrp.usp.br
#  
# $Id$
# 

=head1 NAME 

=head1 SYNOPSIS

=head1 ABSTRACT

=head1 DESCRIPTION

=head1 AUTHOR

Thiago Yukio Kikuchi Oliveira E<lt>stratus@lgmb.fmrp.usp.brE<gt>

Copyright (c) 2005 Regional Blood Center of Ribeirao Preto

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html

=head1 METHODS

=cut

package ProcessBAMStream;
use strict;
use warnings;
#use 5.10.0;
use Carp;
use vars qw($VERSION);
use Moose;
use Bio::SeqIO;
use File::Basename;


$VERSION = '0.10';



=head2 filter_by_primer

 Title   : filter_by_primer
 Usage   : filter_by_primer(<STDIN>)
 Function: 
 Returns : 
 Args    : 

=cut 

sub filter_by_primer {
    my ( $self ) = @_;

    my $print_next = 0;
    while (my $line = <STDIN>) {
        my @f = split( "\t", $line );

        next if $f[6] =~ /\*/;

        if ($print_next) {
            $print_next = 0;
            return $line;
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
                return $line;
            }
        }
    }
}



1;


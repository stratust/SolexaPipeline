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
use Getopt::Long;

my ($input_file);
Usage("Too few arguments") if $#ARGV<0;
GetOptions( 
    "h|?|help" => sub { Usage(); },
    "i=s" => \$input_file,
    
)or Usage();

#Open file

open(my $in, "<",$input_file);


while (<$in>){
    chomp;
    my ($chr,$position, $density, $f, $r) = split "\t", $_;
    my ($start,$end) = split "-", $position;
    say  "$chr\t$start\t$end\t$density\t$f\t$r";
}

close($in);


###############
# Subroutines #
###############
sub Usage {
	my($msg) = @_;
print STDERR "\nERR: $msg\n\n" if $msg;
print STDERR qq[$0  ].q[$Revision$].qq[\n];
print STDERR<<EOU;
Thiago Yukio Kikuchi Oliveira (stratus\@lgmb.fmrp.usp.br) 
(c)2010 Regional Blood Center of Ribeirao Preto

Usage 

    hotspots2bed -i [inputfilename] 

     Options:
       -help            brief help message
       -i         input hotspots file. Required


EOU
exit(1);
}



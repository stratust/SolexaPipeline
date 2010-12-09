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
use Getopt::Long;

Usage("Too few arguments") if $#ARGV<0;

my ($myc_f, $igh_f, $merged_f);

GetOptions(
    "h|?|help"=> sub { Usage(); },
    "myc=s" => \$myc_f,
    "igh=s" => \$igh_f,
    "merged=s" => \$merged_f,

) or Usage();
	   

open(my $in,"<",$merged_f);

my %merged;
foreach (<$in>){
    chomp;
    my @f = split "\t",$_;
    my $id = $f[$#f];
    $id =~ s/a//g;
    $merged{$id}->{name} = $f[($#f - 1)];
    $merged{$id}->{igh} = 0;
    $merged{$id}->{myc} = 0;
}

close($in);

open($in,"<",$igh_f);

foreach (<$in>){
    chomp;
    my @f = split "\t",$_;
    my $id = $f[$#f];
    $id =~ s/a//g;
    $merged{$id}->{igh} += $f[3];
}

close($in);

open($in,"<",$myc_f);

foreach (<$in>){
    chomp;
    my @f = split "\t",$_;
    my $id = $f[$#f];
    $id =~ s/a//g;
    $merged{$id}->{myc} += $f[3];
}

close($in);


print "ID\tLocation\tIgh\tmyc\n";
foreach my $key (sort {$a <=> $b} keys %merged){
    print "a$key\t$merged{$key}{name}\t$merged{$key}->{igh}\t$merged{$key}->{myc}\n";
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



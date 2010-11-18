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
use List::Compare;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Data::Dumper;

# My window size

my $w_size = 10000;


#STEP7 Clustering Translocated pairs
print "--------------------------------\n";
print "STEP9\n";
print "--------------------------------\n";
my $slide_window =
  "STEP9-slide_window";    # Where the tranlocated pairs will be
if ( -e $slide_window ) {

    print "Skipping STEP9...\n\n";

}
else {
    print "Sliding Window...\n\n";
    mkdir $slide_window;
    my @files = 
    open( my $in, '<', 'STEP8-Defining_regions/chr12.txt.sorted' ) or die "cannot open file";

    # Number of reads
#    my @reads = <$in>;
    my @reads;
    my @AoA_windows;
    say "Creting Windows";
    for ( my $i=0; $i <= $#reads; $i++  ){
        for ( my $y = $i; $y < $#reads; $y++){
            chomp $reads[$y];
            my ($chr,$id,$pos) = split ('\t',$reads[$y]);
            unless ($AoA_windows[$i][0]){
                  push(@{$AoA_windows[$i]}, $pos );
             }
             else{
                  push(@{$AoA_windows[$i]}, $pos ) if ($pos <= $AoA_windows[$i][0] + $w_size);
             }

        }
    }
    
    my $first_pos = 0;
    my %frame;
    while (my $line = <$in>){
           chomp $line;
            my ($chr,$id,$pos) = split ('\t',$line);
            if ($first_pos){
                if ($pos <= $first_pos + $w_size){
                    push(@{$frame{$first_pos}},$pos);
                }
                else{
                    $first_pos = $pos;
                    push(@{$frame{$first_pos}},$first_pos);
                    
                }
            }
            else{
                $first_pos = $pos;
                push(@{$frame{$first_pos}},$first_pos);
                next;
            }
            
   
    }
    close($in);
   
#    print Dumper(%frame);

#   exit;
    my @aux;
    say "Removing window contained in other windows";
    foreach my $key (sort {$a <=> $b} keys %frame ){
        unless ($aux[0]){
            push(@aux,$frame{$key});
            next;
        }
        my $lc = List::Compare->new('--unsorted', $aux[$#aux], $frame{$key});
        if ( $lc->is_RsubsetL ){
            next;
        }
        else{
            push(@aux,$frame{$key});
        }
        
    }
    
   
    say "Merging windows with less than window size";
    my @final;
#    print Dumper(@aux);
#    exit;
    if (&is_groups(@aux)){

    foreach my $window (@aux){
        unless (defined $final[0]){
            push(@final,$window);
            next;
        }
        if ($window->[0] < ($final[$#final][$#{$final[$#final]}] + $w_size)){
            my $lc = List::Compare->new('--unsorted',$final[$#final],$window);
            my @a = $lc->get_union;
            $final[$#final] = \@a;
        }
        else{
            push(@final,$window);
        }
    }

    }
    else{
        @final = @aux;
    
    }
    
    
    # Discovering the amount of reads in each window
    my %hash;
    say "Counting amout of reads in each window";
    foreach my $window (@final){
        my $key = @{$window};
        push(@{$hash{$key}}, $window);
    }

    #print Dumper(%hash);

    # Print in descendent order
    foreach my $k ( sort {$b <=> $a} keys %hash){
        last if $k < 3;
        foreach my $window (@{$hash{$k}}){

            say min(@{$window}).' - ' . max(@{$window}) ."\t". $k; 

        }
    }

}

sub is_groups(){
    my (@aux) = @_;
    
    my $is_group = 0;
    foreach (@aux){
       my $scalar = @{$_}; 
       $is_group++ if $scalar > 1;
    }

    return $is_group;
}

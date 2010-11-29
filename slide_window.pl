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
use Data::Merger qw(merger);
# My window size

my $w_size = 2000;


#STEP7 Clustering Translocated pairs
print "--------------------------------\n";
print "STEP9\n";
print "--------------------------------\n";
my $slide_window = "STEP9-slide_window";   # Where the tranlocated pairs will be
if ( -e $slide_window ) {

    print "Skipping STEP9...\n\n";

}
else {
    

    say "Sliding Window...\n";
    mkdir $slide_window;
    my @chrs = (
        qw/
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
          /
    );

    open( my $out, ">", "$slide_window/hotspots.txt" );
    my $y;
    foreach my $chrm (@chrs) {
        next unless ('STEP8-Defining_regions/' . $chrm . '.txt.sorted');
        open( my $in, '<', 'STEP8-Defining_regions/' . $chrm . '.txt.sorted' )
          or die "cannot open file";

        # Number of reads

        my $first_pos = 0;
        my %frame;
        while ( my $line = <$in> ) {
            chomp $line;
            my ( $chr, $id, $pos ) = split( '\t', $line );
            my %h = (
                'chr' => $chr,
                'id'  => $id,
                'pos' => $pos,
            );
            if ($first_pos) {
                if ( $pos <= $first_pos + $w_size ) {
                    push( @{ $frame{$first_pos} }, \%h );
                }
                else {
                    push( @{ $frame{$pos} }, \%h );
                    $first_pos = $pos;
                }
            }
            else {
                $first_pos = $pos;
                push( @{ $frame{$first_pos} }, \%h );
                next;
            }

        }
        close($in);

        foreach ( sort { $a <=> $b } keys %frame ) {

            for ( my $i = 0 ; $i < $#{ $frame{$_} } ; $i++ ) {
                say "Lower value before found "
                  . $frame{$_}->[ ( $i + 1 ) ]->{pos} . "<"
                  . $frame{$_}->[$i]->{pos}
                  if ( $frame{$_}->[ ( $i + 1 ) ]->{pos} <
                    $frame{$_}->[$i]->{pos} );
            }

            #say $_ ." (".$frame{$_}->[0]->{pos}.") ". " => " . $#{$frame{$_}} ;
        }

        #        exit;

        say "Merging windows with less than window size";
        my @final;

        foreach my $key ( sort { $a <=> $b } keys %frame ) {
            
            unless ( defined $final[0] ) {
                push( @final, $frame{$key} );
                next;
            }

                #say $frame{$key}->[0]->{pos} . '<' . ($final[$#final][$#{$final[$#final]}]->{pos} + $w_size ) ;
                if ( $frame{$key}->[0]->{pos} <
                    ( $final[$#final][$#{$final[$#final]}]->{pos} + $w_size ) )
                {

                    my @new = (@{$final[$#final]},@{$frame{$key}});

                    my @sorted =  sort { $a->{pos} <=> $b->{pos} } @new;
                    $final[$#final] = \@sorted;
                    
                    #say "----final-------";
                    #print Dumper($final[$#final]);
                    #say "-------frame----------";
                    #print Dumper($frame{$key});
                    #say "--------Merger---------";
                    #print Dumper($merged_data);
                for (my $i=0; $i < $#{$final[$#final]}; $i++){
                say "Lower value before found "
                  . $final[$#final]->[ ( $i + 1 ) ]->{pos} . "<"
                  . $final[$#final]->[$i]->{pos}
                  if ( $final[$#final]->[ ( $i + 1 ) ]->{pos} < $final[$#final]->[$i]->{pos} );
            }
 
                    
                }
                else {

            push( @final, $frame{$key} );

                            }
        }

    
        #print Dumper(%hash);

        # Print in any order
            foreach my $window ( @final ) {

                my @positions_cluster;
                my $nR = 0;
                my $nF = 0;
                foreach ( @{$window} ) {
                    $y++;
                    push( @positions_cluster, $_->{pos} );
                    if ( $_->{id} =~ m/R$/ ) {
                        $nR++;
                    }
                    else {
                        $nF++;
                    }
                }

                if (   ( $nR / ( $#{$window} + 1 ) ) >= 0.25
                    && ( $nR / ( $#{$window} + 1 ) ) <= 0.75
                    && $#{$window} >= 0 )
                {

                    say $out $chrm . "\t"
                      . min(@positions_cluster) . '-'
                      . max(@positions_cluster) . "\t"
                      . ( $#{$window} + 1 ) . "\t"
                      . $nF . "\t"
                      . $nR;

                }

            }
    }
    say $y;
    close($out);
}




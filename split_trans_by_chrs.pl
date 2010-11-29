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
use Data::Dumper;

# igh or myc
my $primer_used = shift;


#STEP7 Clustering Translocated pairs
print "--------------------------------\n";
print "STEP8\n";
print "--------------------------------\n";

my $translocated_pairs_cluster_output =
  "STEP8-Defining_regions";    # Where the tranlocated pairs will be
if ( -e $translocated_pairs_cluster_output ) {

    print "Skipping STEP8...\n\n";

}
else {

    print "Defining regions in each chromossome...\n\n";
    mkdir $translocated_pairs_cluster_output;
    open( my $in, '<', 'STEP7-clustered_translocated_pairs/cluster_translocated_pairs.sam' );

    my %read;
    my $number_of_pairs;
    while ( my $pair_1 = <$in> ) {
        if ( $pair_1 =~ /^\@/ ) {
            next;
        }
        $number_of_pairs++;

        my $pair_2 = <$in>;

        my ( $id1, $flag1, $chr_1, $pos_read1, $maq1, $cirgar1, $chr_mate1,
            $pos_mate1 )
          = split( "\t", $pair_1 );
        my ( $id2, $flag2, $chr_2, $pos_read2, $maq2, $cirgar2, $chr_mate2,
            $pos_mate2 )
          = split( "\t", $pair_2 );
        
        # Find primer (F or R) and mate who have the primer (1 or 2);
        my $mate_primer = find_mate_primer($pair_1, $primer_used);
        
        my $key;
        my $alignment;
        if ($mate_primer =~ m/1/){
             my %hash;
             $hash{id} = $id1;
             $hash{pos} = $pos_read2;
             push(@{$read{$chr_2}}, \%hash);
        }
        else{
             my %hash;
             $hash{id} = $id1;
             $hash{pos} = $pos_read1;
            push(@{$read{$chr_1}},\%hash);
        }

    }
    close($in);
    
    my %chr_number;
    
    foreach my $key ( keys %read ) {

        open( my $out, '>', $translocated_pairs_cluster_output . "/$key.txt" );

        foreach my $entry ( @{ $read{$key} } ) {
            $chr_number{$key}++;
            print $out "$key\t$entry->{id}\t$entry->{pos}\n";
        }
        close($out);
    }
    
    open( my $info, '>', $translocated_pairs_cluster_output . "/info.stats" );

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
    foreach my $key ( @chr_names) {
        $chr_number{$key} = 0 unless $chr_number{$key};
        print $info "$key\t$chr_number{$key}\n";
    }
    close($info);

    print "Sorting chromossomes per position...\n";
    system('for i in '.$translocated_pairs_cluster_output.'/*.txt; do sort -k3n,3 "$i" > "$i.sorted"; done');


=cut
    foreach my $key ( keys %unique ) {
        my $count = $unique{$key}{count};
        next if $count <= $cutoff;
        $i++;
        my $alignment = $unique{$key}{alignment};
        my $pair_id = $1 if $alignment =~ /^(\S+)\t/;

        my $primer = $unique{$key}{primer};
        my $new_id = "pair" . $i . "_c" . $count . "_" . $primer;
        $alignment =~ s/$pair_id/$new_id/g;

        print $out $alignment;

    }
    close($out);
=cut

}

sub find_primer_region {
    my ( $line, $ref_path_name ) = @_;

    my @f = split( "\t", $line );

    my $read_size;
    my $read = length $f[9];
    if ( $read < 40 ) {
        $read_size = 36;
    }
    else {
        $read_size = 54;
    }

    my ( $primerF_start, $primerR_end, $insert_size, $chr );
    my $igh;
    if ( $ref_path_name !~ m/igH/i ) {

        # Setting myc primers positions
        $chr           = "chr15";
        $primerF_start = 61818182;
        $primerR_end   = 61818392;
        $insert_size   = 114;
        $insert_size   = 0;
        $primerR_end   = $primerR_end + $insert_size;

    }
    else {

# R primer always after the insert doesn't matter the name of the primer (R or F)
# Setting igh primers positions
        $chr           = "chr12";
        $primerF_start = 114664845;
        $primerR_end   = 114665029;
        $insert_size   = 137;
        $insert_size   = 0;
        $primerR_end   = $primerR_end + $insert_size;
        $igh           = 1;
    }

    if (
        (

            # For primer as mate1
            # Primer R
            (
                   $f[2] =~ /$chr/
                && $f[3] >= ( $primerR_end - $read_size - 10 )
                && ( $f[3] + $read_size ) <= ( $primerR_end + 5 )
            )
            ||

            # For primer as mate2
            # Primer R

            (
                ( $f[6] =~ /$chr/ || ( $f[2] =~ /$chr/ && $f[6] eq '=' ) )
                && $f[7] >=
                ( $primerR_end - $read_size - 10 )

                && ( $f[7] + $read_size ) <= ( $primerR_end + 5 )
            )

        )
      )
    {
        unless ($igh) {
            return 'F';
        }
        else {
            return 'R';
        }
    }
    else {
        unless ($igh) {
            return 'R';
        }
        else {
            return 'F';
        }

    }

}

sub find_mate_primer {
    my ( $line, $ref_path_name ) = @_;

    my @f = split( "\t", $line );

    my $read_size;
    my $read = length $f[9];
    if ( $read < 40 ) {
        $read_size = 36;
    }
    else {
        $read_size = 54;
    }

    my ( $primerF_start, $primerR_end, $insert_size, $chr );
    my $igh;
    if ( $ref_path_name !~ m/igH/i ) {

        # Setting myc primers positions
        $chr           = "chr15";
        $primerF_start = 61818182;
        $primerR_end   = 61818392;
        $insert_size   = 114;
        $insert_size   = 0;
        $primerR_end   = $primerR_end + $insert_size;

    }
    else {

# R primer always after the insert doesn't matter the name of the primer (R or F)
# Setting igh primers positions
        $chr           = "chr12";
        $primerF_start = 114664845;
        $primerR_end   = 114665029;
        $insert_size   = 137;
        $insert_size   = 0;
        $primerR_end   = $primerR_end + $insert_size;
        $igh           = 1;
    }

    if (

        # For primer as mate1
        # Primer R
        $f[2] =~ /$chr/
        && $f[3] >= ( $primerR_end - $read_size - 10 )
        && ( $f[3] + $read_size ) <= ( $primerR_end + 5 )
      )
    {
        return 'R1';
    }

    if (

        # For primer as mate2
        # Primer R

        ( $f[6] =~ /$chr/ || ( $f[2] =~ /$chr/ && $f[6] eq '=' ) )
        && $f[7] >=
        ( $primerR_end - $read_size - 10 )

        && ( $f[7] + $read_size ) <= ( $primerR_end + 5 )

      )
    {
        return 'R2';
    }

    if

      # Primer F Mate1
      (    $f[2] =~ /$chr/
        && $f[3] >= ( $primerF_start - 5 )
        && ( $f[3] + $read_size ) <= ( $primerF_start + $read_size + 10 ) )
    {
        return 'F1';
    }

    if (   ( $f[6] =~ /$chr/ || ( $f[2] =~ /$chr/ && $f[6] eq '=' ) )
        && $f[7] >= ( $primerF_start - 5 )
        && ( $f[7] + $read_size ) <= ( $primerF_start + $read_size + 10 ) )
    {
        return 'F2';
    }

}


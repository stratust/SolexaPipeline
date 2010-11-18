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

package ProcessFastq;
use strict;
use warnings;
#use 5.10.0;
use Carp;
use vars qw($VERSION);
use Moose;
use Bio::SeqIO;
use File::Basename;


$VERSION = '0.10';

has "barcode" => (is => 'ro', isa => 'Str');
has "file1" => (is => 'ro', isa => 'Str');
has "file2" => (is => 'ro', isa => 'Str');
has "output_path" => (is => 'ro', isa => 'Str');


=head2 index_barcode

 Title   : index_barcode
 Usage   : index_barcode()
 Function: 
 Returns : 
 Args    : 

=cut 

sub index_barcode {
    my($self) = @_;
    
    my @index;
    my $bc = $self->barcode;
    my $rv_bc = $bc;
    #translate
    $rv_bc =~ tr/[ATCGatcg]/[TAGCtagc]/;
    #reverse
    $rv_bc = reverse $rv_bc;


    
    my @barcode = split('',$bc);
        
    
    # Using Wildcards
    for ( my $i = 0; $i <= $#barcode ; $i++ ) {
        my @aux = split( '', $bc );
        my @aux_reverse = split( '', $rv_bc );
        $aux[$i] = "[ATCGnatcgn\.]";
        $aux_reverse[$i] = "[ATCGNatcgn\.]";
        my $regex = join( '', @aux );
        my $regex_reverse = join( '', @aux_reverse );
       
        
        push( @index, qr/$regex/ );
        push( @index, qr/$regex_reverse/ );
        # A maxium of 1 indels before plus a mismatch
        push( @index, qr/[ATCGNatcgn\.]{1,1}$regex/ );
        push( @index, qr/[ATCGNatcgn\.]{1,1}$regex_reverse/ );
       

    }
   
    # Emulate a maximum of 1 indels after
    substr($bc,0,1) = '';
    push( @index, qr/$bc/ );
    #verse
    substr($rv_bc,0,1) = '';
    push( @index, qr/$rv_bc/ );

    
    return @index;

}


=head2 filter_fastq

 Title   : filter_fastq
 Usage   : filter_fastq()
 Function: 
 Returns : 
 Args    : 

=cut 

sub filter_fastq {
    my ($self) = @_;
    
    # path where generated files will be send
    my $path_out = $self->output_path;

    # grabs the FASTQ parser, specifies the Illumina variant
    my $in1 = Bio::SeqIO->new(
        -format => 'fastq-illumina',
        -file   => $self->file1,
    );

    my $in2 = Bio::SeqIO->new(
        -format => 'fastq-illumina',
        -file   => $self->file2,
    );
 
    mkdir "$path_out" unless ( -e "$path_out");
    # Create object Bio::SeqIO to write new files
    my $out1 = Bio::SeqIO->new(
        -format => 'fastq-illumina',
        -file   => ">$path_out/".basename($self->file1).".barcode_matched",
        -quality_header => 1,
    );
    
    my $out2 = Bio::SeqIO->new(
        -format => 'fastq-illumina',
        -file   => ">$path_out/".basename($self->file2).".barcode_matched",
        -quality_header => 1,

    );

    # Create object Bio::SeqIO to write new umatched files
=cut    
    my $out1_unmatched = Bio::SeqIO->new(
        -format => 'fastq-illumina',
        -file   => ">$path_out/".basename($self->file1).".barcode_unmatched",
        -quality_header => 1,
    );

    my $out2_unmatched = Bio::SeqIO->new(
        -format => 'fastq-illumina',
        -file   => ">$path_out/".basename($self->file2).".barcode_unmatched",
        -quality_header => 1,

    );
=cut
    # Get barcode regex
    my @index = $self->index_barcode();
   
    my $total_seq = 0;
    my $total_match = 0;
    my $total_match_file1 = 0;
    my $total_match_file2 = 0;
    my $total_nomatch = 0;
    my $total_duplicated_barcode = 0;


    while ( my $seq1 = $in1->next_seq ) {
        # Get the second file seq;
        my $seq2 = $in2->next_seq;
       
        # Making IDs became identical
        $seq1->display_id($self->strip_id($seq1->display_id));
        $seq2->display_id($self->strip_id($seq2->display_id));

        # Initialize match counter;
        my $c_match = 0;

        my $seq1_truncated; 
        foreach my $regex (@index){
            if ($seq1->seq =~ m/^$regex/i){
                
                #triming the barcode
                $seq1_truncated = $seq1->trunc($+[0]+1,$seq1->length);

                $c_match++;
                last;
            }
        }
        
        my $seq2_truncated;
        foreach my $regex (@index){
            if ($seq2->seq =~ m/^$regex/i){

                #triming the barcode
                $seq2_truncated = $seq2->trunc($+[0]+1,$seq2->length);               
 
                $c_match++;
                last;
            }
        }
        
        if ($c_match == 1){
            if (defined $seq1_truncated){
                $out1->write_seq($seq1_truncated) ;
                $total_match_file1++;
            }
            else{
                $out1->write_seq($seq1) ;
            }

            if (defined $seq2_truncated){
                $out2->write_seq($seq2_truncated);
                $total_match_file2++;
            }
            else{
                $out2->write_seq($seq2) ;
            }


            #count match
            $total_match++;
        }
        elsif ($c_match < 1){
            $total_nomatch++;
        }
        else{
            $total_duplicated_barcode++; 
        }

        $total_seq++;
    }
   

    # output informations
    open(my $out_info,">","$path_out/barcode_filter_info.txt");

    print $out_info "Filter Information\n";
    print $out_info "-------------------------------------------------------------\n\n";
    print $out_info "Total of sequences:            $total_seq\n";
    print $out_info "Total of matched sequences:    $total_match\n";
    print $out_info "Total of matched seq file1:    $total_match_file1\n";
    print $out_info "Total of matched seq file2:    $total_match_file2\n";
    print $out_info "Total of nomatched sequence:   $total_nomatch\n";
    print $out_info "Total of duplicated barcode:   $total_duplicated_barcode\n";
    print $out_info "-------------------------------------------------------------\n";

    close($out_info);

    # Print in the screen
    print  "Filter Information\n";
    print  "-------------------------------------------------------------\n\n";
    print  "Total of sequences:            $total_seq\n";
    print  "Total of matched sequences:    $total_match\n";
    print  "Total of matched seq file1:    $total_match_file1\n";
    print  "Total of matched seq file2:    $total_match_file2\n";
    print  "Total of nomatched sequence:   $total_nomatch\n";
    print  "Total of duplicated barcode:   $total_duplicated_barcode\n";
    print  "-------------------------------------------------------------\n";

}

=head2 filter_fastq_nobioperl

 Title   : filter_fastq_nobioperl
 Usage   : filter_fastq_nobioperl()
 Function: 
 Returns : 
 Args    : 

=cut 

sub filter_fastq_nobioperl {
    my($self) = @_;
  
    # path where generated files will be send
    my $path_out = $self->output_path;

    # grabs the FASTQ parser, specifies the Illumina variant
    open(my $in1,"<",$self->file1);
    
    open(my $in2,"<",$self->file2);
 
    mkdir "$path_out" unless ( -e "$path_out");

    # Create object Bio::SeqIO to write new files
    open(my $out1,">","$path_out/".basename($self->file1).".barcode_matched");
    open(my $out2,">","$path_out/".basename($self->file2).".barcode_matched");
    

    # Create object Bio::SeqIO to write new umatched files
=cut    
    my $out1_unmatched = Bio::SeqIO->new(
        -format => 'fastq-illumina',
        -file   => ">$path_out/".basename($self->file1).".barcode_unmatched",
        -quality_header => 1,
    );

    my $out2_unmatched = Bio::SeqIO->new(
        -format => 'fastq-illumina',
        -file   => ">$path_out/".basename($self->file2).".barcode_unmatched",
        -quality_header => 1,

    );
=cut
    # Get barcode regex
    my @index = $self->index_barcode();
   
    my $total_seq = 0;
    my $total_match = 0;
    my $total_match_file1 = 0;
    my $total_match_file2 = 0;
    my $total_nomatch = 0;
    my $total_duplicated_barcode = 0;

    $/ = "@";
    while ( my $seq1 = <$in1> ) {
        # Get the second file seq;
        $/ = "@";
        my $seq2 = <$in2>;
        
        next if $seq1 =~ /^\@/;

        $seq1 =~ s/\@$//;
        $seq2 =~ s/\@$//;
        
       
        # Making IDs became identical
        my $id1 = $1  if ($seq1 =~ /^(\S+)\n/);
        my $id1_strip = $self->strip_id($id1);
        $seq1 =~ s/$id1/$id1_strip/g;

        my $id2 = $1  if ($seq2 =~ /^(\S+)\n/);
        my $id2_strip = $self->strip_id($id2);
        $seq2 =~ s/$id2/$id2_strip/g;


        # Initialize match counter;
        my $c_match = 0;
        my $seq1_truncated =0;
        
        my($sid1,$sequence1,$quality1) = ($1,$2,$3) if ($seq1 =~ /^(\S+)\n(\S+)\n\S+\n(\S+)\n/ );
        foreach my $regex (@index){

            if ($sequence1 =~ m/^$regex/i){
                my $seq_truncated = substr $sequence1, $+[0], length $sequence1;
                my $qual_truncated = substr $quality1, $+[0], length $quality1;
                $seq1 = "$sid1\n$seq_truncated\n\+$sid1\n$qual_truncated\n";
                #$seq1 =~ s/\Q$sequence\E/\Q$seq_truncated\E/;
                #$seq1 =~ s/\Q$quality\E/\Q$qual_truncated\E/;

                $seq1_truncated = 1;
                $c_match++;
                last;
            }
        }
        
        my $seq2_truncated =0;

        my ($sid2,$sequence2,$quality2) = ($1,$2,$3) if ($seq2 =~ /^(\S+)\n(\S+)\n\S+\n(\S+)\n/ );
#        exit if $total_seq == 4;
        foreach my $regex (@index){
            if ($sequence2 =~ m/^$regex/i){
                my $seq_truncated = substr $sequence2, $+[0], length $sequence2;
                my $qual_truncated = substr $quality2, $+[0], length $quality2;
 
                $seq2 = "$sid2\n$seq_truncated\n\+$sid2\n$qual_truncated\n";

                #$seq2 =~ s/\Q$sequence\E/\Q$seq_truncated\E/;
                #$seq2 =~ s/\Q$quality\E/\Q$qual_truncated\E/;

                $seq2_truncated = 1;

                $c_match++;
                last;
            }
        }
        
        if ($c_match == 1){
            if ($seq1_truncated){
                $total_match_file1++;
            }
            print $out1 '@'.$seq1;

            if ($seq2_truncated){
                $total_match_file2++;
            }

            print $out2 '@'.$seq2;

            #count match
            $total_match++;
        }
        elsif ($c_match < 1){
            $total_nomatch++;
        }
        else{
            $total_duplicated_barcode++; 
        }

        $total_seq++;
    }
   

    # output informations
    open(my $out_info,">","$path_out/barcode_filter_info.txt");

    print $out_info "Filter Information\n";
    print $out_info "-------------------------------------------------------------\n\n";
    print $out_info "Total of sequences:            $total_seq\n";
    print $out_info "Total of matched sequences:    $total_match\n";
    print $out_info "Total of matched seq file1:    $total_match_file1\n";
    print $out_info "Total of matched seq file2:    $total_match_file2\n";
    print $out_info "Total of nomatched sequence:   $total_nomatch\n";
    print $out_info "Total of duplicated barcode:   $total_duplicated_barcode\n";
    print $out_info "-------------------------------------------------------------\n";

    close($out_info);

    # Print in the screen
    print  "Filter Information\n";
    print  "-------------------------------------------------------------\n\n";
    print  "Total of sequences:            $total_seq\n";
    print  "Total of matched sequences:    $total_match\n";
    print  "Total of matched seq file1:    $total_match_file1\n";
    print  "Total of matched seq file2:    $total_match_file2\n";
    print  "Total of nomatched sequence:   $total_nomatch\n";
    print  "Total of duplicated barcode:   $total_duplicated_barcode\n";
    print  "-------------------------------------------------------------\n";


}


=head2 strip_id

 Title   : strip_id
 Usage   : strip_id()
 Function: 
 Returns : 
 Args    : 

=cut 

sub strip_id {
    my ( $self, $display_id ) = @_;
    $display_id =~ s/_read_\d//g;
    $display_id =~ s/\/\d$/\//g;
    return $display_id;
}


1;


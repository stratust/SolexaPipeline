use Bio::SeqIO;
use strict;
use warnings;


    # create one SeqIO object to read in, and another to write out
    # *STDIN is a 'globbed' filehandle with the contents of Standard In
    my $seqin = Bio::SeqIO->new(-file     => shift,
                                -format => 'fasta');
    my $seqout = Bio::SeqIO->new(-file   => ">".shift,
                                 -format => 'fasta');

    # write each entry in the input file to the output file
    while (my $inseq = $seqin->next_seq) {
         $seqout->write_seq($inseq);
    }
    exit;


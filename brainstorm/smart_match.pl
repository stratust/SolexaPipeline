use strict;
use warnings;
use 5.10.0;
use List::Compare;

my @array1 = (1,2,3,4,5);
my @array2 = (1,2,3,4,);

my  $lc = List::Compare->new('--unsorted', \@array1, \@array2);

 say 'Arra1 is a subset of array2' if $lc->is_LsubsetR;
 say 'Array2 is a subset of array1' if $lc->is_RsubsetL;


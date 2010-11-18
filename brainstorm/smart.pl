my @ar = (0..10);
my @ar2 = @ar[1..4];
for my $i (@ar) {
    # matches for $i==4 only
    print "$i ~~ slice? ", isTrue($i ~~ @ar[1..4]);

    # matches properly
    print "\t$i ~~ ar2? ", isTrue($i ~~ @ar2), "\n";
}

sub isTrue {
    return ($_[0])?"true":"false";
}


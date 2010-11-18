use Algorithm::NeedlemanWunsch;
    sub score_sub {
        if (!@_) {
            return -2; # gap penalty
        }

        return ($_[0] eq $_[1]) ? 1 : -1;
    }

    my $matcher = Algorithm::NeedlemanWunsch->new(\&score_sub);
    my $score = $matcher->align(
               \@a,
               \@b,
               {   align     => \&on_align,
                   shift_a => \&on_shift_a,
                   shift_b => \&on_shift_b,
                   select_align => \&on_select_align
               });

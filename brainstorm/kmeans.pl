use Algorithm::Kmeanspp;
  
  # input documents
  my %documents = (
      Alex => { 'Pop'     => 10, 'R&B'    => 6, 'Rock'   => 4 },
      Bob  => { 'Jazz'    => 8,  'Reggae' => 9                },
      Thiago  => { 'Jazz'    => 8,  'Reggae' => 9                },
      Lilian  => { 'Jazz'    => 8,  'Reggae' => 4                },
      Dave => { 'Classic' => 4,  'World'  => 4                },
      Ted  => { 'Jazz'    => 9,  'Metal'  => 2, 'Reggae' => 6 },
      Fred => { 'Hip-hop' => 3,  'Rock'   => 3, 'Pop'    => 3 },
      Sam  => { 'Classic' => 8,  'Rock'   => 1                },
  );
  
  my $kmp = Algorithm::Kmeanspp->new;
  
  foreach my $id (keys %documents) {
      $kmp->add_document($id, $documents{$id});
  }
  
  my $num_cluster = 5;
  my $num_iter    = 20;
  $kmp->do_clustering($num_cluster, $num_iter);             
  
  # show clustering result
  foreach my $cluster (@{ $kmp->clusters }) {
      print join "\t", @{ $cluster };
      print "\n";
  }
  # show cluster centroids
  foreach my $centroid (@{ $kmp->centroids }) {
      print join "\t", map { sprintf "%s:%.4f", $_, $centroid->{$_} }
          keys %{ $centroid };
      print "\n";
  }

#!/usr/bin/perl -w

@algs=("GREEDY_S1_P1","GREEDY_S1_P2","GREEDY_S1_P3","GREEDY_S1_P4","GREEDY_S1_P5","GREEDY_S1_P6","GREEDY_S1_P7","GREEDY_S2_P1","GREEDY_S2_P2","GREEDY_S2_P3","GREEDY_S2_P4","GREEDY_S2_P5","GREEDY_S2_P6","GREEDY_S2_P7","GREEDY_S3_P1","GREEDY_S3_P2","GREEDY_S3_P3","GREEDY_S3_P4","GREEDY_S3_P5","GREEDY_S3_P6","GREEDY_S3_P7","GREEDY_S4_P1","GREEDY_S4_P2","GREEDY_S4_P3","GREEDY_S4_P4","GREEDY_S4_P5","GREEDY_S4_P6","GREEDY_S4_P7","GREEDY_S5_P1","GREEDY_S5_P2","GREEDY_S5_P3","GREEDY_S5_P4","GREEDY_S5_P5","GREEDY_S5_P6","GREEDY_S5_P7","GREEDY_S6_P1","GREEDY_S6_P2","GREEDY_S6_P3","GREEDY_S6_P4","GREEDY_S6_P5","GREEDY_S6_P6","GREEDY_S6_P7","GREEDY_S7_P1","GREEDY_S7_P2","GREEDY_S7_P3","GREEDY_S7_P4","GREEDY_S7_P5","GREEDY_S7_P6","GREEDY_S7_P7");

print STDOUT "## OVERALL RANKING ##\n";
@NDIMS = ("1","2","3");
@SLAPROBS = ("0.0","0.25","0.5");
@NUMSERVICES= ("100","200","500");
process();

print STDOUT "## FOR NDIM = 1 ##\n";
@NDIMS = ("1");
@SLAPROBS = ("0.0","0.25","0.5");
@NUMSERVICES= ("100","200","500");
process();

print STDOUT "## FOR NDIM = 2 ##\n";
@NDIMS = ("2");
@SLAPROBS = ("0.0","0.25","0.5");
@NUMSERVICES= ("100","200","500");
process();

print STDOUT "## FOR NDIM = 3 ##\n";
@NDIMS = ("3");
@SLAPROBS = ("0.0","0.25","0.5");
@NUMSERVICES= ("100","200","500");
process();

print STDOUT "## FOR SLAPROB = 0.0 ##\n";
@NDIMS = ("1","2","3");
@SLAPROBS = ("0.0");
@NUMSERVICES= ("100","200","500");
process();

print STDOUT "## FOR SLAPROB = 0.25 ##\n";
@NDIMS = ("1","2","3");
@SLAPROBS = ("0.25");
@NUMSERVICES= ("100","200","500");
process();

print STDOUT "## FOR SLAPROB = 0.5 ##\n";
@NDIMS = ("1","2","3");
@SLAPROBS = ("0.5");
@NUMSERVICES= ("100","200","500");
process();

print STDOUT "## FOR N = 100 ##\n";
@NDIMS = ("1","2","3");
@SLAPROBS = ("0.0","0.25","0.5");
@NUMSERVICES= ("100");
process();

print STDOUT "## FOR N = 200 ##\n";
@NDIMS = ("1","2","3");
@SLAPROBS = ("0.0","0.25","0.5");
@NUMSERVICES= ("200");
process();

print STDOUT "## FOR N = 500 ##\n";
@NDIMS = ("1","2","3");
@SLAPROBS = ("0.0","0.25","0.5");
@NUMSERVICES= ("500");
process();


### MAIN SUBROUTINE
sub process
{
  # number of failures when somebody else found a solution
  %num_fail = ();
  %dfbest = ();
  %dfbest_count = ();
  foreach (@algs) {
    $alg = $_;
    $num_fail{$alg} = 0;
    $num_win{$alg} = 0;
    $dfbest{$alg} = 0;
    $dfbest_max{$alg} = 0;
    $dfbest_count{$alg} = 0;
  }
  
  $total_xps=0;
  
  foreach (@NUMSERVICES) {
    $NUMSERVICE=$_;
    foreach (@NDIMS) {
      $NDIM=$_;
      foreach (@SLAPROBS) {
        $SLAPROB=$_;
        $pattern = "LARGE_$NDIM"."_".$NDIM."/problem-$NDIM-$NDIM-64-$NUMSERVICE-*-*-*-$SLAPROB-.5-*-*.result";
        print STDERR "$pattern:";
        @lines = `cat $pattern`;
        print STDERR ($#lines+1)." lines of data...";
        $current_line = 0;
        while ($current_line <= $#lines) {
          %yields = (); 
          ## Gather data for an xp
          while (!($lines[$current_line] eq "##########\n")) {
            @tokens = split(/\|/,$lines[$current_line]);
            if (is_known_alg($tokens[0])) {
              $yields{$tokens[0]} = $tokens[1];
            } 
            if ($tokens[0] eq "LPBOUND") {
              $lpbound = $tokens[1];
            }
            $current_line++;
          }
          $current_line++;
    
          ($num_successes,$best_yield) = process_xp(%yields);
          ## Process the xp
          if (($lpbound >= 0) && ($num_successes > 0)) {
            $total_xps++;
            foreach (@algs) {
              $alg = $_;
              if (($yields{$alg} == -1.0) && ($num_successes > 0)) {
                $num_fail{$alg}++;
              }
              if ($yields{$alg} != -1.0) {
                if ($best_yield == 0) {
                  $deviation = 0;
                } else {
                  $deviation = 100 * ($best_yield - $yields{$alg}) / $best_yield;
                }
                $dfbest{$alg} += $deviation;
                $dfbest_count{$alg}++;
                if ($deviation > $dfbest_max{$alg}) {
                  $dfbest_max{$alg} = $deviation;
                }
                if ($deviation == 0) {
                  $num_win{$alg}++;
                }
              }
            }
          }
        }
        print STDERR "processed!\n";
      }    
    }
  }
  print STDERR "\n";
  
  ### Finalize values
  foreach (@algs) {
    $dfbest{$_} /= $dfbest_count{$_};
    $num_win{$_} = 100*$num_win{$_}/$total_xps;
    $num_fail{$_} = 100*$num_fail{$_}/$total_xps;
  }
  
  # Print output
  print "  NDIMS:     "; foreach (@NDIMS) { print $_." ";} print "\n";
  print "  SLAPROBS:  "; foreach (@SLAPROBS) { print $_." ";} print "\n";
  print "  #SERVICES: "; foreach (@NUMSERVICES) { print $_." ";} print "\n";
  print STDOUT "  NON-DOMINATED ALGORITHMS:\n";
  print STDOUT "  <alg> <\%fail> <ave\%dfb> <max\%dfb> <\%win>\n";
  foreach (@algs) {
    if (is_dominated($_)) {
      print STDOUT "    (";
    }
    printf(STDOUT "  %s %.2f %.2f %.2f %.2f\n",$_,$num_fail{$_},$dfbest{$_},$dfbest_max{$_},$num_win{$_});
  }
}
 
###### HELPER SUBROUTINES #######

## is_dominated()
sub is_dominated {
  my ($alg) = @_;

  foreach (@algs) {
    $other_alg = $_;
    if ( 
         (($dfbest{$other_alg} <= $dfbest{$alg}) &&
          ($num_fail{$other_alg} < $num_fail{$alg})) ||
         (($dfbest{$other_alg} < $dfbest{$alg}) &&
          ($num_fail{$other_alg} <= $num_fail{$alg})) ) {
#         ($dfbest_max{$other_alg} < $dfbest_max{$alg}))  {
#         ($num_win{$other_alg} > $num_win{$alg}) ) {
     return 1;
   }
  }
  return 0;
}

## process_xp
sub process_xp {
  my (%yields) = @_;
  my $success_count = 0;
  my $max_yield = -1;

  foreach (@algs) {
    if ($yields{$_} != -1.0) {
      $success_count++;
    }
    if ($yields{$_} > $max_yield) {
      $max_yield = $yields{$_};
    }
  }
  return ($success_count, $max_yield);
}

## is_known_alg 
sub is_known_alg {
  my ($alg_name) = @_;
  foreach (@algs) {
    if ($_ eq $alg_name) {
      return 1;
    }
  }
  return 0;
}



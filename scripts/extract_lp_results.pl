#!/usr/bin/perl -w

@algs=("RRND","RRNZ","FASTDIVING","SLOWDIVING","METAGREEDY");

print STDOUT "## OVERALL RANKING ##\n";
#@NDIMS = ("1","2","3");
@NDIMS = ("1");
@SLAPROBS = ("0.0","0.25","0.5");
#@NUMSERVICES= ("100","200","500");
@NUMSERVICES= ("100","200");
process();

exit(0);


### MAIN SUBROUTINE
sub process
{
  # number of failures when somebody else found a solution
  %num_fails = ();
  %dfbound = ();
  %numsuccesses = ();
  %num_wins = ();
  %all_dfbound = ();
  foreach (@algs) {
    $alg = $_;
    $num_fails{$alg} = 0;
    $dfbound{$alg} = 0;
    $numsuccesses{$alg} = 0;
    $num_wins{$alg} = 0;
    $num_fails{$alg} = 0;
  }


  for ($i=0; $i<=$#algs; $i++) {
    for ($j=0; $j<=$#algs; $j++) {
	$higher_yield_than[$i][$j] = 0;
	$succeeds_when_other_fails[$i][$j] = 0;
	$num_both_succeed[$i][$j] = 0;
    }
  }

  $num_unfeasible=0;  
  $num_no_success=0;  
  $num_some_successes=0;
  
  foreach (@NUMSERVICES) {
    $NUMSERVICE=$_;
    foreach (@NDIMS) {
      $NDIM=$_;
      foreach (@SLAPROBS) {
        $SLAPROB=$_;
        $pattern = "LARGE_$NDIM"."_".$NDIM."_ALL/problem-$NDIM-$NDIM-64-$NUMSERVICE-*-*-*-$SLAPROB-.5-*-*.result";
        print STDERR "$pattern:";
        @lines = `cat $pattern`;
        print STDERR ($#lines+1)." lines of data...";
        $current_line = 0;
        while ($current_line < $#lines) {
#          print STDERR "Acquiring data for a record ($current_line / $#lines)...";
          %yields = (); 
	  $lpbound_yield = -1;
          ## Gather data for an xp
	  # Look for a record for our algorithms
	  $num_found_algs = 0;
          $num_required_algs = ($#algs+1) + 1; #all + lpbound
          $done = 0;
 	  while($num_found_algs < $num_required_algs) {
#            print STDERR "LINE # $current_line: $lines[$current_line]\n";
            if (index($lines[$current_line],"##") == 0) {
              $current_line++;
              if ($current_line >= $#lines) {
                $done = 1;
                last;
              }
              next;
            }
            @tokens = split(/\|/,$lines[$current_line]);
            if (is_known_alg($tokens[0])) {
  	      if (defined($yields{$tokens[0]})) {
                print STDERR "Incoherent data... \n";
	        print STDERR "   missing result for an alg\n";
		exit(1);
              }
              $yields{$tokens[0]} = $tokens[1];
              $num_found_algs++;
            } elsif ($tokens[0] eq "LPBOUND") {
  	      if ($lpbound_yield != -1) {
                print STDERR "Incoherent data... \n";
	        print STDERR "   missing result for an alg\n";
		exit(1);
              }
              $lpbound_yield = $tokens[1];
              $num_found_algs++;
            }
	
            $current_line++;
          }

          if ($done) {
            last;
          }
#	  print STDERR "acquired\n";

          ##
          #  Do the one-to-one comparison
          #  
          for ($i=0; $i<=$#algs; $i++) {
            for ($j=0; $j<=$#algs; $j++) {
	      $alg1 = $algs[$i];
	      $alg2 = $algs[$j];
	      $yield1 = $yields{$alg1};
	      $yield2 = $yields{$alg2};
   	      if (($yield1 >= 0) &&( $yield2 == -1)) {
	        $succeeds_when_other_fails[$i][$j]++;
              }
              if (($yield1 >= 0) &&( $yield2 >= 0)) {
                $num_both_succeed[$i][$j]++;
              }
              if (($yield1 >= 0) &&( $yield2 >= 0)&&($yield1 > $yield2)) {
	        $higher_yield_than[$i][$j]++;
              }
            }
          }
    
          ## Process the xp
	  # $num_successes: how many algs found a solution
   	  # $best_yield: what's the best yield that was achieved
          ($num_successes,$best_yield) = process_xp(%yields);

	  if ($lpbound_yield <= 0) {
            $num_unfeasible++;
          } else {
            if ($num_successes == 0) {
              $num_no_success++;
            } else {
              $num_some_successes++;
              foreach (@algs) {
                $alg = $_;
                if ($yields{$alg} == -1.0) {
                  $num_fails{$alg}++;
                } else { # Success of the algorithm
    		  $numsuccesses{$alg}++;
                  $dfbound{$alg} += $lpbound_yield - $yields{$alg};
                  $all_dfbound{$alg}[1+$#{ $all_dfbound{$alg} }] = $lpbound_yield - $yields{$alg};
		  if ($yields{$alg} == $best_yield) {
                    $num_wins{$alg}++;
                  }
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

  
  print STDERR "Computing 90th percentile...";
  %percentile = ();
  foreach (@algs) {
    $alg = $_;
    @{$all_dfbound{$alg}} = sort {$a <=> $b} @{$all_dfbound{$alg}};
    $last_element = $#{ $all_dfbound{$alg} };
    $last_element = int(.9 * $last_element);
    $percentile{$alg} = $all_dfbound{$alg}[$last_element];
  }

  print STDERR "conputed\n";

  
  ### Finalize values
  foreach (@algs) {
    if ($numsuccesses{$_} > 0) {
      $dfbound{$_} /= $numsuccesses{$_};
    } else {
      $dfbound{$_} = -1;
    }
    $num_wins{$_} = 100*$num_wins{$_}/$num_some_successes;
    $num_fails{$_} = 100*$num_fails{$_}/$num_some_successes;
  }
  
  # Print output
  print "  NDIMS:     "; foreach (@NDIMS) { print $_." ";} print "\n";
  print "  SLAPROBS:  "; foreach (@SLAPROBS) { print $_." ";} print "\n";
  print "  #SERVICES: "; foreach (@NUMSERVICES) { print $_." ";} print "\n";
  print STDOUT "-------------------------------------------------------\n";
  print STDOUT "Total num xps: ".
     ($num_unfeasible + $num_no_success + $num_some_successes)."\n";
  print STDOUT "Unfeasible xps: $num_unfeasible\n";
  print STDOUT "No success xps: $num_no_success\n";
  print STDOUT "No some successes xps: $num_some_successes\n";
  print STDOUT "-------------------------------------------------------\n";
  print STDOUT "  ALG   ave-d-f-bound  90-percentile % wins  % fails\n";
  foreach (@algs) {
    printf(STDOUT "  %s \t%.2f \t%.2f \t%.2f \t%.2f\n",$_,$dfbound{$_},$percentile{$_},$num_wins{$_},$num_fails{$_});
  }
  print STDOUT "-------------------------------------------------------\n";
  for ($i=0; $i<=$#algs; $i++) {
    $alg1 = $algs[$i];
    for ($j=0; $j<$i; $j++) {
      $alg2 = $algs[$j];
      print STDOUT "$alg1 vs. $alg2:\n";
      print STDOUT "\t$alg1 succeeds when $alg2 fails ".$succeeds_when_other_fails[$i][$j]." times\n";
      print STDOUT "\t$alg2 succeeds when $alg1 fails ".$succeeds_when_other_fails[$j][$i]." times\n";
      print STDOUT "\t$alg1 produces a higher yield than $alg2 when both succeed ".$higher_yield_than[$i][$j]." times\n";
      print STDOUT "\t$alg2 produces a higher yield than $alg1 when both succeed ".$higher_yield_than[$j][$i]." times\n";
    }
  }

  ##### LATEX TABLE GENERATION: Small overall table #####
  print STDOUT "\n";
  print STDOUT "\\begin{tabular}{|l|c|c|c|c|}\n";
  print STDOUT "\\hline\n";
  print STDOUT "Algorithm  & Average & 90th perc. & \\% wins & \\% fails\\\\\n";
  print STDOUT "\\hline\n";
  foreach (@algs) {
    $alg = $_;
    printf( STDOUT "$alg &");
    printf( STDOUT "%.2f &",$dfbound{$alg});
    printf( STDOUT "%.2f &",$percentile{$alg});
    printf( STDOUT "%.2f &",$num_wins{$alg});
    printf( STDOUT "%.2f \\\\\n",$num_fails{$alg});
  }
  print STDOUT "\\hline\n";
  print STDOUT "\\end{tabular}\n";
  print STDOUT "\n";

  ##### LATEX TABLE GENERATION: Large competition table #####
  print STDOUT "\n";
  print STDOUT "\\begin{tabular}{|l|c|c|c|c|c|}\n";
  print STDOUT "\\cline{2-6}\n";
  print STDOUT "\\multicolumn{1}{l|}{} & \\rrnd & \\rrndnz & \\fastdiving & \\slowdiving & \\metagreedy\\\\\n";
  print STDOUT "\\hline\n";
  for ($i=0; $i<=$#algs; $i++) {
    print STDOUT "\t".$algs[$i]." ";
    for ($j=0; $j<=$#algs; $j++) {
      if ($i == $j) {
        print STDOUT "& \\GRAY ";
      } else {
        printf( STDOUT "& %.2f",(100*$succeeds_when_other_fails[$i][$j] / $numsuccesses{$algs[$i]}));
        print STDOUT " /";
        printf( STDOUT " %.2f",(100*$higher_yield_than[$i][$j] / $num_both_succeed[$i][$j]));
      }
    }
    print STDOUT "\\\\\n";
    print STDOUT "\\hline\n";
  }
 print STDOUT "\\end{tabular}\n";
}
  

 
###### HELPER SUBROUTINES #######

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



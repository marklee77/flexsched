#!/usr/bin/perl -w

@slacks=("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9");

@NDIMS = ("1","2","3");
@SLAPROBS = ("0.0","0.25","0.5");
@NUMSERVICES = ("8", "10", "12");

#foreach (@NDIMS) {
  #$NDIM=$_;
  #foreach (@SLAPROBS) {
    #$SLAPROB=$_;
    #foreach (@NUMSERVICES) {
      #$NUMSERVICE=$_;
      #$count = 0;
      #@ave_rel_diffs = ();
      #@max_rel_diffs = ();
      #foreach (@slacks) {
        #$slack = $_;
        #$pattern = "SMALL_$NDIM"."_".$NDIM."/problem-$NDIM-$NDIM-4-$NUMSERVICE-$slack-*-*-$SLAPROB-*-*-*.result";
#

$ave_string=		"ave =           [ ";
$max_string=		"max =           [ ";
$fail_string=		"fail =          [ ";
$failcount_string=	"fail_count =    [ ";
$percentile_75_string=	"percentile_75 = [ ";
$percentile_90_string=	"percentile_90 = [ ";
$percentile_99_string=	"percentile_99 = [ ";
foreach (@slacks) {
  $slack = $_;
  $NDIM = "3";
  $NUMSERVICES="8";
  $SLAPROB="0.5";
  $SIGMA="2.0";
  $pattern = "SMALL_$NDIM"."_".$NDIM."/problem-$NDIM-$NDIM-4-$NUMSERVICES-$slack-*-$SIGMA-$SLAPROB-*-*-*.result";

#  print STDERR $pattern."\n";
($ave_rel_diff, $max_rel_diff, $fail_count, $percentage_fails, $percentile_75, $percentile_90, $percentile_99) = 
          process_data($pattern);
$ave_string.="$ave_rel_diff ";
$max_string.="$max_rel_diff ";
$fail_string.="$percentage_fails ";
$failcount_string.="$fail_count ";
$percentile_75_string.="$percentile_75 ";
$percentile_90_string.="$percentile_90 ";
$percentile_99_string.="$percentile_99 ";

}
$ave_string.="];";
$max_string.="];";
$fail_string.="];";
$failcount_string.="];";
$percentile_75_string.="];";
$percentile_90_string.="];";
$percentile_99_string.="];";

print STDOUT "\n";
print STDOUT $ave_string."\n";
print STDOUT $max_string."\n";
print STDOUT $percentile_75_string."\n";
print STDOUT $percentile_90_string."\n";
print STDOUT $percentile_99_string."\n";
print STDOUT $fail_string."\n";
print STDOUT $failcount_string."\n";

$slack_string = "slack = [ ";
foreach(@slacks) {
  $slack_string .= $_." ";
}
$slack_string .= "];";

print STDOUT $slack_string."\n";

print STDOUT "hold on\n";
print STDOUT "plot(slack, ave);\n";
print STDOUT "plot(slack, max);\n";
print STDOUT "plot(slack, percentile_75);\n";
print STDOUT "plot(slack, percentile_90);\n";
print STDOUT "plot(slack, percentile_99);\n";


sub process_data {
  ($pattern) = @_;

  $fail_count = 0;
  $total_count = 0;
  $ave_rel_diff = 0;
  $max_rel_diff = 0;
  $percentage_fails = 0;
  print STDOUT ".";
  @lines = `cat $pattern`;
  $current_line = 0;
  @rel_diffs = ();
  while ($current_line <= $#lines) {
    $total_count++;
    @tokens = split(/\|/,$lines[$current_line++]);
    $lpbound = -2;
    $milp = -2;
    if ($tokens[0] eq "LPBOUND") {
    $lpbound = $tokens[1];
    } else {
      $milp = $tokens[1];
    }
    @tokens = split(/\|/,$lines[$current_line++]);
    if ($tokens[0] eq "LPBOUND") {
      $lpbound = $tokens[1];
    } else {
      $milp = $tokens[1];
    }
    if (($lpbound == -2) || ($milp == -2)) {
      print STDERR "Missing Results!!!\n";
      exit(1);
    }
    if ($milp < 0) {
      $percentage_fails++;
      $fail_count++;
    }
    if (($milp < 0) || ($lpbound < 0)) {
      next;
    }
    $rel_diff = 100*($lpbound - $milp) / $lpbound;
    $rel_diffs[$#rel_diffs+1] = $rel_diff;
  }

  $percentage_fails = $percentage_fails * 100 / $total_count;

  if ($#rel_diffs == -1) {
    return (-1,-1,$fail_count,$percentage_fails,-1,-1,-1);
  }

  # Compute the average and the max
  $ave_rel_diff = 0;
  $max_rel_diff = 0;
  foreach (@rel_diffs) {
    $value = $_;
    $ave_rel_diff += $value;
    if ($value > $max_rel_diff) {
      $max_rel_diff = $value;
    }
  }
  $ave_rel_diff /= ($#rel_diffs + 1);

  # Compute the XXX percentile
  @rel_diffs = sort {$a <=> $b} @rel_diffs;
  $percentile_75 = $rel_diffs[$#rel_diffs*.75];
  $percentile_90 = $rel_diffs[$#rel_diffs*.90];
  $percentile_99 = $rel_diffs[$#rel_diffs*.99];
 
  $ave_rel_diff = (int($ave_rel_diff * 100) / 100);
  $max_rel_diff = (int($max_rel_diff * 100) / 100);
  $percentage_fails = (int($percentage_fails * 100) / 100);
  $percentile_75 = (int($percentile_75 * 100) / 100);
  $percentile_90 = (int($percentile_90 * 100) / 100);
  $percentile_99 = (int($percentile_99 * 100) / 100);
  return ($ave_rel_diff, $max_rel_diff, $fail_count, $percentage_fails, $percentile_75, $percentile_90, $percentile_99);
}


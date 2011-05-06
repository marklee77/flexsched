#!/usr/bin/perl -w

$NDIM=1;

$prefix = "LARGE_".$NDIM."_".$NDIM;

unlink($prefix."_ALL");
mkdir($prefix."_ALL");

chdir($prefix."_navet/");
@all = glob("*.result");
chdir("..");
print STDERR ($#all+1)."\n";

foreach (@all) {
  $filename = $_;


  $file_new =    $prefix."_ALL/$filename";
  $file_navet =  $prefix."_navet/$filename";
  $file_lani =   $prefix."_lani/$filename";
  $file_breeze = $prefix."_breeze/$filename";

  $command_line = "cat $file_navet ";

  if (-r $file_lani) {
    $command_line .= " $file_lani ";
  } else {
    #print STDOUT "No file $file_lani\n";
  }

  if (-r $file_breeze) {
    $command_line .= " $file_breeze ";
  } else {
    #print STDOUT "No file $file_breeze\n";
  }
 
  # Create the new file
  $command_line.= " >> $file_new";
# print STDOUT $command_line."\n";

  `$command_line`;
  
}

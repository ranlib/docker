#!/usr/local/bin/perl

use TIGR::Foundation;
use AMOS::ParseFasta;

$tf = new TIGR::Foundation;

if (!defined $tf){
    die ("Bad foundation\n");
}

if ($ARGV[0] eq '-h') {
  die "Report the length of the sequences in a FASTA file. If no file is provided, standard input is used\nUSAGE: $0 [fasta_file1] [fasta_file2] ... [fasta_filen] \n";
}

sub report
{
  my $fr = shift;

  die ("Bad reader\n")
    if (!defined $fr);

  while (($head, $body) = $fr->getRecord()){
      $head =~ /(\S+)/;
      $id = $1;
      print "$id ", length($body), "\n";
  }
}


if (scalar @ARGV == 0)
{
  $fr = new AMOS::ParseFasta(\*STDIN);
  report($fr);
}
else
{
  foreach my $file (@ARGV)
  {
    open(IN, $file) || $tf->bail("Cannot open $file: $!\n");
    $fr = new AMOS::ParseFasta(\*IN);
    report($fr);
  }
}


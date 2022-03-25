#!/usr/bin/perl

use TIGR::Foundation;
use AMOS::ParseFasta;

$tf = new TIGR::Foundation;

if (!defined $tf){
    die ("Bad foundation\n");
}

if ((scalar @ARGV == 0) || ($ARGV[0] eq '-h') || ($ARGV[0] =~ /\D/)) {
  die "Select sequences at least the specified length\nUSAGE: $0 min_len [fasta_file1] [fasta_file2] ... [fasta_filen] \n";
}

my $MIN_LEN = shift @ARGV;

sub report
{
  my $fr = shift;

  die ("Bad reader\n")
    if (!defined $fr);

  while (($head, $body) = $fr->getRecord())
  {
      if (length($body) >= $MIN_LEN)
      {
        print ">$head\n$body\n";
      }
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


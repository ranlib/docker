#!/usr/bin/perl

my $USAGE = "fixfastq in.fq > out.fq\n";

if ((scalar @ARGV > 0) && ($ARGV[0] eq "-h"))
{
  print $USAGE;
  print "  Chop the sequence or quality string to be the minimum of the two.\n";
  print "  Also replace any ambiguity codes in the reads with N\n";
  exit 0;
}

my $trims = 0;
my $trimq  = 0;

my $trimsb = 0;
my $trimqb = 0;

while (<>)
{
  print $_;
  my $s = <>;  chomp $s;
  my $h2 = <>; chomp $h2;
  my $q = <>;  chomp $q;

  ## clean up ambiguity codes
  $s =~ s/[^ACGTacgtN]/N/g;

  ## trim sequences or quality strings
  
  my $sl = length($s);
  my $ql = length($q);

  if ($sl < $ql)
  {
    $trimq++;
    $trimqb += $ql-$sl;
    $q = substr($q,0,$sl);
  }
  elsif ($ql < $sl)
  {
    $trims++;
    $trimsb += $sl-$ql;
    $s = substr($s,0,$ql);
  }

  print "$s\n$h2\n$q\n";
}

print STDERR "trimseq: $trims $trimsb bp  trimq: $trimq $trimqb bp\n";

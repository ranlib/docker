#!/usr/bin/perl -w
use strict;

sub printGC
{
  my $seqname = shift;
  return if !defined $seqname;

  my $l = shift;
  my $g = shift;
  my $a = shift;

  my $ratio = sprintf("%.02f", ($g+$a) ? ($g*100 / ($g + $a)) : 0); 

  print "$seqname $l $ratio\n";
}

my $seqname = undef;

my $g = 0;
my $a = 0;
my $l = 0;

while (<>)
{
  if (/^>(\S+)/)
  {
    printGC($seqname, $l, $g, $a);

    $seqname = $1;

    $l = 0;
    $g = 0;
    $a = 0;
  }
  else
  {
    chomp;

    $l += length($_);
    $g += tr/gGcC//;
    $a += tr/aAtT//;
  }
}

printGC($seqname, $l, $g, $a);

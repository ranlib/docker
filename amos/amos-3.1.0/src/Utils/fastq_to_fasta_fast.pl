#!/usr/bin/perl 

## Warning no santity checking, just blindly outputs the name and sequence

while (<>)
{
  my $h1 = $_;
  my $s  = <>;
  my $h2 = <>;
  my $q  = <>;

  die "ERROR: expected '@' but saw $h1" if substr($h1, 0, 1) != '@';

  print ">", substr($h1,1);
  print $s;
}

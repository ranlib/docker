#!/usr/bin/perl
# $Id: revFasta.pl,v 1.2 2004/01/21 14:33:16 shumwaym Exp $ 
#
# revFasta - Reverse complement the fasta file or specified record.
#
# Written by -  Martin Shumway
#
my $HELPTEXT = qq~
Reverse complement the fasta file or specified record in a multifasta file.

  revFasta <fasta file>  [options]
    options:
      -i <id>      Reverse complement the specified id only

~;
#
# (C) Copyright 2001-2003  The Institute for Genomic Research (TIGR)
#     All rights reserved.
#
#

# =============================== Pragmas and imports ======================
use strict;
use TIGR::Foundation;
use TIGR::FASTAreader;
use TIGR::FASTArecord;

# Normally used only in testing
# use warnings;

# ================================ Globals ====================================

my $MY_VERSION = " Version 1.0 (Build " . (qw/$Revision: 1.2 $/ )[1] . ")";
my @MY_DEPENDS =
(
  "TIGR::Foundation",
  "TIGR::FASTAreader",
  "TIGR::FASTArecord",
);

# Reference to TF object
my $tf = new TIGR::Foundation;

# ================================= Procedures ===================================
MAIN:
{
  # == Program Setup ==
  $tf->addDependInfo(@MY_DEPENDS);
  $tf->setHelpInfo($HELPTEXT);
  $tf->setVersionInfo($MY_VERSION);
  $tf->setUsageInfo($HELPTEXT);

  # Handle the input options
  my $id = undef;
  my $result  = $tf->TIGR_GetOptions
                (
                  "i=s",   \$id,
                );
  $tf->bail("Command line parsing failed") if ($result == 0);
  $tf->printUsageInfoAndExit() if ($#ARGV != 0);
  my $fastafile = $ARGV[0] if (defined $ARGV[0]);
  $tf->bail("Could not open \'$fastafile\' ($!)") if (! -r $fastafile);

  # == Gather Inputs ==
  $tf->logLocal("Reading fasta records from $fastafile", 9);
  my @errors = ();
  my $fr = new TIGR::FASTAreader($tf, \@errors, $fastafile);
  $tf->bail("Error creating FASTA Reader object") if (! defined $fr);
  if ($#errors >= 0)
  {
    for my $e (@errors)
    {
      logerr("$e\n");
    }
    $tf->bail("Invalid FASTA on input from $fastafile");
  }

  # == Write outputs ==
  $tf->logLocal("Writing output...", 9);
  while ( $fr->hasNext() )
  {
    my $r = $fr->next();
    my $a = $r->getIdentifier();
    my $h = $r->getHeader();
    next if (defined $id  &&  $a ne $id);
    $tf->logLocal("Found contig $a", 9);
    my $data = $r->reverseComplementData();
    my $rr = new TIGR::FASTArecord($h, $data);
    print $rr->toString();
  }

  exit 0;
}

#! /usr/bin/env perl

# blat2nucmer is a script to convert BLAT's .psl files to the nucmer format
# Author: Nagarjun Vijay

# From http://genome.ucsc.edu/goldenPath/help/blatSpec.html :
# A .psl file describes a series of alignments in a dense easily parsed text
# format.  It begins with a five line header which describes each field. 
# Following this is one line for each alignment with a tab between each field. 
# The fields are describe below in  a format suitable for many relational
# databases.
#
#   matches int unsigned ,      # Number of bases that match that aren't repeats
#   misMatches int unsigned ,   # Number of bases that don't match
#   repMatches int unsigned ,   # Number of bases that match but are part of repeats
#   nCount int unsigned ,       # Number of 'N' bases
#   qNumInsert int unsigned ,   # Number of inserts in query
#   qBaseInsert int unsigned ,  # Number of bases inserted in query
#   tNumInsert int unsigned ,   # Number of inserts in target
#   tBaseInsert int unsigned ,  # Number of bases inserted in target
#   strand char(2) ,            # + or - for query strand, optionally followed by + or – for target strand
#   qName varchar(255) ,        # Query sequence name
#   qSize int unsigned ,        # Query sequence size
#   qStart int unsigned ,       # Alignment start position in query
#   qEnd int unsigned ,         # Alignment end position in query
#   tName varchar(255) ,        # Target sequence name
#   tSize int unsigned ,        # Target sequence size
#   tStart int unsigned ,       # Alignment start position in target
#   tEnd int unsigned ,         # Alignment end position in target
#   blockCount int unsigned ,   # Number of blocks in alignment. A block contains no gaps.
#   blockSizes longblob ,       # Size of each block in a comma separated list
#   qStarts longblob ,          # Start of each block in query in a comma separated list
#   tStarts longblob ,          # Start of each block in target in a comma separated list
#
# In general the coordinates in psl files are “zero based half open.” The first
# base in a sequence is numbered zero rather than one. When representing a range
# the end coordinate is not included in the range. Thus the first 100 bases of a
# sequence are represented as 0-100, and the second 100 bases are represented as
# 100-200. There is a another little unusual feature in the .psl format. It has
# to do with how coordinates are handled on the negative strand. In the qStart/
# qEnd fields the coordinates are where it matches from the point of view of the
# forward strand (even when the match is on the reverse strand). However on the
# qStarts[] list, the coordinates are reversed.
#
# Here's an example of a 30-mer that has 2 blocks that align on the minus strand
# and 2 blocks on the plus strand (this sort of stuff happens in real life in
# response to assembly errors sometimes).
#
# 0         1         2         3 tens position in query
# 0123456789012345678901234567890 ones position in query
#             ++++          +++++ plus strand alignment on query
#     --------    ----------      minus strand alignment on query
#
# Plus strand:
#     qStart 12 qEnd 31 blockSizes 4,5 qStarts 12,26
# Minus strand:
#     qStart 4 qEnd 26 blockSizes 10,8 qStarts 5,19
#
# Essentially the minus strand blockSizes and qStarts are what you would get if
# you reverse complemented the query.However the qStart and qEnd are non-
# reversed. To get from one to the other:
#     qStart = qSize - revQEnd
#     qEnd = qSize - revQStart

use strict;
use warnings;
use TIGR::Foundation;

my $VERSION = '$Revision$ ';
my $HELP = qq~
.USAGE.
  blat2nucmer -i psl_file

.DESCRIPTION.
  psl2nucmer converts BLAT .psl files into a tab-delimited format that nucmer
  can use. The results are printed on screen.

.KEYWORDS.
  converter, nucmer, blat

~;

my $base = new TIGR::Foundation();
if (! defined $base) {
    die("Error: Could not initialize foundation.\n");
}

$base->setVersionInfo($VERSION);
$base->setHelpInfo($HELP);

my $pslfile;

my $err = $base->TIGR_GetOptions(
    "i=s"       => \$pslfile,
);


if (! defined $pslfile){
    die "You must specify an input .psl file with option -i\n";
}

blat2nucmer($pslfile);

exit;




sub blat2nucmer {
  my ($pslfile) = @_;
  open my $in, "<", $pslfile or die $!;
  while ( my $z = <$in> ) {
    my @values = split(/\t/, $z);
    	if($values[8] =~ /\+/){#writing positive strand as is
    print
      ($values[15]+1)."\t".#zero to one based
      $values[16]."\t|\t".
      ($values[11]+1)."\t".#zero to one based
      $values[12]."\t|\t".
      abs($values[16]-$values[15])."\t".
      abs($values[12]-$values[11])."\t|\t".
      (sprintf "%.2f",(($values[0]+$values[2]+$values[3])*100)/abs($values[12]-$values[11]))."\t|\t".
      $values[14]."\t".
      $values[10]."\t|\t".
      (sprintf "%.2f",(abs($values[16]-$values[15])*100/$values[14]))."\t".
      (sprintf "%.2f",(abs($values[12]-$values[11])*100/$values[10]))."\t|\t".
      $values[13]."\t".
      $values[9]."\n";
	}
	elsif($values[8] =~ /\-/){#for negative strand writing qry coordinates in reverse
    print
      ($values[15]+1)."\t".#zero to one based
      $values[16]."\t|\t".
      $values[12]."\t|\t".
      ($values[11]+1)."\t".#zero to one based
      abs($values[16]-$values[15])."\t".
      abs($values[12]-$values[11])."\t|\t".
      (sprintf "%.2f",(($values[0]+$values[2]+$values[3])*100)/abs($values[12]-$values[11]))."\t|\t".
      $values[14]."\t".
      $values[10]."\t|\t".
      (sprintf "%.2f",(abs($values[16]-$values[15])*100/$values[14]))."\t".
      (sprintf "%.2f",(abs($values[12]-$values[11])*100/$values[10]))."\t|\t".
      $values[13]."\t".
      $values[9]."\n";
	}
	else{#few blocks match on positive and few on negative strand. Ignore these hits
	}
  }
  close $in;
  return 1;
}




#!/usr/bin/env perl

use strict;
use Getopt::Long;
use FileHandle;

my $infile;
my $reffile;
my $reference;
my $length=0;
my $alignment;

GetOptions(
   'i=s'      => \$infile,
   'r=s'      => \$reffile,
   'a=f'      => \$alignment,
   'h|help'   => sub{usage()},
);

read_reference($reffile);
read_gaps($infile);

sub read_reference
{
my ($header,@seq,$sequence);
my $fh=FileHandle->new($reffile)|| die "$!";
if ($fh->open("< $reffile")){
   $/=">";
   while(<$fh>){
   $_=~ s/\>//g;
   unless ($_){next;};
   ($header,@seq)=split /\n/,$_;
   $sequence=join "",@seq;
   $reference=length $sequence;
   }
}
}

sub read_gaps
{
my $line;
my $fh=FileHandle->new($infile)|| die "$!";
if ($fh->open("< $infile")){
   $/="\n";
   while (<$fh>){
      unless ($_){next;};
      my ($ref,$start,$end,$len,$query)=split /\t/,$_;
      $length+=$len;
   }
}
my $difference =$length/$reference;
my $aln_threshold = (100-$alignment)/100;
if ($difference < $aln_threshold){print "0";}
if ($difference >= $aln_threshold){print "1";}
}

sub usage
{
print STDERR"
Usage:
   $0 -i <gaps file> -r <reference file> 
";
exit;
}

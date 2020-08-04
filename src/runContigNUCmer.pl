#!/usr/bin/env perl

######################################################
# Written by Sanaa Ahmed
# Nov. 30, 2012

# Given a directory containing fasta files, runs nucmer
# Asks for a reference,
#  If reference not given,
#  picks one from the files randomly.
######################################################

use strict;
use FindBin qw($RealBin);
use Getopt::Long;
use File::Basename;
use Parallel::ForkManager;
use lib "$RealBin/../lib/";
use phame::nucmer;

# require "./runNUCmer.pl";

# set up environments
$ENV{PATH} = "$RealBin:$RealBin/../ext/bin:$ENV{PATH}";

my $breaklen   = 200;
my $mincluster = 65;
my $diagfactor = 0.12;
my $maxgap     = 90;
my $minmatch   = 20;
my $identity   = 0;
my $gap_cutoff = 0;
my ( $working_dir, $reference, $list, $prefix );
my @query_list;
my @reference_list;
my $snp_indel;
my $snp_n;
my $indel_n;
my $gaps;
my $ref_gaps;
my $query_gaps;
my $type;
my $thread = 2;
my $outdir = `pwd`;
$outdir =~ s/\n//;

GetOptions(
    'q|dir=s'       => \$working_dir,
    'r|reference=s' => \$reference,
    'd|outdir=s'    => \$outdir,
    't|thread=i'    => \$thread,
    'l|list=s'      => \$list,
    'y|type=i'      => \$type,
    'h|help'        => sub { usage() }
);

&usage() unless ($working_dir);

my $max_thread =
  ( $^O =~ /darwin/ )
  ? `sysctl hw.ncpu | awk '{print \$2}'`
  : `grep -c ^processor /proc/cpuinfo`;
if ( $thread < 1 || $thread > $max_thread ) {
    die("-thread value must be between than 1 and $max_thread.\n");
}

if ( !-e $outdir )        { mkdir "$outdir"; }
if ( $outdir =~ /.+\/$/ ) { my $temp = chop($outdir); }
chdir $outdir;
if ( $working_dir =~ /.+\/$/ ) { my $temp = chop($working_dir); }
my $query_dir = $working_dir . '/files';

if ( !-e "gaps" )   { mkdir "gaps"; }
if ( !-e "snps" )   { mkdir "snps"; }
if ( !-e "stats" )  { mkdir "stats"; }
if ( !-e "nucmer" ) { mkdir "nucmer"; }

my $options =
"--maxmatch -b $breaklen -c $mincluster -d $diagfactor -g $maxgap -l $minmatch ";
$gap_cutoff = "-l $gap_cutoff";
$identity   = "-i $identity";

my %queries;
open( IN, $list );
while (<IN>) {
    chomp;
    $queries{$_}++;
    $queries{$_}++;
    $queries{$_}++;
}
close IN;

read_directory($query_dir);
if ( $type == 1 ) { run_all_nucmer( @reference_list, @query_list ); }
if ( $type == 2 ) { run_ref_nucmer(@query_list); }

#print "NUCmer on all contig/draft genomes complete.\n\n";
cleanup();

sub read_directory {
    my $dir = shift;
    my ( $name, $path, $suffix ) = fileparse( "$reference", qr/\.[^.]*/ );
    print "\nUsing given reference:  $reference\n";

    opendir( PARENT, $dir ) || die "Directory $dir does not exist!";
    while ( my $files = readdir(PARENT) ) {
        next if ( $files =~ /^..?$/ );
        if ( $files =~ /(.+_contig)/ ) {
            if ( exists $queries{$1} ) {
                my $query = $dir . '/' . $files;
                push( @query_list, $query );
            }
        }
        else {
            my $ref = $dir . '/' . $files;
            push( @reference_list, $ref );
        }
    }
    closedir(PARENT);
}

sub run_ref_nucmer {
    my $pm = new Parallel::ForkManager($thread);
    $pm->run_on_finish(
        sub { my ( $pid, $exit_code, $ident, $exit_signal, $core_dump ) = @_; }
    );

    for ( my $i = 0 ; $i <= $#query_list ; $i++ ) {
        $pm->start($i) and next;
        my %hash;
        my ( $name, $path, $suffix ) = fileparse( "$reference", qr/\.[^.]*/ );
        my ( $qname, $qpath, $qsuffix ) =
          fileparse( "$query_list[$i]", qr/\.[^.]*/ );
        $prefix = "$name\_$qname";
        my $query = $query_list[$i];
        nucmer::run_nucmer(
            {
                outdir            => $outdir,
                prefix            => $prefix,
                ref_fasta         => $reference,
                full_genome_fasta => $query,
                options           => $options
            }
        );

        nucmer::run_delta_filter_snp(
            {
                outdir   => $outdir,
                identity => $identity,
                prefix   => $prefix
            }
        );
        nucmer::run_show_snps(
            {
                outdir => $outdir,
                prefix => $prefix
            }
        );

        my $snp_indel = `SNP_INDEL_count.pl $outdir/'nucmer'/$prefix.snps`;
        $snp_indel =~ s/\n//;
        ( $snp_n, $indel_n ) = split /\t/, $snp_indel;

        nucmer::run_delta_filter_gap(
            {
                outdir   => $outdir,
                identity => $identity,
                prefix   => $prefix
            }

        );

        nucmer::run_show_coords(
            {
                outdir => $outdir,
                prefix => $prefix
            }
        );

        my $gaps1 =
`parseGapsNUCmer.pl $gap_cutoff $outdir/'nucmer'/$prefix.coords 2>/dev/null`;
        ( $ref_gaps, $query_gaps, undef ) = split /\n/, $gaps1;

        ( $ref_gaps, $query_gaps, undef ) = split /\n/, $gaps;

        open( OUT,  ">$outdir/$prefix\_snp_indel.txt" ) || die "$!";
        open( OUT1, ">$outdir/$prefix\_gaps.txt" )      || die "$!";
        print OUT $snp_indel;
        print OUT1 $gaps;
        close OUT;
        close OUT1;

        $pm->finish(0);
    }
    $pm->wait_all_children;
    print "NUCmer on all contig/draft genomes complete.\n\n";
}

sub run_all_nucmer {
    my $pm = new Parallel::ForkManager($thread);
    $pm->run_on_finish(
        sub { my ( $pid, $exit_code, $ident, $exit_signal, $core_dump ) = @_; }
    );

    for ( my $i = 0 ; $i <= $#reference_list ; $i++ ) {
        for ( my $j = 0 ; $j <= $#query_list ; $j++ ) {
            $pm->start() and next;
            my ( $name, $path, $suffix ) =
              fileparse( "$reference_list[$i]", qr/\.[^.]*/ );
            my ( $qname, $qpath, $qsuffix ) =
              fileparse( "$query_list[$j]", qr/\.[^.]*/ );

            $prefix = "$name\_$qname";
            my $query     = $query_list[$j];
            my $reference = $reference_list[$i];
            nucmer::run_nucmer(
                {
                    outdir            => $outdir,
                    prefix            => $prefix,
                    ref_fasta         => $reference,
                    full_genome_fasta => $query,
                    options           => $options
                }
            );
            nucmer::run_delta_filter_snp(
                {
                    outdir   => $outdir,
                    identity => $identity,
                    prefix   => $prefix
                }
            );
            nucmer::run_show_snps(
                {
                    outdir => $outdir,
                    prefix => $prefix
                }
            );

            my $snp_indel = `SNP_INDEL_count.pl $outdir/'nucmer'/$prefix.snps`;
            $snp_indel =~ s/\n//;
            ( $snp_n, $indel_n ) = split /\t/, $snp_indel;

            nucmer::run_delta_filter_gap(
                {
                    outdir   => $outdir,
                    identity => $identity,
                    prefix   => $prefix
                }

            );
            nucmer::run_show_coords(
                {
                    outdir => $outdir,
                    prefix => $prefix
                }
            );

            my $gaps =
`parseGapsNUCmer.pl $gap_cutoff $outdir/'nucmer'/$prefix.coords 2>/dev/null`;
            ( $ref_gaps, $query_gaps, undef ) = split /\n/, $gaps;
            my $check =
`checkNUCmer.pl -i $outdir/'nucmer'/$prefix.gaps -r $reference_list[$i]`;

            if ( $check == 1 ) {
                print
"$query_list[$j] aligned < 25% of the $reference_list[$i] genome\n";
            }

            ( $ref_gaps, $query_gaps, undef ) = split /\n/, $gaps;

            open( OUT,  ">$outdir/$prefix\_snp_indel.txt" ) || die "$!";
            open( OUT1, ">$outdir/$prefix\_gaps.txt" )      || die "$!";
            print OUT $snp_indel;
            print OUT1 $gaps;
            close OUT;
            close OUT1;

            $pm->finish(0);
        }
    }
    $pm->wait_all_children;

    print "NUCmer on all contig/draft genomes complete.\n\n";
}

sub cleanup {
    `mv $outdir/*.snps $outdir/snps`;

    # `mv $outdir/*.gaps $outdir/nucmer`;
    `mv $outdir/*_snp_indel.txt $outdir/nucmer`;
    `mv $outdir/*_gaps.txt $outdir/nucmer`

      # `mv $outdir/*_gaps.txt $outdir/stats`;
      # `mv $outdir/*.coords $outdir/stats`;
}

sub usage {
    print <<Usage;
$0 -r <reference_fasta> -q <query_fasta> -d <output_directory>

Usage
    exit;
}

#!/usr/bin/env perl

######################################################
# Written by Sanaa Ahmed
#        Migun Shakya

# Given a directory containing fasta files, runs nucmer
#  Runs nucmer only on files ending in .fna
#  Runs nucmer on all against all

######################################################

use strict;
use warnings;
use FindBin qw($RealBin);
use Getopt::Long;
use File::Basename;
use Parallel::ForkManager;
use lib "$RealBin/../lib/";
use lib "$RealBin/../ext/lib/perl5";
use phame::PhaME;
use phame::nucmer;
use misc_funcs;
use Cwd 'abs_path';

# set up environments
$ENV{PATH} = "$RealBin:$RealBin/../ext/bin:$ENV{PATH}";

my $ref_genome;
my $breaklen        = 200;
my $mincluster      = 65;
my $diagfactor      = 0.12;
my $maxgap          = 90;
my $minmatch        = 20;
my $identity        = 0;
my $gap_cutoff      = 0;
my $repeat_identity = 97;
my $len_cutoff      = 100;
my $cutoff          = 0;
my ( $query_dir, $thread, $list, $code );
my @query;
my $outdir = `pwd`;
$outdir =~ s/\n//;
my $options = "--maxmatch ";
my $check;
my $gaps2;
my $buildSNPdb;
my $prefix1;
my $first_name;
my $second_name;
my $first_fasta;
my $second_fasta;
my $delta_file1;
my $gfilt_f;

GetOptions(
    'r|ref_genome=s' => \$ref_genome,     # string option to specify ref file
    'b|buildSNPdb=i' => \$buildSNPdb,
    'q|querydir=s'   => \$query_dir,      # string option to specifu query dir
    'd|outdir=s'     => \$outdir,         # output directory
    't|thread=i'     => \$thread,         # integer, number of threads
    'l|list=s'       => \$list,           # string, file with list of fasta
    'c|code=s'       => \$code,
    'f|cutoff=s'     => \$cutoff,         #
    'h|help'         => sub { usage() }
);

&usage() unless ($query_dir);

# Get full path of everything
$query_dir = Cwd::abs_path($query_dir);
$outdir    = Cwd::abs_path($outdir);
if ( !-e $outdir ) { mkdir $outdir; }

my $max_thread =
  ( $^O =~ /darwin/ )
  ? `sysctl hw.ncpu | awk '{print \$2}'`
  : `grep -c ^processor /proc/cpuinfo`;
if ( $thread < 1 || $thread > $max_thread ) {
    die("-thread value must be between 1 and $max_thread.\n");
}

# Change directory.
chdir $outdir;
if ( !-e "gaps" )   { mkdir "gaps"; }
if ( !-e "snps" )   { mkdir "snps"; }
if ( !-e "stats" )  { mkdir "stats"; }
if ( !-e "nucmer" ) { mkdir "nucmer"; }
if ( !-e "temp" )   { mkdir "temp"; }

if ( $code !~ /virus/ ) {
    my $options =
"--maxmatch -b $breaklen -c $mincluster -d $diagfactor -g $maxgap -l $minmatch ";
}

$gap_cutoff = "-l $gap_cutoff";
$identity   = "-i $identity";

my $snp_indel;
my $snp_n;
my $indel_n;
my $gaps;
my $ref_gaps;
my $query_gaps;

#$ENV{PATH}= "$bindir:/usr/local/bin:/usr/bin:/bin:/opt/apps/bin";

########################## Main ################################################

my $sketch_output = File::Spec->catfile( $query_dir, "sketch_output.txt" );

# print "$sketch_output\n";
my @remove_list;

@query = misc_funcs::read_directory( $query_dir, @remove_list );

# run self nucmers on the reference genomes
run_self_nucmer(@query);

# run nucmer, align complete genomes against reference
if ( $buildSNPdb == 0 ) {
    run_ref_nucmer( @query, $ref_genome );
}
elsif ( $buildSNPdb == 1 ) {
    run_all_nucmer(@query);
}
cleanup();

################################################################################

sub run_self_nucmer {

    # function that runs self nucmer on complete genomes
    my $ident      = 95;
    my $len_cutoff = 75;

    my $pm = new Parallel::ForkManager($thread);
    $pm->run_on_finish(
        sub {
            my ( $pid, $exit_code, $ident, $exit_signal, $core_dump ) = @_;
        }
    );

    foreach my $reference (@query) {

        $pm->start and next;
        my ( $name, $path, $suffix ) = fileparse( "$reference", qr/\.[^.]*/ );

        my $coords =
          File::Spec->catpath( $outdir, "stats", $name . '_repeat_coords.txt' );
        my $stat =
          File::Spec->catpath( $outdir, "stats", $name . '_repeat_stats.txt' );
        my $fasta =
          File::Spec->catpath( $outdir, "stats", $name . '_norepeats.fna' );

        if ( -e $coords && -e $stat && -e $fasta ) {
            print "Self-NUCmer already complete for $name.\n";
        }
        else {

            my $command =
"get_repeat_coords.pl -l $len_cutoff -i $repeat_identity -o $coords -s $stat $reference";
            print "[RUNNING:]$command\n";
            if ( system($command) ) { die "Error running $command.\n"; }

            my $remove_repeats =
              "removeRepeats.pl -f $reference -c $coords -o $fasta";
            print "[RUNNING:]$remove_repeats\n";
            if ( system($remove_repeats) ) {
                die "Error running $remove_repeats.\n";
            }
            $pm->finish(0);
            $pm->wait_all_children;
        }
    }
    print "\nRepeat coords for all references found.\n";
    print "Self-NUCmer complete\n";
    print
      "==================================================================\n";
}

################################################################################
sub run_all_nucmer {

    # A function that run all vs. all nucmer
    my $iteration = PhaME::combo( 2, @query );
    my $pm        = new Parallel::ForkManager($thread);

    # creating threads?
    $pm->run_on_finish(
        sub {
            my ( $pid, $exit_code, $ident, $exit_signal, $core_dump ) = @_;
        }
    );

    while ( my @combo = $iteration->() ) {
        $pm->start(@combo) and next;
        my %hash;
        my ( $first_name, $first_path, $first_suffix ) =
          fileparse( "$combo[0]", qr/\.[^.]*/ );
        my ( $second_name, $second_path, $second_suffix ) =
          fileparse( "$combo[1]", qr/\.[^.]*/ );
        my $reference = "$first_path$first_name$first_suffix";
        my $query     = "$second_path$second_name$second_suffix";
        $first_name  =~ s/\.fna//;
        $second_name =~ s/\.fna//;
        my $prefix1 = $first_name . '_' . $second_name;
        my $prefix2 = $second_name . '_' . $first_name;

        my $first_fasta =
          File::Spec->catpath( $outdir, "stats",
            $first_name . '_norepeats.fna' );
        my $second_fasta =
          File::Spec->catpath( $outdir, "stats",
            $second_name . '_norepeats.fna' );

################################################################################
        # run nucmer command, first vs. second genome
        nucmer::run_nucmer(
            {
                outdir            => $outdir,
                prefix            => $prefix1,
                ref_fasta         => $first_fasta,
                full_genome_fasta => $second_fasta,
                options           => $options
            }
        );

        # run nucmer command, second vs. first genome
        nucmer::run_nucmer(
            {
                outdir            => $outdir,
                prefix            => $prefix2,
                ref_fasta         => $second_fasta,
                full_genome_fasta => $first_fasta,
                options           => $options
            }
        );

        nucmer::run_delta_filter_snp(
            {
                outdir   => $outdir,
                identity => $identity,
                prefix   => $prefix1
            }
        );

        nucmer::run_delta_filter_snp(
            {
                outdir   => $outdir,
                identity => $identity,
                prefix   => $prefix2
            }
        );

        nucmer::run_show_snps(
            {
                outdir => $outdir,
                prefix => $prefix1
            }
        );
        nucmer::run_show_snps(
            {
                outdir => $outdir,
                prefix => $prefix2
            }
        );
        $snp_indel = `SNP_INDEL_count.pl $outdir/'nucmer'/$prefix1.snps`;
        $snp_indel =~ s/\n//;
        ( $snp_n, $indel_n ) = split /\t/, $snp_indel;
        $snp_indel = `SNP_INDEL_count.pl $outdir/'nucmer'/$prefix2.snps`;
        $snp_indel =~ s/\n//;
        ( $snp_n, $indel_n ) = split /\t/, $snp_indel;

        # }
################################################################################

        nucmer::run_delta_filter_gap(
            {
                outdir   => $outdir,
                identity => $identity,
                prefix   => $prefix1
            }
        );

        nucmer::run_delta_filter_gap(
            {
                outdir   => $outdir,
                identity => $identity,
                prefix   => $prefix2
            }
        );

################################################################################
        nucmer::run_show_coords(
            {
                outdir => $outdir,
                prefix => $prefix1
            }
        );
        nucmer::run_show_coords(
            {
                outdir => $outdir,
                prefix => $prefix2
            }
        );

        my $gaps1 =
`parseGapsNUCmer.pl $gap_cutoff $outdir/"nucmer"/$prefix1.coords 2>/dev/null`;
        ( $ref_gaps, $query_gaps, undef ) = split /\n/, $gaps1;
        my $gaps2 =
`parseGapsNUCmer.pl $gap_cutoff $outdir/"nucmer"/$prefix2.coords 2>/dev/null`;
        ( $ref_gaps, $query_gaps, undef ) = split /\n/, $gaps2;
################################################################################
        my $check1 =
`checkNUCmer.pl -i $outdir/'nucmer'/$first_name\_$second_name.gaps -r $reference`;

        if ( $check1 eq 1 ) {
            print "$second_name aligned < 25% of the $first_name genome\n";
        }
        my $check2 =
`checkNUCmer.pl -i $outdir/'nucmer'/$second_name\_$first_name.gaps -r $reference`;
        if ( $check1 eq 1 ) {
            print "$second_name aligned < 25% of the $first_name genome\n";
        }
        $pm->finish(0);
    }

    $pm->wait_all_children;

    print "\nNUCmer on all reference genomes complete.\n";
    print
      "==================================================================\n";
}

################################################################################
sub run_ref_nucmer {
    print "Aligning to the reference genome" . "$ref_genome\n\n";
    my $iteration = PhaME::combo( 2, @query );
    my $pm        = new Parallel::ForkManager($thread);

    $pm->run_on_finish(
        sub {
            my ( $pid, $exit_code, $ident, $exit_signal, $core_dump ) = @_;
        }
    );

    foreach my $full_genome (@query) {
        $pm->start and next;
        my ( $full_genome_name, $full_genome_path, $full_genome_suffix ) =
          fileparse( "$full_genome", qr/\.[^.]*/ );
        my ( $ref_name, $ref_path, $ref_suffix ) =
          fileparse( $ref_genome, qr/\.[^.]*/ );
        $ref_name         =~ s/\.fna//;
        $full_genome_name =~ s/\.fna//;
        my $prefix1 = $ref_name . '_' . $full_genome_name;
        my $ref_fasta =
          File::Spec->catpath( $outdir, "stats", $ref_name . '_norepeats.fna' );

        my $full_genome_fasta =
          File::Spec->catpath( $outdir, "stats",
            $full_genome_name . '_norepeats.fna' );

        nucmer::run_nucmer(
            {
                outdir            => $outdir,
                prefix            => $prefix1,
                ref_fasta         => $ref_fasta,
                full_genome_fasta => $full_genome_fasta,
                options           => $options
            }
        );

        nucmer::run_delta_filter_snp(
            {
                outdir   => $outdir,
                identity => $identity,
                prefix   => $prefix1
            }
        );
        nucmer::run_show_snps(
            {
                outdir => $outdir,
                prefix => $prefix1
            }
        );

        $snp_indel = `SNP_INDEL_count.pl $outdir/'nucmer'/$prefix1.snps`;
        $snp_indel =~ s/\n//;
        ( $snp_n, $indel_n ) = split /\t/, $snp_indel;

        nucmer::run_delta_filter_gap(
            {
                outdir   => $outdir,
                identity => $identity,
                prefix   => $prefix1
            }
        );

        nucmer::run_show_coords(
            {
                outdir => $outdir,
                prefix => $prefix1
            }
        );

        my $gaps1 =
`parseGapsNUCmer.pl $gap_cutoff $outdir/'nucmer'/$prefix1.coords 2>/dev/null`;
        ( $ref_gaps, $query_gaps, undef ) = split /\n/, $gaps1;

        $pm->finish(0);
    }
    $pm->wait_all_children;
    print "\nNUCmer on all reference genomes complete.\n";
    print
      "==================================================================\n";
}

sub cleanup {

    if ( -e "$outdir/*.mgaps" ) { unlink "$outdir/*.mgaps" }
    if ( -e "$outdir/*.ntref" ) { unlink "$outdir/*.ntref" }
    `cat $outdir/'stats'/*_repeat_stats.txt > $outdir/'stats'/repeat_stats.txt`;

    # `mv $outdir/*.snps $outdir/snps`;
    # `mv $outdir/*_stats.txt $outdir/stats`;
    # `mv $outdir/*_coords.txt $outdir/stats`;
    # `mv $outdir/*.coords $outdir/stats`;
    `mv $outdir/*norepeats.fna $outdir/temp`;
}

sub usage {
    print <<Usage;
$0 -r <reference genome> -q <query_dir> -d <output_directory> -t <# threads> -l <file containing list of all genomes to run nucmer on> -c <code>

Usage
    exit;
}

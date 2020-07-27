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
use PhaME;
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
if ( !-e "gaps" )  { mkdir "gaps"; }
if ( !-e "snps" )  { mkdir "snps"; }
if ( !-e "stats" ) { mkdir "stats"; }
if ( !-e "temp" )  { mkdir "temp"; }

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

# my @remove_list = PhaME::filter_genomes( $sketch_output, $ref_genome, $cutoff,
# $query_dir );

# print "$ref_genome\n";
# print "$cutoff\n";
# print "$query_dir\n";
# print "Warning: These genomes will be removed as they do not pass given ANI cutoff of $cutoff% with reference:$ref_genome\n";
# print join("\n", @remove_list);

@query = misc_funcs::read_directory( $query_dir, @remove_list );

# my $querySize = @query;
# if ($querySize < 3) {
#     die "There are less than 3 genomes, please add more genomes or relax the cutoff!\n"
# }

# run self nucmers on the reference genomes
run_self_nucmer(@query);

# run nucmer, align complete genomes against reference
if ( $buildSNPdb == 0 ) {
    run_ref_nucmer( @query, $ref_genome );
}
elsif ( $buildSNPdb == 1 ) {
    run_nucmer(@query);
}
cleanup();

################################################################################

sub run_self_nucmer {
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

        if ( -e $coords && -e $stat ) {
            print "\nRepeat coords already calculated.\n";
        }
        else {

            my $command =
"get_repeat_coords.pl -l $len_cutoff -i $repeat_identity -o $coords -s $stat $reference";
            if ( system($command) ) { die "Error running $command.\n"; }

            my $remove_repeats =
              "removeRepeats.pl -f $reference -c $coords -o $fasta";
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
sub run_nucmer {
    my $iteration = PhaME::combo( 2, @query );
    my $pm        = new Parallel::ForkManager($thread);

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

        my $first_fasta  = $outdir . '/' . $first_name . '_norepeats.fna';
        my $second_fasta = $outdir . '/' . $second_name . '_norepeats.fna';

        #need to write a split function here to just get the equivalent of
        # if ( $first_name !~ /$ref_genome/ && $ref_genome ) {
        print "Running nucmer on $prefix1\n";
        my $nucmer_command1 =
          "nucmer $options -p $prefix1 $first_fasta $second_fasta  2>/dev/null";
        print "[RUNNING:] $nucmer_command1\n";
        if ( system($nucmer_command1) ) {
            die "Error running nucmer_command1 $nucmer_command1.\n";
        }

        my $filter_command1 =
"delta-filter -1 $identity $outdir/$prefix1.delta > $outdir/$prefix1.snpfilter";
        if ( system($filter_command1) ) {
            die "Error running filter_command1 $filter_command1.\n";
        }

        my $snp_command1 =
          "show-snps -CT $outdir/$prefix1.snpfilter > $outdir/$prefix1.snps";
        if ( system($snp_command1) ) {
            die "Error running snp_command1 $snp_command1.\n";
        }
        $snp_indel = `SNP_INDEL_count.pl $outdir/$prefix1.snps`;
        $snp_indel =~ s/\n//;
        ( $snp_n, $indel_n ) = split /\t/, $snp_indel;

        my $filter_command2 =
"delta-filter -1 $identity $outdir/$prefix1.delta > $outdir/$prefix1.gapfilter";
        if ( system($filter_command2) ) {
            die "Error running filter_command2 $filter_command2.\n";
        }

        my $coords_command1 =
"show-coords -clTr $outdir/$prefix1.gapfilter > $outdir/$prefix1.coords";
        if ( system($coords_command1) ) {
            die "Error running coords_command1 $coords_command1.\n";
        }
        my $gaps1 =
          `parseGapsNUCmer.pl $gap_cutoff $outdir/$prefix1.coords 2>/dev/null`;
        ( $ref_gaps, $query_gaps, undef ) = split /\n/, $gaps1;

        my $check =
`checkNUCmer.pl -i $outdir/$first_name\_$second_name.gaps -r $reference`;
        if ( $check == 1 ) {
            print "$second_name aligned < 25% of the $first_name genome\n";
        }

        # }

        # if ( $second_name !~ /$ref_genome/ && $ref_genome ) {
        print "Running nucmer on $prefix2\n";
        my $nucmer_command2 =
          "nucmer $options -p $prefix2 $second_fasta $first_fasta  2>/dev/null";
        if ( system($nucmer_command2) ) {
            die "Error running nucmer_command2 $nucmer_command2.\n";
        }

        my $filter_command3 =
"delta-filter -1 $identity $outdir/$prefix2.delta > $outdir/$prefix2.snpfilter";
        if ( system($filter_command3) ) {
            die "Error running filter_command3 $filter_command3.\n";
        }

        my $snp_command2 =
          "show-snps -CT $outdir/$prefix2.snpfilter > $outdir/$prefix2.snps";
        if ( system($snp_command2) ) {
            die "Error running snp_command2 $snp_command2.\n";
        }
        $snp_indel = `SNP_INDEL_count.pl $outdir/$prefix2.snps`;
        $snp_indel =~ s/\n//;
        ( $snp_n, $indel_n ) = split /\t/, $snp_indel;

        my $filter_command4 =
"delta-filter $identity $outdir/$prefix2.delta > $outdir/$prefix2.gapfilter";
        if ( system($filter_command4) ) {
            die "Error running filter_command4 $filter_command4.\n";
        }

        my $coords_command2 =
"show-coords -clTr $outdir/$prefix2.gapfilter > $outdir/$prefix2.coords";
        if ( system($coords_command2) ) {
            die "Error running coords_command2 $coords_command2.\n";
        }

        my $gaps2 =
          `parseGapsNUCmer.pl $gap_cutoff $outdir/$prefix2.coords 2>/dev/null`;
        ( $ref_gaps, $query_gaps, undef ) = split /\n/, $gaps2;
        $check =
          `checkNUCmer.pl -i $outdir/$second_name\_$first_name.gaps -r $query`;
        if ( $check == 1 ) {
            print "$first_name aligned < 25% of the $second_name genome\n";
        }

        # }

        $pm->finish(0);
    }

    $pm->wait_all_children;

    print "\nNUCmer on all reference genomes complete.\n";
    print
      "==================================================================\n";
}

################################################################################
sub run_ref_nucmer {
    print "Aligning to the reference genome " . "$ref_genome\n\n";
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
        my $prefix1   = $ref_name . '_' . $full_genome_name;
        my $ref_fasta = $outdir . '/' . $ref_name . '_norepeats.fna';
        my $full_genome_fasta =
          $outdir . '/' . $full_genome_name . '_norepeats.fna';

        print "Running nucmer on $prefix1\n";
        my $nucmer_command1 =
"nucmer $options -p $prefix1 $ref_fasta $full_genome_fasta  2>/dev/null";
        if ( system($nucmer_command1) ) {
            die "Error running nucmer_command1 $nucmer_command1.\n";
        }

        my $filter_command1 =
"delta-filter -1 $identity $outdir/$prefix1.delta > $outdir/$prefix1.snpfilter";
        if ( system($filter_command1) ) {
            die "Error running filter_command1 $filter_command1.\n";
        }

        my $snp_command1 =
          "show-snps -CT $outdir/$prefix1.snpfilter > $outdir/$prefix1.snps";
        if ( system($snp_command1) ) {
            die "Error running snp_command1 $snp_command1.\n";
        }
        $snp_indel = `SNP_INDEL_count.pl $outdir/$prefix1.snps`;
        $snp_indel =~ s/\n//;
        ( $snp_n, $indel_n ) = split /\t/, $snp_indel;

        my $filter_command2 =
"delta-filter -1 $identity $outdir/$prefix1.delta > $outdir/$prefix1.gapfilter";
        if ( system($filter_command2) ) {
            die "Error running filter_command2 $filter_command2.\n";
        }

        my $coords_command1 =
"show-coords -clTr $outdir/$prefix1.gapfilter > $outdir/$prefix1.coords";
        if ( system($coords_command1) ) {
            die "Error running coords_command1 $coords_command1.\n";
        }
        my $gaps1 =
          `parseGapsNUCmer.pl $gap_cutoff $outdir/$prefix1.coords 2>/dev/null`;
        ( $ref_gaps, $query_gaps, undef ) = split /\n/, $gaps1;

# my $aligned_percentage
#     = `checkNUCmer.pl -i $outdir/$ref_name\_$full_genome_name.gaps -r $ref_fasta`;

# if ( $aligned_percentage < 25 ) {
#     print
#         "$full_genome_name aligned < 25% ($aligned_percentage) of the $ref_name genome\n";
# }
# elsif ( $aligned_percentage >= 25 ) {
#     print
#         "$full_genome_name aligned > 25% ($aligned_percentage) of the $ref_name genome\n";
# }

        $pm->finish(0);
    }
    $pm->wait_all_children;
    print "\nNUCmer on all reference genomes complete.\n";
    print
      "==================================================================\n";
}

sub combo {
    my $by = shift;
    return sub { () }
      if !$by || $by =~ /\D/ || @_ < $by;
    my @list = @_;

    my @position = ( 0 .. $by - 2, $by - 2 );
    my @stop     = @list - $by .. $#list;
    my $end_pos  = $#position;
    my $done     = undef;

    return sub {
        return () if $done;
        my $cur = $end_pos;
        {
            if ( ++$position[$cur] > $stop[$cur] ) {
                $position[ --$cur ]++;
                redo if $position[$cur] > $stop[$cur];
                my $new_pos = $position[$cur];
                @position[ $cur .. $end_pos ] = $new_pos .. $new_pos + $by;
            }
        }
        $done = 1 if $position[0] == $stop[0];
        return @list[@position];
    }
}

sub cleanup {

    #if (-e "$outdir/*.coords"){`rm $outdir/*.coords`};
    if ( -e "$outdir/*.mgaps" ) { unlink "$outdir/*.mgaps" }
    if ( -e "$outdir/*.ntref" ) { unlink "$outdir/*.ntref" }
    `cat $outdir/*_repeat_stats.txt > $outdir/repeat_stats.txt`;
    `mv $outdir/*.snps $outdir/snps`;
    `mv $outdir/*.gaps $outdir/gaps`;
    `mv $outdir/*_stats.txt $outdir/stats`;
    `mv $outdir/*_coords.txt $outdir/stats`;
    `mv $outdir/*.coords $outdir/stats`;
    `mv $outdir/*norepeats.fna $outdir/temp`;
}

sub usage {
    print <<Usage;
$0 -r <reference genome> -q <query_dir> -d <output_directory> -t <# threads> -l <file containing list of all genomes to run nucmer on> -c <code>

Usage
    exit;
}

#!/usr/bin/env perl

######################################################
# Written by Sanaa Ahmed
# Nov. 30, 2012
# Modified by Migun Shakya
# 032718 passed thread variables to bowtie
# Given a directory containing paired read files, runs bowtie
######################################################

use strict;
use FindBin qw($RealBin);
use Getopt::Long;
use File::Basename;
use Parallel::ForkManager;

# set up environments
$ENV{PATH} = "$RealBin:$RealBin/../ext/bin:$ENV{PATH}";

###############Set up variables#################################################
my (
    $indir, $reference, $prefix, $thread,
    $list,  $aligner,   $ploidy, $snp_filter
);
my @command;
my $outdir = `pwd`;
$outdir =~ s/\n//;
###############Set up variables#################################################

GetOptions(
    'q|querydir=s'   => \$indir,
    'r|reference=s'  => \$reference,
    'd|outdir=s'     => \$outdir,
    't|thread=i'     => \$thread,
    'l|list=s'       => \$list,
    'a|aligner=s'    => \$aligner,
    'p|ploidy=s'     => \$ploidy,
    's|snp_filter=f' => \$snp_filter,
    'h|help'         => sub { usage() }
);

&usage() unless ( $indir && $aligner );

####Checking the number of allowed threads######################################
my $max_thread =
  ( $^O =~ /darwin/ )
  ? `sysctl hw.ncpu | awk '{print \$2}'`
  : `grep -c ^processor /proc/cpuinfo`;
if ( $thread < 1 || $thread > $max_thread ) {
    die("-thread value must be between than 1 and $max_thread.\n");
}
################################################################################

####################Setting up directories######################################
if ( !-e $outdir )        { mkdir "$outdir"; }
if ( $outdir =~ /.+\/$/ ) { my $temp = chop($outdir); }
chdir $outdir;

if ( $indir =~ /.+\/$/ ) { my $temp = chop($indir); }
my $querydir = $indir;

################################################################################

read_directory($querydir);

####################Setting up threads##########################################
my $pm = new Parallel::ForkManager(1);

$pm->run_on_finish(    # called BEFORE the first call to start()
    sub {
        my ( $pid, $ident ) = @_;
        opendir( TEMP, "$outdir" );
        my $lines = 0;
        while ( my $file = readdir(TEMP) ) {
            if ( $file =~ /.+\.gaps$/ ) {
                my $gap_file = $outdir . '/' . $file;
                open( FILE, "$gap_file" );
                while (<FILE>) { $lines++; }
                if ( $lines == 1 ) {

                    # `rm -f $gap_file`;
                    print "Removed $gap_file\n";
                    $lines = 0;
                }

                # TODO: will create a folder to keep all outputs from mapping
                else { `mv $gap_file $outdir/gaps`; }
            }
        }
    }
);

################################################################################
foreach my $command (@command) {
    $pm->start and next;
    if ( system($command) ) { die "Error running $command.\n"; }
    $pm->finish;
}
$pm->wait_all_children;

print "Read Mapping complete\n";

################################################################################

sub read_directory {
    my $dir = shift;
    my ( $name, $path, $suffix, $file, $query );
    my $ref = 0;
    my @fastq;
    my %queries;
    if ($reference) {
        ( $name, $path, $suffix ) = fileparse( "$reference", qr/\.[^.]*/ );
    }

    open( IN, $list );
    while (<IN>) {
        chomp;
        $queries{$_}++;
    }
    close IN;

    my $temp = 0;
    my %check;
    opendir( PARENT, $dir );
    while ( my $files = readdir(PARENT) ) {
        next if ( $files =~ /^..?$/ );
        if ( $files !~ /$name/
            && ( $files =~ /(.+)[_.]R1.*\.fa?s?t?q$/ ) )
        {
            $temp = $1 . '_pread';
            if ( exists $queries{$temp} && !exists $check{$temp} ) {
                $check{$temp}++;
                $query  = $dir . '/' . $files;
                $prefix = "$name";
                mapping_function( $query, $prefix, $temp, $thread,
                    $snp_filter );
            }
            $temp = $1 . '_read';
            if ( exists $queries{$temp} && !exists $check{$temp} ) {
                $check{$temp}++;
                $query  = $dir . '/' . $files;
                $prefix = "$name";
                mapping_function( $query, $prefix, $temp, $thread,
                    $snp_filter );
            }
            $temp = $1 . '_sread';
            if ( exists $queries{$temp} && !exists $check{$temp} ) {
                $check{$temp}++;
                $query  = $dir . '/' . $files;
                $prefix = "$name";
                mapping_function( $query, $prefix, $temp, $thread,
                    $snp_filter );
            }
        }
        if ( $files !~ /$name/ && $files =~ /(.+)\.fa?s?t?q$/ ) {
            $temp = $1 . '_sread';
            if ( exists $queries{$temp} && !exists $check{$temp} ) {
                $check{$temp}++;
                $query  = $dir . '/' . $files;
                $prefix = "$name";
                mapping_function( $query, $prefix, $temp, $thread,
                    $snp_filter );
            }
        }
    }
    closedir(PARENT);
}

################################################################################
sub mapping_function {
    my $read       = shift;
    my $prefix     = shift;
    my $temp       = shift;
    my $thread     = shift;
    my $snp_filter = shift;
    my $read1;
    my $read2;
    my $readu;
    my $map_cmd;
    my ( $name, $path, $suffix ) = fileparse( "$read", qr/\.[^.]*/ );
    my $bowtie_options = '-p $thread';

    if ( $temp =~ /pread/i ) {

        if ( $name =~ /(.+)([_.]R)(\d)(.*)$/ ) {
            $prefix .= "\_$1";
            $read1 = $path . $1 . $2 . $3 . $4 . $suffix;
            $read2 = $path . $1 . $2 . '2' . $4 . $suffix;
        }

        my $map_plot = File::Spec->catpath( $outdir, $prefix . '_plots.pdf' );
        if ( -e $map_plot ) {
            print "Paired reads, $read1 and $read2, are already mapped";
            $map_cmd =
"runReadsToGenome.pl -snp_filter $snp_filter -ploidy $ploidy -p '$read1 $read2' -ref $reference -pre $prefix -d $outdir -aligner $aligner -cpu $thread -consensus 0 -skip_aln";
        }
        else {
            $map_cmd =
"runReadsToGenome.pl -snp_filter $snp_filter -ploidy $ploidy -p '$read1 $read2' -ref $reference -pre $prefix -d $outdir -aligner $aligner -cpu $thread -consensus 0";
        }
        print "\n[RUNNING:] $map_cmd\n";
        push( @command, $map_cmd );

    }
    elsif ( $temp =~ /sread/i ) {
        $prefix .= "\_$name";
        my $map_plot = File::Spec->catpath( $outdir, $prefix . '_plots.pdf' );

        if ( -e $map_plot ) {
            print "Single read $read is already mapped";
            $map_cmd =
"runReadsToGenome.pl -snp_filter $snp_filter -ploidy $ploidy -u $read -ref $reference -pre $prefix -d $outdir -aligner $aligner -cpu $thread -consensus 0 -skip_aln";
        }
        else {
            $map_cmd =
"runReadsToGenome.pl -snp_filter $snp_filter -ploidy $ploidy -u $read -ref $reference -pre $prefix -d $outdir -aligner $aligner -cpu $thread -consensus 0";
        }
        print "\n[RUNNING:] $map_cmd\n";
        push( @command, $map_cmd );
    }

    elsif ( $temp =~ /_read/i ) {

        if ( $name =~ /(.+)([_.]R)(\d)(.*)$/ ) {
            $prefix .= "\_$1";
            $read1 = $path . $1 . $2 . $3 . $4 . $suffix;
            $read2 = $path . $1 . $2 . '2' . $4 . $suffix;
            $readu = $path . $1 . $suffix;
        }
        my $map_plot = File::Spec->catpath( "maps", $prefix . '.pdf' );
        if ( -e $map_plot ) {
            print "Mapping step already finished";
            my $map_cmd =
"runReadsToGenome.pl -snp_filter $snp_filter -ploidy $ploidy -p $read1,$read2 -u $readu -ref $reference -pre $prefix -d $outdir -aligner $aligner -cpu $thread `-bowtie_options '-p $thread'";
        }
        else {

            my $map_cmd =
"runReadsToGenome.pl -snp_filter $snp_filter -ploidy $ploidy -p $read1,$read2 -u $readu -ref $reference -pre $prefix -d $outdir -aligner $aligner -cpu $thread `-bowtie_options '-p $thread' -skip_aln";
        }

        print "\n[RUNNING:] $map_cmd\n";

        push( @command, $map_cmd );
    }
}

sub usage {
    print <<"END";
    Usage:perl $0 
        -r          <reference_fasta> 
        -q          <query_dir> , contains reads files in the format *R[1-2].fastq, if paired, elase *.fastq
        -d          <output_directory>, if not provided will write files in current directory
        -t          <# of threads> 
        -l          <list file>
        -a          <aligner bwa|bowtie>

Synopsis:
    perl $0 -r <reference_fasta> -q <query_dir> -d <output_directory> -t <# of threads> -a <aligner bwa|bowtie>
END

    exit;
}


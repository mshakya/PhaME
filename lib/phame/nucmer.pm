#!/usr/bin/env perl

use strict;
use warnings;
use FileHandle;
use File::Basename;
use Cwd;

package nucmer;

sub run_nucmer {

    # Runs nucmer, returns delta file, moves it to nucmer folder
    my ($args) = @_;

    my $prefix            = $args->{prefix};
    my $ref_fasta         = $args->{ref_fasta};
    my $full_genome_fasta = $args->{full_genome_fasta};
    my $opt               = $args->{options};

    if ( !-e "nucmer" ) {
        mkdir "nucmer";
    }
    my $delta_file =
      File::Spec->catpath( $args->{outdir}, "nucmer",
        $args->{prefix} . '.delta' );
    if ( -e $delta_file ) {
        print
"nucmer complete for $args->{ref_fasta} and $args->{full_genome_fasta}.\n";
    }
    else {
        print "Running nucmer on $prefix\n";
        my $nc1 =
          "nucmer $opt -p $prefix $ref_fasta $full_genome_fasta  2>/dev/null";
        print "[RUNNING:] $nc1\n";
        if ( system($nc1) ) {
            die "Error running nc1 $nc1.\n";
        }
        `mv $args->{prefix}.delta nucmer/`;
    }
}

sub run_delta_filter_snp {

    # Runs delta filter, and writes the output in .snpfilter
    my ($args) = @_;

    my $outdir   = $args->{outdir};
    my $identity = $args->{identity};
    my $prefix   = $args->{prefix};

    my $snpfilter_file =
      File::Spec->catpath( $outdir, "nucmer", $prefix . '.snpfilter' );
    my $delta_file =
      File::Spec->catpath( $outdir, "nucmer", $prefix . '.delta' );
    if ( -e $snpfilter_file ) {
        print "nucmer already complete for $prefix.\n";
    }
    else {
        my $filter_command1 =
          "delta-filter -1 $identity $delta_file > $snpfilter_file";

        if ( system($filter_command1) ) {
            die "Error running filter_command1 $filter_command1.\n";
        }
    }
}

sub run_delta_filter_gap {

    # Runs delta filter, and writes the output in .gapfilter
    my ($args) = @_;

    my $out_dir  = $args->{outdir};
    my $identity = $args->{identity};
    my $prefix   = $args->{prefix};

    my $gapfilter_file =
      File::Spec->catpath( $out_dir, "nucmer", $prefix . '.gapfilter' );
    my $delta_file =
      File::Spec->catpath( $out_dir, "nucmer", $prefix . '.delta' );
    if ( -e $gapfilter_file ) {
        print "nucmer already complete for $prefix.\n";
    }
    else {
        my $filter_command1 =
          "delta-filter $identity $delta_file > $gapfilter_file";

        if ( system($filter_command1) ) {
            die "Error running filter_command1 $filter_command1.\n";
        }
    }
}

sub run_show_snps {
    my ($args) = @_;

    my $outdir = $args->{outdir};
    my $prefix = $args->{prefix};

    my $snpfilter_file =
      File::Spec->catpath( $outdir, "nucmer", $prefix . '.snpfilter' );
    my $snp_file = File::Spec->catpath( $outdir, "nucmer", $prefix . '.snps' );
    if ( -e $snp_file ) {
        print "show-snps already complete for $prefix.\n";
    }
    else {
        my $snp_command1 = "show-snps -CT $snpfilter_file > $snp_file";
        if ( system($snp_command1) ) {
            die "Error running snp_command1 $snp_command1.\n";
        }
    }
}

sub run_show_coords {
    my ($args) = @_;
    my $outdir = $args->{outdir};
    my $prefix = $args->{prefix};

    my $gfilt_file =
      File::Spec->catpath( $outdir, "nucmer", $prefix . '.gapfilter' );

    my $coord_file =
      File::Spec->catpath( $outdir, "nucmer", $prefix . '.coords' );
    if ( -e $coord_file ) {
        print "show-coords already complete for $prefix.\n";
    }
    else {
        my $coords_command1 = "show-coords -clTr $gfilt_file > $coord_file";
        if ( system($coords_command1) ) {
            die "Error running coords_command1 $coords_command1.\n";
        }
    }
}

1;

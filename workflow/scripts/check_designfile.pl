#!/usr/bin/perl -w
#check_designfile.pl

use strict;
use warnings;

my $pe = shift @ARGV;
my $dfile = shift @ARGV;
open OUT, ">design.valid.txt" or die $!;
open DFILE, "<$dfile" or die $!;
my $head = <DFILE>;
chomp($head);
$head =~ s/FullPathTo//g;
my @colnames = split(/\t/,$head);
my %newcols = map {$_=> 1} @colnames;

unless (grep(/FqR1/,@colnames)) {
    die "Missing Sequence Files in Design File: FqR1\n";
}
unless (grep(/SampleID/,@colnames)) {
    die "Missing SampleID in Design File\n";
}

if ($pe eq 'pe') {
    unless (grep(/FqR2/,@colnames)) {
	die "Missing Sequence Files in Design File: FqR2\n";
    }
}else {
    delete $newcols{FqR2};
}

my @cols = sort {$a cmp $b} keys %newcols;
print OUT join("\t",@cols),"\n";
my @grp = ('a','b');
my $lnct = 0;
while (my $line = <DFILE>) {
    chomp($line);
    $line =~ s/FullPathTo//g;
    my @row = split(/\t/,$line);
    my %hash;
    foreach my $i (0..$#row) {
	next unless ($newcols{$colnames[$i]});
	$row[$i] =~ s/-//g unless ($colnames[$i] =~ m/Fq/);
	$hash{$colnames[$i]} = $row[$i];
    }
    if ($hash{SampleID} =~ m/^\d/) {
	$hash{SampleID} =~ s/^/S/;
    }
    $hash{SampleName} = $hash{SampleID} unless ($hash{SampleName});
    $hash{SubjectID} = $hash{SampleID} unless ($hash{SubjectID});
    $hash{FinalBam} = $hash{SampleID}.".final.bam";
    $hash{Bam} = $hash{SampleID}.".bam";
    unless ($hash{SampleGroup}) {
	my $j = $lnct %% 2;
	$hash{SampleGroup} = $grp[$j];
    }
    $hash{SampleGroup} =~ s/_//g;
    unless ($hash{FqR1} =~ m/.fastq.gz/) {
        my $name = $hash{FqR1};
        $name =~ s/.f.*/.fastq.gz/;
        unless ($hash{FqR1} eq $name) {
	    $hash{FqR1} = $name;
            unless (-e ($name)) {
                    die "Unable to find fastq read 1\n${name}\n";
            }
	}
    }
    $hash{FqR2} = 'na' unless ($hash{FqR2});
    unless ($hash{FqR2} eq 'na') {
        my $name = $hash{FqR2};
        $name =~ s/.f.*/.fastq.gz/;
        unless ($hash{FqR2} eq $name) {
	    $hash{FqR2} = $name;
            unless (-e ($name)) {
                    die "Unable to find fastq read 2\n${name}\n";
            }
	}
    }
    my @line;
    foreach my $f (@cols) {
	push @line, $hash{$f};
    }
    print OUT join("\t",@line),"\n";
    print join(",",$hash{SampleID},$hash{FqR1},$hash{FqR2}),"\n";
    $lnct ++;
}

#!/usr/bin/env perl
use strict;
use warnings;

#This perl script will read in the multi-fasta file to verify validity of the ORF.
#It will detect in-frame stop codons, DNA "ambiguity" characters other than N, lack of stop codon at the end of the sequence and lack of start codon at the beginning of start codon.

my $orffile = $ARGV[0];
my $counter=0;
my $orfcount=0;
my $header='';
my $name;
my $seq='';
my $len;
my $newseq;
my %sequence;
my $i;
my $start = $ARGV[1];

open(ORF, "$orffile") || die "File not found: '$orffile'.\n";
while (<ORF>) {
        chomp($_);
        if ($_=~/^>(.*)/) {
                if ($seq ne ""){$sequence{$header}=$seq;}
                $header = $1;$counter++;$seq = '';
        } else {$seq.=$_;}
}
$sequence{$header}=$seq;
foreach $name (sort {lc $a cmp lc $b} keys %sequence) {
        $len = length($sequence{$name});
                print "$name\n";
                $newseq=substr($sequence{$name}, 0, $start);
                print "$newseq\t";
                for ( $i=$start; $i<$len; $i=$i+3 ) {
                $newseq=substr($sequence{$name}, $i, 3);
                print "$newseq\t";
                }
        print "\n";
}

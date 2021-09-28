#!/usr/bin/perl/

open IN, "<read1_match_read2_alignment.txt" or die;
open OUT, ">alignment_visualization.txt" or die;

my $line, $id, $cigar, $mdtag, $temp_mut_loci, $genomic_align_loci, $read1_seq, $read1_mut_loci, $strand, $genomic_seq, $read1_mut, $genomic_mut, $mut_count, $mut_rate;
my $soft_clip, $read1_mut_loci_vis, $deletion_length, $insertion_length, $match_length, $number, $letter, $seq, $loci, $seq_del, $read1_seq_vis, $genomic_seq_vis, $cigar_loci;
my @eles =(), @match_length =(), @deletion_length=(), @insertion_length=(), @read1=(), @genomic=();

print OUT "<IN>";

while ($line=<IN>) {
	chomp($line);
	@eles = split ("\t", $line);
	$id = $eles[0];
	$cigar = $eles[1];
	$mdtag = substr($eles[2],5);
	$temp_mut_loci = $eles[3];
	$genomic_align_loci = $eles[4];
	$read1_seq = $eles[5];
	$read1_mut_loci = $eles[6];
	$strand = $eles[7];
	$genomic_seq = $eles[8];
	$read1_mut = $eles[9];
	$genomic_mut = $eles[10];
	$mut_count = $eles[11];
	$mut_rate = $eles[12];
	$read1_seq_vis = undef;
    $genomic_seq_vis = undef;
    $seq = undef;
	$soft_clip = 0;
	
	
	if ($cigar =~ /^(\d+)(S).*/) {
        $soft_clip = $1;
        $cigar =~ s/\d+S//s;
    }

 	
	$genomic_seq = substr($genomic_seq, $genomic_align_loci);
	$read1_seq = substr($read1_seq, $soft_clip);
	
	$loci = 0;
	$cigar_loci = 0;
	while ($cigar) {
	    $cigar =~ /^(\d+)([\D])(.*)/;
	    $number = $1;
	    $letter = $2;
	    if ($letter eq "M") {
	        $seq = substr($read1_seq, $loci, $number);
	    }
	    if ($letter eq "D") {
            $seq = "-" x $number;
	    }
	    if ($letter eq "I") {
	        $seq_del = substr($read1_seq, $loci, $number);
	    }
	    $read1_seq_vis = $read1_seq_vis . $seq;
	    $read1_seq = substr($read1_seq, $number);
	    $cigar_loci = length($number) + 1;
	    $cigar = substr($cigar, $cigar_loci);	
	}
	
	$cigar = $eles[1];
	if ($cigar =~ /^(\d+)(S).*/) {
        $soft_clip = $1;
        $cigar =~ s/\d+S//s;
    }
    
    $loci = 0;	
    $cigar_loci = 0;
	while ($cigar) {
	    $cigar =~ /^(\d+)([\D])(.*)/;
	    $number = $1;
	    $letter = $2;
	    if ($letter eq "M") {
	        $seq = substr($genomic_seq, $loci, $number);
	    }
	    if ($letter eq "I") {
            $seq = "-" x $number;
	    }
	    if ($letter eq "D") {
	        $seq_del = substr($genomic_seq, $loci, $number);
	    }
	    $genomic_seq_vis = $genomic_seq_vis . $seq;
	    $genomic_seq = substr($genomic_seq, $number);
	    $cigar_loci = length($number) + 1;
	    $cigar = substr($cigar, $cigar_loci);	
	}
	$genomic_seq_vis = substr($genomic_seq_vis,0,(length($read1_seq_vis)));
		
    print OUT "$line\n";
    print OUT "$read1_seq_vis\n";
    @read1 = split (//, $read1_seq_vis);
    @genomic = split (//, $genomic_seq_vis);
    foreach (0..(length($read1_seq_vis))) {
      if($read1[$_] and $genomic[$_]){
        if ($read1[$_] eq $genomic[$_]) {            
            print OUT "|";
        }
        else {
            if ($read1[$_] eq "A" and $genomic[$_] eq "G") {
                print OUT "*";
            }
            else {
                print OUT "-";
            }
        }
      }
    }
    print OUT "\n$genomic_seq_vis\n\n";

}

close IN, OUT;
        	
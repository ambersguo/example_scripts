#!/usr/bin/perl/

open IN, "<read1_match_read2_bwa_mdtag_filtered.sam" or die;
open OUT, ">read1_match_read2_bwa_mdtag_filtered_mdtag_parsing_2.txt" or die;

my $count =0, $count1 = 0, $count2 =0;
while ($line=<IN>) {
    chomp($line);
    @eles=split(/\t/, $line);
    $CIGAR = $eles[5];
    $soft_clip = 0;
    if ($eles[5] =~ /^(\d+)(S).*/) {
        $soft_clip = $1;
    }
    $mdtag = substr ($eles[12], 5);
    if ($mdtag =~ /(.*\D)(\d+)(G)(.*)/) {            
        while ($mdtag =~ m/(.*\D)(\d+)([G])(.*)/g) {           
            $head1 = $1;
            $head2 = $2;
            @deletion = ();
            $total_deletion_length = 0;
            if ($head1 =~ /(.*)\^(.*)/) {
                while ($head1 =~ m/([\^])(\D+)/g) {
                     $deletion = $2;
                     push (@deletion, $deletion);
                 }                 
                foreach (@deletion) {
                     $deletion_length = length($_);
                     $total_deletion_length = $total_deletion_length + $deletion_length;
                 } 

                $head1 =~ tr/^//d;
                @d_D_split = split(/(?<=\d)(?=\D)|(?<=\D)(?=\d)/, $head1);
                foreach $d_D_split (@d_D_split) {
                    if ($d_D_split =~ /\d+/) {              
                        $count = $count + $d_D_split;
                    }
                    else {
           
                        $count = $count + (length($d_D_split));
        
                    }               
                }       
            
                $count = $count - $total_deletion_length;    
      
            }
            
            else {          
                while ($head1 =~ m/(\d+)\D/g) {
                    $count = $count + $1 +1;
                } 
            }
            $count = $count + $head2 + $soft_clip;
            
            $count_constant = $count;
            while ($CIGAR =~ m/(.*\D)(\d+)([I])(.*)/g ) {
                $CIGAR_head = $1;
                $CIGAR_count2 = $2;
                @CIGAR_split = split(/(?<=\d)(?=\D)|(?<=\D)(?=\d)/, $CIGAR_head);
                foreach $CIGAR_split (@CIGAR_split) {
                    if ($CIGAR_split =~ /\d+/) {              
                    $CIGAR_count1 = $CIGAR_count1 + $CIGAR_split;
                    }    
                }
                if (($CIGAR_count1 < $count_constant) or ($CIGAR_count1 == $count_constant)) {

                    $count = $count + $CIGAR_count2;
                }
                $CIGAR_count1 = 0;
            }           
            push (@count_array, $count);
            $count = 0;                                     
        }
        
        if ($mdtag =~ /^(\d+)([G])(.*)/) {           
            $count = $1 + $soft_clip;
            
            $count_constant = $count;
            while ($CIGAR =~ m/(.*\D)(\d+)([I])(.*)/g ) {
                $CIGAR_head = $1;
                $CIGAR_count2 = $2;
                @CIGAR_split = split(/(?<=\d)(?=\D)|(?<=\D)(?=\d)/, $CIGAR_head);
                foreach $CIGAR_split (@CIGAR_split) {
                    if ($CIGAR_split =~ /\d+/) {              
                    $CIGAR_count1 = $CIGAR_count1 + $CIGAR_split;
                    }    
                }
                if (($CIGAR_count1 < $count_constant) or ($CIGAR_count1 == $count_constant)) {

                    $count = $count + $CIGAR_count2;
                }
                $CIGAR_count1 = 0;
            } 
            push (@count_array, $count);
            $count = 0;
        }      
                                      
        foreach (@count_array) {
            $newline = join("\t", $line, $_);
            push (@newline, $newline);
        }
        @count_array = ();           
    }
}

@sorted_line = sort {
@a_fields = split /\t/, $a;
@b_fields = split /\t/, $b;
$a_fields[0] cmp $b_fields[0]  
} @newline;

foreach $sorted_line (@sorted_line) {
	print OUT "$sorted_line\n";
}

close IN, OUT;

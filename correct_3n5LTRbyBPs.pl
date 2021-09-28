#!/usr/bin/perl

use File::Copy;

opendir DIR_bp, "/Volumes/Seagate/Data/GeneTherapyProject/PT2716/3n5LTRbyBPs" or die "Cannot open directory:$!";
@files = readdir DIR_bp;
@months = ();
foreach $file (@files) {
    if ($file =~ /^\d+m$/) {
        push (@months, $file);
    }
}

@cells = qw(
CD19
CD3
NK
PMN
CD14
CD-34
);

foreach $month (@months) {
    opendir DIR_bp, "/Volumes/Seagate/Data/GeneTherapyProject/PT2716/3n5LTRbyBPs/$month" or die "Cannot open directory:$!";
    opendir DIR_reads, "/Volumes/Seagate/Data/GeneTherapyProject/PT2716/3n5LTRbyreads/$month" or die "Cannot open directory:$!";
    
    @files_bp = readdir DIR_bp;
    @files_reads = readdir DIR_reads;
    
    foreach $cell (@cells) {
    
        $check0 = $check1 = 0;
        foreach $file_bp (@files_bp) {
            if ($file_bp =~ /$cell/) {
                open IN_bp, "</Volumes/Seagate/Data/GeneTherapyProject/PT2716/3n5LTRbyBPs/$month/$file_bp" or die;
                $file2out = $file_bp;
                $check0 = 1;
            }
        }
        foreach $file_reads (@files_reads) {
            if ($file_reads =~ /$cell/) {
                open IN_reads, "</Volumes/Seagate/Data/GeneTherapyProject/PT2716/3n5LTRbyreads/$month/$file_reads" or die;
                $check1 = 1;
            }
        }
        
        if($check0 and $check1) {
            @cors = @rests = ();
            $n = $m = $total_BP = $total_BP1 = $total_BP2 = $total_reads = $total_reads1 = $total_reads2 = $check2 = $check3 = 0;
            %cors_hash = %cor2line = %cor2read = %cor2percent = ();
        
            open OUT, ">/Volumes/Seagate/Data/GeneTherapyProject/PT2716/3n5LTRbyBPscorrected/$month/$file2out" or die;
            while ($line = <IN_bp>) {
                $n++;
                if ($n > 10) {
                    chomp($line);
                    @eles = split(/\t/,$line);
                    $cor = $eles[1].$eles[2].$eles[3];
                    $BPcount = $eles[4];
                    if ($BPcount > 20) {
                        push (@cors, $cor);
                        $total_BP1 = $total_BP1 + $BPcount;
                        $check2 = 1;
                    }
                    else {
                        push (@rests, $cor);
                        $total_BP2 = $total_BP2 + $BPcount;
                    }
                    $cor2line{$cor} = $line;
                }
                else {
                    print OUT "$line";
                }
            }
            %cors_hash = map { $_ => 1 } @cors;

            if ($check1 == 1 and $check2 == 1) {
                while ($line = <IN_reads>) {
                    $m++;
                    if ($m > 10) {
                        chomp($line);
                        @eles = split(/\t/,$line);
                        $cor = $eles[1].$eles[2].$eles[3];
                        $readcount = $eles[4];
                        $percent = $eles[5];
                        $total_reads = $total_reads + $readcount;
                        if(exists($cors_hash{$cor})) {
                            $cor2read{$cor} = $readcount;
                            $cor2percent{$cor} = $percent;
                            $total_reads1 = $total_reads1 + $readcount;
                        }
                    }
                }
                $total_reads2 = $total_reads - $total_reads1;
                $total_BP1_correct = int(($total_reads1/$total_reads2) * $total_BP2);
                $total_BP = $total_BP1_correct + $total_BP2;
    
                foreach $cor (@cors) {
                    $newline = "";
                    $percent_correct = $cor2percent{$cor};
                    $BPcount_correct = int($total_BP * $percent_correct);
                    $line = $cor2line{$cor};
                    @eles = split(/\t/,$line);
                    $length = @eles;
                    foreach $i (0..($length -1)) {
                        if ($i eq 4) {
                            $newline = $newline . "$BPcount_correct\t";
                        }
                        elsif ($i eq 5) {
                            $newline = $newline . "$percent_correct\t";
                        }
                        else {
                            $newline = $newline . "$eles[$i]\t";
                        }
                    }
                    $newline = $newline . $eles[-1];
                    print OUT "$newline\n";
                    $check3 = $check3 + $BPcount_correct;
                }
            }
        
            if ($check2 == 1 and $check3 == 0) { 
                print "Files by breakpoints and by reads do not match, need to check: $month\t$cell\n";
            }
            foreach $rest (@rests) {
                print OUT "$cor2line{$rest}\n";
            }
            close OUT;
        }
        else {
            print "File not created since at least one of corresponding file does not exist: $month\t$cell\n";
        }
        close IN_bp, IN_reads;
    }
}

#!/usr/bin/perl

use File::Copy;

@months = qw(
1m
2m
3m
6m
9m
12m
15m
18m
21m
24m
30m
36m
42m
);

@cells = qw(
CD19
CD3
NK
PMN
CD14
);

@ISlist = ();
foreach $cell (@cells) {
    opendir DIR, "/Volumes/Seagate/Data/GeneTherapyProject/PT2716/$cell" or die "Cannot open directory:$!";
    @files = readdir DIR;
    foreach $file (@files) {
        @filename = split (/\_/,$file);
        $month = $filename[0];
        open IN, "</Volumes/Seagate/Data/GeneTherapyProject/PT2716/$cell/$file" or die;
        
        $total = "total_" . $cell . "_" . $month;
        $$total = 0;
        $cell_month = $cell . "_" . $month;
        $n = 0;
        while ($line=<IN>) {
            $n++;
            if ($n > 10) {
                chomp($line);
                @eles = split(/\t/,$line);
                $IS = $eles[1].$eles[2].$eles[3];
                $IS_cell_month = $IS . "&" . $cell . "&" . $month;
                $count = $eles[4];
                $$IS_cell_month = $$IS_cell_month + $count;
                $$total = $$total + $count;
                unless ($seen{$IS}++) {
                    push (@ISlist, $IS);
                }
                unless ($cell_month{$IS_cell_month}++) {
                    push (@$cell_month, $IS_cell_month);
                }
            }            
        }
        
        foreach $IS_cell_month (@$cell_month) {
            $IS2ratio{$IS_cell_month} = ($$IS_cell_month/$$total)*100;
        }
        close IN;
    }
    closedir DIR;
}

open OUT, ">all_IS_frequencies_cross_all_samples.txt";
print OUT "IS\t";
foreach $cell (@cells) {
    foreach $month (@months) {
        print OUT "$cell\_$month\t";
    }
}
print OUT "sum\n";

foreach $IS (@ISlist) {
    $sum = 0;
    $newline = "$IS\t";
    foreach $cell (@cells) {
        foreach $month (@months) {
            $IS_cell_month = $IS . "&" . $cell . "&" . $month;
            if ($IS2ratio{$IS_cell_month}) {
                $newline = $newline . "$IS2ratio{$IS_cell_month}\t";
            }
            else {
                $newline = $newline .  "0\t";
            }
            $sum = $sum + $IS2ratio{$IS_cell_month};
        }
    }   
    $newline = $newline .  "$sum\n";
    push (@newlines, $newline);
}

@sorted_lines = sort {
	@a_fields = split /\t/, $a;
	@b_fields = split /\t/, $b;	
	
	$b_fields[-1] <=> $a_fields[-1]	
				
} @newlines;

foreach $sorted_line (@sorted_lines) {
    print OUT "$sorted_line";
}
    
close OUT;


open OUTR, ">all_IS_frequencies_cross_all_samples_for_line_plot.txt";
print OUTR "Lineage\tMonth\tIS\tFrequency\n";

foreach $IS (@ISlist) {
    foreach $cell (@cells) {
        foreach $month (@months) {
            $month_new = substr($month,0,-1);
            print OUTR "$cell\t$month_new\t$IS\t";
            $IS_cell_month = $IS . "&" . $cell . "&" . $month;

            if ($IS2ratio{$IS_cell_month}) {
            $freq = $IS2ratio{$IS_cell_month};
                print OUTR "$freq\n";
            }
            else {
                print OUTR  "0\n";
            }
        }
    }   
}

close OUTR;

@last_ISlist = ();
foreach $cell (@cells) {
    opendir DIR, "/Volumes/Seagate/Data/GeneTherapyProject/PT2716/$cell" or die "Cannot open directory:$!";
    @files = readdir DIR;
    foreach $file (@files) {
        @filename = split (/\_/,$file);
        $month = $filename[0];
        
        if ($month eq $months[-1]) {
            open IN, "</Volumes/Seagate/Data/GeneTherapyProject/PT2716/$cell/$file" or die;
            
            $last_total = "last_total_" . $cell . "_" . $month;
            $$last_total = 0;
            $last_cell_month = "last_" . $cell . "_" . $month;
            $n = 0;
            while ($line=<IN>) {
                $n++;
                if ($n > 10) {
                    chomp($line);
                    @eles = split(/\t/,$line);
                    $IS = $eles[1].$eles[2].$eles[3];;
                    $last_IS_cell_month = "last_" . $IS . "&" . $cell . "&" . $month;
                    $count = $eles[4];
                    $$last_IS_cell_month = $$last_IS_cell_month + $count;
                    $$last_total = $$last_total + $count;
                    unless ($last_seen{$IS}++) {
                        push (@last_ISlist, $IS);
                    }
                    unless ($last_cell_month{$last_IS_cell_month}++) {
                        push (@$last_cell_month, $last_IS_cell_month);
                    }
                }            
            }
            foreach $last_IS_cell_month (@$last_cell_month) {
                $lastIS2ratio{$last_IS_cell_month} = ($$last_IS_cell_month/$$last_total)*100;
            }
            close IN;
        }
    }
    closedir DIR;
}

open LAST, ">all_IS_frequencies_cross_last_time_point.txt";
print LAST "IS\t";
foreach $cell (@cells) {
    if ($cell =~ /CD19/) {
        print LAST "B\_$months[-1]\t";
    }
    elsif ($cell =~ /CD3/) {
        print LAST "T\_$months[-1]\t";
    }
    else {
        print LAST "$cell\_$months[-1]\t";
    }
}
print LAST "\n";

@last_new_lines = ();
foreach $IS (@last_ISlist) {

    $last_newline = "$IS\t";
    foreach $cell (@cells) {
        $last_IS_cell_month = "last_" . $IS . "&" . $cell . "&" . $months[-1];
        if ($lastIS2ratio{$last_IS_cell_month}) {
            $last_newline = $last_newline . "$lastIS2ratio{$last_IS_cell_month}\t";
        }
        else {
            $last_newline = $last_newline .  "0\t";
        }
    }   
    push (@last_newlines, $last_newline);
}

@last_sorted_lines = sort {
	@a_fields = split /\t/, $a;
	@b_fields = split /\t/, $b;	
	
	$a_fields[3] <=> $b_fields[3]	
				
} @last_newlines;


foreach $last_sorted_line (@last_sorted_lines) {
    print LAST "$last_sorted_line\n";
}
    
close LAST;


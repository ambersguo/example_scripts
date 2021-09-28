#!/usr/bin/perl/

@lib = qw(
extraction10-11-17
extraction6-12-17
all
);

foreach $cell (@lib) {
    $file2in = "PBMC_CD4_final_extractions_". $cell . ".txt";
    open IN, "<$file2in" or die;
    <IN>;
    while ($line = <IN>) {
        @eles = split (/\t/,$line);
        $count = $eles[11];
        $cor = $eles[22];
        $gene_ID = $eles[23];
        $cell_cor = $cell ."_".$cor;
        $$cell_cor{$cor} = $count;
        unless ($seen{$cor}++) {
            $cor2gene{$cor} = $gene_ID;
        }
        if ($eles[11] > 2) {       
            push (@$cell,$cor);
        } 
        if ($eles[11] == 2) {
            $sample = $eles[1];
            $LTR = $eles[19];
            @id1 = split (/\-/,$eles[0]);
            $celltype = $id1[1];
            
            unless ($$cell{$cor}++) {
                $cor2sample{$cor} = $eles[1];
                $cor2LTR{$cor} = $eles[19];
                @id2 = split (/\-/,$eles[0]);
                $cor2celltype{$cor} = $id2[1];
            }
            else {
                if ($sample eq $cor2sample{$cor}) {
                    if ($LTR eq $cor2LTR{$cor}) {
                        push (@$cell,$cor);
                    }
                    else {
                        if ($celltype ne $cor2celltype{$cor}) {
                            push (@$cell,$cor);
                        }
                    }
                }
                else {
                    push (@$cell,$cor);
                }
            }
        }
    }
    close IN;
}
                
open IN2, "<all_IS_unique_PBMC_CD4_all_extractions.txt" or die;
open OUT1, ">clone_info_PBMC_CD4_extractions.txt" or die;
open OUT2, ">top_gene_info_PBMC_CD4_extractions.txt" or die;
open OUT3, ">clone_summary.txt" or die;

<IN2>;
print OUT1 "gene\tIS\textraction10-11-17_count\textraction6-12-17_count\tall_count\textraction10-11-17_mark\textraction6-12-17_mark\tall_mark\n";

%allhash = map { $_ => 1 } @all;
while ($line = <IN2>) {
    chomp ($line);
    @eles = split (/\t/,$line);
    $cor = $eles[0];
    $ex1_count = $eles[1];
    $ex2_count = $eles[2];
    $all_count = $eles[3];
    $gene = $cor2gene{$cor};
    $newline = $gene."\t".$line;
    
    foreach $cell (@lib) {
        $mark = $cell . "_mark";
        $cell_cor = $cell ."_".$cor;
        %cellhash = map { $_ => 1 } @$cell;
        if(exists($cellhash{$cor})) {
            $$mark = "clone";
        }
        else {
            if(exists($allhash{$cor}) and $$cell_cor{$cor}) {
                $$mark = "clone";
            }
            else {
                $$mark = "nonclone";
            }
        }
        $newline = $newline . "\t$$mark";
    }
    print OUT1 "$newline\n";
    
    unless ($seen2{$gene}++) {
        push (@genelist,$gene);
    }
    
    $gene2ex1count{$gene}= $gene2ex1count{$gene} + $ex1_count;
    $gene2ex2count{$gene}= $gene2ex2count{$gene} + $ex2_count;
    $gene2allcount{$gene}= $gene2allcount{$gene} + $all_count;
    
    push (@$gene,$cor);   
    
}

print OUT2 "gene\textraction10-11-17_gene_count\textraction6-12-17_gene_count\tall_gene_count\n";
foreach $gene (@genelist) {
    print OUT2 "$gene\t$gene2ex1count{$gene}\t$gene2ex2count{$gene}\t$gene2allcount{$gene}\t@$gene\n";
    foreach $cor1 (@$gene) {
        if ($cor1 =~ /(.*)([+-])(.*)/) {
            $chr1 = $1;
            $strand1 = $2;
            $IS1 = $3;
        }    
        foreach $cor2 (@$gene) {
            if ($cor2 =~ /(.*)([+-])(.*)/) {
                $chr2 = $1;
                $strand2 = $2;
                $IS2 = $3;
            }

            if (($chr1 eq $chr2) and ($strand1 eq $strand2) and (!($IS1 == $IS2)) and (abs($IS1-$IS2)<10)) {
                unless ($seen3{$cor2}++) {
                    push (@check, $cor2);
                }
            }
        }
    }
}      
        
close IN2, OUT1, OUT2;
        
foreach $cell (@lib) {
    $file2in = "PBMC_CD4_final_extractions_". $cell . ".txt";
    $file2out = "PBMC_CD4_final_extractions_". $cell . "_check.txt";
    open IN, "<$file2in" or die;
    open OUT, ">$file2out" or die;
    $header = <IN>;
    chomp ($header);
    print OUT "$header\tmark\n";
    while ($line = <IN>) {
        chomp ($line);
        @eles = split (/\t/,$line);
        $cor = $eles[22];
        %checkhash = map { $_ => 1 } @check;
        if(exists($checkhash{$cor})) {
            $mark = "need_to_check";
        }
        else {
            $mark = "";
        }
        print OUT "$line\t$mark\n";
    }
    close IN, OUT;
}

open IN3, "<clone_info_PBMC_CD4_extractions.txt" or die;
open OUT3, ">clone_summary.txt" or die;

<IN3>;
$ex1_clone_count=$ex2_clone_count=$all_clone_count=$ex1_IS_count=$ex2_IS_count=$all_IS_count=0;
while ($line = <IN3>) {
    chomp ($line);
    @eles = split (/\t/,$line); 
    $ex1_mark = $eles[5];
    $ex2_mark = $eles[6];
    $all_mark = $eles[7];
    if ($ex1_mark eq "clone") {
        $ex1_clone_count++;
        $ex1_IS_count = $ex1_IS_count + $eles[2];
    }
    if ($ex2_mark eq "clone") {
        $ex2_clone_count++;
        $ex2_IS_count = $ex2_IS_count + $eles[3];
    }    
    if ($all_mark eq "clone") {
        $all_clone_count++;
        $all_IS_count = $all_IS_count + $eles[4];        
    }
}

print OUT3 "extraction\ttotal_clone_count\ttotal_IS_count_in_clone\n";
print OUT3 "extraction10-11-17\t$ex1_clone_count\t$ex1_IS_count\n";
print OUT3 "extraction6-12-17\t$ex2_clone_count\t$ex2_IS_count\n";
print OUT3 "all\t$all_clone_count\t$all_IS_count\n";

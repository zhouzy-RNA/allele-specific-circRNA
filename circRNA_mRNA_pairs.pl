use strict;
use warnings;


my $input_linear_gtf;
my $input_circRNAs;
my $output_circRNA_gtf;
my $output_circRNAmRNA_pair;
my $minus_small2large;

my @a=@ARGV;
if(@a==10){
	for(my $i=0; $i < @a; $i+=2){
		if($a[$i] eq "-input_linear_gtf"){
			$input_linear_gtf=$a[$i+1];
		}
		if($a[$i] eq "-input_circRNAs"){
			$input_circRNAs=$a[$i+1];
		}
		if($a[$i] eq "-output_circRNA_gtf"){
			$output_circRNA_gtf=$a[$i+1];
		}
		if($a[$i] eq "-output_circRNAmRNA_pair"){
			$output_circRNAmRNA_pair=$a[$i+1];
		}
		if($a[$i] eq "-minus_small2large"){
			$minus_small2large=$a[$i+1];
		}
	}
}
else{
	print "perl circRNA_mRNA_pairs.pl -input_linear_gtf  test.protein_coding.exons.gtf  -input_circRNAs   test.merged.all.circ.temp4 -output_circRNA_gtf test.circRNAs.gtf -output_circRNAmRNA_pair  test.circRNAmRNA_pair -minus_small2large 0 or 1\n";
	exit();
}
my %gtf_trans;
my $total_tran_number;
load_gtf_file_with_details();
print "load_gtf_file_with_details done\ntotal_tran_number: $total_tran_number\n";

my %circRNA_nodes;
my $valid_circRNA_number=0;

open INPUT, "<", $input_circRNAs or die "cannot open $input_circRNAs!\n";
while(<INPUT>){
	chomp;
	my @a=split(/\s+/);
	$circRNA_nodes{$valid_circRNA_number}{"circRNA_id"}=$a[0];
	$circRNA_nodes{$valid_circRNA_number}{"chrom"}=$a[1];
	$circRNA_nodes{$valid_circRNA_number}{"tStart"}=$a[2];
	$circRNA_nodes{$valid_circRNA_number}{"tEnd"}=$a[3];
	$circRNA_nodes{$valid_circRNA_number}{"strand"}=$a[6];
	$valid_circRNA_number++;
}
close(INPUT);

print "circRNA_nodes done\nvalid_circRNA_number: $valid_circRNA_number\n";

get_circ_struct();

print "get_circ_struct done\n";

for(my $ii=0; $ii<$valid_circRNA_number; $ii++){
	$circRNA_nodes{$ii}{"length"}=0;
	if($circRNA_nodes{$ii}{"paired_mRNA_length"} > 1){
		for(my $jj=0; $jj<$circRNA_nodes{$ii}{"block_count"}; $jj++){
		
			$circRNA_nodes{$ii}{"length"} += $circRNA_nodes{$ii}{"startend"}{$jj}{2};
		}
	}
}

open OUTPUT1, ">", $output_circRNA_gtf or die "cannot open $output_circRNA_gtf!\n";
open OUTPUT2, ">", $output_circRNAmRNA_pair or die "cannot open $output_circRNAmRNA_pair!\n";

for (my $kk=0; $kk<$valid_circRNA_number; $kk++){
	if($circRNA_nodes{$kk}{"paired_mRNA_length"} > 1){
		my $temp = $circRNA_nodes{$kk}{"junct_start"} + $circRNA_nodes{$kk}{"length"} -1;
		print OUTPUT2 "$circRNA_nodes{$kk}{\"circRNA_id\"}\t$circRNA_nodes{$kk}{\"paired_mRNA_transID\"}\t$circRNA_nodes{$kk}{\"junct_start\"}\t$temp\t$circRNA_nodes{$kk}{\"strand\"}\t$circRNA_nodes{$kk}{\"paired_mRNA_start\"}\t$circRNA_nodes{$kk}{\"paired_mRNA_end\"}\t$circRNA_nodes{$kk}{\"paired_mRNA_length\"}\t$circRNA_nodes{$kk}{\"tStart\"}\t$circRNA_nodes{$kk}{\"tEnd\"}\t$circRNA_nodes{$kk}{\"length\"}\t$circRNA_nodes{$kk}{\"chrom\"}\n";
		if($circRNA_nodes{$kk}{"strand"} eq "-" && $circRNA_nodes{$kk}{"block_count"} >1){
			for(my $block_index=$circRNA_nodes{$kk}{"block_count"}-1; $block_index>=0; $block_index--){
				my $temp =$circRNA_nodes{$kk}{"block_count"} - $block_index;
				print OUTPUT1 "$circRNA_nodes{$kk}{\"chrom\"}\tcircRNA\texon\t$circRNA_nodes{$kk}{\"startend\"}{$block_index}{0}\t $circRNA_nodes{$kk}{\"startend\"}{$block_index}{1}\t.\t$circRNA_nodes{$kk}{\"strand\"}\t.\tgene_id \"$circRNA_nodes{$kk}{\"paired_mRNA_ID\"}\"; transcript_id \"$circRNA_nodes{$kk}{\"circRNA_id\"}\"; exon_number \"$temp\";\n";
			}
		}else{
			for(my $block_index=0; $block_index < $circRNA_nodes{$kk}{"block_count"}; $block_index++){
				my $temp = $block_index+1;
				
				print OUTPUT1 "$circRNA_nodes{$kk}{\"chrom\"}\tcircRNA\texon\t$circRNA_nodes{$kk}{\"startend\"}{$block_index}{0}\t $circRNA_nodes{$kk}{\"startend\"}{$block_index}{1}\t.\t$circRNA_nodes{$kk}{\"strand\"}\t.\tgene_id \"$circRNA_nodes{$kk}{\"paired_mRNA_ID\"}\"; transcript_id \"$circRNA_nodes{$kk}{\"circRNA_id\"}\"; exon_number \"$temp\";\n";			
			}
		}
	}
}


close(OUTPUT1);
close(OUTPUT2);



sub get_circ_struct{
	
	my $firstexon;
	my $lastexon;
	my $mRNA_start;
	my $mRNA_end;
	my $pairedcount=0;
	my $wellpaired=0;
	for(my $cc = 0; $cc < $valid_circRNA_number; $cc++){
		my $foundmatched = 0;
		my $maxlength = 1;
		my $tranidmaxlen = "";
		my $geneidmaxlen = "";
		for(my $mm = 0; $mm < $total_tran_number; $mm++){
			if($gtf_trans{$mm}{"tStart"}<=$circRNA_nodes{$cc}{"tStart"} && $gtf_trans{$mm}{"tEnd"}>=$circRNA_nodes{$cc}{"tEnd"} && $gtf_trans{$mm}{"strand"} eq $circRNA_nodes{$cc}{"strand"} && $gtf_trans{$mm}{"chrom"} eq $circRNA_nodes{$cc}{"chrom"}){
				$firstexon=-1;
				$lastexon=-1;
				for(my $ee=0; $ee < $gtf_trans{$mm}{"block_count"};$ee++){
					if($gtf_trans{$mm}{"startend"}{$ee}{0}<=$circRNA_nodes{$cc}{"tStart"} &&  $gtf_trans{$mm}{"startend"}{$ee}{1}>=$circRNA_nodes{$cc}{"tStart"}){
						$firstexon=$ee;
						last;
					}
				}
				for(my $ee=0; $ee < $gtf_trans{$mm}{"block_count"};$ee++){
					if($gtf_trans{$mm}{"startend"}{$ee}{0}<=$circRNA_nodes{$cc}{"tEnd"} &&  $gtf_trans{$mm}{"startend"}{$ee}{1}>=$circRNA_nodes{$cc}{"tEnd"}){
						$lastexon=$ee;
						last;
					}
				}	
				if($firstexon != -1 && $lastexon != -1){
					$foundmatched=1;
					if($gtf_trans{$mm}{"length"} > $maxlength){
						$maxlength=$gtf_trans{$mm}{"length"};
						$tranidmaxlen=$gtf_trans{$mm}{"tranid"};
						$mRNA_start=$gtf_trans{$mm}{"tStart"};
						$mRNA_end=$gtf_trans{$mm}{"tEnd"};
						$geneidmaxlen=$gtf_trans{$mm}{"gene_id"};
					}				
				}
				
			}
		}
		if($foundmatched==1){
			$circRNA_nodes{$cc}{"paired_mRNA_length"}=$maxlength;
			$circRNA_nodes{$cc}{"paired_mRNA_transID"}=$tranidmaxlen;
			$circRNA_nodes{$cc}{"paired_mRNA_ID"}=$geneidmaxlen;
			$circRNA_nodes{$cc}{"paired_mRNA_start"}=$mRNA_start;
			$circRNA_nodes{$cc}{"paired_mRNA_end"}=$mRNA_end;
			$pairedcount++;
		}else{
			$circRNA_nodes{$cc}{"paired_mRNA_length"}=0;
			$circRNA_nodes{$cc}{"paired_mRNA_transID"}="";			
		}
	}
	print "$valid_circRNA_number\t$pairedcount\nseacher max done\n";
	for(my $cc = 0; $cc < $valid_circRNA_number; $cc++){
		if($circRNA_nodes{$cc}{"paired_mRNA_length"} > 1){
			for(my $mm = 0; $mm < $total_tran_number; $mm++){
				if($circRNA_nodes{$cc}{"paired_mRNA_transID"} eq $gtf_trans{$mm}{"tranid"}){
					$firstexon=-1;
					$lastexon=-1;
					for(my $ee=0; $ee < $gtf_trans{$mm}{"block_count"};$ee++){
						if($gtf_trans{$mm}{"startend"}{$ee}{0}<=$circRNA_nodes{$cc}{"tStart"} &&  $gtf_trans{$mm}{"startend"}{$ee}{1}>=$circRNA_nodes{$cc}{"tStart"}){
							$firstexon=$ee;
							last;
						}
					}
					for(my $ee=0; $ee < $gtf_trans{$mm}{"block_count"};$ee++){
						if($gtf_trans{$mm}{"startend"}{$ee}{0}<=$circRNA_nodes{$cc}{"tEnd"} &&  $gtf_trans{$mm}{"startend"}{$ee}{1}>=$circRNA_nodes{$cc}{"tEnd"}){
							$lastexon=$ee;
							last;
						}
					}
					if($firstexon != -1 && $lastexon != -1){
						$wellpaired++;
						if($firstexon == $lastexon){
							$circRNA_nodes{$cc}{"block_count"}=1;
							$circRNA_nodes{$cc}{"startend"}{0}{0}=$circRNA_nodes{$cc}{"tStart"};
							$circRNA_nodes{$cc}{"startend"}{0}{1}=$circRNA_nodes{$cc}{"tEnd"};
							$circRNA_nodes{$cc}{"startend"}{0}{2}=$circRNA_nodes{$cc}{"tEnd"}-$circRNA_nodes{$cc}{"tStart"}+1;
						}else{
							
							$circRNA_nodes{$cc}{"block_count"}=1;
							$circRNA_nodes{$cc}{"startend"}{$circRNA_nodes{$cc}{"block_count"}-1}{0}=$circRNA_nodes{$cc}{"tStart"};
							$circRNA_nodes{$cc}{"startend"}{$circRNA_nodes{$cc}{"block_count"}-1}{1}=$gtf_trans{$mm}{"startend"}{$firstexon}{1};
							$circRNA_nodes{$cc}{"startend"}{$circRNA_nodes{$cc}{"block_count"}-1}{2}=$gtf_trans{$mm}{"startend"}{$firstexon}{1} - $circRNA_nodes{$cc}{"tStart"} + 1;
								
							if($circRNA_nodes{$cc}{"strand"} eq "+" || ($circRNA_nodes{$cc}{"strand"} eq "-" && $minus_small2large==1)){
								for(my $kk=1; $kk<($lastexon-$firstexon); $kk++){
									$circRNA_nodes{$cc}{"block_count"}++;
									$circRNA_nodes{$cc}{"startend"}{$circRNA_nodes{$cc}{"block_count"}-1}{0}=$gtf_trans{$mm}{"startend"}{$firstexon+$kk}{0};
									$circRNA_nodes{$cc}{"startend"}{$circRNA_nodes{$cc}{"block_count"}-1}{1}=$gtf_trans{$mm}{"startend"}{$firstexon+$kk}{1};
									$circRNA_nodes{$cc}{"startend"}{$circRNA_nodes{$cc}{"block_count"}-1}{2}=$gtf_trans{$mm}{"startend"}{$firstexon+$kk}{1} - $gtf_trans{$mm}{"startend"}{$firstexon+$kk}{0} + 1;
								}				
							}else{
								for(my $kk=$firstexon-1; $kk>$lastexon; $kk--){
									$circRNA_nodes{$cc}{"block_count"}++;
									$circRNA_nodes{$cc}{"startend"}{$circRNA_nodes{$cc}{"block_count"}-1}{0}=$gtf_trans{$mm}{"startend"}{$kk}{0};
									$circRNA_nodes{$cc}{"startend"}{$circRNA_nodes{$cc}{"block_count"}-1}{1}=$gtf_trans{$mm}{"startend"}{$kk}{1};
									$circRNA_nodes{$cc}{"startend"}{$circRNA_nodes{$cc}{"block_count"}-1}{2}=$gtf_trans{$mm}{"startend"}{$kk}{1} - $gtf_trans{$mm}{"startend"}{$kk}{0} + 1;
								}				
							}
								
							$circRNA_nodes{$cc}{"block_count"}++;
							$circRNA_nodes{$cc}{"startend"}{$circRNA_nodes{$cc}{"block_count"}-1}{0}=$gtf_trans{$mm}{"startend"}{$lastexon}{0};
							$circRNA_nodes{$cc}{"startend"}{$circRNA_nodes{$cc}{"block_count"}-1}{1}=$circRNA_nodes{$cc}{"tEnd"};
							$circRNA_nodes{$cc}{"startend"}{$circRNA_nodes{$cc}{"block_count"}-1}{2}=$circRNA_nodes{$cc}{"tEnd"} - $gtf_trans{$mm}{"startend"}{$lastexon}{0} +1;
						}
						if($circRNA_nodes{$cc}{"strand"} eq "+"){
								my $junct_start=0; 
								for(my $jjj=0; $jjj<$firstexon; $jjj++){
									$junct_start += $gtf_trans{$mm}{"startend"}{$jjj}{2};
								}
								$junct_start = $junct_start + $circRNA_nodes{$cc}{"tStart"} - $gtf_trans{$mm}{"startend"}{$firstexon}{0};
								$circRNA_nodes{$cc}{"junct_start"} = $junct_start;
						}elsif(($circRNA_nodes{$cc}{"strand"} eq "-" && $minus_small2large==1)){
								my $junct_start=0;
								for(my $jjj=$lastexon+1; $jjj<$gtf_trans{$mm}{"block_count"}; $jjj++){
									$junct_start += $gtf_trans{$mm}{"startend"}{$jjj}{2};
								}
								$junct_start = $junct_start + $gtf_trans{$mm}{"startend"}{$firstexon}{1} -  $circRNA_nodes{$cc}{"tEnd"};
								$circRNA_nodes{$cc}{"junct_start"} = $junct_start;	
						}else{
								my $junct_start=0;
								for(my $jjj=$lastexon-1; $jjj>=0; $jjj--){
									$junct_start += $gtf_trans{$mm}{"startend"}{$jjj}{2};
								}
								$junct_start = $junct_start + $gtf_trans{$mm}{"startend"}{$lastexon}{1} -  $circRNA_nodes{$cc}{"tEnd"};
								$circRNA_nodes{$cc}{"junct_start"} = $junct_start;				
						}					
					}
				}					
			}
		}
	}
}		

sub load_gtf_file_with_details{
	open INPUT, "<", $input_linear_gtf or die "cannot open $input_linear_gtf!\n";
	my $trans_index=-1;
	while(<INPUT>){
		chomp;
		last if(/^a/);
		my @a=split(/\t/);
		my $gene_id;
		my $transcript_id;
		my $exon_number;  
		#gene_id "ENSG00000186092"; transcript_id "ENST00000335137"; exon_number "1";
		if($a[8]=~/gene_id "(.+?)";/){
			$gene_id=$1;
		}
		if($a[8]=~/transcript_id "(.+?)";/){
			$transcript_id=$1;
		}
		if($a[8]=~/exon_number "(.+?)";/){
			$exon_number=$1;
		}
		my $gene_name;
		my $trans_name;
		my $biotype_name;
		#//gene_name "CLIC6"; gene_biotype "protein_coding"; transcript_name "CLIC6-201";
		if($a[8]=~/gene_name "(.+?)";/){
			$gene_name = $1;
		}else{
			$gene_name = "";
		}
		if($a[8]=~/gene_biotype "(.+?)";/){
			$biotype_name = $1;
		}else{
			$biotype_name = "";
		}
		if($a[8]=~/transcript_name "(.+?)";/){
			$trans_name = $1;
		}else{
			$trans_name = "";
		}
		if($exon_number==1){
			$trans_index += 1;
			$gtf_trans{$trans_index}{"tranid"}=$transcript_id;
			$gtf_trans{$trans_index}{"chrom"}=$a[0];
			$gtf_trans{$trans_index}{"gene_id"}=$gene_id;
			$gtf_trans{$trans_index}{"gene_name"}=$gene_name;
			$gtf_trans{$trans_index}{"tran_name"}=$trans_name;
			$gtf_trans{$trans_index}{"biotype_name"}=$biotype_name;
			$gtf_trans{$trans_index}{"strand"}=$a[6];
			$gtf_trans{$trans_index}{"block_count"}=1;
			$gtf_trans{$trans_index}{"startend"}{0}{0}=$a[3];
			$gtf_trans{$trans_index}{"startend"}{0}{1}=$a[4];
			$gtf_trans{$trans_index}{"startend"}{0}{2}=$a[4]-$a[3]+1;

		}else{
			$gtf_trans{$trans_index}{"block_count"}++;
			$gtf_trans{$trans_index}{"startend"}{$gtf_trans{$trans_index}{"block_count"}-1}{0}=$a[3];
			$gtf_trans{$trans_index}{"startend"}{$gtf_trans{$trans_index}{"block_count"}-1}{1}=$a[4];
			$gtf_trans{$trans_index}{"startend"}{$gtf_trans{$trans_index}{"block_count"}-1}{2}=$a[4]-$a[3]+1;
		}
		
	}
	close(INPUT);
	$total_tran_number=$trans_index;
	for(my $j = 0; $j < $total_tran_number; $j++){
		if($gtf_trans{$j}{"strand"} eq "+"){
			$gtf_trans{$j}{"tStart"} = $gtf_trans{$j}{"startend"}{0}{0};
			$gtf_trans{$j}{"tEnd"} = $gtf_trans{$j}{"startend"}{$gtf_trans{$j}{"block_count"}-1}{1};
		}else{
			if($minus_small2large == 1){ #负链由小到大
				$gtf_trans{$j}{"tStart"} = $gtf_trans{$j}{"startend"}{0}{0};
				$gtf_trans{$j}{"tEnd"} = $gtf_trans{$j}{"startend"}{$gtf_trans{$j}{"block_count"}-1}{1};
			}else{#负链由大到小
				$gtf_trans{$j}{"tStart"} = $gtf_trans{$j}{"startend"}{$gtf_trans{$j}{"block_count"}-1}{0};
				$gtf_trans{$j}{"tEnd"} = $gtf_trans{$j}{"startend"}{0}{1};
			}
		}
	}
	for(my $j = 0; $j < $total_tran_number; $j++){
		for(my $i = 0; $i < $gtf_trans{$j}{"block_count"}; $i++){
			$gtf_trans{$j}{"length"} += $gtf_trans{$j}{"startend"}{$i}{2};
		}	
	}	
}















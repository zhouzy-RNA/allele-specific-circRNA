use strict;
use warnings;


my $input_pair;
my $input_sam;
my $input_snp;
my $input_mRNA_simple_bed;
my $output_reformed;
my $chr;

my @a=@ARGV;
my $num=@a;
print "$num\n";
if(@a==12){
	for(my $i=0; $i < @a; $i+=2){
		if($a[$i] eq "-input_pair"){
			$input_pair=$a[$i+1];
		}
		if($a[$i] eq "-input_sam"){
			$input_sam=$a[$i+1];
		}
		if($a[$i] eq "-input_snp"){
			$input_snp=$a[$i+1];
		}
		if($a[$i] eq "-input_mRNA_simple_bed"){
			$input_mRNA_simple_bed=$a[$i+1];
		}
		if($a[$i] eq "-output_reformed"){
			$output_reformed=$a[$i+1];
		}
		if($a[$i] eq "-chr"){
			$chr=$a[$i+1];
		}
	}
}
else{
	print "perl statisticscomplex.pl -input_pair test.circRNAmRNA_pair.numberid.txt -input_sam  test.sam.reformed -input_snp chr1_overlapping.chr_posit.vcf -input_mRNA_simple_bed    kept_longest_mRNAs.1_22.NewNumberID.bed.txt -output_reformed test.readscoverjuncsnp -chr 1\n";
	exit();
}





my $max_block_number=1000;

open(A,"$input_pair") or die "cannot open $input_pair!\n";
my %pair;
my $pairnum=0;
while(<A>){
	chomp;
	my @a = split(/\t/);
	next if($a[11] != $chr);
	$pair{$pairnum}{"circRNAid"}= $a[0];
	$pair{$pairnum}{"mRNAid"}= $a[1];
	$pair{$pairnum}{"junct_start"}= $a[2] + 1;
	$pair{$pairnum}{"junct_end"}= $a[3] + 1;
	$pair{$pairnum}{"strand"}= $a[4];
	$pair{$pairnum}{"mRNA_start"}= $a[5];
	$pair{$pairnum}{"mRNA_end"}= $a[6];
	$pair{$pairnum}{"mRNA_length"}= $a[7];
	$pair{$pairnum}{"cirRNA_start"}= $a[8];
	$pair{$pairnum}{"cirRNA_end"}= $a[9];
	$pair{$pairnum}{"cirRNA_length"}= $a[10];
	$pair{$pairnum}{"chromosome"}= $a[11];
	$pairnum++;	
}
close(A);

open(A,"$input_snp") or die "cannot open $input_snp!\n";
my %snp;
my $snpnum=0;
while(<A>){
	chomp;
	my @a = split(/\s+/);
	$snp{$snpnum}{"chrom"}= $a[0];
	$snp{$snpnum}{"posit"}= $a[1];
	$snp{$snpnum}{"rs"}= $a[2];
	$snp{$snpnum}{"ref"}= $a[3];
	$snp{$snpnum}{"alt"}= $a[4];
	$snpnum++;	
}
close(A);



my %mRNA_simple_bed;
(my $mRNA_simple_bednum, *mRNA_simple_bed)=load_simple_bed($input_mRNA_simple_bed, \%mRNA_simple_bed);

*mRNA_simple_bed = calculate_relativeposition_of_eachexon_simple_bed($mRNA_simple_bednum, \%mRNA_simple_bed);


open(A,"$input_sam") or die "cannot open $input_sam!\n";
open(C,">$output_reformed") or die "cannot open $output_reformed!\n";
my %sam;
my $read_start_relative=0;
my $read_end_relative=0;

my $strand;
my $junction_relative;

my $read_absolute_start=0;
my $read_absolute_end=0;

while(<A>){
	chomp;
	next if(/^[@#]/);
	my @a=split(/\t/);
	$sam{"qname"}=$a[0];
	$sam{"rname"}=$a[1];
	$sam{"pos"}=$a[2];
	$read_start_relative=$a[2];
	$sam{"cigar"}=$a[3];
	$sam{"rnext"}=$a[4];
	$sam{"pnext"}=$a[5];
	$sam{"tlen"}=$a[6];			
	$sam{"mergedread"}=$a[7];
	$sam{"end"}=$a[8];	
	$read_end_relative=$a[8];
 	my $iscircRNA=-1;
 	my $iscoverjunction=0;
	my $tempvalue;
	my %out_node;
	my $cover_junction_absolute;
	my %outsnp;
	my $matchedsnpcount;
	my $msms;
	my $snp_relative_position;
	my $read_1bp;
	my $read_left_length;
	my $read_right_length;
	my $read_right_seq;
	my $read_left_seq;
	my $out_nucleotide_position;
	my $junction_absolute_start=0;
	my $junction_absolute_end=0;
	my $temp_value;
	for(my $pair_idx=0; $pair_idx< $pairnum; $pair_idx++){
		if($pair{$pair_idx}{"mRNAid"} == $sam{"rname"}){
			$strand = $pair{$pair_idx}{"strand"};
			$iscircRNA=0;
			$iscoverjunction=0;
			if($sam{"pos"} <= $pair{$pair_idx}{"junct_start"} && $sam{"end"} >= $pair{$pair_idx}{"junct_start"}){
				$iscoverjunction = 1;
				$junction_relative = $pair{$pair_idx}{"junct_start"};				
			} 
			if($sam{"pos"} <= $pair{$pair_idx}{"junct_end"} && $sam{"end"} >= $pair{$pair_idx}{"junct_end"}){
				$iscoverjunction = 1;
				$junction_relative = $pair{$pair_idx}{"junct_end"};	
			}			
			if($iscoverjunction == 1){

				($read_absolute_start, $strand) = get_absolutepositon_of_relativepostion($pair{$pair_idx}{"mRNAid"}, $sam{"pos"}, $mRNA_simple_bednum, \%mRNA_simple_bed);
				($read_absolute_end, $strand) = get_absolutepositon_of_relativepostion($pair{$pair_idx}{"mRNAid"}, $sam{"end"}, $mRNA_simple_bednum, \%mRNA_simple_bed);
				if($read_absolute_start>$read_absolute_end) {
					$temp_value=$read_absolute_start;
					$read_absolute_start=$read_absolute_end;
					$read_absolute_end=$temp_value;
				}			
				*out_node = get_exons_between_two_absolute_positions($pair{$pair_idx}{"mRNAid"}, $read_absolute_start, $read_absolute_end, $mRNA_simple_bednum, \%mRNA_simple_bed, \%out_node);
				($junction_absolute_start, $strand) = get_absolutepositon_of_relativepostion($pair{$pair_idx}{"mRNAid"}, $pair{$pair_idx}{"junct_start"}, $mRNA_simple_bednum, \%mRNA_simple_bed);
				($junction_absolute_end, $strand) = get_absolutepositon_of_relativepostion($pair{$pair_idx}{"mRNAid"}, $pair{$pair_idx}{"junct_end"}, $mRNA_simple_bednum, \%mRNA_simple_bed);				
				if($sam{"pos"} <= $pair{$pair_idx}{"junct_start"} && $sam{"end"} >= $pair{$pair_idx}{"junct_start"}){
					$cover_junction_absolute = $junction_absolute_start;
				}
				if($sam{"pos"} <= $pair{$pair_idx}{"junct_end"} && $sam{"end"} >= $pair{$pair_idx}{"junct_end"}){
					$cover_junction_absolute = $junction_absolute_end;
				}
				if($junction_absolute_start>$junction_absolute_end) { 
					$temp_value=$junction_absolute_start;
					$junction_absolute_start=$junction_absolute_end;
					$junction_absolute_end=$temp_value;
				}			
				($matchedsnpcount, *outsnp) = are_snps_on_exons(\%out_node, \%snp, $snpnum, \%outsnp, $junction_absolute_start, $junction_absolute_end);
				if($matchedsnpcount >= 1){
					for($msms=0; $msms<$matchedsnpcount; $msms++){
						$snp_relative_position = get_relativepositon_of_absolutepostion(\%mRNA_simple_bed, $mRNA_simple_bednum, $pair{$pair_idx}{"mRNAid"}, $outsnp{$msms}{"posit"});
						$read_1bp = substr($sam{"mergedread"}, $snp_relative_position - $sam{"pos"}, 1);
						if($strand eq "-"){
							if($read_1bp eq "A"){
								$read_1bp = "T";
							}elsif($read_1bp eq "T"){
								$read_1bp = "A";
							}elsif($read_1bp eq "G"){
								$read_1bp = "C";
							}elsif($read_1bp eq "C"){
								$read_1bp = "G";
							}
						}
						if($read_1bp eq $outsnp{$msms}{"ref"} || $read_1bp eq $outsnp{$msms}{"alt"}){
							print C "hitmRNA\t$sam{\"rname\"}\t$sam{\"qname\"}\t$pair{$pair_idx}{\"circRNAid\"}\t$pair{$pair_idx}{\"mRNAid\"}\t$outsnp{$msms}{\"posit\"}\t$snp_relative_position\t$outsnp{$msms}{\"ref\"}\t$outsnp{$msms}{\"alt\"}\t$outsnp{$msms}{\"rs\"}\t$read_1bp\t";
							print C "$junction_relative\t$read_absolute_start\t$read_absolute_end\t$read_start_relative\t$read_end_relative\t$strand\t$cover_junction_absolute\t$sam{\"mergedread\"}\n";
						}
					}
				}	
			}
			next if($pair_idx + 1 == $pairnum);
			if($pair{$pair_idx + 1}{"mRNAid"} == $pair{$pair_idx}{"mRNAid"}){
				next;
			}else{
				last;
			}
		}

		if($pair{$pair_idx}{"circRNAid"} == $sam{"rname"}){
			$strand = $pair{$pair_idx}{"strand"};
			$iscircRNA=0;	
			if($pair{$pair_idx}{"cirRNA_length"} % 2 ==0){
				$junction_relative = $pair{$pair_idx}{"cirRNA_length"} / 2 + 1;
				if($sam{"pos"} <= $junction_relative && $sam{"end"} >= $junction_relative){
					$iscoverjunction = 1;
					$read_left_length = $junction_relative - $sam{"pos"};
					$read_right_length = $sam{"end"} - $sam{"pos"} + 1 -$read_left_length;
				}
			}else{
				$junction_relative = int($pair{$pair_idx}{"cirRNA_length"} / 2) + 2;
				if($sam{"pos"} <= $junction_relative && $sam{"end"} >= $junction_relative){
					$iscoverjunction = 1;
					$read_left_length = $junction_relative - $sam{"pos"};
					$read_right_length = $sam{"end"} - $sam{"pos"} + 1 -$read_left_length;
				}				
			}
			
			if($iscoverjunction==1){
				$read_right_seq = substr($sam{"mergedread"}, $read_left_length, $read_right_length);
				$read_left_seq = substr($sam{"mergedread"}, 0, $read_left_length);				
			
				($read_absolute_start, $strand) = get_absolutepositon_of_relativepostion($pair{$pair_idx}{"mRNAid"}, $pair{$pair_idx}{"junct_start"}, $mRNA_simple_bednum, \%mRNA_simple_bed);
				($read_absolute_end, $strand) = get_absolutepositon_of_relativepostion($pair{$pair_idx}{"mRNAid"}, $pair{$pair_idx}{"junct_start"} + $read_right_length -1, $mRNA_simple_bednum, \%mRNA_simple_bed);
				if($read_absolute_start>$read_absolute_end) {
					$temp_value=$read_absolute_start;
					$read_absolute_start=$read_absolute_end;
					$read_absolute_end=$temp_value;
				}
				*out_node = get_exons_between_two_absolute_positions($pair{$pair_idx}{"mRNAid"}, $read_absolute_start, $read_absolute_end, $mRNA_simple_bednum, \%mRNA_simple_bed, \%out_node);			
							
				($junction_absolute_start, $strand) = get_absolutepositon_of_relativepostion($pair{$pair_idx}{"mRNAid"}, $pair{$pair_idx}{"junct_start"}, $mRNA_simple_bednum, \%mRNA_simple_bed);
				($junction_absolute_end, $strand) = get_absolutepositon_of_relativepostion($pair{$pair_idx}{"mRNAid"}, $pair{$pair_idx}{"junct_end"}, $mRNA_simple_bednum, \%mRNA_simple_bed);
				$cover_junction_absolute=$junction_absolute_start;

				if($junction_absolute_start>$junction_absolute_end) {
					$temp_value=$junction_absolute_start;
					$junction_absolute_start=$junction_absolute_end;
					$junction_absolute_end=$temp_value;
				}	

				($matchedsnpcount, *outsnp) = are_snps_on_exons(\%out_node, \%snp, $snpnum, \%outsnp, $junction_absolute_start, $junction_absolute_end);

				if($matchedsnpcount >= 1){
					for($msms=0; $msms<$matchedsnpcount; $msms++){
						$snp_relative_position = get_relativepositon_of_absolutepostion(\%mRNA_simple_bed, $mRNA_simple_bednum, $pair{$pair_idx}{"mRNAid"}, $outsnp{$msms}{"posit"});
						$out_nucleotide_position = $snp_relative_position - $pair{$pair_idx}{"junct_start"};
						$read_1bp = substr($read_right_seq, $out_nucleotide_position, 1);
						if($strand eq "-"){
							if($read_1bp eq "A"){
								$read_1bp = "T";
							}elsif($read_1bp eq "T"){
								$read_1bp = "A";
							}elsif($read_1bp eq "G"){
								$read_1bp = "C";
							}elsif($read_1bp eq "C"){
								$read_1bp = "G";
							}
						}
						if($read_1bp eq $outsnp{$msms}{"ref"} || $read_1bp eq $outsnp{$msms}{"alt"}){
							print C "hitcircRNA\t$sam{\"rname\"}\t$sam{\"qname\"}\t$pair{$pair_idx}{\"circRNAid\"}\t$pair{$pair_idx}{\"mRNAid\"}\t$outsnp{$msms}{\"posit\"}\t$snp_relative_position\t$outsnp{$msms}{\"ref\"}\t$outsnp{$msms}{\"alt\"}\t$outsnp{$msms}{\"rs\"}\t$read_1bp\t";
							print C "$pair{$pair_idx}{\"junct_start\"}\t$read_absolute_start\t$read_absolute_end\t$read_start_relative\t$read_end_relative\t$strand\t$cover_junction_absolute\t$read_right_seq\n";
						}
					}
				}
				
				($read_absolute_start, $strand) = get_absolutepositon_of_relativepostion($pair{$pair_idx}{"mRNAid"}, $pair{$pair_idx}{"junct_end"} - $read_left_length + 1, $mRNA_simple_bednum, \%mRNA_simple_bed);
				($read_absolute_end, $strand) = get_absolutepositon_of_relativepostion($pair{$pair_idx}{"mRNAid"}, $pair{$pair_idx}{"junct_end"}, $mRNA_simple_bednum, \%mRNA_simple_bed);
				if($read_absolute_start>$read_absolute_end) {
					$temp_value=$read_absolute_start;
					$read_absolute_start=$read_absolute_end;
					$read_absolute_end=$temp_value;
				}
				
				
				*out_node = get_exons_between_two_absolute_positions($pair{$pair_idx}{"mRNAid"}, $read_absolute_start, $read_absolute_end, $mRNA_simple_bednum, \%mRNA_simple_bed, \%out_node);
				if($strand eq "-"){
					$cover_junction_absolute=$junction_absolute_start;
				}else{
					$cover_junction_absolute=$junction_absolute_end;
				}
	
				($matchedsnpcount, *outsnp) = are_snps_on_exons(\%out_node, \%snp, $snpnum, \%outsnp, $junction_absolute_start, $junction_absolute_end);
				if($matchedsnpcount >= 1){
					for($msms=0; $msms<$matchedsnpcount; $msms++){
						$snp_relative_position = get_relativepositon_of_absolutepostion(\%mRNA_simple_bed, $mRNA_simple_bednum, $pair{$pair_idx}{"mRNAid"}, $outsnp{$msms}{"posit"});
						$out_nucleotide_position = $read_left_length - ($pair{$pair_idx}{"junct_end"} - $snp_relative_position);
						$read_1bp = substr($read_left_seq, $out_nucleotide_position - 1, 1);
						if($strand eq "-"){
							if($read_1bp eq "A"){
								$read_1bp = "T";
							}elsif($read_1bp eq "T"){
								$read_1bp = "A";
							}elsif($read_1bp eq "G"){
								$read_1bp = "C";
							}elsif($read_1bp eq "C"){
								$read_1bp = "G";
							}
						}
						if($read_1bp eq $outsnp{$msms}{"ref"} || $read_1bp eq $outsnp{$msms}{"alt"}){
							print C "hitcircRNA\t$sam{\"rname\"}\t$sam{\"qname\"}\t$pair{$pair_idx}{\"circRNAid\"}\t$pair{$pair_idx}{\"mRNAid\"}\t$outsnp{$msms}{\"posit\"}\t$snp_relative_position\t$outsnp{$msms}{\"ref\"}\t$outsnp{$msms}{\"alt\"}\t$outsnp{$msms}{\"rs\"}\t$read_1bp\t";
							print C "$pair{$pair_idx}{\"junct_start\"}\t$read_absolute_start\t$read_absolute_end\t$read_start_relative\t$read_end_relative\t$strand\t$cover_junction_absolute\t$read_left_seq\n";
						}
					}
				}				
			}
			last;
		}
	}	
	
}
close(A);
close(B);
close(C);

sub are_snps_on_exons{
	my ($outnode, $snp1, $snp1num, $out_snp, $start, $end)=@_;
	my $matchedsnpcount1=0;
	my $snp_position;
	my $minsite = findclosesite($start,\%$snp1,$snp1num);
	my $maxsite = findclosesite($end,\%$snp1,$snp1num) + 1;
	for(my $ss=$minsite; $ss < $maxsite; $ss++){
		$snp_position = $snp1->{$ss}->{"posit"};
		next if(!$outnode->{"block_count"});
		for(my $pp=0; $pp < $outnode->{"block_count"}; $pp++){
			
			next if($start > $snp_position || $end < $snp_position);
			if($outnode->{"startend"}->{$pp}->{0} <= $snp_position && $outnode->{"startend"}->{$pp}->{1} >= $snp_position){
				$out_snp->{$matchedsnpcount1}->{"posit"}=$snp_position;
				$out_snp->{$matchedsnpcount1}->{"ref"}=$snp1->{$ss}->{"ref"};
				$out_snp->{$matchedsnpcount1}->{"alt"}=$snp1->{$ss}->{"alt"};
				$out_snp->{$matchedsnpcount1}->{"rs"}=$snp1->{$ss}->{"rs"};
				$matchedsnpcount1++;
				last;
			}			
		}
	}	
	return ($matchedsnpcount1, $out_snp);
}

sub findclosesite{
	my ($site, $snp2, $snp2num)=@_;
#	print "$site\t$$snp2num\t$snp2->{1}->{"posit"}";
#	sleep(1);
	my $closesite=0;
	my $left=0;
	my $right=$snp2num;
	my $i;
	if($site <= $snp2->{0}->{"posit"}){
		$i=0;
	}elsif($site >= $snp2->{$snp2num - 1}->{"posit"}){
		$i=$snp2num - 1;
	}else{
		$i = int(($left+$right)/2);
		while($snp2->{$i}->{"posit"}>$site || $snp2->{$i+1}->{"posit"}<$site){
			if($snp2->{$i}->{"posit"}>$site){
				$right=$i;
			}
			if($snp2->{$i+1}->{"posit"}<$site){
				$left=$i;
			}
			$i=int(($left+$right)/2);
		}		
	}
	return $i;
}

sub get_exons_between_two_absolute_positions{
	my ($pairRNAid, $start, $tEnd,  $RNAnum, $RNA, $outnode) = @_;
	my $tran_idx;
	my $exon_idx;
	my $firstexon;
	my $lastexon;
	my $wellpaired=0;
	my $kk;
	my $jjj;
	my $junct_start=0;
	my $ee;
	for($tran_idx=0; $tran_idx<$RNAnum; $tran_idx++){
		if($RNA->{$tran_idx}->{"transcript_numberId"} == $pairRNAid){
					$firstexon = -1;
					$lastexon = -1;
					for($ee = 0; $ee < $RNA->{$tran_idx}->{"block_count"}; $ee++){
						
						if($RNA->{$tran_idx}->{"startend"}->{$ee}->{0} <= $start && $RNA->{$tran_idx}->{"startend"}->{$ee}->{1} >= $start){
							$firstexon = $ee;
							last;
						}
					}
					for($ee = 0; $ee < $RNA->{$tran_idx}->{"block_count"}; $ee++){
						if($RNA->{$tran_idx}->{"startend"}->{$ee}->{0} <= $tEnd && $RNA->{$tran_idx}->{"startend"}->{$ee}->{1} >= $tEnd){
							$lastexon = $ee;
							last;
						}
					}
					if($firstexon != -1 && $lastexon != -1){
						if($RNA->{$tran_idx}->{"strand"} eq "+"){
							if($firstexon == $lastexon){
								$outnode->{"block_count"} = 1;
								$outnode->{"startend"}->{0}->{0} = $start;
								$outnode->{"startend"}->{0}->{1} = $tEnd;
								$outnode->{"startend"}->{0}->{2} = $tEnd - $start +1;
							}else{
								$outnode->{"block_count"} = 1;
								$outnode->{"startend"}->{0}->{0} = $start;
								$outnode->{"startend"}->{0}->{1} = $RNA->{$tran_idx}->{"startend"}->{$firstexon}->{1};
								$outnode->{"startend"}->{0}->{2} = $RNA->{$tran_idx}->{"startend"}->{$firstexon}->{1} - $start +1;
								for($kk=1;$kk<($lastexon - $firstexon); $kk++){
									$outnode->{"block_count"}++;
									$outnode->{"startend"}->{$kk}->{0} = $RNA->{$tran_idx}->{"startend"}->{$firstexon+$kk}->{0};
									$outnode->{"startend"}->{$kk}->{1} = $RNA->{$tran_idx}->{"startend"}->{$firstexon+$kk}->{1};
									$outnode->{"startend"}->{$kk}->{2} = $RNA->{$tran_idx}->{"startend"}->{$firstexon+$kk}->{1} - $RNA->{$tran_idx}->{"startend"}->{$firstexon+$kk}->{0} + 1;
								}
								$outnode->{"block_count"}++;
								$outnode->{"startend"}->{$kk}->{0} = $RNA->{$tran_idx}->{"startend"}->{$lastexon}->{0};
								$outnode->{"startend"}->{$kk}->{1} = $tEnd;
								$outnode->{"startend"}->{$kk}->{2} = $tEnd - $RNA->{$tran_idx}->{"startend"}->{$lastexon}->{0} + 1;
							}					
						}else{
							if($firstexon == $lastexon){
								$outnode->{"block_count"} = 1;
								$outnode->{"startend"}->{0}->{0} = $start;
								$outnode->{"startend"}->{0}->{1} = $tEnd;
								$outnode->{"startend"}->{0}->{2} = $tEnd - $start +1;
							}else{
								$outnode->{"block_count"} = 1;
								$outnode->{"startend"}->{0}->{0} = $RNA->{$tran_idx}->{"startend"}->{$lastexon}->{0};
								$outnode->{"startend"}->{0}->{1} = $tEnd;								
								$outnode->{"startend"}->{0}->{2} = $tEnd - $RNA->{$tran_idx}->{"startend"}->{$lastexon}->{0} +1;
								for($kk=1;$kk <$firstexon - $lastexon; $kk++){
									$outnode->{"block_count"}++;
									$outnode->{"startend"}->{$kk}->{0} = $RNA->{$tran_idx}->{"startend"}->{$lastexon+$kk}->{0};
									$outnode->{"startend"}->{$kk}->{1} = $RNA->{$tran_idx}->{"startend"}->{$lastexon+$kk}->{1};
									$outnode->{"startend"}->{$kk}->{2} = $RNA->{$tran_idx}->{"startend"}->{$lastexon+$kk}->{1} - $RNA->{$tran_idx}->{"startend"}->{$lastexon+$kk}->{0} + 1;
								}
								$outnode->{"block_count"}++;							
								$outnode->{"startend"}->{$kk}->{0} = $start;
								$outnode->{"startend"}->{$kk}->{1} = $RNA->{$tran_idx}->{"startend"}->{$firstexon}->{1};
								$outnode->{"startend"}->{$kk}->{2} = $RNA->{$tran_idx}->{"startend"}->{$firstexon}->{1} - $start + 1;
							}									
						}								
					}	
		}
	}
	return $outnode;
}

sub get_relativepositon_of_absolutepostion{
	my ($RNA, $RNAnum, $pairRNAid, $posit)=@_;
	my $tran_idx;
	my $exon_idx;
	my $length = 0;
	my $out_position=0;
	for($tran_idx=0; $tran_idx<$RNAnum; $tran_idx++){
		if($RNA->{$tran_idx}->{"transcript_numberId"} == $pairRNAid){
			if($RNA->{$tran_idx}->{"strand"} eq "+"){
				for($exon_idx = 0; $exon_idx < $RNA->{$tran_idx}->{"block_count"}; $exon_idx++){
					if($RNA->{$tran_idx}->{"startend"}->{$exon_idx}->{0} <= $posit && $RNA->{$tran_idx}->{"startend"}->{$exon_idx}->{1} >= $posit){
						$out_position = $RNA->{$tran_idx}->{"relative_startend"}->{$exon_idx}->{0} + ($posit - $RNA->{$tran_idx}->{"startend"}->{$exon_idx}->{0});
						last;
					}
				}
			}else{
				for($exon_idx = 0; $exon_idx < $RNA->{$tran_idx}->{"block_count"}; $exon_idx++){
					if($RNA->{$tran_idx}->{"startend"}->{$exon_idx}->{0} <= $posit && $RNA->{$tran_idx}->{"startend"}->{$exon_idx}->{1} >= $posit){
						$out_position = $RNA->{$tran_idx}->{"relative_startend"}->{$exon_idx}->{0} + ($RNA->{$tran_idx}->{"startend"}->{$exon_idx}->{1} - $posit);
						last;
					}						
				}
			}
		}
	}
	return $out_position;
}

sub get_absolutepositon_of_relativepostion{
	my ($pairRNAid, $sampos, $RNAnum, $RNA) = @_;
	my $tran_idx;
	my $exon_idx; 
	my $length=0;
	my $out_position=0;
	my $strandtemp;
	for($tran_idx=0; $tran_idx<$RNAnum; $tran_idx++){
		if($RNA->{$tran_idx}->{"transcript_numberId"} == $pairRNAid){
			$strandtemp = $RNA->{$tran_idx}->{"strand"};
			if($RNA->{$tran_idx}->{"strand"} eq "+"){
				for($exon_idx = 0; $exon_idx < $RNA->{$tran_idx}->{"block_count"}; $exon_idx++){
					$length += $RNA->{$tran_idx}->{"startend"}->{$exon_idx}->{2};
					if($length >= $sampos){
						$out_position = $RNA->{$tran_idx}->{"startend"}->{$exon_idx}->{0} + ($sampos - $RNA->{$tran_idx}->{"relative_startend"}->{$exon_idx}->{0});
						last;
					}
				}
			}else{
				for($exon_idx = 0; $exon_idx < $RNA->{$tran_idx}->{"block_count"}; $exon_idx++){
					$length +=  $RNA->{$tran_idx}->{"startend"}->{$exon_idx}->{2};
					if($length >= $sampos){
					$out_position = $RNA->{$tran_idx}->{"startend"}->{$exon_idx}->{0} + ($RNA->{$tran_idx}->{"relative_startend"}->{$exon_idx}->{1} - $sampos);
						last;
					}						
				}
			}
		}
	}
	return ($out_position, $strandtemp);
}

sub calculate_relativeposition_of_eachexon_simple_bed{
	my ($simple_bednum, $simple_bed)=@_;
	my $tran_idx;
	my $exon_idx;
	for($tran_idx=0; $tran_idx < $simple_bednum; $tran_idx++){
		$simple_bed->{$tran_idx}->{"relative_startend"}->{0}->{0} = 1;
		$simple_bed->{$tran_idx}->{"relative_startend"}->{0}->{1} = $simple_bed->{$tran_idx}->{"startend"}->{0}->{2};
		for($exon_idx=1;$exon_idx<$simple_bed->{$tran_idx}->{"block_count"};$exon_idx++){
			$simple_bed->{$tran_idx}->{"relative_startend"}->{$exon_idx}->{0}=$simple_bed->{$tran_idx}->{"relative_startend"}->{$exon_idx-1}->{1} + 1;
			$simple_bed->{$tran_idx}->{"relative_startend"}->{$exon_idx}->{1}=$simple_bed->{$tran_idx}->{"relative_startend"}->{$exon_idx}->{0} + $simple_bed->{$tran_idx}->{"startend"}->{$exon_idx}->{2} - 1;	
		}
	}
}	



sub load_simple_bed{
	my ($inputfile,$simple_bed)=@_;
	open(A,"$inputfile") or die "cannot open $inputfile!\n";
	my $simple_bednum=0;
	my $exon_total_number=0;
	my $rowid=0;
	while(<A>){
		chomp;
		my @a=split(/\t/);
		next if($a[0] != $chr);
		$simple_bed->{$simple_bednum}->{"start"}= $a[1];
		$simple_bed->{$simple_bednum}->{"end"}= $a[2];
		$simple_bed->{$simple_bednum}->{"exon_number_id"}= $a[5];		
		if($a[5]==1){
			if($rowid>0){
				$simple_bednum++;
			}
			$simple_bed->{$simple_bednum}->{"transcript_numberId"}= $a[4];
			$simple_bed->{$simple_bednum}->{"chromosome"}= $a[0];
			$simple_bed->{$simple_bednum}->{"strand"}= $a[3];
			$simple_bed->{$simple_bednum}->{"block_count"}= 1;
			$simple_bed->{$simple_bednum}->{"startend"}->{0}->{0}= $a[1];
			$simple_bed->{$simple_bednum}->{"startend"}->{0}->{1}= $a[2];
			$simple_bed->{$simple_bednum}->{"startend"}->{0}->{2}= $a[2] - $a[1] + 1;
			$simple_bed->{$simple_bednum}->{"relative_startend"}->{0}->{0}=-1;
			$simple_bed->{$simple_bednum}->{"relative_startend"}->{0}->{1}=-1;
			
			$exon_total_number=1;
		}else{
			$exon_total_number++;
			if($exon_total_number > $max_block_number){
				print "ERROR, exon_total_number>1000, $exon_total_number\n";
			}
			$simple_bed->{$simple_bednum}->{"block_count"}++;
			$simple_bed->{$simple_bednum}->{"startend"}->{$exon_total_number-1}->{0}= $a[1];
			$simple_bed->{$simple_bednum}->{"startend"}->{$exon_total_number-1}->{1}= $a[2];
			$simple_bed->{$simple_bednum}->{"startend"}->{$exon_total_number-1}->{2}= $a[2] - $a[1] + 1;
			$simple_bed->{$simple_bednum}->{"relative_startend"}->{$exon_total_number-1}->{0}=-1;
			$simple_bed->{$simple_bednum}->{"relative_startend"}->{$exon_total_number-1}->{1}=-1;			
		}
		$rowid++;
	}
	close(A);
	$simple_bednum++;
	return  ($simple_bednum, $simple_bed);
}

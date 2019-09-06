use strict;
use warnings;

my ($input_hitnode, $input_simple, $input_chromosome, $output_statistic);
my @a=@ARGV;
if(@a==8){
	for(my $i=0; $i < @a; $i+=2){
		if($a[$i] eq "-input_hitnode"){
			$input_hitnode=$a[$i+1];
		}
		if($a[$i] eq "-input_simple"){
			$input_simple=$a[$i+1];
		}
		if($a[$i] eq "-input_chromosome"){
			$input_chromosome=$a[$i+1];
		}
		if($a[$i] eq "-output_statistic"){
			$output_statistic=$a[$i+1];
		}
	}
}
else{
	print "perl statisticby_pairAndSnp.pl -input_hitnode test.readscoverjuncsnp -input_simple XXX  -input_chromosome  1  -output_statistic test.statistic.out\n";
	exit();
}

open(A,"$input_hitnode") or die "cannot open $input_hitnode!\n";
my %hitnode;
my $hitnodenum=0;

while(<A>){
	chomp;
	my @a = split(/\s+/);
	$hitnode{$hitnodenum}{"hittype"} = $a[0];
	$hitnode{$hitnodenum}{"circRNAid"} = $a[3];
	$hitnode{$hitnodenum}{"mRNAid"} = $a[4];
	$hitnode{$hitnodenum}{"snp_position"} = $a[5];
	$hitnode{$hitnodenum}{"snp_ref"} = $a[7];
	$hitnode{$hitnodenum}{"snp_alt"} = $a[8];
	$hitnode{$hitnodenum}{"snp_rs"} = $a[9];
	$hitnode{$hitnodenum}{"nucleotide"} = $a[10];
	$hitnode{$hitnodenum}{"junction"} = $a[17];	

	$hitnodenum++;
}
close(A);

my $match = 0;
my %pairandsnp;
my $pairandsnpnum=0;
for(my $i = 0; $i < $hitnodenum; $i++){
	$match =0;
	for(my $j = 0; $j <$pairandsnpnum; $j++){
	
		if($hitnode{$i}{"mRNAid"} == $pairandsnp{$j}{"mRNAid"} && $hitnode{$i}{"circRNAid"} == $pairandsnp{$j}{"circRNAid"} && $hitnode{$i}{"snp_position"} == $pairandsnp{$j}{"snp_pos"}){
			$match = 1;
			if($hitnode{$i}{"hittype"} eq "hitmRNA"){
				$pairandsnp{$j}{"mRNA_Acount"}++ if($hitnode{$i}{"nucleotide"} eq "A");
				$pairandsnp{$j}{"mRNA_Tcount"}++ if($hitnode{$i}{"nucleotide"} eq "T");
				$pairandsnp{$j}{"mRNA_Ccount"}++ if($hitnode{$i}{"nucleotide"} eq "C");
				$pairandsnp{$j}{"mRNA_Gcount"}++ if($hitnode{$i}{"nucleotide"} eq "G");
			}else{
				$pairandsnp{$j}{"circRNA_Acount"}++ if($hitnode{$i}{"nucleotide"} eq "A");
				$pairandsnp{$j}{"circRNA_Tcount"}++ if($hitnode{$i}{"nucleotide"} eq "T");
				$pairandsnp{$j}{"circRNA_Ccount"}++ if($hitnode{$i}{"nucleotide"} eq "C");
				$pairandsnp{$j}{"circRNA_Gcount"}++ if($hitnode{$i}{"nucleotide"} eq "G");				
			}
		}
	}
	if($match == 0 && $hitnode{$i}{"nucleotide"}=~/[ATCG]/){
		$pairandsnp{$pairandsnpnum}{"circRNAid"}=$hitnode{$i}{"circRNAid"};  			
		$pairandsnp{$pairandsnpnum}{"mRNAid"}=$hitnode{$i}{"mRNAid"}; 
		$pairandsnp{$pairandsnpnum}{"junction"}=$hitnode{$i}{"junction"}; 
		$pairandsnp{$pairandsnpnum}{"snp_pos"}=$hitnode{$i}{"snp_position"};  
		$pairandsnp{$pairandsnpnum}{"snp_rs"}=$hitnode{$i}{"snp_rs"};
		$pairandsnp{$pairandsnpnum}{"snp_ref"}=$hitnode{$i}{"snp_ref"};
		$pairandsnp{$pairandsnpnum}{"snp_alt"}=$hitnode{$i}{"snp_alt"};	

		$pairandsnp{$pairandsnpnum}{"mRNA_Acount"}=0;
		$pairandsnp{$pairandsnpnum}{"mRNA_Ccount"}=0;
		$pairandsnp{$pairandsnpnum}{"mRNA_Gcount"}=0;
		$pairandsnp{$pairandsnpnum}{"mRNA_Tcount"}=0;


		$pairandsnp{$pairandsnpnum}{"circRNA_Acount"}=0;
		$pairandsnp{$pairandsnpnum}{"circRNA_Ccount"}=0;
		$pairandsnp{$pairandsnpnum}{"circRNA_Gcount"}=0;
		$pairandsnp{$pairandsnpnum}{"circRNA_Tcount"}=0;

			
		if($hitnode{$i}{"hittype"} eq "hitmRNA"){
			$pairandsnp{$pairandsnpnum}{"mRNA_Acount"}++ if($hitnode{$i}{"nucleotide"} eq "A");
			$pairandsnp{$pairandsnpnum}{"mRNA_Tcount"}++ if($hitnode{$i}{"nucleotide"} eq "T");
			$pairandsnp{$pairandsnpnum}{"mRNA_Ccount"}++ if($hitnode{$i}{"nucleotide"} eq "C");
			$pairandsnp{$pairandsnpnum}{"mRNA_Gcount"}++ if($hitnode{$i}{"nucleotide"} eq "G");
		}else{
			$pairandsnp{$pairandsnpnum}{"circRNA_Acount"}++ if($hitnode{$i}{"nucleotide"} eq "A");
			$pairandsnp{$pairandsnpnum}{"circRNA_Tcount"}++ if($hitnode{$i}{"nucleotide"} eq "T");
			$pairandsnp{$pairandsnpnum}{"circRNA_Ccount"}++ if($hitnode{$i}{"nucleotide"} eq "C");
			$pairandsnp{$pairandsnpnum}{"circRNA_Gcount"}++ if($hitnode{$i}{"nucleotide"} eq "G");				
		}	
		$pairandsnpnum++;
	}	
}

open(A,">$output_statistic") or die "cannot open $output_statistic!\n";
for(my $x=0; $x < $pairandsnpnum; $x++){
	my $n1=0;
	my $n2=0;
	my $n3=0;
	my $n4=0;
	$n1 = $pairandsnp{$x}{"mRNA_Acount"} if($pairandsnp{$x}{"snp_ref"} eq "A");
	$n1 = $pairandsnp{$x}{"mRNA_Tcount"} if($pairandsnp{$x}{"snp_ref"} eq "T");
	$n1 = $pairandsnp{$x}{"mRNA_Ccount"} if($pairandsnp{$x}{"snp_ref"} eq "C");
	$n1 = $pairandsnp{$x}{"mRNA_Gcount"} if($pairandsnp{$x}{"snp_ref"} eq "G");
	
	$n2 = $pairandsnp{$x}{"mRNA_Acount"} if($pairandsnp{$x}{"snp_alt"} eq "A");
	$n2 = $pairandsnp{$x}{"mRNA_Tcount"} if($pairandsnp{$x}{"snp_alt"} eq "T");
	$n2 = $pairandsnp{$x}{"mRNA_Ccount"} if($pairandsnp{$x}{"snp_alt"} eq "C");
	$n2 = $pairandsnp{$x}{"mRNA_Gcount"} if($pairandsnp{$x}{"snp_alt"} eq "G");	
	
	$n3 = $pairandsnp{$x}{"circRNA_Acount"} if($pairandsnp{$x}{"snp_ref"} eq "A");
	$n3 = $pairandsnp{$x}{"circRNA_Tcount"} if($pairandsnp{$x}{"snp_ref"} eq "T");
	$n3 = $pairandsnp{$x}{"circRNA_Ccount"} if($pairandsnp{$x}{"snp_ref"} eq "C");
	$n3 = $pairandsnp{$x}{"circRNA_Gcount"} if($pairandsnp{$x}{"snp_ref"} eq "G");
	
	$n4 = $pairandsnp{$x}{"circRNA_Acount"} if($pairandsnp{$x}{"snp_alt"} eq "A");
	$n4 = $pairandsnp{$x}{"circRNA_Tcount"} if($pairandsnp{$x}{"snp_alt"} eq "T");
	$n4 = $pairandsnp{$x}{"circRNA_Ccount"} if($pairandsnp{$x}{"snp_alt"} eq "C");
	$n4 = $pairandsnp{$x}{"circRNA_Gcount"} if($pairandsnp{$x}{"snp_alt"} eq "G");	
	
	my $t1=$n1+$n2;
	my $t2=$n3+$n4;
	if($n1+$n2>0 && $n3+$n4>0 && $n1+$n3>0 && $n2+$n4>0){	
		print A "$input_simple\t$input_chromosome\t$pairandsnp{$x}{\"snp_pos\"}\t$pairandsnp{$x}{\"snp_ref\"}\t$pairandsnp{$x}{\"snp_alt\"}\t$pairandsnp{$x}{\"junction\"}\t$pairandsnp{$x}{\"mRNAid\"}\t$pairandsnp{$x}{\"circRNAid\"}\t$n1\t$n2\t$n3\t$n4\t$t1\t$t2\tVALID\n";
	}else{
		print A "$input_simple\t$input_chromosome\t$pairandsnp{$x}{\"snp_pos\"}\t$pairandsnp{$x}{\"snp_ref\"}\t$pairandsnp{$x}{\"snp_alt\"}\t$pairandsnp{$x}{\"junction\"}\t$pairandsnp{$x}{\"mRNAid\"}\t$pairandsnp{$x}{\"circRNAid\"}\t$n1\t$n2\t$n3\t$n4\t$t1\t$t2\tINVALID\n";
	}
}

close(A);





















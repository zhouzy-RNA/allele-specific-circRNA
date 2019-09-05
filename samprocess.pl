use strict;
use warnings;


my $input_sam;
my $output_reformed;
my $output_invalidreads;

my @a=@ARGV;
if(@a==6){
	for(my $i=0; $i < @a; $i+=2){
		if($a[$i] eq "-input_sam"){
			$input_sam=$a[$i+1];
		}
		if($a[$i] eq "-output_reformed"){
			$output_reformed=$a[$i+1];
		}
		if($a[$i] eq "-output_invalidreads"){
			$output_invalidreads=$a[$i+1];
		}
	}
}
else{
	print "perl samprocess.pl -input_sam  test.sam -output_reformed  test.sam.reformed -output_invalidreads test.messages\n";
	exit();
}

open(A,"$input_sam") or die "cannot open $input_sam!\n";
open(B,">$output_reformed") or die "cannot open $output_reformed!\n";
open(C,">$output_invalidreads") or die "cannot open $output_invalidreads!\n";	
while(<A>){
	chomp;
	my $row1=$_;
	next if($row1=~/^[@#]/);
	chomp(my $row2=<A>);		
#	next if($row2=~/^[@#]/);
	my @a=split(/\t/,$row1);
	my @b=split(/\t/,$row2);
	next if($a[2] eq '\*' || $b[2] eq '\*');
	if($a[0] ne $b[0]){
		print C "unpaired reads: $a[0]\t$b[0]\n";
		next;
	}
	if($a[4]<20 || $b[4]<20){
		print C "mapping quality: $a[0]\t$a[4]\t$b[4]\n";
		next;
	}		
	if($a[5]=~/[IDSHP=X]/ || $b[5]=~/[IDSHP=X]/){
		print C "invalid cigar: $a[0]\t$a[6]\t$b[6]\n";
		next;
	}		
	if(abs($a[8]) > 3*length($a[9]) || abs($b[8]) > 3*length($b[9]) || abs($a[8]) <= length($a[9])){
		print C "invalid length: $a[0]\t$a[8]\t$b[8]\n";
		next;
	}
		
	if($a[8]<0){	
		my @c=@a;
		@a=@b;
		@b=@c;
	}	
	my $read1part;
	my $read2part;
	my $merged_seq;
	my $end_position;
	if(length($a[9]) + length($b[9]) > $a[8]){			
		my $read_1_start = $a[8] - length($b[9]);
		my $read_2_end = length($a[9]) + length($b[9]) - $a[8];	
		my $read1_overlap=substr($a[9], $read_1_start, length($a[9]) - $read_1_start);		
		my $read2_overlap = substr($b[9], 0, $read_2_end);
		my $new_overlap="";
		$read1part = substr($a[9], 0, $read_1_start);
		if($read1_overlap eq $read2_overlap){			
			$merged_seq = $read1part.$b[9];
		}else{
			for(my $i=0; $i < length($read1_overlap); $i++){
				if(substr($read1_overlap, $i, 1) ne substr($read2_overlap, $i, 1)){
					$new_overlap .= "N";
				}else{
					$new_overlap .= substr($read1_overlap, $i, 1);
				}
			}
			$read2part = substr($b[9], $read_2_end, length($b[9]) - $read_2_end);
			$merged_seq = $read1part . $new_overlap . $read2part;
								
		}
	}elsif(length($a[9]) + length($b[9]) == $a[8]){
		$merged_seq = $a[9] . $b[9];
	}else{
		my $gaplength = $a[8] - (length($a[9]) + length($b[9]));
		for(my $i=0; $i < $gaplength; $i++){
			$a[9] .= "-";
		}
		$merged_seq = $a[9] . $b[9];
	}
	$end_position = $a[3] + length($merged_seq) - 1;
	print B "$a[0]\t$a[2]\t$a[3]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$merged_seq\t$end_position\n";
}	
close(A);
close(B);
close(C); 

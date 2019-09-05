use strict;
use warnings;

my $input_circRNA_seq;
my $output_linear_seq;

my @a=@ARGV;
if(@a==4){
	for(my $i=0; $i < @a; $i+=2){
		if($a[$i] eq "-input_circRNA_seq"){
			$input_circRNA_seq=$a[$i+1];
		}
		if($a[$i] eq "-output_linear_seq"){
			$output_linear_seq=$a[$i+1];
		}
	}
}
else{
	print "perl splice_linear_from_circRNA.pl  -input_circRNA_seq  test.transcript.fa -output_linear_seq  test.linear.transcript.fa\n";
	exit();
}

open INPUT, "<", $input_circRNA_seq or die "cannot open $input_circRNA_seq!\n";
open OUTPUT, ">", $output_linear_seq or die "cannot open $output_linear_seq!\n";
while(<INPUT>){
	chomp;
	print OUTPUT "$_\n";
	my $fa=<INPUT>;
	if($fa){
		chomp($fa);
	}else{
		next;
	}
	my $len=length($fa);
	if($len % 2 ==0){
		my $left = substr($fa, 0, $len/2);
		my $right = substr($fa, $len/2, $len/2);
		$fa = $right . $left;
	}else{
		my $left = substr($fa, 0, int($len/2));
		my $right = substr($fa, $len/2, int($len/2)+1);
		$fa = $right . $left;		
	}
	print OUTPUT "$fa\n";
}
close(INPUT);
close(OUTPUT);

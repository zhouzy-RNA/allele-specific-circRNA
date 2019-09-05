use strict;
use warnings;

my $input_sam;
my $output_validreads;
my $output_invalidreads;

my @a=@ARGV;
if(@a==6){
	for(my $i=0; $i < @a; $i+=2){
		if($a[$i] eq "-input_sam"){
			$input_sam=$a[$i+1];
		}
		if($a[$i] eq "-output_validreads"){
			$output_validreads=$a[$i+1];
		}
		if($a[$i] eq "-output_invalidreads"){
			$output_invalidreads=$a[$i+1];
		}
	}
}
else{
	print "perl sam_pairedend_isvalid.pl  -input_sam   test.NumberID.sorted.sam.temp_noheader -output_validreads   test.NumberID.sorted.validpairs.sam -output_invalidreads test.NumberID.sorted.invalidpairs\n";
	exit();
}

open INPUT, "<", $input_sam or die "cannot open $input_sam!\n";
open OUTPUT1, ">", $output_validreads or die "cannot open $output_validreads!\n";
open OUTPUT2, ">", $output_invalidreads or die "cannot open $output_invalidreads!\n";
while(<INPUT>){
	chomp;
	my $reads1=$_;
	my @a=split(/\t/,$reads1);
	my $reads2=<INPUT>;
	chomp $reads2;
	my @b=split(/\t/,$reads2);
	next if($a[2] eq "*" || $b[2] eq "*");
	if($a[0] eq $b[0]){
		for(my $i=0;$i<10;$i++){
			print OUTPUT1 "$a[$i]\t";
		}
		print OUTPUT1 "\n";
		for(my $i=0;$i<10;$i++){
			print OUTPUT1 "$b[$i]\t";
		}
		print OUTPUT1 "\n";
	}else{
		for(my $i=0;$i<10;$i++){
			print OUTPUT2 "$a[$i]\t";
		}
		print OUTPUT2 "\n";
		for(my $i=0;$i<10;$i++){
			print OUTPUT2 "$b[$i]\t";
		}
		print OUTPUT2 "\n";		
	}
}
close(INPUT);
close(OUTPUT1);
close(OUTPUT2);

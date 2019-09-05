use strict;
use warnings;

my $input_fasta;
my $output;

my @a=@ARGV;
if(@a==4){
	for(my $i=0; $i < @a; $i+=2){
		if($a[$i] eq "-input_fasta"){
			$input_fasta=$a[$i+1];
		}
		if($a[$i] eq "-output"){
			$output=$a[$i+1];
		}
	}
}
else{
	print "perl ConcatExons.pl -input_fasta output.fa.temp -output  test.transcript.fa \n";
	exit();
}

open INPUT, "<", $input_fasta or die "cannot open $input_fasta!\n";
my %seqhashtemp;
my %strand;
my %tranname;
my $trannum=0;
while(<INPUT>){
	chomp;
	my @a=split(/\t/);
	chomp($seqhashtemp{$a[5]}{$a[2]}=<INPUT>);
	$strand{$a[5]}=$a[4];
	if($a[6] == 1){
		$tranname{$trannum}=$a[5];
		$trannum++;
	}
	
}
close(INPUT);

open OUTPUT, ">", $output or die "cannot open $output!\n";
for(my $i=0; $i<$trannum; $i++){
	my $tempseq;
	if($strand{$tranname{$i}} eq "+"){
		foreach my $m(sort {$a<=>$b} keys %{$seqhashtemp{$tranname{$i}}}){
			$tempseq .= $seqhashtemp{$tranname{$i}}{$m};
		}
	}else{
		foreach my $m(sort {$b<=>$a} keys %{$seqhashtemp{$tranname{$i}}}){
			$tempseq .= $seqhashtemp{$tranname{$i}}{$m};		
		}
	}
	print OUTPUT ">$tranname{$i}\n$tempseq\n";
}
close(OUTPUT);


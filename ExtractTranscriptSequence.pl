use strict;
use warnings;


my $genome;
my $transcript;
my $output;

my @a=@ARGV;
if(@a==6){
	for(my $i=0; $i < @a; $i+=2){
		if($a[$i] eq "-genome"){
			$genome=$a[$i+1];
		}
		if($a[$i] eq "-transcript"){
			$transcript=$a[$i+1];
		}
		if($a[$i] eq "-output"){
			$output=$a[$i+1];
		}
	}
}
else{
	print "perl ExtractTranscriptSequence.pl -genome test.masked.fa.1row -transcript test.transcript.bed -output output.fa\n";
	exit();
}

open INPUT, "<", $genome or die "cannot open $genome!\n";
my %genome_hash;
while(<INPUT>){
	chomp;
	my $tempchrname=$1 if(/^>(.+?)$/);
	chomp(my $fa=<INPUT>);
	$genome_hash{$tempchrname}=$fa;
}
close(INPUT);

open INPUT, "<", $transcript or die "cannot open $transcript!\n";
open OUTPUT, ">", $output or die "cannot open $output!\n";
while(<INPUT>){
	chomp;
	my @a=split(/\t/);
	my $tempgenomelen=length($genome_hash{$a[0]});
	my $tempseq;
	if($a[2]>$tempgenomelen){
		print "ERROR, stop_pos > chr length, chr_lenth: $tempgenomelen, stop_pos: $a[2]\n";
	}else{
		$tempseq=substr($genome_hash{$a[0]},$a[1]-1,$a[2]-$a[1]+1);
	}
	if($a[3] eq "+"){
		print OUTPUT ">$a[4]@<$a[0]\t$_\n$tempseq\n";
	}else{
		my $tranlen=length($tempseq);
		my $reverseseq;
		for(my $i = $tranlen - 1; $i >= 0; $i--){
			my $base=substr($tempseq,$i,1);
			if($base eq "A"){
				$base = "T";
			}elsif($base eq "T"){
				$base = "A";
			}elsif($base eq "C"){
				$base = "G";
			}elsif($base eq "G"){
				$base = "C";
			}elsif($base eq "N"){
				$base = "N";
			}else{
				print "ERROR: INVALIDE CHARACTER  $base\n";
			}
			$reverseseq .= $base;
		}
		print OUTPUT ">$a[4]@<$a[0]\t$_\n$reverseseq\n";
	}
}
close(INPUT);
close(OUTPUT);




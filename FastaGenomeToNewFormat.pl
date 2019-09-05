use strict;
use warnings;


my $origin_fasta;
my $seq1row_fasta;


my @a=@ARGV;
if(@a==4){
	for(my $i=0; $i < @a; $i+=2){
		if($a[$i] eq "-origin_fasta"){
			$origin_fasta=$a[$i+1];
		}
		if($a[$i] eq "-seq1row_fasta"){
			$seq1row_fasta=$a[$i+1];
		}
	}
}
else{
	print "perl FastaGenomeToNewFormat.pl -origin_fasta xxx.fa -seq1row_fasta xxx.fa.1row\n";
	exit();
}
open INPUT, "<", $origin_fasta or die "cannot open $origin_fasta!\n";
open OUTPUT, ">", $seq1row_fasta or die "cannot open $seq1row_fasta!\n";
my $rowid=0;
while(<INPUT>){
	chomp;
	if(/^>/){
		my @a=split(/\s+/);
		if($rowid==0){
			$rowid=1;
			print OUTPUT "$a[0]\n";
		}else{
			print OUTPUT "\n$a[0]\n";
		}
	}else{
		print OUTPUT "$_";
	}
}
close(INPUT);
close(OUTPUT);

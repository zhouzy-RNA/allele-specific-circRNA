use strict;
use warnings;

my $input_genome;
my $input_snp;
my $output_maskedgenome;
my @a=@ARGV;
if(@a==6){
	for(my $i=0; $i < @a; $i+=2){
		if($a[$i] eq "-input_genome"){
			$input_genome=$a[$i+1];
		}
		if($a[$i] eq "-input_snp"){
			$input_snp=$a[$i+1];
		}
		if($a[$i] eq "-output_maskedgenome"){
			$output_maskedgenome=$a[$i+1];
		}
	}
}
else{
	print "perl GenomeSequenceMaskN.pl -input_genome xxx.fa.1row -input_snp xxx.vcf -output_maskedgenome xxx.masked.fa.1row\n";
	exit();
}
my %vcf_hash;
open INPUTVCF, "<", $input_snp or die "cannot open $input_snp!\n";
while(<INPUTVCF>){
	chomp;
	next if(/^#/);
	my @a=split(/[(\s+)(\t)]/);
	$vcf_hash{$a[0]}{$a[1]}++;
}
close(INPUTVCF);

my %genome_hash;
my %chrname_list;
my $chrnum=1;
open INPUTFA, "<", $input_genome or die "cannot open $input_genome!\n";
open OUTPUT, ">", $output_maskedgenome or die "cannot open $output_maskedgenome!\n";
while(<INPUTFA>){
	chomp;
	print OUTPUT "$_\n";
	my $tempchrname=$1 if(/^>(.+?)$/);
	chomp(my $fa=<INPUTFA>);
	my $len=length($fa);
	my $temp;
	for(my $i=0; $i<$len; $i++){
		if($vcf_hash{$tempchrname}{$i+1}){
			$temp .= "N";
		}else{
			$temp .= substr($fa,$i,1);
		}
	}
	print OUTPUT "$temp\n";
}
close(INPUTFA);
close(OUTPUT);

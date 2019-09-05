# allele-specific-circRNA

This repository is used to authenticate allele specific ciriRNA in organism.

the total ten scripts are related. The latter script may use the output of the foregoing, when you use this repository, be care for the order.

usage

step 1. FastaGenomeToNewFormat.pl    ### transform each chromsome of the reference genome to one row

perl FastaGenomeToNewFormat.pl -origin_fasta xxx.fa -seq1row_fasta

	-origin_fasta  ###   the reference genome
	-seq1row_fasta   ### output file
    
  
step 2. GenomeSequenceMaskN.pl   ### this step need the vcf files you called or downlod from the website to masked the snp with "N". 

perl GenomeSequenceMaskN.pl -input_genome xxx.fa.1row -input_snp xxx.vcf -output_maskedgenome xxx.masked.fa.1row

	-input_genome    ### the fisrt step output
	-input_snp     ### vcf file recording chromosome number, snp position, rs number, ref snp, alt snp, like "1       629241  rs10458597      C       T"
	-output_maskedgenome	### output file name
 
 

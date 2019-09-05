# allele-specific-circRNA

This repository is used to authenticate allele specific ciriRNA in organism.

the total ten scripts are related. The latter script may use the output of the foregoing, when you use this repository, be care for the order.

usage

step one. FastaGenomeToNewFormat.pl    ### transform each chromsome of the reference genome to one row

perl FastaGenomeToNewFormat.pl -origin_fasta xxx.fa -seq1row_fasta xxx.fa.1row

	-origin_fasta  ###   the reference genome
	-seq1row_fasta   ### output file
    
  
step two. GenomeSequenceMaskN.pl   ### this step need the vcf files you called or downlod from the website to masked the snp with "N". 

perl GenomeSequenceMaskN.pl -input_genome xxx.fa.1row -input_snp xxx.vcf -output_maskedgenome xxx.masked.fa.1row

	-input_genome    ### the step one output
	-input_snp     ### this vcf file recording chromosome number, snp position, rs number, ref snp, alt snp.
	-output_maskedgenome	### output file name


step three. circRNA_mRNA_pairs.pl	### utilize the gtf file and the output of CIRI2.pl, which is used to  authenticate ciriRNA, to generate the gtf file of circRNA and the file about the infpmation of circRNA and the most length mRNA piared with it. 
 
perl circRNA_mRNA_pairs.pl -input_linear_gtf  xxx.protein_coding.exons.gtf -input_circRNAs   xxx.circ -output_circRNA_gtf circRNAs.gtf -output_circRNAmRNA_pair  circRNAmRNA_pair -minus_small2large $minus_small2large

	-input_linear_gtf	### only contain the line that contain "exon" and "protein_coding"
	-input_circRNAs		### this file contains circRNA_ID, chr, circRNA_start, circRNA_end, circRNA_type, gene_id and strand.
	-output_circRNA_gtf	### this file is simalar to the input_linear_gtf
	-output_circRNAmRNA_pair	### this file contains circRNA name, mRNA name, junct_start, junct_end, strand, mRNA_start, mRNA_end, mRNA_length, circRNA_start, circRNA_end, circRNA_length, chromosome number.
	-minus_small2large 	### 0 or 1, when minus strand is large to samll, like gtf for ensembl, minus_small2large is 1, otherwise minus_small2large is 0
	
step four. ExtractTranscriptSequence.pl ConcatExons.pl splice_linear_from_circRNA.pl	###






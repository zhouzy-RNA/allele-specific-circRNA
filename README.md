# allele-specific-circRNA

### This repository is used to authenticate allele specific ciriRNA in organism.

### the total ten scripts are related. The latter script may use the output of the foregoing, when you use this repository, be care for the order.

##usage

###step one. FastaGenomeToNewFormat.pl    ### transform each chromsome of the reference genome to one row

perl FastaGenomeToNewFormat.pl -origin_fasta xxx.fa -seq1row_fasta xxx.fa.1row

	-origin_fasta  ###   the reference genome
	-seq1row_fasta   ### output file
    
  
step two. GenomeSequenceMaskN.pl   ### this step need the vcf files you called or downlod from the website to masked the snp with "N". 

# mask snp
perl GenomeSequenceMaskN.pl -input_genome xxx.fa.1row -input_snp xxx.vcf -output_maskedgenome xxx.masked.fa.1row

	-input_genome    ### the step one output
	-input_snp     ### this vcf file recording chromosome number, snp position, rs number, ref snp, alt snp.
	-output_maskedgenome	### output file name

## the precondition of step two

# split chromosome
for((i=1;i<=${chr_num};i++));
do
	sed -n "$[i*2-1],$[i*2]p" xxx.fa.1row >xxx.chrx.fa.1row;
done

# extract the overlop of the reference vcf and the vcf you called and filtered

vcf-isec -a -f -n +2 xxx.vcf.gz xxx.new.vcf.gz chrx_overlapping.vcf







step three. circRNA_mRNA_pairs.pl	### utilize the gtf file and the output of CIRI2.pl, which is used to  authenticate ciriRNA, to generate the gtf file of circRNA and the file about the information of circRNA and the most length mRNA piared with it. 
 
perl circRNA_mRNA_pairs.pl -input_linear_gtf  xxx.protein_coding.exons.gtf -input_circRNAs   xxx.circ -output_circRNA_gtf circRNAs.gtf -output_circRNAmRNA_pair  circRNAmRNA_pair -minus_small2large $minus_small2large

	-input_linear_gtf	### only contain the line that contain "exon" and "protein_coding"
	-input_circRNAs		### this file contains circRNA_ID, chr, circRNA_start, circRNA_end, circRNA_type, gene_id and strand.
	-output_circRNA_gtf	### this file is simalar to the input_linear_gtf
	-output_circRNAmRNA_pair	### this file contains circRNA name, mRNA name, junct_start, junct_end, strand, mRNA_start, mRNA_end, mRNA_length, circRNA_start, circRNA_end, circRNA_length, chromosome number.
	-minus_small2large 	### 0 or 1, when minus strand of the gtf is large to samll, like gtf for ensembl, minus_small2large is 1, otherwise minus_small2large is 0
	
	
step four. ExtractTranscriptSequence.pl ConcatExons.pl splice_linear_from_circRNA.pl	### this three scripts are used to  generate the transcript fasta according to the masked genome. Run the third sciript only if the transcript is circRNA. Before this you need to rename the transcript according to the chromosome and the order of the transcript position, both the mRNA and circRNA.

perl ExtractTranscriptSequence.pl -genome xxx.masked.fa.1row -transcript xxx.NewNumberID.bed.sorted.txt -output xxx.transcript.temp

	-genome 	### masked genome according to step one
	-transcript 	### this circRNA's file can generate from the gtf of circRNA, the mRNA's file can generate from the gtf and the circRNAmRNA_pair, containing chromosome number, RNA_start, RNA_end, strand, RNA name, exon number, like "1 100906852 100906965 + 1000001 1"
	-output xxx.transcript.temp	### output
	
perl ConcatExons.pl -input_fasta xxx.transcript.temp -output  xxx.NewNumberID.transcript.fa 

	-input_fasta 	### the output of the ExtractTranscriptSequence.pl 
	-output  	### output
	
perl splice_linear_from_circRNA.pl  -input_circRNA_seq  xxx.NewNumberID.transcript.fa -output_linear_seq  xxx.linear.NewNumberID.transcript.fa

	-input_circRNA_seq	### the output of splice_linear_from_circRNA.pl
	-output_linear_seq  	### output


step five. sam_pairedend_isvalid.pl	### extract the paired reads from the bam generated by the raw data and the transcript fasta.

perl sam_pairedend_isvalid.pl  -input_sam   xxx.NumberID.sorted.sam.temp_noheader -output_validreads   xxx.NumberID.sorted.validpairs.sam -output_invalidreads xxx.NumberID.sorted.invalidpairs


	-input_sam   ### sorted sam
	-output_validreads   ### valid reads
	-output_invalidreads 	### invalid reads


step six. samprocess.pl		### merge the paired reads to form one reads

perl samprocess.pl -input_sam xxx.chr.x.sam -output_reformed xxx.chr.x.sam.reformed -output_invalidreads xxx.chr.x.sam.invalidreads

	-input_sam 	### the output of step five split according to chromosome
	-output_reformed 	### the merged reads
	-output_invalidreads 	### output


step seven. statisticscomplex.pl	### generate the information of snp in reads about its position to genome (absolute pos) or gene (relative pos) and the back splicing site position

perl statisticscomplex.pl -input_pair circRNAmRNA_pair.numberid.txt -input_sam xxx.chr.x.sam.reformed -input_snp xxx.vcf -input_mRNA_simple_bed kept_longest_mRNAs.1_x.NewNumberID.bed -output_reformed xxx.chr.x.readscoverjuncsnp -chr x

	-input_pair	### the one of output files of step three (circRNAmRNA_pair). (renamed the transcript)
	-input_sam 	### the output of the step six. 
	-input_snp  	### this vcf is same as the second step 
	-input_mRNA_simple_bed 	### the mRNA bed is same as the forth step
	-output_reformed 	### this output contains hitmRNA or hitcircRNA, hitRNA name, reads name, circRNA id, mRNA id, snp position, snp relative position, ref anp, slt snp, re number, current bp,   junction relative, read absolute start, read absolute end, read start relative, read end relative, strand, cover junction absolute, merged read
	-chr 		### chromosome number


step eight. statisticby_pairAndSnp.pl	### the statistic result of ref or alt and the number of linear RNA and circRNA of 
each.

perl statisticby_pairAndSnp.pl -input_hitnode xxx.chr.x.readscoverjuncsnp -input_sample xxx -input_chromosome x -output_statistic xxx.x.statistic

		-input_hitnode		### the output of step eight
		-input_simple		### simple name
		-input_chromosome 	### chromosome number
		-output_statistic 	### output contains simple name, chromosome number, snp pos, ref, alt, junction pos, mRNA id, circRNA id, linear mRNA ref number, linear mRNA alt number, circRNA ref number, circRNA alt number, valid or invalid

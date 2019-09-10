# allele-specific-circRNA

#### This repository is used to authenticate allele specific ciriRNA in organism.

#### the total ten scripts are related. The latter script may use the output of the foregoing, when you use this repository, be care for the order.


## usage

#### note: XXX represents file name. X in chrx or the sole x represent the numnber of chromosome (1 to 22). Before you run the script, please read the precondition of each step.

### step one. FastaGenomeToNewFormat.pl

#### transform each chromsome of the reference genome to one row

	`perl FastaGenomeToNewFormat.pl -origin_fasta xxx.fa -seq1row_fasta xxx.fa.1row

		-origin_fasta  ###   the reference genome
		-seq1row_fasta   ### output file`
    
  
### step two. GenomeSequenceMaskN.pl   

#### this step need the vcf files you called or downlod from the website to masked the snp with "N". 

#### the precondition of step two
	
##### sort the reference vcf

	vcf-sort -c xxx.vcf >xxx.sort.vcf
	
##### convert the sorted reference vcf

	cat xxx.sort.vcf | vcf-convert -v 4.1 -r xxx.fa >xxx.new.vcf
	
##### compress the comverted vcf

	bgzip -c xxx.new.vcf >xxx.new.vcf.gz
	
##### index

	tabix -p vcf xxx.new.vcf.gz
	
##### extract the overlop of the compressed reference vcf and the vcf you called by chromose, which also need to compress and index

	vcf-isec -a -f -n +2 hard_filter_chrx.vcf.gz xxx.new.vcf.gz >chrx_overlapping.vcf

##### genrate the required vcf

	grep -v "^#"  chrx_rlapping.vcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5;}' | awk '{if(length($4)==1 && length($5)==1) {print $0;}}' > chrx_overlapping.chr_posit.vcf
	

	perl GenomeSequenceMaskN.pl -input_genome xxx.fa.1row -input_snp chrx_overlapping.chr_posit.vcf -output_maskedgenome xxx.masked.fa.1row

		-input_genome    ### the step one output
		-input_snp     ### this vcf file recording chromosome number, snp position, rs number, ref snp, alt snp.
		-output_maskedgenome	### output file name


### step three. circRNA_mRNA_pairs.pl	

#### utilize the gtf file and the output of CIRI2.pl, which is used to  authenticate ciriRNA, to generate the gtf file of circRNA and the file about the information of circRNA and the most length mRNA piared with it. 

#### the precondition of step three

##### generate the required reference gtf

	sed -n '6,$p' xxx.f | awk '$3 == "exon"{print $0}' |   \
	sed 's/transcript_version.*exon_number/exon_number/g' |  \
 	sed 's/gene_version.*transcript_id/transcript_id/g' | grep "protein_coding" \
 	> xxx.protein_coding.exons.gtf

##### generate the required input_circRNAs
	
	cat human.ciri2.xxx.out >> merged.all.circ.temp		### human.ciri2.xxx.out is the output of CIRI2.pl
	
	awk 'BEGIN { OFS="\t"} { print $1"\t"$2"\t"$3"\t"$4"\t"$9"\t"$10"\t"$11;  }' merged.all.circ.temp \
 	| awk '$5 == "exon" {print $0;}' | sed 's/"\t"//g' | sed 's/:/_/g'| sed 's/|/_/g' \
	 | awk '{if ($7 == "+") {print $1"_plus""\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;} \
	 else {print $1"_minus""\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;}}' \
	 > merged.all.circ.temp2

	cat ${output_file_path}/output/merged.all.circ.temp2 | LC_ALL=C  sort -k1,1  | \
	uniq > merged.all.circ.temp3

	awk ' $2 ~ /^[123456789]/' ${output_file_path}/output/merged.all.circ.temp3 \
	> merged.all.circ.temp4
 
 
	perl circRNA_mRNA_pairs.pl -input_linear_gtf  xxx.protein_coding.exons.gtf -input_circRNAs   merged.all.circ.temp4 -output_circRNA_gtf circRNAs.gtf -output_circRNAmRNA_pair  circRNAmRNA_pair -minus_small2large $minus_small2large

		-input_linear_gtf	### only contain the line that contain "exon" and "protein_coding"
		-input_circRNAs		### this file contains circRNA_ID, chr, circRNA_start, circRNA_end, circRNA_type, gene_id and strand.
		-output_circRNA_gtf	### this file is simalar to the input_linear_gtf
		-output_circRNAmRNA_pair	### this file contains circRNA name, mRNA name, junct_start, junct_end, strand, mRNA_start, mRNA_end, mRNA_length, circRNA_start, circRNA_end, circRNA_length, chromosome number.
		-minus_small2large 	### 0 or 1, when minus strand of the gtf is large to samll, like gtf for ensembl, minus_small2large is 1, otherwise minus_small2large is 0

	
### step four. ExtractTranscriptSequence.pl ConcatExons.pl splice_linear_from_circRNA.pl

#### this three scripts are used to  generate the transcript fasta according to the masked genome. Run the third sciript only if the transcript is circRNA. 

#### the precondition of step four

##### split the masked reference genome

	for((i=1;i<=${chr_num};i++));
	do
		sed -n "$[i*2-1],$[i*2]p" xxx.masked.fa.1row >xxx.chrx.masked.fa.1row;
	done
##### generate the circRNA's bed

	grep -v "^#"  ${output_file_path}/output/circRNAs.gtf | \
 	awk 'BEGIN {FS="\t";OFS="\t";} $3 == "exon" {print $1"\t"$4"\t"$5"\t"$7"\t"$9;}'  \
	| sed 's/ //g'  \
	| sed 's/gene_id.*ranscript_id//g' \
	| sed 's/transcript_version.*exon_number//g' \
	| sed 's/;gene_name.*//g' \
	| sed 's/\"//g' \
	| sed 's/;/\t/g' | sed 's/exon_number//g' \
	> circRNAs.bed
	
	awk 'BEGIN { OFS="\t"}  $1 ~ /^[123456789]/{print}' circRNAs.bed > circRNAs.1_22.bed
	
##### generate the mRNA's bed

	awk '{print $2;}' circRNAmRNA_pair | sort | uniq > final.mRNA.IDs
	
	awk 'BEGIN { OFS="\t"} NR==FNR{a[$1]++};NR>FNR&&a[$5]{print}' final.mRNA.IDs \
	xxx.protein_coding.exons.bed >kept_longest_mRNAs.1_22.bed
	
##### generate the new transcript ID
	
	cat ${output_file_path}/output/kept_longest_mRNAs.1_22.bed \
	circRNAs.1_22.bed  \
	| awk '{ print $1"\t"$5; }' | LANG=C sort | uniq > merged.circ.mRNA.temp001

	awk 'BEGIN { OFS = "\t" } {if ($1 ~ /^[123456789]/)  print  $1*1000000+NR"\t"$1"\t"$2;}' merged.circ.mRNA.temp001 \
	> merged.circ.mRNA.new.transcriptID

##### generate new transcript ID bed

	awk 'BEGIN { OFS="\t"} NR==FNR{a[$3]=$1};NR>FNR{$5=a[$5];print}' \
	merged.circ.mRNA.new.transcriptID kept_longest_mRNAs.1_22.bed >kept_longest_mRNAs.1_22.NewNumberID.bed


	awk 'BEGIN { OFS="\t"} NR==FNR{a[$3]=$1};NR>FNR{$5=a[$5];print}' \
	merged.circ.mRNA.new.transcriptID circRNAs.1_22bed >circRNAs.1_22.NewNumberID.bed

	sort -k5,5n -k6,6n circRNAs.1_22.NewNumberID.bed >circRNAs.1_22.NewNumberID.bed.sorted.txt
	
	mv circRNAs.1_22.NewNumberID.bed.sorted.txt circRNAs.1_22.NewNumberID.bed
	

	perl ExtractTranscriptSequence.pl -genome xxx.chrx.masked.fa.1row -transcript xxx.NewNumberID.bed -output xxx.transcript.temp

		-genome 	### generate from the precondition of step four
		-transcript 	### this circRNA's file can generate from the gtf of circRNA, the mRNA's file can generate from the gtf and the circRNAmRNA_pair, containing chromosome number, RNA_start, RNA_end, strand, RNA name, exon number, like "1 100906852 100906965 + 1000001 1"
		-output xxx.transcript.temp	### output
	
	perl ConcatExons.pl -input_fasta xxx.transcript.temp -output  xxx.NewNumberID.transcript.fa 

		-input_fasta 	### the output of the ExtractTranscriptSequence.pl 
		-output  	### output
	
	perl splice_linear_from_circRNA.pl  -input_circRNA_seq  xxx.NewNumberID.transcript.fa -output_linear_seq  circRNA.linear.NewNumberID.transcript.fa

		-input_circRNA_seq	### the output of splice_linear_from_circRNA.pl
		-output_linear_seq  	### output
	

### step five. sam_pairedend_isvalid.pl	

#### extract the paired reads from the bam generated by the raw data and the transcript fasta.

#### the precondition of step five

##### generate the reference transcript fasta
	
	cat circRNAs.linear.NewNumberID.transcript.fa longest.mRNAs.NewNumberID.transcript.fa \
	> merged.circRNAs.protein_coding.NewNumberID.fa	
	
##### generated bam according to the raw data and the transcript fasta.
	
	${bowtie2_path}/bowtie2 -p 4 -q -x merged.circRNAs.protein_coding.NewNumberID  \
  	-1 xxx_1.fastq \
  	-2 xxx_2.fastq \
 	| samtools view -bS - > bowtie2_xxx.NewNumberID.bam
  
	samtools sort -n -m 10000000000 bowtie2_xxx.NewNumberID.bam bowtie2_xxx.NumberID.sorted 

	samtools view -h -o bowtie2_xxx.NumberID.sorted.sam bowtie2_xxx.NumberID.sorted.bam

	grep -v "^@" bowtie2_xxx.NumberID.sorted.sam  > bowtie2_xxx.NumberID.sorted.sam.temp_noheader
	

	perl sam_pairedend_isvalid.pl  -input_sam   xxx.NumberID.sorted.sam.temp_noheader -output_validreads   xxx.NumberID.sorted.validpairs.sam -output_invalidreads xxx.NumberID.sorted.invalidpairs


		-input_sam   ### sam file
		-output_validreads   ### valid reads
		-output_invalidreads 	### invalid reads


### step six. samprocess.pl	

#### merge the paired reads to form one reads

#### the precondition of step six

#####
	for i in  $(seq 1 22)
	do 
		awk ' $3 >= '${i}'000000 && $3 <= '${i}'999999 { print $0; } '  \
		bowtie2.xxx.NumberID.sorted.validpairs.sam > bowtie2.xxx.chr.x.sam
	done
	

	perl samprocess.pl -input_sam xxx.chr.x.sam -output_reformed xxx.chr.x.sam.reformed -output_invalidreads xxx.chr.x.sam.invalidreads

		-input_sam 	### the sam splited according to chromosome
		-output_reformed 	### the merged reads
		-output_invalidreads 	### output


### step seven. statisticscomplex.pl	

#### generate the information of snp in reads about its position to genome (absolute pos) or gene (relative pos) and the back splicing site position

#### the precondition of step seven

##### generate the input_pair according to the 

	awk 'BEGIN { OFS="\t"} NR==FNR{a[$3]=$1};NR>FNR{$1=a[$1];$2=a[$2];print}' \
	merged.circ.mRNA.new.transcriptID circRNAmRNA_pair >circRNAmRNA_pair.numberid.txt


### step eight. statisticby_pairAndSnp.pl	

#### the statistic result of ref or alt and the number of linear RNA and circRNA of each.
	
	perl statisticby_pairAndSnp.pl -input_hitnode xxx.chr.x.readscoverjuncsnp -input_sample xxx -input_chromosome x -output_statistic xxx.x.statistic

		-input_hitnode		### the output of step eight
		-input_simple		### simple name
		-input_chromosome 	### chromosome number
		-output_statistic 	### output contains simple name, chromosome number, snp pos, ref, alt, junction pos, mRNA id, circRNA id, linear mRNA ref number, linear mRNA alt number, circRNA ref number, circRNA alt number, valid or invalid
		

	perl statisticscomplex.pl -input_pair circRNAmRNA_pair.numberid.txt -input_sam xxx.chr.x.sam.reformed -input_snp xxx.vcf -input_mRNA_simple_bed kept_longest_mRNAs.1_x.NewNumberID.bed -output_reformed xxx.chr.x.readscoverjuncsnp -chr x

		-input_pair	### the one of output files of step three (circRNAmRNA_pair). (renamed the transcript)
		-input_sam 	### the output of the step six. 
		-input_snp  	### this vcf is same as the second step 
		-input_mRNA_simple_bed 	### the mRNA bed is same as the forth step
		-output_reformed 	### this output contains hitmRNA or hitcircRNA, hitRNA name, reads name, circRNA id, mRNA id, snp position, snp relative position, ref anp, slt snp, re number, current bp,   junction relative, read absolute start, read absolute end, read start relative, read end relative, strand, cover junction absolute, merged read
		-chr 		### chromosome number
		

### step nine. oddsratio.R

#### generate the valid oddsrate of ecah snp

#### the precondition of step seven

##### merge all statistic file of each simple
	
	cat xxx.*.statistic > xxx.statistic
	
	
	Rscript ${perl_script_path}/oddsratio.R xxx output
	
	grep -v "sample" xxx.statistics.oddsratio.final | awk '$13>1 && $14>1 && $15>1{print $2"_"$3"_"$6"_"$7"_"$8}'  \
	| LANG=C sort  | awk -F, '{arr[$1]++}END{for (a in arr) print a, arr[a]}' \
	| sort  -k1,1n -k2,2n >xxx.valid.oddsratio.final
	

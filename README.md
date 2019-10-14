#Tutorial of a Suite to Identify Allele Specific circRNAs

#### This repository can be used to identify allele specific circRNAs in organisms.
#### This repository consists of ten script files. We recommend users to use these script files consequently following our guide.


## Usage
## Note that ‘xxx’ represents a file name. ‘X’ in ‘chrx’ or a singular ‘x’ stands for the number of chromosome (1 to 22).

###===============================
### Step 1. 
### Transform each chromosome sequence of the reference genome into a row.
	perl FastaGenomeToNewFormat.pl -origin_fasta xxx.fa -seq1row_fasta xxx.fa.1row
		-origin_fasta   the reference genome.
		-seq1row_fasta  the output file.  
  
###===============================
### Step 2. 
### Use the VCF files you called as input and mask SNPs as "N". 
	
### Step 2.1 Sort the reference vcf.
	vcf-sort -c xxx.vcf >xxx.sort.vcf                                            
	
### Step 2.2 Convert the sorted reference vcf.
	cat xxx.sort.vcf | vcf-convert -v 4.1 -r xxx.fa >xxx.new.vcf                      
	
### Step 2.3 Compress the converted vcf.
	bgzip -c xxx.new.vcf >xxx.new.vcf.gz                                       
	
### Step 2.4 Index vcf.
	tabix -p vcf xxx.new.vcf.gz                                                
	
### Step 2.5 Extract the overlap of the compressed reference vcf and the vcf you called by chromosome, which also need to be compressed and indexed.
	vcf-isec -a -f -n +2 hard_filter_chrx.vcf.gz xxx.new.vcf.gz >chrx_overlapping.vcf

### Step 2.6 Generate the required vcf.
	grep -v "^#"  chrx_rlapping.vcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5;}' | awk '{if(length($4)==1 && length($5)==1) {print $0;}}' > chrx_overlapping.chr_posit.vcf
	
 ### Step 2.7 Mask SNPs as ‘N’.
	perl GenomeSequenceMaskN.pl -input_genome xxx.fa.1row -input_snp chrx_overlapping.chr_posit.vcf -output_maskedgenome xxx.masked.fa.1row
	-input_genome    the output of step 1.
	-input_snp     this vcf file includes chromosome number, snp position, rs number, ref snp, and alt snp.
	-output_maskedgenome	the output file name.

###===============================
### Step 3. 
### Utilize the gtf file and the output of CIRI2.pl, which is a tool to identify circRNAs, to generate the gtf file of circRNAs and the file about the information of circRNAs and the longest mRNAs matched with corresponding circRNAs. 

### Step 3.1 Generate the required reference gtf.
	sed -n '6,$p' xxx.f | awk '$3 == "exon"{print $0}' |   \
	sed 's/transcript_version.*exon_number/exon_number/g' |  \
 	sed 's/gene_version.*transcript_id/transcript_id/g' | grep "protein_coding" \
 	> xxx.protein_coding.exons.gtf

### Step 3.2 Generate the required input circRNAs.	
	cat human.ciri2.xxx.out >> merged.all.circ.temp 
# human.ciri2.xxx.out is the output of CIRI2.pl
	
	awk 'BEGIN { OFS="\t"} { print $1"\t"$2"\t"$3"\t"$4"\t"$9"\t"$10"\t"$11;  }' merged.all.circ.temp \
 	| awk '$5 == "exon" {print $0;}' | sed 's/"\t"//g' | sed 's/:/_/g'| sed 's/|/_/g' \
	 | awk '{if ($7 == "+") {print $1"_plus""\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;} \
	 else {print $1"_minus""\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;}}' \
	 > merged.all.circ.temp2

	cat ${output_file_path}/output/merged.all.circ.temp2 | LC_ALL=C  sort -k1,1  | \
	uniq > merged.all.circ.temp3

	awk ' $2 ~ /^[123456789]/' ${output_file_path}/output/merged.all.circ.temp3 \
	> merged.all.circ.temp4
 
### Step 3.3 Get matched mRNAs of circRNAs.  
	perl circRNA_mRNA_pairs.pl -input_linear_gtf  xxx.protein_coding.exons.gtf -input_circRNAs   merged.all.circ.temp4 -output_circRNA_gtf circRNAs.gtf -output_circRNAmRNA_pair  circRNAmRNA_pair -minus_small2large 0
	-input_linear_gtf	only contain the line that contain "exon" and "protein_coding"
	-input_circRNAs		this file contains circRNA_ID, chr, circRNA_start, circRNA_end, circRNA_type, gene_id and strand.
	-output_circRNA_gtf	this file is simalar to the input_linear_gtf
	-output_circRNAmRNA_pair	this file contains circRNA name, mRNA name, junct_start, junct_end, strand, mRNA_start, mRNA_end, mRNA_length, circRNA_start, circRNA_end, circRNA_length, chromosome number.
	-minus_small2large 	0 or 1, when minus strand of the gtf is large to samll, like gtf for ensembl, minus_small2large is 1, otherwise minus_small2large is 0
	
###===============================
### Step 4. 
### Generate the transcript fasta according to the masked genome. 

### Step 4.1 Split the masked reference genome. Please notice the order of chromosome
	for((i=1;i<=${chr_num};i++));
	do
		sed -n "$[i*2-1],$[i*2]p" xxx.masked.fa.1row >xxx.chrx.masked.fa.1row;
	done
	
### Step 4.2 Generate the circRNA's bed file.
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
	
### Step 4.3 Generate the mRNA's bed file.
	awk '{print $2;}' circRNAmRNA_pair | sort | uniq > final.mRNA.IDs
	
	awk 'BEGIN { OFS="\t"} NR==FNR{a[$1]++};NR>FNR&&a[$5]{print}' final.mRNA.IDs \
	xxx.protein_coding.exons.bed >kept_longest_mRNAs.1_22.bed
	
### Step 4.4 Generate the new transcript ID.	
	cat ${output_file_path}/output/kept_longest_mRNAs.1_22.bed \
	circRNAs.1_22.bed  \
	| awk '{ print $1"\t"$5; }' | LANG=C sort | uniq > merged.circ.mRNA.temp001

	awk 'BEGIN { OFS = "\t" } {if ($1 ~ /^[123456789]/)  print  $1*1000000+NR"\t"$1"\t"$2;}' merged.circ.mRNA.temp001 \
	> merged.circ.mRNA.new.transcriptID

### Step 4.5 Generate new transcript ID bed.
	awk 'BEGIN { OFS="\t"} NR==FNR{a[$3]=$1};NR>FNR{$5=a[$5];print}' \
	merged.circ.mRNA.new.transcriptID kept_longest_mRNAs.1_22.bed >kept_longest_mRNAs.1_22.NewNumberID.bed

	awk 'BEGIN { OFS="\t"} NR==FNR{a[$3]=$1};NR>FNR{$5=a[$5];print}' \
	merged.circ.mRNA.new.transcriptID circRNAs.1_22bed >circRNAs.1_22.NewNumberID.bed

	sort -k5,5n -k6,6n circRNAs.1_22.NewNumberID.bed >circRNAs.1_22.NewNumberID.bed.sorted.txt
	
	mv circRNAs.1_22.NewNumberID.bed.sorted.txt circRNAs.1_22.NewNumberID.bed
	
### Step 4.6 Generate the transcript fasta according to the masked genome.
	perl ExtractTranscriptSequence.pl -genome xxx.chrx.masked.fa.1row -transcript xxx.NewNumberID.bed -output xxx.transcript.temp
	-genome 	generated from Step 4.1.
	-transcript 	this circRNA's file can be obtained from the gtf of circRNAs, the mRNA's file can be obtained from the gtf and the circRNAmRNA_pair, containing chromosome number, RNA_start, RNA_end, strand, RNA name, exon number. i.e. "1 100906852 100906965 + 1000001 1".
	-output the output file.
	
	perl ConcatExons.pl -input_fasta xxx.transcript.temp -output  xxx.NewNumberID.transcript.fa 
	-input_fasta 	the output of the ExtractTranscriptSequence.pl script file.
	-output  	the output file.
	
	perl splice_linear_from_circRNA.pl  -input_circRNA_seq  xxx.NewNumberID.transcript.fa -output_linear_seq  circRNA.linear.NewNumberID.transcript.fa
	-input_circRNA_seq	the output of splice_linear_from_circRNA.pl script file.
	-output_linear_seq  	the output file.	

###===============================
### Step 5. 
### Extract the paired reads from the bam generated by the raw data and the transcript fasta.

### Step 5.1 Generate the reference transcript fasta.	
	cat circRNAs.linear.NewNumberID.transcript.fa longest.mRNAs.NewNumberID.transcript.fa \
	> merged.circRNAs.protein_coding.NewNumberID.fa	
	
### Step 5.2 Generated bam according to the raw data and the transcript fasta.	
	${bowtie2_path}/bowtie2 -p 4 -q -x merged.circRNAs.protein_coding.NewNumberID  \
  	-1 xxx_1.fastq \
  	-2 xxx_2.fastq \
 	| samtools view -bS - > bowtie2_xxx.NewNumberID.bam
  
	samtools sort -n -m 10000000000 bowtie2_xxx.NewNumberID.bam bowtie2_xxx.NumberID.sorted 

	samtools view -h -o bowtie2_xxx.NumberID.sorted.sam bowtie2_xxx.NumberID.sorted.bam

	grep -v "^@" bowtie2_xxx.NumberID.sorted.sam  > bowtie2_xxx.NumberID.sorted.sam.temp_noheader
	
### Step 5.3 Extract the paired reads from the bam generated by the raw data and the transcript fasta.
	perl sam_pairedend_isvalid.pl  -input_sam   xxx.NumberID.sorted.sam.temp_noheader -output_validreads   xxx.NumberID.sorted.validpairs.sam -output_invalidreads xxx.NumberID.sorted.invalidpairs
	-input_sam   the sam file.
	-output_validreads   the valid reads.
	-output_invalidreads 	the invalid reads.

###===============================
### Step 6. 
### Merge the paired reads to form one reads.

### Step 6.1 Filtering.
	for i in  $(seq 1 22)
	do 
		awk ' $3 >= '${i}'000000 && $3 <= '${i}'999999 { print $0; } '  \
		bowtie2.xxx.NumberID.sorted.validpairs.sam > bowtie2.xxx.chr.x.sam
	done
	
### Step 6.2 Merge the paired reads to form one reads.
	perl samprocess.pl -input_sam xxx.chr.x.sam -output_reformed xxx.chr.x.sam.reformed -output_invalidreads xxx.chr.x.sam.invalidreads
	-input_sam 	the sam splited according to chromosome.
	-output_reformed 	the merged reads.
	-output_invalidreads 	the output file.


###===============================
### Step 7. 
### Generate the information of snp in reads about its position to genome (absolute pos) or gene (relative pos) and the back splicing site position.

### Step 7.1 Generate the input pairs. 
	awk 'BEGIN { OFS="\t"} NR==FNR{a[$3]=$1};NR>FNR{$1=a[$1];$2=a[$2];print}' \
	merged.circ.mRNA.new.transcriptID circRNAmRNA_pair >circRNAmRNA_pair.numberid.txt
	
### Step 7.2 Generate the information of snp in reads about its position to genome (absolute pos) or gene (relative pos) and the back splicing site position. 
	perl statisticscomplex.pl -input_pair circRNAmRNA_pair.numberid.txt -input_sam xxx.chr.x.sam.reformed -input_snp xxx.vcf -input_mRNA_simple_bed kept_longest_mRNAs.1_x.NewNumberID.bed -output_reformed xxx.chr.x.readscoverjuncsnp -chr x
	-input_pair	the one of output files of step three (circRNAmRNA_pair). (renamed the transcript).
	-input_sam 	the output of the step six. 
	-input_snp  	this vcf is same as the second step. 
	-input_mRNA_simple_bed 	the mRNA bed is same as the forth step.
	-output_reformed 	this output contains hitmRNA or hitcircRNA, hitRNA name, reads name, circRNA id, mRNA id, snp position, snp relative position, ref anp, slt snp, re number, current bp,   junction relative, read absolute start, read absolute end, read start relative, read end relative, strand, cover junction absolute, merged read.
	-chr 	chromosome number.

###===============================
### Step 8. 
### The statistic result of ref or alt and the number of linear RNA and circRNA of each.
	perl statisticby_pairAndSnp.pl -input_hitnode xxx.chr.x.readscoverjuncsnp -input_sample xxx -input_chromosome x -output_statistic xxx.x.statistic
	-input_hitnode		the output of step eight.
	-input_sample		sample name.
	-input_chromosome 	chromosome number.
	-output_statistic 	output contains simple name, chromosome number, snp pos, ref, alt, junction pos, mRNA id, circRNA id, linear mRNA ref number, linear mRNA alt number, circRNA ref number, circRNA alt number, valid or invalid.		

###===============================
### Step 9. 
### Generate the valid odds rate of each snp.

### Step 9.1 Merge all statistic file of each sample.	
	cat xxx.*.statistic > xxx.statistic
	
### Step 9.2 Generate the valid odds rate of each snp. 
	Rscript ${perl_script_path}/oddsratio.R xxx output
	
	grep -v "sample" xxx.statistics.oddsratio.final | awk '$13>1 && $14>1 && $15>1{print $2"_"$3"_"$6"_"$7"_"$8}'  \
	| LANG=C sort  | awk -F, '{arr[$1]++}END{for (a in arr) print a, arr[a]}' \
	| sort  -k1,1n -k2,2n >xxx.valid.oddsratio.final
	

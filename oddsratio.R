args<-commandArgs(T)
folder<-args[2]
file <-paste0(args[1],".statistic") 

newnold_id<-"merged.circ.mRNA.new.transcriptID" 

setwd(folder)
my_data = read.table(file=paste0(folder,args[1],"/",file), as.is=T,  header = FALSE )
my_data
new_data<-my_data[my_data$V15 == "VALID",]
new_data<-new_data[,c(1:12)]
new_data

id_data = read.table(file=paste0(folder,newnold_id), as.is=T,  header = FALSE )
id_count<-dim(id_data)[1]

id_count

mat_1<-match(new_data[, 7], id_data[, 1])
new_data[, 7]<-id_data[mat_1, 3]

mat_2<-match(new_data[, 8], id_data[, 1])
new_data[, 8]<-id_data[mat_2, 3]

oddsratioWald.proc <- function(n00, n01, n10, n11, alpha = 0.05){	 
	myOR <- (n00 * n11)/(n01 * n10)
	siglog <- sqrt((1/n00) + (1/n01) + (1/n10) + (1/n11))
	logOR<-log(myOR)
	loglo<-exp(logOR-1.96*siglog)
	loghi<-exp(logOR+1.96*siglog)
	oframe<-data.frame(LowerCI = loglo, OR = myOR, UpperCI = loghi, alpha = alpha)
	oframe
}

rowcount<-dim(new_data)[1]

new_data$oddsratio    <- rep(0, rowcount) 
new_data$ci_low      <- rep(0, rowcount) 
new_data$ci_high     <- rep(0, rowcount) 


col_names<-c ("sample","chromosome","SNP","ref","alt","junction","mRNA","circRNA",
"n1","n2","n3","n4","oddsratio","ci_low","ci_high")
colnames(new_data)<-col_names



for(XX in 1:rowcount){ 
	AA <- mRNA_ref     <- new_data[XX,9]
	BB <- mRNA_alt     <- new_data[XX,10]
	CC <- circRNA_ref  <- new_data[XX,11]
	DD <- circRNA_alt  <- new_data[XX,12]
	
	if(AA+BB==0 || AA+CC==0 || BB+DD==0 || CC+DD==0)
	{ next	}
	
	if(AA==0 || BB==0 || CC==0 || DD==0)
	{ AA=AA+0.5;BB=BB+0.5;CC=CC+0.5;DD=DD+0.5 }

	if(AA>0 && BB>0 && CC>0 && DD>0) {
	
		newdf<- oddsratioWald.proc(AA,BB,CC,DD)
		if(newdf$OR<1){ 

			newdf <- oddsratioWald.proc(BB,AA,DD,CC)
		} 		
	
		new_data[XX,13]  <- newdf$OR   
		new_data[XX,14]  <-	newdf$LowerCI   
		new_data[XX,15]  <-	newdf$UpperCI  
		 
	}
	 

} 
select<-c ("sample","chromosome","SNP","ref","alt","junction","mRNA","circRNA",
"n1","n2","n3","n4","oddsratio","ci_low","ci_high" )

save_data <- new_data[,select]
 
write.table(save_data, file = paste0(folder,args[1],"/",args[1],".statistics.oddsratio.final"), row.names = F, quote = F, sep="\t")
 







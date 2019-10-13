##################################################################################################################################
##### Viral and fungal pathogens associated with Pneumolaelaps niutirani (Acari: Laelapidae),
##### a mite found in diseased nests of Vespula wasps (Felden et al. 2019)
##################################################################################################################################

##### Working directory, data files, libraries and stuff #########################################################################

rm(list = ls())
DataDir = "/Users/antoinefelden/Documents/Research/Manuscripts/X-Mites_Wasps/Analysis/01_data"
FigDir = "/Users/antoinefelden/Documents/Research/Manuscripts/X-Mites_Wasps/Analysis/03_figures"

#save(list=ls(all.names=TRUE),file=paste(DataDir,"/mites_wasps_InSoc.RData",sep=""))
#load(paste(DataDir,"/mites_wasps_InSoc.RData",sep=""))

##### Libraries ##################################################################################################################

library("limma")
library("edgeR")
library("RColorBrewer")
library("gplots")
library("rtracklayer")
library("ggplot2")
library("gridExtra")
library("reshape")
library("multcomp")
library("mixOmics")
library("WGCNA")
library("pgirmess")

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

my_palette <- colorRampPalette(c("blue", "white", "red3"))(n = 24)

##### Load data ##################################################################################################################

expression_matrix <- read.table(paste(DataDir,"/genes.TMM.EXPR.matrix",sep=""))

data_rnaseq = read.csv(paste(DataDir,"/gene_count_matrix.csv",sep=""))
x <- DGEList(as.matrix(data_rnaseq[,c(2:length(data_rnaseq))]),remove.zeros=F)

##### BLASTx  ##################################################################################################################

file=paste0(DataDir,"/blastx_viruses.outfmt6")
raw_blast_data = read.table(file,h=TRUE)
colnames(raw_blast_data) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

# Discard isoform data (for Trinity-assembled transcripts)
raw_blast_data$qseqid <- sub("\\_i.*","",raw_blast_data$qseqid)

# Loop to keep best hits for assembled transcripts
data_all_hits_viruses = NULL
p.ident=50
leng=100
for (gene in unique(raw_blast_data$qseqid)) {
  subset_blast_gene <- subset(raw_blast_data,qseqid == gene)
  subset_blast_gene <- subset_blast_gene[subset_blast_gene$length>leng,] # Filter out short sequences
  subset_blast_gene <- subset_blast_gene[subset_blast_gene$pident>p.ident,] # Filter out poorly matching hits
  data_all_hits_viruses <- rbind(data_all_hits_viruses,subset_blast_gene)} # all good hits are in data_all_hits_viruses (including multiple hits per assembled transcript)

# Select the best hit per transcript
data_unique_blast_viruses = NULL
for (gene_loop in unique(data_all_hits_viruses$qseqid)) {
  subset_allhits_gene <- subset(data_all_hits_viruses,qseqid == gene_loop)
  if (nrow(subset_allhits_gene) == 1) { # This if statement selects the best hit based on bit score if there are more than one hit per transcript
    best_match <- subset_allhits_gene } else {
      max_bit <- which.max(subset_allhits_gene$bitscore)
      best_match <- subset_allhits_gene[c(max_bit),]}
  data_unique_blast_viruses <- rbind(data_unique_blast_viruses,best_match)
}

virus_expression_matrix <- expression_matrix[match(data_unique_blast_viruses$qseqid,rownames(expression_matrix)),1:7]
virus_expression_matrix_wasps <- virus_expression_matrix[,1:6]

# Use library sizes of reads mapped to the Vespula genome to get host-library-size-standardised TMM-normalised TPMs
mapped_reads_lib_size <- data.frame("sample" = rownames(x$samples), "alias" = colnames(virus_expression_matrix_wasps), "lib.size" = x$samples$lib.size)
for(n in seq(1:ncol(virus_expression_matrix_wasps))){
  virus_expression_matrix_wasps[,n] <- (virus_expression_matrix_wasps[,n]/mapped_reads_lib_size[n,3])*(mean(mapped_reads_lib_size[,3]))
}

norm_virus_expression_matrix <- cbind(virus_expression_matrix_wasps,"Mites"=virus_expression_matrix[,ncol(virus_expression_matrix)])
norm_virus_expression_matrix$transcript_id <- rownames(norm_virus_expression_matrix)

data_unique_blast_expr_virus <- cbind(data_unique_blast_viruses,norm_virus_expression_matrix)

write.csv(data_unique_blast_expr_virus,file=paste(file,"_",p.ident,"_",leng,"_full.csv",sep=""),row.names = FALSE)
write.table(data_unique_blast_expr_virus$sseqid,file=paste(file,"_",p.ident,"_",leng,"_accession.txt",sep=""),row.names = FALSE,col.names = FALSE, quote = FALSE)

# Input blastx_viruses_50_100_accession.txt into BATCH ENTREZ (https://www.ncbi.nlm.nih.gov/sites/batchentrez) to retrieve taxon IDs
# In BASH, run the below sed commands to modify the BATCH ENTREZ output (i.e. download "Summary") and save file as csv (i.e. output is protein_result_viruses.csv)
#sed -i '' '/aa protein/d' /Users/antoinefelden/Downloads/protein_result.txt
#sed -i '' '/^$/d' /Users/antoinefelden/Downloads/protein_result.txt
#sed -i '' '/]\n/]/s' /Users/antoinefelden/Downloads/protein_result.txt # DOESN'T WORK, FIGURE OUT WHY

hit_names_viruses <- read.csv(paste(DataDir,"/protein_result_viruses.csv", sep=""),h=FALSE)
data_unique_blast_expr_virus$gene <- hit_names_viruses[match(data_unique_blast_expr_virus$sseqid,hit_names_viruses$V2),1]
data_unique_blast_expr_virus$organism <- sub("\\].*", "", sub(".*\\[", "", data_unique_blast_expr_virus$gene))

### High confidence transcripts

HC_transcripts <- data_unique_blast_expr_virus[which(data_unique_blast_expr_virus$pident > 95),]

HC_data_org_hit = NULL
for (org_hit in unique(HC_transcripts$organism)) {
  HC_subset_org_hit <- subset(HC_transcripts, HC_transcripts$organism == org_hit)
  HC_data_org_hit_line <- data.frame("organism" = org_hit, "transcripts" = length(unique(HC_subset_org_hit$gene)), 
                                     "length" = ifelse(nrow(HC_subset_org_hit) == 1, as.character(HC_subset_org_hit$length),
                                                       paste(min(HC_subset_org_hit$length)," - ", max(HC_subset_org_hit$length),sep="")),
                                     "av_length" = ifelse(nrow(HC_subset_org_hit) == 1, as.character(HC_subset_org_hit$length),
                                                          as.character(mean(HC_subset_org_hit$length))),
                                     "pident" = ifelse(nrow(HC_subset_org_hit) == 1, as.character(HC_subset_org_hit$pident),
                                                       paste(min(HC_subset_org_hit$pident)," - ", max(HC_subset_org_hit$pident),sep="")),
                                     "av_pident" = ifelse(nrow(HC_subset_org_hit) == 1, as.character(HC_subset_org_hit$pident),
                                                          as.character(mean(HC_subset_org_hit$pident))),
                                     "D_wasps" = ifelse(sum(HC_subset_org_hit$Diseased_1+HC_subset_org_hit$Diseased_2+HC_subset_org_hit$Diseased_3) > 0,1,0),
                                     "H_wasps"=  ifelse(sum(HC_subset_org_hit$Healthy_1+HC_subset_org_hit$Healthy_2+HC_subset_org_hit$Healthy_3) > 0,1,0),
                                     "Mites" =  ifelse(sum(HC_subset_org_hit$Mites) > 0,1,0))
  HC_data_org_hit <- rbind(HC_data_org_hit,HC_data_org_hit_line)
}

write.csv(HC_data_org_hit,file=paste(DataDir,"/HC_transcripts.csv",sep=""))

sum_viruses = NULL
for (virus_loop in unique(HC_transcripts$organism)){
  viral_subset <- subset(HC_transcripts,HC_transcripts$organism == virus_loop,select=c(1,13:19))
  virus_subset <- melt(viral_subset); virus_subset$condition <- ifelse(virus_subset$variable == "Mites", "Mites", substr(as.character(virus_subset$variable),1,nchar(as.character(virus_subset$variable))-2))
  # Detailed plot of viral transcripts
  virus_subset_no_mites <- subset(virus_subset, virus_subset$condition != "Mites") #remove mites from plot data
  plot_viral_transcripts <- ggplot(data = virus_subset_no_mites, aes(x = virus_subset_no_mites$variable, y = virus_subset_no_mites$value,fill=virus_subset_no_mites$qseqid)) + scale_fill_discrete(name = "Transcript") +
    geom_bar(stat = "identity") +
    theme_bw() + theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
    labs(title=paste(virus_loop,"",sep=""),x="",y="Overall viral load\n(summed host-normalised TMM-normalised TPMs)") +
    theme(plot.title=element_text(size=18, vjust=2),legend.position="right", legend.text=element_text(size=12),
          axis.text.x = element_text(size = 14, colour = "black",angle=90,vjust=0.5),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title.y=element_text(size = 16, colour = "black",vjust=1),
          axis.title.x=element_text(size = 16, colour = "black"))
  print(plot_viral_transcripts)
  dev.copy2pdf(file=paste(FigDir,"/",virus_loop,".pdf",sep=""), width=8, height=6)
  # Sum viral transcripts
  sum_virus_subset <- aggregate(virus_subset[,3],by=list(sample=virus_subset$variable,condition=virus_subset$condition),sum)
  sum_virus_subset$virus <- virus_loop
  sum_viruses <- rbind(sum_viruses,sum_virus_subset)
}

# Filter out zero values to calculate average load when virus present
filt_sum_viruses <- sum_viruses[-which(sum_viruses$x == 0),]
mean_sum_viruses <- aggregate(filt_sum_viruses[,3],by=c(list(filt_sum_viruses$condition,filt_sum_viruses$virus)),mean)

### Low confidence transcripts

LC_transcripts <- data_unique_blast_expr_virus[which(data_unique_blast_expr_virus$pident < 95),]

LC_bitscore_hits <- aggregate(LC_transcripts$bitscore,by=list(LC_transcripts$organism),mean)
colnames(LC_bitscore_hits) <- c("organism","bitscore")
LC_bitscore_hits <- LC_bitscore_hits[order(LC_bitscore_hits$bitscore,decreasing=TRUE),]

LC_data_org_hit = NULL
for (org_hit in unique(LC_bitscore_hits$organism)) {
  LC_subset_org_hit <- subset(LC_transcripts, LC_transcripts$organism == org_hit)
  LC_data_org_hit_line <- data.frame("organism" = org_hit, "transcripts" = length(unique(LC_subset_org_hit$gene)), 
                                     "length" = ifelse(nrow(LC_subset_org_hit) == 1, as.character(LC_subset_org_hit$length),
                                                       paste(min(LC_subset_org_hit$length)," - ", max(LC_subset_org_hit$length),sep="")),
                                     "av_length" = ifelse(nrow(LC_subset_org_hit) == 1, as.character(LC_subset_org_hit$length),
                                                          as.character(mean(LC_subset_org_hit$length))),
                                     "pident" = ifelse(nrow(LC_subset_org_hit) == 1, as.character(LC_subset_org_hit$pident),
                                                       paste(min(LC_subset_org_hit$pident)," - ", max(LC_subset_org_hit$pident),sep="")),
                                     "av_pident" = ifelse(nrow(LC_subset_org_hit) == 1, as.character(LC_subset_org_hit$pident),
                                                          as.character(mean(LC_subset_org_hit$pident))),
                                     "D_wasps" = ifelse(sum(LC_subset_org_hit$Diseased_1+LC_subset_org_hit$Diseased_2+LC_subset_org_hit$Diseased_3) > 0,1,0),
                                     "H_wasps"=  ifelse(sum(LC_subset_org_hit$Healthy_1+LC_subset_org_hit$Healthy_2+LC_subset_org_hit$Healthy_3) > 0,1,0),
                                     "Mites" =  ifelse(sum(LC_subset_org_hit$Mites) > 0,1,0))
  LC_data_org_hit <- rbind(LC_data_org_hit,LC_data_org_hit_line)
}

LC_data_org_hit <- LC_data_org_hit[-which(LC_data_org_hit$H_wasps==0&LC_data_org_hit$D_wasps==0&LC_data_org_hit$Mites==0),]

write.csv(LC_data_org_hit,file=paste(DataDir,"/LC_transcripts.csv",sep=""))

# Manually addeed virus family, virus classification and host range to LC_transcripts.csv > LC_transcripts_edited.csv
LC_data_org_hit_families <- read.csv(paste(DataDir,"/LC_transcripts_edited.csv",sep=""))

arthropod_viruses <- LC_data_org_hit_families[c(which(LC_data_org_hit_families$host == "arthropods")),]
write.csv(arthropod_viruses,file=paste(DataDir,"/LC_viruses_arthropods.csv",sep=""))

other_viruses <- LC_data_org_hit_families[c(which(LC_data_org_hit_families$host != "arthropods")),]
write.csv(other_viruses,file=paste(DataDir,"/LC_viruses_others.csv",sep=""))

## Fig. 2: Venn diagrams for viruses
par(mfrow=c(1,2))

# Arthropod viruses
arthropod_viruses_for_venn <- arthropod_viruses[-which(arthropod_viruses$organism=="Kashmir bee virus"),c("organism","D_wasps","H_wasps","Mites")] # remove low confidence KBV record
HC_records <- data.frame("organism"=c("Moku virus","Kashmir bee virus"),"D_wasps"= c(1,1), "H_wasps"=c(1,1),"Mites"=c(1,0)) # add high confidence records for KBV and Moku
arthropod_viruses_for_venn <- rbind(HC_records,arthropod_viruses_for_venn)
vennDiagram(arthropod_viruses_for_venn[,2:4], circle.col=c("blue","green","red"), font.main=1,
            main="\n\na) Putative viruses\ninfecting arthropods", names = c("Diseased\nwasps\n","Healthy\nwasps\n","Mites"), cex.main=2.5,cex = 2)

# Other viruses
other_viruses_for_venn <- other_viruses[,c("organism","D_wasps","H_wasps","Mites")]
vennDiagram(other_viruses_for_venn[,2:4], circle.col=c("blue","green","red"), font.main=1,
            main="\n\nb) Putative viruses\ninfecting other life forms", names = c("Diseased\nwasps\n","Healthy\nwasps\n","Mites"), cex.main=2.5,cex = 2)

dev.copy2pdf(file=paste(FigDir,"/venn_diagrams_viruses.pdf",sep=""), width=16, height=8)

## Viral quantification

sum_viruses_wasps <- subset(sum_viruses,condition != "Mites")

moku_wasps=subset(sum_viruses_wasps,virus %in% "Moku virus")[-6,]
summary_moku_wasps <- data.frame("condition" = c("Diseased","Healthy"),
                                 "x" = c(mean(subset(moku_wasps$x,moku_wasps$condition=="Diseased")),mean(subset(moku_wasps$x,moku_wasps$condition=="Healthy"))),
                                 "se" = c(sd(subset(moku_wasps$x,moku_wasps$condition=="Diseased"))/3,sd(subset(moku_wasps$x,moku_wasps$condition=="Healthy"))/length(subset(moku_wasps$x,moku_wasps$condition=="Healthy"))))
kw_moku <- kruskal.test(moku_wasps$x,as.factor(moku_wasps$condition))
test_moku_wasps = paste("X2 = ",round(kw_moku$statistic,3),", df = 1, p = ",round(kw_moku$p.value,3),sep="")

kbv_wasps=subset(sum_viruses_wasps,virus %in% "Kashmir bee virus")
summary_kbv_wasps <- data.frame("condition" = c("Diseased","Healthy"),
                                 "x" = c(mean(subset(kbv_wasps$x,kbv_wasps$condition=="Diseased")),mean(subset(kbv_wasps$x,kbv_wasps$condition=="Healthy"))),
                                 "se" = c(sd(subset(kbv_wasps$x,kbv_wasps$condition=="Diseased"))/3,sd(subset(kbv_wasps$x,kbv_wasps$condition=="Healthy"))/length(subset(kbv_wasps$x,kbv_wasps$condition=="Healthy"))))
kw_kbv <- kruskal.test(kbv_wasps$x,as.factor(kbv_wasps$condition))
test_kbv_wasps = paste("X2 = ",round(kw_kbv$statistic,3),", df = 1, p = ",round(kw_kbv$p.value,3),sep="")

# plot of means
plot_moku_mean <- ggplot(data = summary_moku_wasps, aes(x = condition, y = x)) +
  geom_bar(stat = "identity", position = position_dodge(0.90),fill=c("grey40","grey40")) +
  geom_errorbar(aes(ymax = summary_moku_wasps$x + summary_moku_wasps$se, ymin = summary_moku_wasps$x - summary_moku_wasps$se),
                position = position_dodge(0.90), width = 0.25) +
  theme_bw() + theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  labs(title=expression("(a)"~italic("Moku virus")),
       subtitle=test_moku_wasps,
        x="Wasp larvae condition",y="Viral load\n(summed host-normalised TMM-normalised TPMs)") +
  theme(plot.title=element_text(size=18, vjust=2, face = "bold"),legend.position="", legend.text=element_text(size=14), plot.subtitle = element_text(size=16),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.y=element_text(size = 16, colour = "black",vjust=1),
        axis.title.x=element_text(size = 16, colour = "black"))
plot_moku_mean

plot_kbv_mean <- ggplot(data = summary_kbv_wasps, aes(x = condition, y = x)) +
  geom_bar(stat = "identity", position = position_dodge(0.90),fill=c("grey40","grey40")) +
  geom_errorbar(aes(ymax = summary_kbv_wasps$x + summary_kbv_wasps$se, ymin = summary_kbv_wasps$x - summary_kbv_wasps$se),
                position = position_dodge(0.90), width = 0.25) +
  theme_bw() + theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  labs(title=expression("(a)"~italic("Kashmir bee virus")),
       subtitle=test_kbv_wasps,
       x="Wasp larvae condition",y="\n") +
  theme(plot.title=element_text(size=18, vjust=2, face = "bold"),legend.position="", legend.text=element_text(size=14), plot.subtitle = element_text(size=16),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.y=element_text(size = 16, colour = "black",vjust=1),
        axis.title.x=element_text(size = 16, colour = "black"))
plot_kbv_mean

grid.arrange(plot_moku_mean,plot_kbv_mean,ncol=2,nrow=1)

dev.copy2pdf(file=paste(FigDir,"/viruses_transcript_quant.pdf",sep=""), width=12, height=6)

### Fungi ###########

file=paste0(DataDir,"/blastx_fungi.outfmt6")
raw_blast_data = read.table(file,h=TRUE)
colnames(raw_blast_data) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

# Discard isoform data (for Trinity-assembled transcripts)
raw_blast_data$qseqid <- sub("\\_i.*","",raw_blast_data$qseqid)

# Loop to keep best match per gene
data_all_hits_fungi = NULL
p.ident=95
leng=100
for (gene in unique(raw_blast_data$qseqid)) {
  subset_blast_gene <- subset(raw_blast_data,qseqid == gene)
  subset_blast_gene <- subset_blast_gene[subset_blast_gene$length>leng,] # Filter out short sequences
  subset_blast_gene <- subset_blast_gene[subset_blast_gene$pident>p.ident,] # Filter out poorly matching hits
  data_all_hits_fungi <- rbind(data_all_hits_fungi,subset_blast_gene)} # all good hits are in data_all_hits_fungi (including multiple hits per assembled transcript)

# Select the best hit per transcript
data_unique_blast_fungi = NULL
for (gene_loop in unique(data_all_hits_fungi$qseqid)) {
  subset_allhits_gene <- subset(data_all_hits_fungi,qseqid == gene_loop)
  if (nrow(subset_allhits_gene) == 1) { # This if statement selects the best hit based on bit score if there are more than one hit per transcript
    best_match <- subset_allhits_gene } else {
      max_bit <- which.max(subset_allhits_gene$bitscore)
      best_match <- subset_allhits_gene[c(max_bit),]}
  data_unique_blast_fungi <- rbind(data_unique_blast_fungi,best_match)
}

fungi_expression_matrix <- expression_matrix[match(data_unique_blast_fungi$qseqid,rownames(expression_matrix)),1:7]
fungi_expression_matrix_wasps <- fungi_expression_matrix[,1:6]

# Use library sizes of reads mapped to the Vespula genome to get host-library-size-standardised TMM-normalised TPMs
mapped_reads_lib_size <- data.frame("sample" = rownames(x$samples), "alias" = colnames(fungi_expression_matrix_wasps), "lib.size" = x$samples$lib.size)
for(n in seq(1:ncol(fungi_expression_matrix_wasps))){
  fungi_expression_matrix_wasps[,n] <- (fungi_expression_matrix_wasps[,n]/mapped_reads_lib_size[n,3])*(mean(mapped_reads_lib_size[,3]))
}

norm_fungi_expression_matrix <- cbind(fungi_expression_matrix_wasps,"Mites"=fungi_expression_matrix[,ncol(fungi_expression_matrix)])
norm_fungi_expression_matrix$transcript_id <- rownames(norm_fungi_expression_matrix)

data_unique_blast_expr_fungi <- cbind(data_unique_blast_fungi,norm_fungi_expression_matrix)

write.csv(data_unique_blast_expr_fungi,file=paste(file,"_",p.ident,"_",leng,"_full.csv",sep=""),row.names = FALSE)
write.table(data_unique_blast_expr_fungi$sseqid,file=paste(file,"_",p.ident,"_",leng,"_accession.txt",sep=""),row.names = FALSE,col.names = FALSE, quote=FALSE)

# Input blastx_fungi_50_100_accession.txt into BATCH ENTREZ (https://www.ncbi.nlm.nih.gov/sites/batchentrez) to retrieve taxon IDs
# In BASH, run the below sed commands to modify the BATCH ENTREZ output (i.e. download "Summary") and save file as csv (i.e. output is protein_result_fungi.csv)
#sed -i '' '/aa protein/d' /Users/antoinefelden/Downloads/protein_result.txt
#sed -i '' '/^$/d' /Users/antoinefelden/Downloads/protein_result.txt
#sed -i '' '/]\n//' /Users/antoinefelden/Downloads/protein_result.txt # DOESN'T WORK, FIGURE OUT WHY

hit_names_fungi <- read.csv(paste0(DataDir,"/protein_result_fungi.csv"),h=FALSE)
data_unique_blast_expr_fungi$gene <- hit_names_fungi[match(data_unique_blast_expr_fungi$sseqid,hit_names_fungi$V2),1]
data_unique_blast_expr_fungi$organism <- sub("\\].*", "", sub(".*\\[", "", data_unique_blast_expr_fungi$gene))

fungi_transcripts <- data_unique_blast_expr_fungi[which(data_unique_blast_expr_fungi$pident > 95),]

sum_fungi = NULL
for (fungus_loop in unique(fungi_transcripts$organism)){
  fungal_subset <- subset(fungi_transcripts,fungi_transcripts$organism == fungus_loop,select=c(1,13:19))
  fungus_subset <- melt(fungal_subset); fungus_subset$condition <- ifelse(fungus_subset$variable == "Mites", "Mites", substr(as.character(fungus_subset$variable),1,nchar(as.character(fungus_subset$variable))-2))
  # Detailed plot of viral transcripts
  fungus_subset_no_mites <- subset(fungus_subset, fungus_subset$condition != "Mites") #remove mites from plot data
  plot_fungal_transcripts <- ggplot(data = fungus_subset_no_mites, aes(x = fungus_subset_no_mites$variable, y = fungus_subset_no_mites$value,fill=fungus_subset_no_mites$qseqid)) + scale_fill_discrete(name = "Transcript") +
    geom_bar(stat = "identity") +
    theme_bw() + theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
    labs(title=paste(fungus_loop,"",sep=""),x="",y="Overall fungi load\n(summed host-normalised TMM-normalised TPMs)") +
    theme(plot.title=element_text(size=18, vjust=2),legend.position="none", legend.text=element_text(size=12),
          axis.text.x = element_text(size = 14, colour = "black",angle=90,vjust=0.5),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title.y=element_text(size = 16, colour = "black",vjust=1),
          axis.title.x=element_text(size = 16, colour = "black"))
  print(plot_fungal_transcripts)
  dev.copy2pdf(file=paste(FigDir,"/fungi/",fungus_loop,".pdf",sep=""), width=8, height=6)
  # Sum fungal transcripts
  sum_fungus_subset <- aggregate(fungus_subset[,3],by=list(sample=fungus_subset$variable,condition=fungus_subset$condition),sum)
  sum_fungus_subset$fungus <- fungus_loop
  sum_fungi <- rbind(sum_fungi,sum_fungus_subset)
}

# Filter out zero values
filt_sum_fungi <- sum_fungi[-which(sum_fungi$x == 0),]
mean_sum_fungi <- aggregate(filt_sum_fungi[,3],by=c(list(filt_sum_fungi$condition,filt_sum_fungi$fungus)),mean)

sum_fungi_formatted = NULL
for (org_hit in unique(filt_sum_fungi$fungus)) {
  sum_fungi_formatted_subset <- subset(filt_sum_fungi, filt_sum_fungi$fungus == org_hit)
  sum_fungi_formatted_line <- data.frame("organism" = org_hit,
                                         "mean_dis" = round(mean(sum_fungi_formatted_subset[sum_fungi_formatted_subset$condition == "Diseased","x"]),3),
                                         "range_dis" = ifelse(nrow(sum_fungi_formatted_subset[sum_fungi_formatted_subset$condition == "Diseased",]) == 1,
                                                              "-",
                                                          paste("(",round(min(sum_fungi_formatted_subset[sum_fungi_formatted_subset$condition == "Diseased","x"]),3),
                                                                " - ", round(max(sum_fungi_formatted_subset[sum_fungi_formatted_subset$condition == "Diseased","x"]),3),")",sep="")),
                                         "n_dis" = length(unique(sum_fungi_formatted_subset[sum_fungi_formatted_subset$condition == "Diseased","x"])),
                                         "mean_heal" = round(mean(sum_fungi_formatted_subset[sum_fungi_formatted_subset$condition == "Healthy","x"]),3),
                                         "range_heal" = ifelse(nrow(sum_fungi_formatted_subset[sum_fungi_formatted_subset$condition == "Healthy",]) == 1,
                                                               "-",
                                                          paste("(",round(min(sum_fungi_formatted_subset[sum_fungi_formatted_subset$condition == "Healthy","x"]),2),
                                                                " - ", round(max(sum_fungi_formatted_subset[sum_fungi_formatted_subset$condition == "Healthy","x"]),3),")",sep="")),
                                         "n_heal" = length(unique(sum_fungi_formatted_subset[sum_fungi_formatted_subset$condition == "Healthy","x"])),
                                         "mites" = ifelse(isEmpty(sum_fungi_formatted_subset[sum_fungi_formatted_subset$condition == "Mites","x"])==FALSE,
                                                          round(sum_fungi_formatted_subset[sum_fungi_formatted_subset$condition == "Mites","x"],3),"-"))
sum_fungi_formatted <- rbind(sum_fungi_formatted,sum_fungi_formatted_line)
}
sum_fungi_formatted <- sum_fungi_formatted[order(as.character(sum_fungi_formatted$organism)),]
write.csv(sum_fungi_formatted,paste(DataDir,"/fungi_quant.csv",sep=""))

# Extract Aspergillus data

aspergillus_summed <- sum_fungi[grep("Aspergillus",sum_fungi$fungus),]
aspergillus_summed <- subset(aspergillus_summed, !(aspergillus_summed$condition %in% "Mites"))
aspergillus_summed$fungus <- ifelse(aspergillus_summed$fungus == "Aspergillus novofumigatus IBT 16806", "Aspergillus novofumigatus", "Aspergillus spp")

aspergillus_all_summed <- aggregate(aspergillus_summed$x~aspergillus_summed$sample,data=aspergillus_summed[,c(1,3)],FUN=sum)
colnames(aspergillus_all_summed) <- c("sample","x")
aspergillus_all_summed$condition <- c(rep("Diseased",3),rep("Healthy",3))

summary_aspergillus_wasps <- data.frame("condition" = c("Diseased","Healthy"),
                                "x" = c(mean(subset(aspergillus_all_summed$x,aspergillus_all_summed$condition=="Diseased")),mean(subset(aspergillus_all_summed$x,aspergillus_all_summed$condition=="Healthy"))),
                                "se" = c(sd(subset(aspergillus_all_summed$x,aspergillus_all_summed$condition=="Diseased"))/3,sd(subset(aspergillus_all_summed$x,aspergillus_all_summed$condition=="Healthy"))/length(subset(aspergillus_all_summed$x,aspergillus_all_summed$condition=="Healthy"))))

kw_aspergillus <- kruskal.test(aspergillus_all_summed$x,as.factor(aspergillus_all_summed$condition))
test_aspergillus_wasps = paste("X2 = ",round(kw_aspergillus$statistic,3),", df = 1, p < 0.001",sep="")

plot_aspergillus_mean <- ggplot(data = summary_aspergillus_wasps, aes(x = condition, y = x)) +
  geom_bar(stat = "identity", position = position_dodge(0.90),fill=c("grey40","grey40")) +
  geom_errorbar(aes(ymax = summary_aspergillus_wasps$x + summary_aspergillus_wasps$se, ymin = summary_aspergillus_wasps$x - summary_aspergillus_wasps$se),
                position = position_dodge(0.90), width = 0.25) +
  theme_bw() + theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  labs(title=expression("(b) Overall"~italic("Aspergillus")~"load"),
       subtitle=test_aspergillus_wasps,
       x="Wasp larvae condition",y="Summed host-normalised\nTMM-normalised TPM") +
  theme(plot.title=element_text(size=18, vjust=2, face = "bold"),legend.position="", legend.text=element_text(size=14), plot.subtitle = element_text(size=16),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.y=element_text(size = 16, colour = "black",vjust=1),
        axis.title.x=element_text(size = 16, colour = "black"))
plot_aspergillus_mean

dev.copy2pdf(file=paste(FigDir,"/Aspergillus_transcript_quant.pdf",sep=""), width=5, height=5)

### Venn Diagram

fungi_bitscore_hits <- aggregate(fungi_transcripts$bitscore,by=list(fungi_transcripts$organism),mean)
colnames(fungi_bitscore_hits) <- c("organism","bitscore")
fungi_bitscore_hits <- fungi_bitscore_hits[order(fungi_bitscore_hits$bitscore,decreasing=TRUE),]

fungi_data_org_hit = NULL
for (org_hit in unique(fungi_bitscore_hits$organism)) {
  fungi_subset_org_hit <- subset(fungi_transcripts, fungi_transcripts$organism == org_hit)
  fungi_data_org_hit_line <- data.frame("organism" = org_hit, "transcripts" = length(unique(fungi_subset_org_hit$gene)), 
                                     "length" = ifelse(nrow(fungi_subset_org_hit) == 1, as.character(fungi_subset_org_hit$length),
                                                       paste(min(fungi_subset_org_hit$length)," - ", max(fungi_subset_org_hit$length),sep="")),
                                     "av_length" = ifelse(nrow(fungi_subset_org_hit) == 1, as.character(fungi_subset_org_hit$length),
                                                          as.character(mean(fungi_subset_org_hit$length))),
                                     "pident" = ifelse(nrow(fungi_subset_org_hit) == 1, as.character(fungi_subset_org_hit$pident),
                                                       paste(min(fungi_subset_org_hit$pident)," - ", max(fungi_subset_org_hit$pident),sep="")),
                                     "av_pident" = ifelse(nrow(fungi_subset_org_hit) == 1, as.character(fungi_subset_org_hit$pident),
                                                          as.character(mean(fungi_subset_org_hit$pident))),
                                     "D_wasps" = ifelse(sum(fungi_subset_org_hit$Diseased_1+fungi_subset_org_hit$Diseased_2+fungi_subset_org_hit$Diseased_3) > 0,1,0),
                                     "H_wasps"=  ifelse(sum(fungi_subset_org_hit$Healthy_1+fungi_subset_org_hit$Healthy_2+fungi_subset_org_hit$Healthy_3) > 0,1,0),
                                     "Mites" =  ifelse(sum(fungi_subset_org_hit$Mites) > 0,1,0))
  fungi_data_org_hit <- rbind(fungi_data_org_hit,fungi_data_org_hit_line)
}

write.csv(fungi_data_org_hit,file=paste(DataDir,"/fungi_transcripts.csv",sep=""))

par(mfrow=c(1,1),cex=1)
vennDiagram(fungi_data_org_hit[,7:9], circle.col=c("blue","green","red"), font.main=1,
            names=c("Diseased\nwasps\n","Healthy\nwasps\n","Mites"),main="\n\n(a) Fungi shared between wasps\nand mites",cex.main=2.5,cex = 2)
dev.copy2pdf(file=paste(FigDir,"/fungi_venn_diagram.pdf",sep=""), width=8, height=8)

### qPCR confirmation ###########

qPCR <- read.csv(paste(DataDir,"/qPCR_taqman_array_KBV_Aspergillus.csv",sep=""))
qPCR_KBV <- qPCR[grep("Kashmir bee virus",qPCR$pathogen),]
qPCR_aspergillus <- qPCR[grep("Aspergillus",qPCR$pathogen),]

plot_KBV_qPCR <- ggplot(data = qPCR_KBV) + scale_fill_discrete(name = "") +
  geom_segment(aes(y = 0, yend = qPCR_KBV$x, x = qPCR_KBV$sample, xend = qPCR_KBV$sample),size=20) + scale_y_log10() +
  theme_bw() + theme(plot.title = element_text(face="bold", size=18, hjust=0)) + scale_x_discrete(labels=c("Diseased\n1","Diseased\n2","Diseased\n3","Healthy\n1","Healthy\n2","Healthy\n3"))+
  labs(title=expression("(a)"~italic("Kashmir bee virus")~"load (RT-qPCR)"),x="",y=bquote('Pathogen load (Pathogen 2'^'-Ct'~'/Ref 2'^'-Ct'~')')) +
  theme(plot.title=element_text(size=18, vjust=2),legend.position="none", legend.text=element_text(size=12),
        axis.text.x = element_text(size = 14, colour = "black",vjust=0.5),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 16, colour = "black",vjust=1),
        axis.title.x=element_text(size = 16, colour = "black"))
print(plot_KBV_qPCR)

plot_aspergillus_qPCR <- ggplot(data = qPCR_aspergillus) + scale_fill_discrete(name = "") +
  geom_segment(aes(y = 0, yend = qPCR_aspergillus$x, x = qPCR_aspergillus$sample, xend = qPCR_aspergillus$sample),size=20) + scale_y_log10() +
  theme_bw() + theme(plot.title = element_text(face="bold", size=18, hjust=0)) +  scale_x_discrete(labels=c("Diseased\n1","Diseased\n2","Diseased\n3","Healthy\n1","Healthy\n2","Healthy\n3"))+
  labs(title=expression("(b)"~italic("Aspergillus")~"load (RT-qPCR)"),x="",y="\n") +
  theme(plot.title=element_text(size=18, vjust=2),legend.position="none", legend.text=element_text(size=12),
        axis.text.x = element_text(size = 14, colour = "black",vjust=0.5),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 16, colour = "black",vjust=1),
        axis.title.x=element_text(size = 16, colour = "black"))
print(plot_aspergillus_qPCR)

# Detailed plots from RNA-Seq

plot_kbv <- ggplot(kbv_wasps) + scale_fill_discrete(name = "") +
  geom_segment(aes(y = 0, yend = kbv_wasps$x, x = kbv_wasps$sample, xend = kbv_wasps$sample),size=20) + scale_y_log10() +
  scale_x_discrete(labels=c("Diseased\n1","Diseased\n2","Diseased\n3","Healthy\n1","Healthy\n2","Healthy\n3")) +
  labs(title=expression("(c)"~italic("Kashmir bee virus")~"load (RNA-Seq)"),x="",y="\n\n") +
  theme_bw() +   theme(plot.title=element_text(size=18, vjust=2),legend.position="none", legend.text=element_text(size=12),
                       axis.text.x = element_text(size = 14, colour = "black",vjust=0.5),
                       axis.text.y = element_text(size = 14, colour = "black"),
                       axis.title.y=element_text(size = 16, colour = "black",vjust=1),
                       axis.title.x=element_text(size = 16, colour = "black"))

plot_aspergillus <- ggplot(aspergillus_all_summed) + scale_fill_discrete(name = "") +
  geom_segment(aes(y = 0, yend = aspergillus_all_summed$x, x = aspergillus_all_summed$sample, xend = aspergillus_all_summed$sample),size=20) + scale_y_log10() +
  scale_x_discrete(labels=c("Diseased\n1","Diseased\n2","Diseased\n3","Healthy\n1","Healthy\n2","Healthy\n3")) +
  labs(title=expression("(d)"~italic("Aspergillus")~"load (RNA-Seq)"),x="",y="\n\n") +
  theme_bw() +   theme(plot.title=element_text(size=18, vjust=2),legend.position="none", legend.text=element_text(size=12),
                       axis.text.x = element_text(size = 14, colour = "black",vjust=0.5),
                       axis.text.y = element_text(size = 14, colour = "black"),
                       axis.title.y=element_text(size = 16, colour = "black",vjust=1),
                       axis.title.x=element_text(size = 16, colour = "black"))

grid.arrange(plot_KBV_qPCR,plot_aspergillus_qPCR,plot_kbv,plot_aspergillus, ncol=2,nrow=2)
dev.copy2pdf(file=paste(DataDir,"/qPCR_confirmation.pdf",sep=""), width=14, height=12)
## Maya Willey 
#8/10/2023

library(data.table)
library(dplyr)
library(stringr)
library("readxl")
## Setup
#/////////////////////////////////////////////////

  fastacount<-read.csv("/data/CARD_AA/users/willeyml/ViralRNA_Blast/Files/FastaCount.csv", header=FALSE)
  colnames(fastacount)<-c("Count","SampleID")
  fastacount$SeqCount<-as.numeric(fastacount$Count)/4
  bfiles<- list.files(path= "/data/CARD_AA/users/willeyml/ViralRNA_Blast/BLAST_OUT", pattern = ".blast.gz")
  afiles<- list.files(path= "/data/CARD_AA/data/2023_05_16_ROSMAP_sample_FASTQs/",pattern = '.fastq.gz')
  fl<-data.frame(rna_id=bfiles)
  seqid<- read.csv("/data/CARD_AA/users/willeyml/ViralRNA_Blast/Files/SeqID_newfasta.csv",header=T)
  cohort_meta<-fread("/data/CARD_AA/users/willeyml/ViralRNA_Blast/Files/threeCohorts.csv")
  cohort_xtra<-fread("/data/CARD_AA/users/willeyml/dataset_1364_long_03-12-2024.csv")
  cohort_meta<-cohort_meta[cohort_meta$cohort == 'ROSMAP']
  rosmap_meta<-cohort_meta[cohort_meta$rna_id != '']
  fl<-str_remove(fl$rna_id,".blast.gz")
  fl<-data.frame(rna_id=fl)

##Merged Meta
#/////////////////////////////////////////////////
  harmmeta<- fread("/data/CARD_AA/users/willeyml/ViralRNA_Blast/Files/RNAseq_Harmonization_ROSMAP_combined_metadata.csv")
  sl<-fl
  colnames(sl)<-"specimenID"
  hm<- harmmeta[harmmeta$assay == 'rnaSeq']

  nm<- sl %>% left_join(hm)
  ri<- nm[,c("specimenID","individualID","projid")]
  colnames(nm)[colnames(nm) == 'specimenID'] <-"rna_id"
  newmeta<-nm %>% left_join(rosmap_meta)

##new Meta
#/////////////////////////////////////////////////
  cohort_xtra<-fread("/data/CARD_AA/users/willeyml/dataset_1364_long_03-12-2024.csv")  
  tmp_meta<- ri %>% full_join(cohort_xtra, by = join_by(projid))
  tmp_meta$libraryPrep == "polyAselection"
  sl<-fl
  colnames(sl)<-"specimenID"

##labeling samples by metadata: using both
#/////////////////////////////////////////////////
  colnames(newmeta)[colnames(newmeta) == 'specimenID'] <-"rna_id"
  setDT(fl)
  setDT(newmeta)
  tmp_m<- newmeta[fl, mult = "first", on ="rna_id"]
  tmp_meta<-sl %>% left_join(rosmap_meta)
  meta<-data.frame(rna_id=c(paste(tmp_meta$rna_id,"_",tmp_meta$clinAD,"AD","_",tmp_meta$NIA_RI,"_",tmp_meta$libraryBatch,"LB",sep='')))
  #meta<-data.frame(rna_id=c(paste(tmp_meta$rna_id,"_",tmp_meta$pmi,sep='')))
  #meta<-data.frame(rna_id=c(paste(tmp_meta$rna_id,"_",tmp_meta$braaksc,"Bsc",sep='')))

 

##labeling samples by kmeans: run own clustering before using
#/////////////////////////////////////////////////

  #kmean<-read_excel("INSERT_OWN", sheet = "Sheet1")
  #kmean<-kmean[,c(1,4)]
  #kmeans_grp<-fl %>% left_join(kmean)

##labeling viruses by their host or genus
#/////////////////////////////////////////////////

  #seqid_host<-read.csv("/data/CARD_AA/users/willeyml/ViralRNA_Blast/Files/sequences-3.csv", header=TRUE)
  #seqid_host$Host[seqid_host$Host == '']=NA

  #seqid_genus<-read.csv("/data/CARD_AA/users/willeyml/ViralRNA_Blast/Files/sequences-5.csv", header=TRUE)
  #seqid_genus$Genus[seqid_genus$Genus == '']=NA

## Creating datafame, extracting results
#/////////////////////////////////////////////////

indx<-nrow(seqid)-1
Viraldf<- data.frame(index=(c(0:indx)))

for (i in 1:length(bfiles)){
    tmp<- fread(paste("/data/CARD_AA/users/willeyml/ViralRNA_Blast/BLAST_OUT/",bfiles[i],sep=""))
    seq_cnt<- fastacount[fastacount$SampleID == paste(str_remove(bfiles[i],pattern= ".blast.gz")),]
    seq_cnt<- seq_cnt$SeqCount
    ref<-tmp[,2]
    unq<- unique(ref)
    unq$x<-1
        for (k in 1:sum(unq$x)) {
            reflist<-as.character(ref)
            utmp<-as.vector(paste0(unq[k,1]))
            filt<- str_count(reflist, utmp)
            unq[k,2]<-filt/seq_cnt
        }
    
    colnames(unq)<-c("index","Binary")
    #as.numeric(unq$index)
    try(lj<-Viraldf %>% full_join(unq))
    
  if (nrow(tmp)>1) {
        Viraldf[,ncol(Viraldf)+1]<- lj$Binary
    } else { 
    Viraldf[,ncol(Viraldf)+1]<- NA}
   
} 
head(Viraldf)


#count fasta lines with all rna (including human) /4 then divide viral count by total fasta line count 
#add 1 to all then log transform

colnames(Viraldf)= c("Index",meta$rna_id)
Viraldf$Index<- seqid$Accession
rownames(Viraldf)<- Viraldf$Index
head(Viraldf)

Viraldf_full<- Viraldf[which(rowSums(is.na(Viraldf))<639),]
vdf_full<- Viraldf_full[colnames(Viraldf_full) != 'Index']
vdf_full[is.na(vdf_full)]=0
vdf_bil_full<- (vdf_full*1000000000)+1
vdf_log_full<- log2(vdf_bil_full)
vdf_rmv78<- vdf_log_full[,grep(pattern = "7LB|8LB",colnames(vdf_log_full), invert = TRUE)]
colnames(vdf_rmv78)
write.xlsx(vdf_rmv78, file.path("vdf_full_log2_rm78.xlsx"), append = TRUE)
vdf_rmv78 <- read_excel("vdf_full_log2_rm78.xlsx")
vdf_log_full <- read_excel("vdf_full_log2.xlsx")
##t-test and logfc: ClinAD
#/////////////////////////////////////////////////

AD<- vdf_rmv78[,grep(pattern= "_1AD_1",colnames(vdf_rmv78))]
Cont<- vdf_rmv78[,grep(pattern= "_0AD_0",colnames(vdf_rmv78))]
colnames(AD)
vs_df<-data.frame(rownames(AD))
vs_df$pval<-NA
vs_df$logfc<-NA

for (j in 1:561) {
  viralstat<- t.test(AD[j,],Cont[j,]) 
  vs_df[j,2]<- viralstat$p.value
  ADavg<-rowMeans(AD[j,])
  Contavg<-rowMeans(Cont[j,])
  L2FC<- ADavg - Contavg
  vs_df[j,3]<- L2FC
}
 padj<- p.adjust(vs_df$pval, method = "BH")
 vs_df[,4]<-padj
write.xlsx(vs_df, file.path("ViralSeq_RNA_log2_1AD1_0AD0_newseq_missing_rmv78_padjBH.xlsx"), append = TRUE)

###### t-test and logfc: ClinAD and NIA_RI
AD_1<- vdf_log_full[,grep(pattern= "_1AD_1",colnames(vdf_log_full))]
Cont_1<- vdf_log_full[,grep(pattern= "_0AD_0",colnames(vdf_log_full))]

vs_df<-data.frame(rownames(AD_1))
vs_df$pval<-NA
vs_df$logfc<-NA
vs_df$padj<-NA

for (j in 1:788) {
  viralstat<- t.test(AD_1[j,],Cont_1[j,]) 
  vs_df[j,2]<- viralstat$p.value
  ADavg<-rowMeans(AD_1[j,])
  Contavg<-rowMeans(Cont_1[j,])
  L2FC<- ADavg - Contavg
  vs_df[j,3]<- L2FC
 
}
 padj<- p.adjust(vs_df$pval, method = "BH")
 vs_df[,4]<-padj
write.xlsx(vs_df, file.path("ViralSeq_RNA_Stats_log2_1AD1_0AD0_newseq_missing_rmv78_padjBH_788.xlsx"), append = TRUE)

###### t-test and logfc: cohort
ROS<- vdf_log_full[,grep(pattern= "_ROS",colnames(vdf_log_full))]
MAP<- vdf_log_full[,grep(pattern= "_MAP",colnames(vdf_log_full))]

vs_df<-data.frame(rownames(ROS))
vs_df$pval<-NA
vs_df$logfc<-NA
vs_df$padj<-NA

for (j in 1:561) {
  viralstat<- t.test(ROS[j,],MAP[j,]) 
  vs_df[j,2]<- viralstat$p.value
  ROSavg<-rowMeans(ROS[j,])
  MAPavg<-rowMeans(MAP[j,])
  L2FC<- ROSavg - MAPavg
  vs_df[j,3]<- L2FC
 
}
 padj<- p.adjust(vs_df$pval, method = "BH")
 vs_df[,4]<-padj
write.xlsx(vs_df, file.path("ViralSeq_RNA_Stats_log1_ROS_MAP_newseq_missing_padjBH.xlsx"), append = TRUE)

###### t-test and logfc: cohort and AD
ROS_AD<- vdf_log_full[,grep(pattern= "_1AD_1_ROS",colnames(vdf_log_full))]
ROS_CNT<- vdf_log_full[,grep(pattern= "_0AD_0_ROS",colnames(vdf_log_full))]

vs_df<-data.frame(rownames(ROS_AD))
vs_df$pval<-NA
vs_df$logfc<-NA
vs_df$padj<-NA

for (j in 1:788) {
  viralstat<- t.test(ROS_AD[j,],ROS_CNT[j,]) 
  vs_df[j,2]<- viralstat$p.value
  ADavg<-rowMeans(ROS_AD[j,])
  CNTavg<-rowMeans(ROS_CNT[j,])
  L2FC<- ADavg - CNTavg
  vs_df[j,3]<- L2FC
 
}
 padj<- p.adjust(vs_df$pval, method = "BH")
 vs_df[,4]<-padj
write.xlsx(vs_df, file.path("ViralSeq_RNA_Stats_log1_ROS_AD_RI_nf_m_corrREF_padjBH.xlsx"), append = TRUE)

###### t-test and logfc: cohort and AD
MAP_AD<- vdf_log_full[,grep(pattern= "_1AD_1_MAP",colnames(vdf_log_full))]
MAP_CNT<- vdf_log_full[,grep(pattern= "_0AD_0_MAP",colnames(vdf_log_full))]

vs_df<-data.frame(rownames(MAP_AD))
vs_df$pval<-NA
vs_df$logfc<-NA
vs_df$padj<-NA

for (j in 1:788) {
  viralstat<- t.test(MAP_AD[j,],MAP_CNT[j,]) 
  vs_df[j,2]<- viralstat$p.value
  ADavg<-rowMeans(MAP_AD[j,])
  CNTavg<-rowMeans(MAP_CNT[j,])
  L2FC<- ADavg - CNTavg
  vs_df[j,3]<- L2FC
 
}
 padj<- p.adjust(vs_df$pval, method = "BH")
 vs_df[,4]<-padj
write.xlsx(vs_df, file.path("ViralSeq_RNA_Stats_log1_MAP_AD_RI_nf_m_corREF_padjBH.xlsx"), append = TRUE)
###########

blast792_13053<- fread("/data/CARD_AA/users/willeyml/ViralRNA_Blast/Output/792_130530.blast")
proteus<- blast792_13053[grepl(pattern= "9849",blast792_13053$"reference acc.")]
write.xlsx(proteus, file.path("ViralSeq_RNA_Parainfluenza5Virus_792_130530.xlsx"), append = TRUE)
######### Viral tracks

######### Host

  virus_color<-rainbow(25)
  seqhost_list<-unique(seqid_host$Host)

  #clinAD
  anot<- data.frame(ClinAD=as.character(tmp_meta$clinAD))
  anotcolors<- list(ClinAD= c('1' = 'chartreuse3', '0' = 'cornflowerblue'),ViralHost= virus_color)
  rownames(anot)<-colnames(vdf_log_half)
  vdfhalf<-data.frame(Accession=rownames(vdf_log_half))
  host<- vdfhalf %>% left_join(seqid_host)
  anotrow<- data.frame(ViralHost=as.character(host$Host))
  #anotcolors<- list(ViralHost=virus_color)
  rownames(anotrow)<- vdfhalf$Accession

  pdf(("Heatmap_viralRNA_ROSMAP_clinAD_logHalf1_Host.pdf"),onefile=T, height = 20, width = 25)
  pheatmap(mat=vdf_log_half,color=colorRampPalette(c("blue","white","red"))(100), annotation_col = anot, annotation_colors = anot_viruscol, annotation_row = anotrow , scale='row',show_rowname=T,show_colname=F, border_color= "grey",na_col="grey", clustering_distance_cols = "euclidean")
  dev.off()
  anot_viruscol<- list(ViralHost= c(paste(shQuote(seqhost_list),"=" ,virus_color)))
  noqt_col<-noquote(anot_viruscol)
  virus_color<-shQuote(virus_color)
  anot_viruscol<-list(ViralHost=c('activated sludge metagenome' = '#FF0000' ,
  'Paenibacillus larvae' = '#FF3D00' ,
  'Escherichia coli' = '#FF7A00' ))



####Genus

  seqgenus_list<-unique(seqid_genus$Genus)
  virus_color<-rainbow(43)


  #clinAD
  anot<- data.frame(ClinAD=as.character(tmp_meta$clinAD))
  rownames(anot)<-colnames(vdf_log_half)
  vdfhalf<-data.frame(Accession=rownames(vdf_log_half))
  genus<- vdfhalf %>% left_join(seqid_genus)
  anotrow<- data.frame(ViralGenus=as.character(genus$Genus))
  #anotcolors<- list(ViralHost=virus_color)
  rownames(anotrow)<- vdfhalf$Accession

  pdf(("Heatmap_viralRNA_ROSMAP_clinAD_logHalf1_Genus.pdf"),onefile=T, height = 20, width = 25)
  pheatmap(mat=vdf_log_half,color=colorRampPalette(c("blue","white","red"))(100), annotation_col = anot, annotation_colors = anot_viruscol, annotation_row = anotrow , scale='row',show_rowname=T,show_colname=F, border_color= "grey",na_col="grey", clustering_distance_cols = "euclidean")
  dev.off()
  virus_color<-shQuote(virus_color)
  anot_viruscol<- list(ViralGenus= c(paste(shQuote(seqgenus_list),"=" ,virus_color)))
  noqt_col<-noquote(anot_viruscol)

  anot_viruscol<-list(ViralGenus=c( 'Decadevirus' = '#FF0000',
  'Nehpavirus' = '#FF2400',
  'Halcyonevirus' = '#FF4700',
  'Bievrevirus' = '#FF6B00',
  'Jouyvirus' = '#FF8E00'))

##############################

##t-test and logfc: braak score
#/////////////////////////////////////////////////

B0<- vdf_log_full[,grep(pattern= "_0Bsc|_1Bsc|_2Bsc",colnames(vdf_log_full))]
B6<- vdf_log_full[,grep(pattern= "_6Bsc|_5Bsc|_4Bsc",colnames(vdf_log_full))]

vs_df<-data.frame(rownames(B0))
vs_df$pval<-NA
vs_df$logfc<-NA

for (j in 1:561) {
  viralstat<- t.test(B6[j,],B0[j,]) 
  vs_df[j,2]<- viralstat$p.value
  ADavg<-rowMeans(B6[j,])
  Contavg<-rowMeans(B0[j,])
  L2FC<- ADavg - Contavg
  vs_df[j,3]<- L2FC
}
 padj<- p.adjust(vs_df$pval, method = "BH")
 vs_df[,4]<-padj

write.xlsx(vs_df, file.path("ViralSeq_RNA_Stats_log1_B456_B012_newseq.xlsx"), append = TRUE)

##Volcano plot
#/////////////////////////////////////////////////


pdf(("volcanoplot_0011_count_newfasta_adRI_78rem_lfc0.5.pdf"),onefile=T, height = 12, width = 10)
EnhancedVolcano(vs_df,lab= vs_df$rownames.AD_1.,x= "logfc",y="padj",pCutoff = 0.05, FCcutoff = 0.5,xlim= c(-1,1),ylim=c(0,3.5),pointSize=5,drawConnectors = TRUE)
dev.off()


#ignore
#/////////////////////////////////////////////////

tempDF <- melt(v_meta_sub, gene.vars = "braaksc")
colnames(tempDF) <- c("Protein","Sample","Abundance")
head(tempDF)

tpos<-t(v_meta_sub)
tpos<-data.frame(tpos)
  colnames(rosmap_meta)[colnames(rosmap_meta) == 'rna_id'] <-"specimenID"
colnames(tmp_meta)
tm<-sl %>% left_join(tmp_meta)
t<- rosmap_meta %>% left_join(tm)
h <- tmp_meta[complete.cases(tmp_meta$specimenID), ]

AD<- t[grep(pattern= "1",t$clinAD),]
Cont<- t[grep(pattern= "0",t$clinAD),]
p<- ggplot(AD, aes(x = fu_year> 5, y = bisq_tbi, color = specimenID)) +
  geom_line()  + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)

p + theme(legend.position = "none")

AD_year<- AD[AD$fu_year < 6,]
Cont_year<- Cont[Cont$fu_year < 6,]
c<- ggplot(AD_year, aes(x = fu_year, y = cogn_ep, color = specimenID)) +
  geom_line()  + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)

c + theme(legend.position = "none")
Cont_year$fu_year

r<-tmp_meta[complete.cases(tmp_meta$bisq_tbi),]
Cont_year$cogn_ep

#setwd('F:/loqsproject_compile/GSE11086/output/Counts/')
#count_data_raw=read.csv('count.txt',sep = '\t',header = FALSE)
#Unmapped_reads=read.csv('unmapped_count.txt',sep = '\t',header = FALSE)

do_something <- function(mapped_count, unmapped_count, out_path, next_image, another_image) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  
  count_data_raw=read.csv(mapped_count,sep = '\t',header = FALSE)
  Unmapped_reads=read.csv(unmapped_count,sep = '\t',header = FALSE)
  
  colnames(count_data_raw) <- c("File","counts")
  count_data <- count_data_raw %>%
  separate(File, c("Lib", "Index"), "_")
  count_data$Index <-gsub("[^0-9]","",count_data$Index)

  attach(count_data)
  count_data$sRNA_Type[as.integer(Index)==0]<-'Total_reads_to_dm6'
  count_data$sRNA_Type[as.integer(Index)>=1 & as.integer(Index)<=7]<-'Background'
  count_data$sRNA_Type[as.integer(Index)==8]<-'miRNA'
  count_data$sRNA_Type[as.integer(Index)==9]<-'hpRNA'
  count_data$sRNA_Type[as.integer(Index)>=10 & as.integer(Index)<=11]<-'Retrotransposon'
  count_data$sRNA_Type[as.integer(Index)==12]<-'DNA Transposon'
  count_data$sRNA_Type[as.integer(Index)>=13 & as.integer(Index)<=19] <- 'Others'
  count_data$sRNA_Type[as.integer(Index)>=21 & as.integer(Index)<=25]<-'OthersiRNA'
  count_data$sRNA_Type[as.integer(Index)>=26]<-'NewCisNat'
  detach(count_data)
  #count_data

  count_data=aggregate(counts~Lib+sRNA_Type,count_data,sum)
  #head(count_data)

  mapped_reads<-subset(count_data,count_data$sRNA_Type=='Total_reads_to_dm6')

  sRNA_count_data<-subset(count_data,count_data$sRNA_Type!='Total_reads_to_dm6')

  #Add_contitions (to customise)
  #attach(sRNA_count_data)
  #sRNA_count_data$Condition[as.character(sRNA_count_data$Lib)=='B133'|as.character(sRNA_count_data$Lib)=='B137'|as.character(sRNA_count_data$Lib)=='B141']<-'Gfp-control'
  #sRNA_count_data$Condition[as.character(sRNA_count_data$Lib)=='B134'|as.character(sRNA_count_data$Lib)=='B138'|as.character(sRNA_count_data$Lib)=='B142']<-'Loqs_PA'
  #sRNA_count_data$Condition[as.character(sRNA_count_data$Lib)=='B135'|as.character(sRNA_count_data$Lib)=='B139'|as.character(sRNA_count_data$Lib)=='B143']<-'Loqs_PB'
  #sRNA_count_data$Condition[as.character(sRNA_count_data$Lib)=='B136'|as.character(sRNA_count_data$Lib)=='B140'|as.character(sRNA_count_data$Lib)=='B144']<-'Loqs_PD'
  #detach(sRNA_count_data)

  #1st plot

  graph_count_data <- plyr::ddply(sRNA_count_data, c("Lib","sRNA_Type"), summarise,count=counts)

  sRNA_Type_order=list("Background","miRNA","hpRNA","Retrotransposon","DNA Transposon","Others","OthersiRNA","NewCisNat")
  graph_count_data$sRNA_Type <- factor(graph_count_data$sRNA_Type, levels = sRNA_Type_order)
  count_1 <- ggplot(graph_count_data, aes(x=Lib, y=count, fill=Lib)) + 
    geom_bar(position=position_dodge(), stat="identity")+
    facet_grid(sRNA_Type ~ ., scales = 'free_y', shrink = TRUE)+ 
    theme(strip.text.y = element_text(angle = 0))
  count_1
  ggsave(out_path,count_1)

  #2nd plot

  colnames(Unmapped_reads) <- c("File","unmap counts")
  unmap_count_data <- Unmapped_reads %>%
    separate(File, c("Lib"), "_")

  new_stuff <- merge(mapped_reads,unmap_count_data, by = 'Lib')

  new_stuff$map_ratio=new_stuff$counts/(new_stuff$counts+new_stuff$`unmap counts`)
  new_stuff$map_ratio
  png(filename = next_image)
  plot<-par(mfrow=c(1,1))
  MappedPercentage_plot=barplot(new_stuff$map_ratio, names.arg = new_stuff$Lib, main = 'Percentage of reads mapped',las=2,col=c(4,2,7,'darkviolet'))
  dev.off()

  png(filename = another_image)
  lib_list <- unique(count_data$Lib)

  some_dim_1=plyr::round_any(sqrt(length(lib_list)),1, f=ceiling)
  plot<-par(mfrow=c(some_dim_1,some_dim_1))

  head(count_data)
  for (i in lib_list){
    percent=subset(count_data,Lib==i)
    reads_for_lib=percent[which(percent$sRNA_Type=="Total_reads_to_dm6"),3]
    percent[4]=percent[3]/reads_for_lib
    plot <- barplot(percent$counts.1,names.arg = percent$sRNA_Type, main =paste(i,"% mapped reads"),las=2)
  }
  dev.off()
}
do_something(snakemake@input[[1]],snakemake@input[[2]], snakemake@output[[1]],snakemake@output[[2]],snakemake@output[[3]])

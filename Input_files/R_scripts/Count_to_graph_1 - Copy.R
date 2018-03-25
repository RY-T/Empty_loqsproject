#setwd('F:/loqsproject_compile/RNA_LIB_2/output/Counts/')
#count_data_raw=read.csv('count.txt',sep = '\t',header = FALSE)
#spikein_data_raw=read.csv('../../../Spike_in_folder/Output/Spike_in_count.txt',sep = '\t',header = FALSE)
#Unmapped_reads=read.csv('unmapped_count.txt',sep = '\t',header = FALSE)

#plot 1 adapted from 
#http://www.cookbook-r.com/Graphs/Bar_and_line_graphs_(ggplot2)/
#file1='count.txt'
#file2='../../../Spike_in_folder/Output/Spike_in_count.txt'
#file3='unmapped_count.txt'
#out1='count2.png'
#out2='count1.png'
#out3='read_coverage.png'
#out4='type_percentage.png'

#do_something(file1,file2,file3,out1,out2,out3,out4)

do_something <- function(mapped_count,spike_in_norm, unmapped_count, out_path, out_path_1, next_image, another_image) {
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  
  count_data_raw=read.csv(mapped_count,sep = '\t',header = FALSE)
  spikein_data_raw=read.csv(spike_in_norm,sep = '\t',header = FALSE)
  Unmapped_reads=read.csv(unmapped_count,sep = '\t',header = FALSE)
  
  colnames(count_data_raw) <- c("File","counts")
  count_data <- count_data_raw %>%
    separate(File, c("Lib", "Index"), "_")
  count_data$Index <-gsub("[^0-9]","",count_data$Index)
  
  #attaching sRNA Type details
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
  count_data
  
  #removing B145
  count_data<-subset(count_data,count_data$Lib!='B145')
  
  #collasping reads based off sRNA_types
  count_data=aggregate(counts~Lib+sRNA_Type,count_data,sum)
  #head(count_data)
  
  #Processing Spike in data 
  colnames(spikein_data_raw) <- c("Lib","Spike_counts")
  spikein_data<- subset(spikein_data_raw,spikein_data_raw$Lib!="B145.spike.bowtie.txt")
  spikein_data
  spikein_data$Lib <-gsub(".spike.bowtie.txt","",spikein_data$Lib)
  
  #merge count with spikein count
  merged_count_data=merge(count_data,spikein_data,by='Lib')
  
  #calculate % read mapped
  mapped_reads<-subset(merged_count_data,merged_count_data$sRNA_Type=='Total_reads_to_dm6')
  merged_count_data<-subset(merged_count_data,merged_count_data$sRNA_Type!='Total_reads_to_dm6')
  
  #performing spikein normalisation per 1000 reads
  merged_count_data$Norm_counts <- merged_count_data$counts/merged_count_data$Spike_counts*1000
  
  norm_count_data <- subset(merged_count_data, select= -c(counts,Spike_counts) )
  #norm_count_data
  
  #labelling conditions based on lib type
  attach(norm_count_data)
  norm_count_data$Condition[as.character(norm_count_data$Lib)=='B133'|as.character(norm_count_data$Lib)=='B137'|as.character(norm_count_data$Lib)=='B141']<-'Gfp_control'
  norm_count_data$Condition[as.character(norm_count_data$Lib)=='B134'|as.character(norm_count_data$Lib)=='B138'|as.character(norm_count_data$Lib)=='B142']<-'Loqs_PA'
  norm_count_data$Condition[as.character(norm_count_data$Lib)=='B135'|as.character(norm_count_data$Lib)=='B139'|as.character(norm_count_data$Lib)=='B143']<-'Loqs_PB'
  norm_count_data$Condition[as.character(norm_count_data$Lib)=='B136'|as.character(norm_count_data$Lib)=='B140'|as.character(norm_count_data$Lib)=='B144']<-'Loqs_PD'
  detach(norm_count_data)
  head(norm_count_data)
  
  
  norm_count_data
  detach(package:plyr)
  cdata <- norm_count_data %>%
    group_by(Condition, sRNA_Type)%>%
      summarise('mean_norm_count'= mean(Norm_counts))
  library(plyr)
  #?group_by()
  gfp_count_data <- subset(cdata,cdata$Condition=='Gfp_control')
  
  norm_count_graph1 <- inner_join(norm_count_data, gfp_count_data, by = c("sRNA_Type"="sRNA_Type"),copy=TRUE)
  #norm_count_graph1
  norm_count_graph1$fold_change=norm_count_graph1$Norm_counts/norm_count_graph1$mean
  
  #http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
  
  colnames(norm_count_graph1) <- c('Lib','sRNA_Type','A','Condition','B','C','Fold_change')
  fc_count_data <- subset(norm_count_graph1, select= -c(A,B,C,Lib) )
  fc_count_data
  
  
  fc_data <- ddply(fc_count_data, c("Condition", "sRNA_Type"), summarise,
                   N    = length(Fold_change),
                   mean = mean(Fold_change),
                   sd   = sd(Fold_change),
                   se   = sd / sqrt(N))
  
  #head(fc_data)
  #fc_data
  
  count_2 <- ggplot(fc_data, aes(x=sRNA_Type, y=mean, fill=Condition)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                  width=.3,                    
                  position=position_dodge(.9))+
    theme(axis.text.x=element_text(angle=45, hjust = 1))+
    labs(title="Fold change of sRNA classes across loqs-isoform rescue", x='sRNA class',y='Mean Fold Change')
  #count_2
  ggsave(out_path_1, count_2)
  
  #1st plot
  
  #?ddply
  graph_count_data <- plyr::ddply(norm_count_data, c("Lib","sRNA_Type","Condition"), summarise,count=Norm_counts)
  
  #adding factor for the graph output
  graph_count_data$Lib <- factor(graph_count_data$Lib, levels = graph_count_data$Lib[order(graph_count_data$Condition)])
  #changing the orddering
  sRNA_Type_order=list("Background","miRNA","hpRNA","Retrotransposon","DNA Transposon","Others","OthersiRNA","NewCisNat")
  #changing the orddering
  graph_count_data$sRNA_Type <- factor(graph_count_data$sRNA_Type, levels = sRNA_Type_order)
  
  count_1 <- ggplot(graph_count_data, aes(x=Lib, y=count, fill=Condition)) + 
    geom_bar(position=position_dodge(), stat="identity")+
    facet_grid(sRNA_Type ~ ., scales = 'free_y', shrink = TRUE)+ 
    theme(strip.text.y = element_text(angle = 0))+
    labs(title="Read counts for sRNA classes across loqs-isoform rescue", x='sRNA Library',y='Counts per thousand spike-in')
  count_1
  ggsave(out_path,count_1)
  
  
  #?ggplot
  #2nd plot
  Unmapped_reads$V1<-substring(Unmapped_reads$V1,0,4)
  Unmapped_reads<-subset(Unmapped_reads,Unmapped_reads$V1!='B145')
  Unmapped_reads
  
  mapped_reads
  Unmapped_reads
  new_stuff <- cbind(mapped_reads,Unmapped_reads[2])
  new_stuff
  colnames(new_stuff)[5] <- 'unmap'
  new_stuff$ratio=new_stuff$counts/(new_stuff$counts+new_stuff$unmap)
  new_stuff$ratio
  png(filename = next_image)
  plot<-par(mfrow=c(1,1))
  MappedPercentage_plot=barplot(new_stuff$ratio, names.arg = new_stuff$Lib, main = 'Percentage of reads mapped',las=2,col=c(4,2,7,'darkviolet'))
  dev.off()
  
  #3rd plot
  png(filename = another_image)
  plot<-par(mfrow=c(3,4))
  
  lib_list <- unique(count_data$Lib)
  lib_list
  head(count_data)
  for (i in lib_list){
    percent=subset(count_data,Lib==i)
    reads_for_lib=percent[which(percent$sRNA_Type=="Total_reads_to_dm6"),3]
    percent[4]=percent[3]/reads_for_lib
    plot <- barplot(percent$counts.1,names.arg = percent$sRNA_Type, main =paste(i,"% mapped reads"),las=2)
  }
}
do_something(snakemake@input[[1]],snakemake@input[[2]], snakemake@input[[3]], snakemake@output[[1]],snakemake@output[[2]],snakemake@output[[3]],snakemake@output[[4]])

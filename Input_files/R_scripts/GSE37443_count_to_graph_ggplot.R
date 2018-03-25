#setwd('F:/loqsproject_compile/GSE37443/output/Counts/')
#count_data_raw=read.csv('count.txt',sep = '\t',header = FALSE)
#Unmapped_reads=read.csv('unmapped_count.txt',sep = '\t',header = FALSE)

do_something <- function(mapped_count, unmapped_count, out_path, next_image, another_image) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(cowplot)
  
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
  head(mapped_reads)
  head(sRNA_count_data)
  data 
  total_reads_per_lib <- data.frame(mapped_reads[1],mapped_reads[3])
  total_reads_per_lib
  colnames(total_reads_per_lib) <- c('Lib','total_counts')
  data <- full_join(sRNA_count_data,total_reads_per_lib, by = 'Lib')
  sRNA_count_data <- data %>%
    mutate(cpm=counts/total_counts*1000000)
  

  graph_count_data <- plyr::ddply(sRNA_count_data, c("Lib","sRNA_Type"), summarise,count=cpm)

  sRNA_Type_order=list("Background","miRNA","hpRNA","Retrotransposon","DNA Transposon","Others","OthersiRNA","NewCisNat")
  graph_count_data$sRNA_Type <- factor(graph_count_data$sRNA_Type, levels = sRNA_Type_order)
  
  
  graph_count_data$Lib_no <-gsub("GSM","",graph_count_data$Lib)
  graph_count_data_wIR <- subset.data.frame(graph_count_data,as.numeric(graph_count_data$Lib_no)<919410)
  sRNA_library_types_wIR <- unique(graph_count_data_wIR$Lib)
  graph_count_data_w1118 <- subset.data.frame(graph_count_data,as.numeric(graph_count_data$Lib_no)>=919410)
  sRNA_library_types_w118 <- unique(graph_count_data_w1118$Lib)
  count_1 <- ggplot(graph_count_data_wIR, aes(x=Lib, y=count, fill=Lib)) + 
    geom_bar(position=position_dodge(), stat="identity")+
    facet_grid(sRNA_Type ~ ., scales = 'free_y', shrink = TRUE)+ 
    theme(legend.position="none",strip.text.y = element_text(size = 8, angle=0), axis.text.x=element_text(size = 8,angle=45, hjust = 1)) +
    labs(title='Fly head', x='',y='cpm')+
    scale_fill_discrete(name="loqs-isoform rescue",breaks=sRNA_library_types_wIR,
                        labels=c("wt","loqs-PA","loqs-PB","loqs-PA+loqs-PB","loqs-PA+loqs-PD","loqs-PB+loqs-PD",
                                 "loqs-PA+loqs-PB+loqs-PD","loqs-null control","loqs-ko/loqs-ko","lqos-ko/het" ))
  
  count_2 <- ggplot(graph_count_data_w1118, aes(x=Lib, y=count, fill=Lib)) + 
    geom_bar(position=position_dodge(), stat="identity")+
    facet_grid(sRNA_Type ~ ., scales = 'free_y', shrink = TRUE)+ 
    theme(strip.text.y = element_text(angle = 0, size = 8), axis.text.x=element_text(size = 8,angle=45, hjust = 1)) +
    labs(title='Fly Ovary', x='',y='cpm')+
    scale_fill_discrete(name="loqs-isoform rescue",breaks=sRNA_library_types_w118,
                        labels=c("wt","loqs-PA","loqs-PB","loqs-PA+loqs-PB","loqs-PA+loqs-PD","loqs-PB+loqs-PD",
                                 "loqs-PA+loqs-PB+loqs-PD","loqs-null control","loqs-ko/loqs-ko","lqos-ko/het" ))
  count_2
  colnames(Unmapped_reads) <- c("File","unmap counts")
  unmap_count_data <- Unmapped_reads %>%
    separate(File, c("Lib"), "_")
  
  new_stuff <- merge(mapped_reads,unmap_count_data, by = 'Lib')
  
  new_stuff$map_ratio=new_stuff$counts/(new_stuff$counts+new_stuff$`unmap counts`)
  
  total_count_plot <- ggplot(new_stuff, aes(x=Lib, y=counts, fill=rep(sRNA_library_types_wIR,2))) + 
    geom_bar(position=position_dodge(), stat="identity")+
    theme(legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title=element_text(size=10))+
    labs(title = "Total read counts mapped to dm6", x='',y='')
  
  percentage_mapped_plot <- ggplot(new_stuff, aes(x=Lib, y=map_ratio, fill=rep(sRNA_library_types_wIR,2))) + 
    geom_bar(position=position_dodge(), stat="identity")+
    theme(legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title=element_text(size=10))+
    labs(title = "Percentage of reads mapped to dm6", x='',y='' )
  
  
  multi_count <- ggdraw()+
    draw_plot(count_1, 0, 0.20,0.38,0.80)+
    draw_plot(count_2, 0.38, 0.20,0.62,0.80)+
    draw_plot(total_count_plot,0,0,0.4,0.20)+
    draw_plot(percentage_mapped_plot,0.5,0,0.4,0.20)+
    draw_plot_label(c("A", "B","C","D"),x=c(0,0.38,0,0.5),y=c(1, 1, 0.21,0.21) )
  
  multi_count
  ?plot_grid
  save_plot(out_path,multi_count)
  

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

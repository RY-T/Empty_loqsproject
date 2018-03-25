#setwd('F:/loqsproject_compile/GSE17171/output/Counts/')
#count_data_raw=read.csv('count.txt',sep = '\t',header = FALSE)
#xUnmapped_reads=read.csv('unmapped_count.txt',sep = '\t',header = FALSE)

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
  
  total_reads_per_lib <- data.frame(mapped_reads[1],mapped_reads[3])
  total_reads_per_lib
  colnames(total_reads_per_lib) <- c('Lib','total_counts')
  data <- full_join(sRNA_count_data,total_reads_per_lib, by = 'Lib')
  sRNA_count_data <- data %>%
    mutate(cpm=counts/total_counts*1000000)
  
  graph_count_data <- plyr::ddply(sRNA_count_data, c("Lib","sRNA_Type"), summarise,count=cpm)
  sRNA_library_types <- unique(graph_count_data$Lib)
  sRNA_Type_order=list("Background","miRNA","hpRNA","Retrotransposon","DNA Transposon","Others","OthersiRNA","NewCisNat")
  graph_count_data$sRNA_Type <- factor(graph_count_data$sRNA_Type, levels = sRNA_Type_order)
  count_1 <- ggplot(graph_count_data, aes(x=Lib, y=count, fill=Lib)) + 
    geom_bar(position=position_dodge(), stat="identity")+
    facet_grid(sRNA_Type ~ ., scales = 'free_y', shrink = TRUE)+ 
    theme(strip.text.y = element_text(angle = 0), axis.text.x=element_text(size = 8,angle=45, hjust = 1)) +
    labs(x='',y='cpm')+
    scale_fill_discrete(name="Knockdown Condition",breaks=sRNA_library_types,
                        labels=c("Untreated","Dcr-1_knockdown","Dcr-2_knockdown","Loqs-ORF_knockdown","r2d2_knockdown","Lac Z_knockdown"))
  
  count_1
  
  colnames(Unmapped_reads) <- c("File","unmap counts")
  unmap_count_data <- Unmapped_reads %>%
    separate(File, c("Lib"), "_")
  
  new_stuff <- merge(mapped_reads,unmap_count_data, by = 'Lib')
  
  new_stuff$map_ratio=new_stuff$counts/(new_stuff$counts+new_stuff$`unmap counts`)
  
  total_count_plot <- ggplot(new_stuff, aes(x=Lib, y=counts, fill=Lib)) + 
    geom_bar(position=position_dodge(), stat="identity")+
    theme(legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title=element_text(size=10))+
    labs(title = "Total read counts mapped to dm6", x='',y='')
  
  percentage_mapped_plot <- ggplot(new_stuff, aes(x=Lib, y=map_ratio, fill=Lib)) + 
    geom_bar(position=position_dodge(), stat="identity")+
    theme(legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title=element_text(size=10))+
    labs(title = "Percentage of reads mapped to dm6", x='',y='' )
  
  
  multi_count <- ggdraw()+
    draw_plot(count_1, 0, 0.25,1,0.75)+
    draw_plot(total_count_plot,0,0,0.4,0.25)+
    draw_plot(percentage_mapped_plot,0.5,0,0.4,0.25)+
    draw_plot_label(c("A", "B","C"),x=c(0, 0,0.5),y=c(1, 0.29,0.29) )
  #multi_count
  
  ggsave(out_path,multi_count)
  
  #2nd plot
  
  png(filename = next_image)
  plot<-par(mfrow=c(1,1))
  MappedPercentage_plot=barplot(new_stuff$map_ratio, names.arg = new_stuff$Lib, main = 'Percentage of reads mapped',las=2,col=c(4,2,7,'darkviolet'))
  new_stuff
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

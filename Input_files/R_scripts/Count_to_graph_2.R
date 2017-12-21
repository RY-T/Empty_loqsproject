do_something <- function(data_path, out_path, spike_in_norm) {
  d=read.csv(file =data_path,sep = '\t',header = FALSE)
  d$V3<-d$V1
  d$V3<-substring(d$V3,0,4)
  d$V4<-d$V1
  d$V4<-substring(d$V4,5)
  d$V4<-gsub("[^0-9]", "",d$V4)
  d <- d[-c(1)]
  d
  attach(d)
  d$V5[as.integer(V4)>=1 & as.integer(V4)<=7]<-'background'
  d$V5[as.integer(V4)==8]<-'miRNA'
  d$V5[as.integer(V4)==9]<-'hpRNA'
  d$V5[as.integer(V4)>=10 & as.integer(V4)<=11]<-'retrotransposon'
  d$V5[as.integer(V4)==12]<-'DNA Transposon'
  d$V5[as.integer(V4)>=13 & as.integer(V4)<=19] <- 'Others'
  d$V5[as.integer(V4)>=21 & as.integer(V4)<=25]<-'OthersiRNA'
  d$V5[as.integer(V4)>=26]<-'newCisNat'
  detach(d)
  
  #smaller <- subset(d,d$V5!='retrotransposon')
  #smaller<-subset(smaller,smaller$V3!='B145')
  #smaller
  #collapse=aggregate(V2~V5+V3,smaller,sum)
  #plot<-barplot(collapse$V2, names.arg = collapse$V3, col=c(2:5), legend=unique(collapse$V5))
  
  d
  smaller<-subset(d,d$V3!='B145')
  collapse=aggregate(V2~V5+V3,smaller,sum)
  collapse
  e=read.csv(file =spike_in_norm,sep = '\t',header = FALSE)
  e$V3<-e$V1
  e$V3<-substring(e$V3,0,4)
  e$V4<-e$V1
  e$V4<-substring(e$V4,5)
  e$V4<-gsub("[^0-9]", "",e$V4)
  e <- e[-c(1)]
  e
  
  f=merge(collapse,e,by='V3')
  f[5]=f[3]/f[4]
  f
  collapse=cbind(f[1],f[2],f[5])
  gfp_lib=c('B133','B137','B141')
  gfp_lib
  gfp<-subset(collapse, V3=='B133'|V3=='B137'|V3=='B141')
  gfp<-aggregate(V4~V5,gfp,sum)
  gfp
  
  loqs_pa<-subset(collapse,V3=='B134' | V3=='B138'| V3=='B142')
  loqs_pa<-aggregate(V4~V5,loqs_pa,sum)
  loqs_pa<-merge(loqs_pa,gfp,by='V5')
  loqs_pa$V4=loqs_pa$V4.x/loqs_pa$V4.y
  loqs_pa=cbind(loqs_pa[1],loqs_pa[4])
  colnames(loqs_pa)[2]<-'loqs_pa'
  loqs_pa
  
  loqs_pb<-subset(collapse,V3=='B135' | V3=='B139'| V3=='B143')
  loqs_pb<-aggregate(V4~V5,loqs_pb,sum)
  loqs_pb<-merge(loqs_pb,gfp,by='V5')
  loqs_pb$V4=loqs_pb$V4.x/loqs_pb$V4.y
  loqs_pb=cbind(loqs_pb[1],loqs_pb[4])
  colnames(loqs_pb)[2]<-'loqs_pb'
  loqs_pb
  
  loqs_pd<-subset(collapse,V3=='B136'| V3=='B140'| V3=='B144')
  loqs_pd<-aggregate(V4~V5,loqs_pd,sum)
  loqs_pd<-merge(loqs_pd,gfp,by='V5')
  loqs_pd$V4=loqs_pd$V4.x/loqs_pd$V4.y
  loqs_pd=cbind(loqs_pd[1],loqs_pd[4])
  colnames(loqs_pd)[2]<-'loqs_pd'
  loqs_pd
  
  final_plot<-merge(loqs_pa,loqs_pb, by='V5')
  final_plot<-merge(final_plot,loqs_pd, by='V5')
  row.names(final_plot)<-final_plot$V5
  final_plot<-final_plot[-1]
  final_plot<-t(final_plot)
  final_plot
  png(filename = out_path)
  plot<-par(mfrow=c(1,1))
  par(mfrow=c(1,1))
  barplot(final_plot,names.arg = colnames(final_plot),col=c(4,2,3),beside=TRUE)
  legend("topright", fill=c(4,2,3), legend=row.names(final_plot))
  dev.off()
}

do_something(snakemake@input[[1]], snakemake@output[[1]], snakemake@input[[2]])
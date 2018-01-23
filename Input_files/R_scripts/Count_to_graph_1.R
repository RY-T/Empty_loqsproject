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
  
  retrotransposon=subset(collapse,collapse$V5=='retrotransposon')
  mirna=subset(collapse,collapse$V5=='miRNA')
  hprna=subset(collapse,collapse$V5=='hpRNA')
  Dna_trans=subset(collapse,collapse$V5=='DNA Transposon')
  Others=subset(collapse,collapse$V5=='Others')
  background=subset(collapse,collapse$V5=='background')
  OthersiRNA=subset(collapse,collapse$V5=='OthersiRNA')
  NewCisNat=subset(collapse,collapse$V5=='newCisNat')
  
  png(filename = out_path)
  plot<-par(mfrow=c(4,2))
  backgroundplot<-barplot(background$V4, names.arg = background$V3, main = 'background',las=2,col=c(1,2,3,7))
  mirnaplot<-barplot(mirna$V4, names.arg = mirna$V3, main = 'mirna',las=2,col=c(1,2,3,7))
  hprnaplot<-barplot(hprna$V4, names.arg = hprna$V3, main = 'hprna',las=2,col=c(1,2,3,7))
  retroplot<-barplot(retrotransposon$V4, names.arg = retrotransposon$V3, main = 'retrotransposon',las=2,col=c(1,2,3,7))
  Dna_trans_plot=barplot(Dna_trans$V4, names.arg = Dna_trans$V3, main = 'DNA Transposon',las=2,col=c(1,2,3,7))
  Others_plot=barplot(Others$V4, names.arg = Others$V3, main = 'Others',las=2,col=c(1,2,3,7))
  OthersiRNA_plot=barplot(OthersiRNA$V4, names.arg = OthersiRNA$V3, main = 'OthersiRNA',las=2,col=c(1,2,3,7))
  NewCisNatplot<-barplot(NewCisNat$V4, names.arg = NewCisNat$V3, main = 'NewCisNat',las=2,col=c(1,2,3,7))
  dev.off()
}

do_something(snakemake@input[[1]], snakemake@output[[1]], snakemake@input[[2]])
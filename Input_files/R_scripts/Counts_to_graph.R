setwd("E:/loqsproject_compile/Seq_map_directory_M/")

d=read.csv(file ='count.txt',sep = '\t',header = FALSE)
d$V3<-d$V1
d$V3<-substring(d$V3,0,4)
d$V4<-d$V1
d$V4<-substring(d$V4,5)
d$V4<-gsub("[^0-9]", "",d$V4)
d <- d[-c(1)]
d
attach(d)
d$V5[V4>=1 & V4<=7]<-'background'
d$V5[V4==8]<-'miRNA'
d$V5[V4==9]<-'hpRNA'
d$V5[V4>=10 & V4<=11]<-'retrotransposon'
d$V5[V4==12]<-'DNA Transposon'
d$V5[is.element(V4,13:19)] <- 'Others'
detach(d)
d

smaller<-subset(d,d$V3!='B145')
collapse=aggregate(V2~V5+V3,smaller,sum)
collapse

e=read.csv(file ='./../Spike_in_folder/Output/Spike_in_count.txt',sep = '\t',header = FALSE)
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
collapse
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
par(mfrow=c(1,1))
barplot(final_plot,names.arg = colnames(final_plot),col=c(4,2,3),beside=TRUE)
legend("topright", fill=c(4,2,3), legend=row.names(final_plot))

colnames(final_plot)
retrotransposon=subset(collapse,collapse$V5=='retrotransposon')
mirna=subset(collapse,collapse$V5=='miRNA')
hprna=subset(collapse,collapse$V5=='hpRNA')
Dna_trans=subset(collapse,collapse$V5=='DNA Transposon')
Others=subset(collapse,collapse$V5=='Others')
background=subset(collapse,collapse$V5=='background')
collapse
png(file= 'something.png')
par(mfrow=c(3,2))
retroplot<-barplot(retrotransposon$V4, names.arg = retrotransposon$V3, main = 'retrotransposon')
mirnaplot<-barplot(mirna$V4, names.arg = mirna$V3, main = 'mirna')
hprnaplot<-barplot(hprna$V4, names.arg = hprna$V3, main = 'hprna')
Dna_trans_plot=barplot(Dna_trans$V4, names.arg = Dna_trans$V3, main = 'DNA Transposon')
Others_plot=barplot(Others$V4, names.arg = Others$V3, main = 'Others')
background_plot=barplot(background$V4, names.arg = background$V3, main = 'background')
#smaller <- subset(d,d$V5!='retrotransposon')
#smaller<-subset(smaller,smaller$V3!='B145')
#smaller
#png('bar_chart.png')
plot<-barplot(collapse$V2, names.arg = collapse$V3, col=c(2:5), legend=unique(collapse$V5))
dev.off()
collapse

barplot(collapse$V2, names.arg = collapse$V3, col=c(2:6), legend=unique(collapse$V5))

help(png)
lib=unique(d$V3)
for(i in 1:length(lib)){
  assign(lib[i],subset(d, d$V3==lib[i]))
}
B133

help("barplot")
library(ggplot2)
library(dplyr)
libs=group_by(smaller, V3,V5)
#https://stackoverflow.com/questions/1660124/how-to-sum-a-variable-by-group
summary=summarise(libs, counts=sum(V2))
summary
ggplot(smaller)

help(ggplot)

help(barplot)
counts <- table(mtcars$vs, unique(d$V3))
barplot(counts, main="Car Distribution by Gears and VS",
        xlab="Number of Gears", col=c("darkblue","red"),
        legend = rownames(counts), beside=TRUE)

subset(B133,)
sum(B133$V2)

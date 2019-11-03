###################################################################################################
###########                                                                             ###########
###########             ABC model selection with abc package (Csillery et al., 2010)    ###########
###########                                                                             ###########
###################################################################################################

library(abc)
library(abc.data)

#Vector of col.names
statcols=c('gr_K_1','gr_mean_K','gr_H_1','gr_mean_H','gr_mean_D',
           'gr_sd_D','mean_FS','sd_FS','gr_Pi_1','gr_mean_Pi',
           'gr_sd_Pi','K_1','K_2','mean_K','H_1','H_2','mean_H',
           'mean_D','sd_D','mean_FS','sd_FS','Pi_1','Pi_2','mean_Pi','sd_Pi')

##Calculating mean of simulated data
##If means are alread calculated, skip this code
ss=list()
ss_mean=list()
n.scen=7
for (i in 1:n.scen) {
  ss[[i]]=list()
  files=list.files(paste('simulations_arlsumstat_results/scen',i,sep=''))
  for (j in 1:3) {
    ss[[i]][[j]]=read.table(paste('simulations_arlsumstat_results/scen',i,'/',files[j],sep=''),header=T,sep='\t')
  }
  value=list()
  for (k in 1:NCOL(ss[[i]][[j]])) {
    var=cbind(ss[[i]][[1]][,k],ss[[i]][[2]][,k],ss[[i]][[3]][,k])
    value[[k]]=apply(var,1,mean)
  }
  data=c()
  for (k in 1:length(value)) {
    data=cbind(data,value[[k]])
  }
  ss_mean[[i]]=data[,-26]
  colnames(ss_mean[[i]])=statcols
  write.table(ss_mean[[i]],paste('simulations_arlsumstat_results/scen',i,'_ss_mean',sep=''),row.names=F,col.names = T,quote=F,sep='\t')
}

##If means are already calculated, read them into the environment
for (i in 1:n.scen) {
  ss_mean[[i]]=read.table(paste('simulations_arlsumstat_results/scen',i,'_ss_mean',sep=''),header=T,sep='\t')
}

##Converting observed arp from DnaSP into DNA format
##I don't remember exactly why I needed to do this, but it seems this was necessary to use arlsumstat to calculate the SS on the observed data.
##In any case, if you reach the point where you have the mean of each SS across all three markers for the observed data, it doesn't matter how
##you got there (as long as you got arlsumstat to work on your files correctly and the files had the correct data).

hap=read.table('ruf_cytb.hap')
hap_n=read.table('hap_n_cytb')
sequence=c()

for (i in 1:NROW(hap_cytb)) {
  a=which(as.character(hap_cytb_n[,1])==as.character(hap_cytb[i,1]))
  sequence=c(sequence,rep(as.character(hap_cytb[i,2]),hap_cytb_n[a,2]))
}

data=data.frame(c(paste('2_',1:42,sep=''),paste('1_',1:16,sep='')),rep(1,NROW(sequence)),sequence)
write.table(data,'observed/data_cytb',row.names=F,col.names=F,quote=F,sep='\t')

hap=read.table('observed/ruf_g3pdh.hap',sep='\t')
hap_n=read.table('observed/hap_n_g3')
sequence=c()

for (i in 1:NROW(hap_n)) {
  a=which(as.character(hap[,1])==as.character(hap_n[i,1]))
  sequence=c(sequence,rep(as.character(hap[a,2]),hap_n[i,2]))
}

data=data.frame(c(paste('2_',1:56,sep=''),paste('1_',1:22,sep='')),rep(1,NROW(sequence)),sequence)
write.table(data,'observed/data_g3pdh',row.names=F,col.names=F,quote=F,sep='\t')

hap=read.table('observed/ruf_tgfb2.hap')
hap_n=read.table('observed/hap_n_tgf')
sequence=c()

for (i in 1:NROW(hap_n)) {
  a=which(as.character(hap[,1])==as.character(hap_n[i,1]))
  sequence=c(sequence,rep(as.character(hap[a,2]),hap_n[i,2]))
}

data=data.frame(c(paste('2_',1:76,sep=''),paste('1_',1:28,sep='')),rep(1,NROW(sequence)),sequence)
write.table(data,'observed/data_tgfb2',row.names=F,col.names=F,quote=F,sep='\t')

##Retrieving observed data and creating observed vector (mean across markers)
##If mean of observed values is already calculated, skip this.
obs=list()
obs_mean=list()
files=list.files(paste('observed',sep=''),pattern='ss.txt')
for (i in 1:3) {
  obs[[i]]=read.table(paste('observed/',files[j],sep=''),header=T,sep='\t')
}
value=list()
for (k in 1:NCOL(obs[[1]])) {
  var=cbind(obs[[1]][,k],obs[[2]][,k],obs[[3]][,k])
  value[[k]]=apply(var,1,mean)
}
data=c()
for (k in 1:length(value)) {
  data=cbind(data,value[[k]])
}
obs_mean=data[,-26]
write.table(obs_mean,'observed/obs_mean',row.names=F,col.names = T,quote=F,sep='\t')

##Retrieving already calculated observed data
obs_mean=read.table('observed/obs_mean',sep='\t')
obs_mean=obs_mean$V1

##Combining all simulations
ssmean=rbind(ss_mean[[1]],ss_mean[[2]],ss_mean[[3]],ss_mean[[4]],ss_mean[[5]],ss_mean[[6]],ss_mean[[7]])
colnames(ssmean)=statcols

#Creating models labels vector
models=c(rep('Scenario 1',1000000),rep('Scenario 2',1000000),rep('Scenario 3',1000000),
         rep('Scenario 4',1000000),rep('Scenario 5',1000000),rep('Scenario 6',1000000),
         rep('Scenario 7',1000000))

##Creating summary statistics label vector
statsname=c("Mean number of alleles over loci output for each population (at group level)",
            "Mean number of alleles over loci and population (at group level)",
            "Mean heterozygosity over loci output for each population (at group level)",
            "Mean heterozygosity over loci and population (at group level)",
            "Tajima's D reported for each population (at group level)",
            "s.d. Tajima's D over all populations (at group level)",
            "Mean Fu's F over all pops",
            "s.d. Fu's F over all pops",
            "Mean number of pairwise differences for northern population",
            "Mean number of pairwise differences for each pop",
            "s.d. of the mean number of pairwise differences over pops",
            "Mean number of alleles over loci output for northern population",
            "Mean number of alleles over loci output for southern population",
            "Mean number of alleles",
            "Mean heterozygosity over loci output for northen population",
            "Mean heterozygosity over loci output for southern population",
            "Mean heterozygosity",
            "Mean Tajima's D",
            "s.d. Tajima's D",
            "Mean Fu's F",
            "s.d. Fu's F",
            "Mean number of pairwise differences of northern population",
            "Mean number of pairwise differences of southern population",
            "Mean number of pairwise differences",
            "s.d. Mean number of pairwise differences")

#Boxplots
library(ggplot2)
for (i in 1:NCOL(ssmean)) {
  data=data.frame(var=ssmean[,i],models)
  ggplot(data=data,aes(x=models,y=var))+geom_boxplot()+
    ggtitle(paste(statsname[i],sep=''))+
    ggsave(paste('results/boxplot_',colnames(ssmean)[i],'.pdf',sep=''))
}

##PCA
pca=list()
for (i in 1:length(ss_mean)) {
  data=rbind(ss_mean[[i]],obs_mean)
  pca[[i]]=prcomp(data)
}
colors=c('blue','darkorange1','red','lightblue','darkviolet','green3','lightpink4','black')
for (i in 1:length(pca)) {
  data=data.frame(group=c(rep('Simulated',1000000),'Observed'),pca[[i]]$x)
  ggplot(data=data,aes(x=PC1,y=PC2,color=group))+geom_point()+
    scale_color_manual(values=c(colors[8],colors[i]))+guides(color=FALSE)+
    ggtitle(paste('PCA for Scenario ',i,sep=''))+
    ggsave(paste('results/pca_scenario',i,'.tiff',sep=''))
}

del=c()
##Cross-validation
##Checking variables variance
for (i in 1:NCOL(ssmean)) {
  value=var(ssmean[1:1000000,i])
  if (value==0) {
    del=rbind(del,c('scen1',i))
  }
  value=var(ssmean[1000001:2000000,i])
  if (value==0) {
    del=rbind(del,c('scen2',i))
  }
  value=var(ssmean[2000001:3000000,i])
  if (value==0) {
    del=rbind(del,c('scen3',i))
  }
  value=var(ssmean[3000001:4000000,i])
  if (value==0) {
    del=rbind(del,c('scen4',i))
  }
  value=var(ssmean[4000001:5000000,i])
  if (value==0) {
    del=rbind(del,c('scen5',i))
  }
  value=var(ssmean[5000001:6000000,i])
  if (value==0) {
    del=rbind(del,c('scen6',i))
  }
  value=var(ssmean[6000001:7000000,i])
  if (value==0) {
    del=rbind(del,c('scen7',i))
  }
}
remove=as.numeric(unique(del[,2]))

ssmean=ssmean[,-remove]
obs_mean=obs_mean[-remove]

cv.modsel <- cv4postpr(models, ssmean, nval=100, tol=.05, method="mnlogistic")
s <- summary(cv.modsel)
plot(cv.modsel, names.arg=c(paste('scen',1:7,sep='')))

modsel <- postpr(obs_mean, models, ssmean, 
                    tol=.05, method="mnlogistic")

summary(modsel)

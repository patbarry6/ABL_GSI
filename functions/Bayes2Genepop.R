####Bayes2Genepop####
library(tidyverse)
library(plyr)

Bayes2Genepop_mixture<-function(infile,outfile,LociFile,digits,SampPrefix="Mix",Project=paste(infile,"Bayes2Genepop",sep="")){

if(file.exists(file.path(infile))){
  dat<-readLines(infile)
  }else{cat("Infle not found! Check to make sure the path to the file is correct.")
    }  

if(file.exists(file.path(LociFile))){
  Loci<-read.csv(LociFile,stringsAsFactors = F,header=F)
  } else{cat("LociFile not found! Check to make sure the path to the file is correct.")
    }

  
nloci<-nrow(Loci)

Genos<-lapply(1:length(dat), function(x) str_split(str_trim(dat[x]),"\\s+"))
GenoMat<-matrix(data=unlist(Genos),nrow=length(dat),ncol=nloci,byrow=T)

GenoMatAll<-matrix(data=NA,nrow=nrow(GenoMat),ncol=nloci*2)
#each locus could have a different number of alleles (if uSat)
#so we should loop over the loci to get each individuals genotype
for (L in 1:nloci){
  temp<-GenoMat[,L]
  alleles<-length(unlist(str_split(temp[1],"")))
  GenoAll<-lapply(1:length(temp),function(Li) rep(seq(1,alleles),as.numeric(unlist(str_split(temp[Li],"")))))
  Zeros<-which(sapply(1:length(GenoAll),function(x) length(GenoAll[[x]])<1))
    if(length(Zeros)>0){
      for (z in 1:length(Zeros)){
      GenoAll[[Zeros[z]]]<-c(0,0)
      }
    }#if we have any ungenotyped individuals
  GenoMatAll[,((L*2-1):(L*2))]<-matrix(data=unlist(GenoAll),nrow=length(GenoAll),ncol=2,byrow=T)
  }

dat2<-lapply(1:nrow(GenoMatAll),function(x) paste(formatC(as.numeric(GenoMatAll[x,seq(1,nloci*2,2)]),width=digits,flag="0"),formatC(as.numeric(GenoMatAll[x,seq(2,nloci*2+1,2)]),width=digits,flag="0"),sep=""))
dat2<-matrix(unlist(dat2),ncol=nloci,nrow=nrow(GenoMatAll),byrow=T)
dat2<-cbind(paste(SampPrefix,seq(1,nrow(GenoMat)),sep="."),",",dat2)

write.table(Project,outfile,quote=F,row.names=F,col.names=F)
write.table(Loci,outfile,quote=F,row.names=F,col.names=F,append=T)
write.table("POP",outfile,quote=F,row.names=F,col.names=F,append=T)
write.table(dat2,outfile,quote=F,row.names=F,col.names=F,append=T)
}

##################################

Bayes2Genepop_baseline<-function(infile,outfile,LociFile,digits,SampPrefix="Mix",Project=paste(infile,"Bayes2Genepop",sep="")){
  #Need a baseline file from Chuck
}

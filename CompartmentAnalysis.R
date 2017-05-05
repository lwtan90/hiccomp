library(wavethresh)

helper1<-function(data){
	a<-c(sign(data),-1,1)
	b1<-(100*(1-(length(which(a[2:length(a)]==a[1]))-1)/(length(a)-3)))
	b2<-table(a[2:length(a)])-1
  return(c(b1,b2))	
}

helper2<-function(data){
	100*min(data)/sum(data)
}


# x: the data table with
    # col1 = an ID
    # col2 = chromosome (e.g. chr1)
    # col3 = start of region
    # col4 = end of region
    # col5 = strand
    # col6 = PC1
# chrom: the chromosome that is analyzed
# B: the number of subsamples to generate to assess data noise
# cut_disagreement: the cut-off for the disagreement method
    # how many subsample-based signs are different from the original's
# cut_instability: the cut-off for the instability method
    # how often the subsample-based signs differ
# by: disagreement or instability
# thresh_type / threshPolicy / threshValue / return.threshold / scale.threshold: parameters for the wavelet decomposition
# plot_trend: whether the wavelet estimated trend will be included in the plot
# runs: if TRUE compartment triplets of the form A0A or B0B are adjusted by converting the 0 to the value of the neighbors

ChangeIndex<-function(x,chrom,B=999,cut_disagreement=5,cut_instability=5,by="Disagreement",threshType="soft",threshPolicy="BayesThresh",threshValue=0,return.threshold = FALSE,scale.threshold=1.5,plot_trend=FALSE,runs=TRUE){

    Alldata<-x[which(x[,2]==chrom),]
	if(nrow(Alldata)==0){
		stop("No data for this chromosome")
	}
    Alldata<-Alldata[sort.list(as.numeric(Alldata[,3])),]
	data<-Alldata[,6]

	pos<-apply(cbind(apply(cbind(as.character(Alldata[,2]),as.numeric(Alldata[,3])),1,paste,collapse=":"),
		as.numeric(Alldata[,4])),1,paste,collapse="-")
	J<-ceiling(log(length(data),2))
	wavdata<-c(data,data[length(data):1])[1:2^J]
	wy <- wd(wavdata)
	if(return.threshold==TRUE){
		thresh <- threshold(wy,type=threshType,policy="universal",value=threshValue,return.threshold=return.threshold)
		thresh <- threshold(wy,type=threshType,policy="manual",value=scale.threshold*thresh)
	} else {
		thresh <- threshold(wy,type=threshType,policy=threshPolicy)
	}
	yr <- wr(thresh)[1:length(data)]
	resids<-data-yr

	test<-yr
	for(b in 1:B){
		bresids<-sample(resids,length(resids),replace=FALSE)
		randomData<-data+bresids
		randomWavData<-c(randomData,randomData[length(randomData):1])[1:2^J]
		wy <- wd(randomWavData)
		thresh <- threshold(wy, type="soft",policy = "BayesThresh")
		test <- cbind(test,wr(thresh)[1:length(randomData)])
	}

	res<-apply(test,1,helper1)
	index<-cbind(res[1,],apply(res[2:3,],2,helper2))

	cols<-rep(2,length(data))
	if(by=="Disagreement"){
		cut<-cut_disagreement
		ii<-index[,1]
		cols[which(ii>cut)]<-1
	} else {
		cut<-cut_instability
		ii<-index[,2]
		cols[which(ii>cut)]<-1
	}
	fpos<-rep("",nrow(index))
	fpos[which(cols==1)]<-pos[which(cols==1)]

	compartments<-rep("A",nrow(index))
	compartments[which(Alldata[,6]<0)]<-"B"
	compartments[which(ii>cut)]<-"0"

    if(runs==TRUE){
        adj1<-paste(compartments,collapse="")
        w<-gregexpr("A0A",adj1)
        if(w[[1]][1]!= -1){
            compartments[as.numeric(w[[1]])+1]<-"A"
        }
        w<-gregexpr("B0B",adj1)
        if(w[[1]][1]!= -1){
            compartments[as.numeric(w[[1]])+1]<-"B"
        }
    }
    
	cols<-rep(2,nrow(index))
	cols[which(Alldata[,6]<0)]<-3
	cols[which(compartments=="0")]<-1

	par(mfrow=c(2,2))
	mainLeg<-unlist(strsplit(chrom,"r"))[2]
	plot(Alldata[,3],data,cex=0.2,col=cols,xlab="sorted HiC regions by genomic (start) location",ylab="PC1",
		main=paste("Analysis of Chromosome ",mainLeg,sep=""),cex.axis=0.8)
	if(plot_trend==TRUE){
		lines(Alldata[,3],yr,col=4,lwd=0.5)
	}
	hist(resids,xlab="Model residuals",sub=paste("KS-test P-value = ",
		round(ks.test(resids,"pnorm",0,sqrt(var(resids)))$p.value,3),sep=""),main=paste("Analysis of Chromosome ",mainLeg,sep=""))
	
	barplot(index[,1],names=fpos,las=2,cex.names=0.1,ylab="Disagreement P",xlab="sorted HiC regions by genomic (start) location",
			main=paste("Analysis of Chromosome ",mainLeg,sep=""))
	barplot(index[,2],names=fpos,las=2,cex.names=0.1,ylab="Instability Index",xlab="sorted HiC regions by genomic (start) location",
			main=paste("Analysis of Chromosome ",mainLeg,sep=""))


	Result<-cbind(Alldata,index,compartments)
	colnames(Result)<-c("peakID","chr","start","end","strand","PC1","Disagreement P","Instability Index","Compartment")

  return(Result)
}







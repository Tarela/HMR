#data<-read.table("~/Desktop/mESC_GSM1562337_CBX7_HMsig.bed",row.names=4,header=T)
data<-read.table("ABL_GSM969562_EZH2_HMsig.bed",row.names=4,header=T)

#rawsig <- data[,"mESC_GSM1954953_H3K9me3"]
rawsig <- data[,"GSM969563_H3K27me3"]
sudo_count <- min(rawsig[which(rawsig>0)])/2
sig <- log10(rawsig + sudo_count)

maxG<-function(cutoff,usedata){
	# function of estimating G given cutoff
	# G : inter-class variance
	# method refer to "Otus' method"
	w0 <- length(which(usedata < cutoff))/length(usedata)
	w1 <- length(which(usedata >= cutoff))/length(usedata)
    u0 <- mean(usedata[which(usedata < cutoff)])
    u1 <- mean(usedata[which(usedata >= cutoff)])
    g <- w0*w1*(u0-u1)**2
    return(g)
}
#print(1)

signal2cutoff <- function(rawsig){
	# main function
	# input data: a vector of signal (HMsignal on CR peak, linear scale)
	# method: go through all possible cutoff, find the cutoff corresponding to maximum G
	# output: the cutoff, a vector of selected cutoff candidates (Gbins) and a vector of G value for each cutoff candidates (G) 
	
	# log transform, sudo_count = min(linear signal)/2
	sudo_count <- min(rawsig[which(rawsig>0)])/2
	sig <- log10(rawsig + sudo_count)

	# separate the section of (log) signal to N cutoffs	
	Gbins <- seq(min(sig)+0.1,max(sig)-0.1,0.01)
	
	# estimate G for each cutoff candidates
	G <- unlist(lapply(Gbins,maxG, sig))
	
	# select cutoff at the first time G meats its maximum value
	NCcut<-seq(min(sig)+0.1,max(sig)-0.1,0.01)[which(G==max(G))]

	group_detail <- rep(0, length(sig))
	group_detail[which(NCcut <= sig)] <- 1
	
	# output the grouping result: list contains 3 items 
	# item1: NC cutoff
	# item2: 2 column for cutoff candidates and corresponded G
	# item3: Nrow = peak number,Ncolumn = 2, c1 for logsignal, c2 for group number, 0 for lowHM group (solo, non-canonical), 1 for highHM group (ensemble, canonical)
	
	return(list(NCcut, cbind(Gbins,G), cbind(sig, group_detail) ))
}

testdata <- data[,"GSM969563_H3K27me3"]
testout <- signal2cutoff(testdata)

### 
hist(testout[[3]][,1],n=200,xlim=c(min(testout[[3]][,1])+0.1,max(testout[[3]][,1])-0.1),main="Otus' method on ABL EZH2",xlab="log10 H3K27me3 signal",ylab="frequency")
mtext("inter-class variance",4)
par(new=T)
plot(testout[[2]],type="l",col="red",xlim=c(min(testout[[3]][,1])+0.1,max(testout[[3]][,1])-0.1),main="",xlab="",ylab="",axes=F)
axis(side=4)
abline(v= testout[[1]],col="blue")
legend("topleft",legend=c("inter-class var","NC cutoff"),lwd=3,bty="n",col=c("red","blue"))



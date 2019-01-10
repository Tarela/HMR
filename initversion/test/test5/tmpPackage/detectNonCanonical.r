
# read parameter

a<-commandArgs(T)

outname <- a[1]
signalname <- a[2]
cutoff <- as.numeric(a[3])
Alpha <- as.numeric(a[4])
LambdaChoice <- (a[5])
topN <- a[6]
R2cutoff <- 0.1
tmp_rpackage_dir <- a[7]

if("foreach" %in% installed.packages()[,"Package"]){
	library(foreach)
}else{
	install.packages("foreach",dependencies=TRUE,lib=tmp_rpackage_dir,repos="https://mirrors.tongji.edu.cn/CRAN/")
	library(foreach,lib.loc=tmp_rpackage_dir)
}

if("glmnet" %in% installed.packages()[,"Package"]){
	library(glmnet)
}else{
	install.packages("glmnet",dependencies=TRUE,lib=tmp_rpackage_dir,repos="https://mirrors.tongji.edu.cn/CRAN/")
	library(glmnet,lib.loc=tmp_rpackage_dir)
}

library(methods)

### function

unilinear_permute <- function(TFov,HMsig){
    lmresult <- summary(lm(HMsig ~ TFov))
    fgR2 <- lmresult$adj.r.squared
    #print(lmresult)
    lmcoeff <- lmresult$coefficients['TFov','Estimate']
    lmP <- lmresult$coefficients['TFov',c('Pr(>|t|)')]
    #if(fgR2 >= R2cutoff){
    bgR2s <- c()
    for(i in 1:999){
        bgTFov <- rep(0,length(TFov))
        bgTFov[sample(length(TFov),length(which(TFov > 0)))] <- 1
        bgR2 <- summary(lm(HMsig ~ bgTFov))$adj.r.squared
        bgR2s <- c(bgR2s,bgR2)
    }
    permuteP <- (length(which(fgR2 < bgR2s))+1)/(length(bgR2s)+1)
    return(c(permuteP,fgR2,lmcoeff))
    #}else{
    #    return(c("NA",fgR2,lmcoeff))
    #}
}

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

signal2cutoff <- function(rawsig){
    # function for separate NC peaks considering HMsignal
    # input data: a vector of signal (HMsignal on CR peak, linear scale)
    # method: go through all possible cutoff, find the cutoff corresponding to maximum G
    # output: the cutoff, a vector of selected cutoff candidates (Gbins) and a vector of G value for each cutoff candidates (G) 
    
    sig <- rawsig
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
    # item3: Nrow = peak number,Ncolumn = 2, c1 for signal, c2 for group number, 0 for lowHM group (solo, non-canonical), 1 for highHM group (ensemble, canonical)
    
    return(list(NCcut, cbind(Gbins,G), cbind(sig, group_detail) ))
}

# step1 data preprocess
peakov <- read.table(paste0(outname,"_peakov.bed"))
signal <- read.table(paste0(outname,"_HMsig.bed"))

candidate_list <- as.vector(read.table(paste0(outname,"_cofactor_candidate_list.txt"))[,1])
bedncol <- ncol(peakov) - length(candidate_list)
colnames(peakov) <- c(paste0("c",seq(1:bedncol)),candidate_list)
colnames(signal) <- c(paste0("c",seq(1:bedncol)),signalname)
rownames(peakov) <- paste0("r",seq(1,nrow(peakov)))
rownames(signal) <- paste0("r",seq(1,nrow(peakov)))

TFov <- peakov[,candidate_list]
TFov[TFov > 1] <- 1
HMsig <- as.matrix(signal[,signalname])
colnames(HMsig) <- signalname
rownames(HMsig) <- rownames(signal)


### step2, prepare X,Y for model selection 
lindata <- cbind(HMsig,TFov)
bind_sum <- apply(TFov,1,sum)
## sites with 90% factors overlapped are excluded
use1_lindata <- lindata[names(bind_sum[which(bind_sum < ncol(TFov)*0.9)]),]
## factors with < 100 cobinding events are excluded
TF_sum <- apply(use1_lindata[,colnames(TFov)],2,sum)
use_lindata <- lindata[names(bind_sum[which(bind_sum < ncol(TFov)*0.9)]),c(colnames(HMsig),names(TF_sum)[which(TF_sum>=100)])]
## form X, Y
Y <- as.matrix(scale(use_lindata[,colnames(HMsig)]))
colnames(Y) <- colnames(HMsig)
X <- as.matrix(use_lindata[,(ncol(HMsig)+1):ncol(use_lindata)])
peakX <- peakov[rownames(X),1:bedncol]
### elastic net model
set.seed(1007)
## model selection, alpha=0.5
#if(ncol(Y) == 1){
        ## uni-substrate case
    #print(Alpha)
    #glmmod = glmnet(x=X, y=Y, family="gaussian", alpha=Alpha, nlambda=200)
    #print("A")
cv.glmmod <- cv.glmnet(x=X,y=Y,alpha=Alpha,family="gaussian")
if (LambdaChoice == "1se"){
    useS <- max(1,which(cv.glmmod$glmnet.fit$lambda == cv.glmmod$lambda.1se)[1])
}else{
    useS <- max(1,which(cv.glmmod$glmnet.fit$lambda == cv.glmmod$lambda.min)[1])
}

coeff_this <- coef(cv.glmmod) 
coTF_name <- names(coeff_this[which(coeff_this[,min(useS,ncol(coeff_this))] != 0),min(useS,ncol(coeff_this))])
coTF_name_use <-  coTF_name[which(coTF_name!="(Intercept)" & coTF_name!="")]
all_use_coTF_names <- coTF_name_use

## selected co-factors
coTF_binding_matrix <- as.matrix(X[,unique(all_use_coTF_names)])

## separate NC/C peaks based on HMsignal using otsu's method
NCotsu <- signal2cutoff(Y)
peakgroup <- NCotsu[[3]][,2]


### elastic net CV plot
pdf(file=paste0(outname,"_elnet_lambdaSelection.pdf"))
plot(cv.glmmod)
if (LambdaChoice == "1se"){
    abline(v=log(cv.glmmod$lambda.1se),col="blue",lwd=2)
    legend("topleft",legend="lambda.1se",lwd=3,bty="n",col='blue')
}else{
    abline(v=log(cv.glmmod$lambda.min),col="blue",lwd=2)
    legend("topleft",legend="lambda.min",lwd=3,bty="n",col='blue')
}
dev.off()


### re-fit by a uni-variate linear model, output pvalue(empirical), r-square and beta coeff for each selected co-factors
noTFdetected = 0

if(ncol(coTF_binding_matrix)==0){
    print("no significant candidates")
    noTFdetected = 1
}else{
    summary_table_raw <- as.matrix(t(apply(coTF_binding_matrix,2,unilinear_permute,Y)))
    colnames(summary_table_raw) <- c("Pvalue","R2","coeff")
    summary_table_tmp <- suppressWarnings(summary_table_raw[which(as.numeric(summary_table_raw[,"Pvalue"]) <= cutoff & as.numeric(summary_table_raw[,"coeff"]) < 0 & as.numeric(summary_table_raw[,"R2"] >= R2cutoff)),])
    summary_table <- summary_table_tmp[order(summary_table_tmp[,"R2"],decreasing=TRUE),]
    if(nrow(summary_table) ==0 ){
        print("no significant candidates")
        noTFdetected = 1
    }
}

if(noTFdetected==1){  
    write.table("no non-classic function detected",file=paste0(outname,"_elnetNC.txt"),quote=F,sep="	",row.names=F,col.names=F)
    write.table("no non-classic function detected",file=paste0(outname,"_filterNC.txt"),quote=F,sep="	",row.names=F,col.names=F)
    write.table("no non-classic function detected",file=paste0(outname,'_NCsummary.txt'),quote=F,sep="	",row.names=F,col.names=F)
}else{
    ### summarize output NC table according to topN, if topN==all, output all coTF that pass cutoff
    if(topN == "all"){
        topN <- nrow(summary_table)
    }else{
        topN <- as.numeric(topN)
    }
    
    out_table_raw <- as.matrix(summary_table[1:topN,])
    
    ### summary and output steps
    ### for each predicted coTF, output the cobinding NC sites as a bed file
    if(!file.exists("nonClassicPeaks/")){dir.create("nonClassicPeaks")}
    
    num_NCsites <- c()
    pdf(file=paste0(outname,"_coTF_HMsignal.pdf"))
    if(topN == 1){
        par(mar=c(4,4,2,2))    
    }else if(topN == 2){
        par(mfrow=c(1,2),mar=c(4,4,2,2))
    }else if(topN == 3){
        par(mfrow=c(2,2),mar=c(4,4,2,2))
    }else if(topN == 4){
        par(mfrow=c(2,2),mar=c(4,4,2,2))
    }else{
        par(mfrow=c(3,2),mar=c(4,4,2,2))   
    }
    
    for(topnum in 1:nrow(out_table_raw)){
        coTF = rownames(out_table_raw)[topnum]
        coTF_cobinding <- X[,coTF]
        cobinding_NC_peak <- peakX[which(coTF_cobinding>0 & peakgroup == 0),]
        write.table(cobinding_NC_peak, file=paste0("nonClassicPeaks/",outname,"_",coTF,"_top",topnum,"nonclassic_peaks.bed"),quote=F,sep="	",row.names=F,col.names=F)
        num_NCsites <- c(num_NCsites, nrow(cobinding_NC_peak))
        if(topnum %in% 1:5){
            boxplot(Y[which(coTF_cobinding>0 & peakgroup == 0)],Y[which(coTF_cobinding==0 | peakgroup > 0)],
                names=c("non-classic peak","classic peak"),ylab=paste0("HM signal"),main=paste0(coTF),
                outline=T,cex.main=1)
            legend("topleft",legend=paste0("#NCpeak = ",nrow(cobinding_NC_peak)),bty="n")
        }
    }
    dev.off()
    
    out_table <- cbind(out_table_raw, num_NCsites)
    
    write.table(summary_table_raw,file=paste0(outname,"_elnetNC.txt"),quote=F,sep="	",row.names=T,col.names=T)
    write.table(summary_table,file=paste0(outname,"_filterNC.txt"),quote=F,sep="	",row.names=T,col.names=T)
    write.table(out_table,file=paste0(outname,'_NCsummary.txt'),quote=F,sep="	",row.names=T,col.names=T)
} 
  

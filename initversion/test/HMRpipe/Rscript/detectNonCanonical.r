library(glmnet)


# read parameter

a<-commandArgs(T)

outname <- a[1]
signalname <- a[2]
usepq <- a[3]
cutoff <- as.numeric(a[4])
Alpha <- as.numeric(a[5])
LambdaChoice <- (a[6])
topN <- a[7]
R2cutoff <- 0.1

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
print(summary(X))
print(summary(Y))
### elastic net model
set.seed(1007)
## model selection, alpha=0.5
if(ncol(Y) == 1){
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
    coeff.glmmod <- coef(glmmod) 
    coeff_this <- coeff.glmmod
    coTF_name <- names(coeff_this[which(coeff_this[,min(useS,ncol(coeff_this))] != 0),min(useS,ncol(coeff_this))])
    coTF_name_use <-  coTF_name[which(coTF_name!="(Intercept)" & coTF_name!="")]
    all_use_coTF_names <- coTF_name_use
}else{
    ## multi-substrate case
    glmmod = glmnet(x=X, y=Y, family="mgaussian", alpha=Alpha, nlambda=200)
    cv.glmmod <- cv.glmnet(x=X,y=Y,alpha=Alpha,family="mgaussian")
    if (LambdaChoice == "1se"){
        useS <- max(1,which(cv.glmmod$glmnet.fit$lambda == cv.glmmod$lambda.1se)[1])
    }else{
        useS <- max(1,which(cv.glmmod$glmnet.fit$lambda == cv.glmmod$lambda.min)[1])
    }
    coeff.glmmod <- coef(glmmod) 
    all_use_coTF_names <- c()
    for(y in names(coeff.glmmod)){
        coeff_this <- coeff.glmmod[[y,exact=F]]
        coTF_name <- names(coeff_this[which(coeff_this[,min(useS,ncol(coeff_this))] != 0),min(useS,ncol(coeff_this))])
        coTF_name_use <-  coTF_name[which(coTF_name!="(Intercept)" & coTF_name!="")]
        all_use_coTF_names <- c(all_use_coTF_names, coTF_name_use)
    }    
}
## selected co-factors
coTF_binding_matrix <- as.matrix(X[,unique(all_use_coTF_names)])

### re-fit by a uni-variate linear model, output beta, pvalue, r-square for each selected co-factors
unilinear_real <- function(singleX,useY){
    singleLM <- summary(lm(useY~singleX))
    return(c(singleLM$coefficients['singleX',c('Estimate','Pr(>|t|)')],singleLM$adj.r.squared))
}

unilinear_permute <- function(TFov,HMsig){
    lmresult <- summary(lm(HMsig ~ TFov))
    fgR2 <- lmresult$adj.r.squared
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
    return(c(permuteP,fgR2))
    #}else{
    #    return(c("NA",fgR2))
    #}
}

if(ncol(coTF_binding_matrix)==0){
    print("no significant candidates")
}else{
#summary_table <- c()
summary_table_raw <- as.matrix(t(apply(coTF_binding_matrix,2,single_round_linear,Y)))
colnames(summary_table_raw) <- c("Pvalue","R2")

if(usepq == "Q"){
    Qval <- p.adjust(summary_table_raw[,"Pvalue"])
    tmp_sumtable <- cbind(Qval,summary_table_raw[,"R2"])
    colnames(tmp_sumtable) <- c("Qvalue","R2")
}else{
    tmp_sumtable <- summary_table_raw
}

summary_table <- tmp_sumtable[order(tmp_sumtable[,"R2"],decreasing=TRUE),]

if(topN == "all"){
    topN <- nrow(summary_table)
}else{
    topN <- as.numeric(topN)
}

out_table <- c()
selected_coTF <- c()
N <- 0
for(coTF in rownames(summary_table)){
    if(summary_table[coTF,1] < cutoff){
        N <- N+1
        if(N <= topN){
            out_table <- rbind(out_table,summary_table[coTF,])
            selected_coTF <- c(selected_coTF,coTF)
        }
    }
}
out_table <- as.matrix(out_table)
rownames(out_table) <- selected_coTF
colnames(out_table) <- colnames(summary_table)
write.table(out_table,file=paste0(outname,'_NCsummary.txt'),quote=F,sep="\t",row.names=T,col.names=T)


#if(ncol(coTF_binding_matrix)==0){
#
#}
#
#### output
#allout <-c()
#    for(yname in colnames(Y)){
#        allout <- cbind(allout,matrix(rep("NA",3),nrow=1))
#    }
#}else{
#    for(yname in colnames(Y)){
#    single_linear_coeff <- as.matrix(t(apply(coTF_binding_matrix,2,single_round_linear,Y[,yname])))
#    allout <- cbind(allout,single_linear_coeff)
#    }
#}
#colnames(allout) <- paste(rep(colnames(Y),each=3),c('Estimate','Pval','adjRsq'),sep="_")
#write.table(allout,file=paste0('elastic_mgaussian/',crid,'_outinfo.txt'),quote=F,sep="\t",row.names=T,col.names=T)
#
#
#
#











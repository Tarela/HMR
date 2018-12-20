library(glmnet)

### parameters
a<-commandArgs(T)
crid <- a[1] ### e.g. mESC_GSM1562337_CBX7
### read data
crid <- "mESC_GSM656523_RNF2"

HMsig_raw <- read.table(paste0('/Volumes/Tarela/S1bk/Project/CR/Data/overlap_signal_matrix/',crid,'_HMsig.bed'),row.names=4,header=T)
TFov_raw <- read.table(paste0('/Volumes/Tarela/S1bk/Project/CR/Data/overlap_signal_matrix/',crid,'_TFov.bed'),row.names=4,header=T)

### fetch details
gsmid <- strsplit(crid,"_")[[1]][2]
crname <- strsplit(crid,"_")[[1]][3]
hmname <- colnames(HMsig_raw)[10:ncol(HMsig_raw)]

### step1, data preprocess
use_rownames <- intersect(rownames(HMsig_raw),rownames(TFov_raw))
## prepare HM signal
HMsig <- as.matrix(HMsig_raw[use_rownames,hmname])
colnames(HMsig) <- hmname
rownames(HMsig) <- use_rownames

## prepare TF overlap
get_factor_name <- function(full_name){
return(strsplit(full_name,"_")[[1]][3])
}
TFnames_raw <- colnames(TFov_raw)[10:ncol(TFov_raw)]
# samples for same factors are excluded
TFov <- TFov_raw[use_rownames,TFnames_raw[which(unlist(lapply(TFnames_raw,get_factor_name))!=crname)]]
TFov[TFov > 1] <- 1

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

### elastic net model
set.seed(1007)
## model selection, alpha=0.5
if(ncol(Y) == 1){
    ## uni-substrate case
    glmmod = glmnet(x=X, y=Y, family="gaussian", alpha=0.5, nlambda=200)
    cv.glmmod <- cv.glmnet(x=X,y=Y,alpha=0.5,family="gaussian")
    useS <- max(1,which(cv.glmmod$glmnet.fit$lambda == cv.glmmod$lambda.1se)[1])
    coeff.glmmod <- coef(glmmod) 
    coeff_this <- coeff.glmmod
    coTF_name <- names(coeff_this[which(coeff_this[,min(useS,ncol(coeff_this))] != 0),min(useS,ncol(coeff_this))])
    coTF_name_use <-  coTF_name[which(coTF_name!="(Intercept)" & coTF_name!="")]
    all_use_coTF_names <- coTF_name_use
}else{
    ## multi-substrate case
    glmmod = glmnet(x=X, y=Y, family="mgaussian", alpha=0.5, nlambda=200)
    cv.glmmod <- cv.glmnet(x=X,y=Y,alpha=0.5,family="mgaussian")
    useS <- max(1,which(cv.glmmod$glmnet.fit$lambda == cv.glmmod$lambda.1se)[1])
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
single_round_linear <- function(singleX,useY){
    singleLM <- summary(lm(useY~singleX))
    return(c(singleLM$coefficients['singleX',c('Estimate','Pr(>|t|)')],singleLM$adj.r.squared))
}

### output
allout <-c()
if(ncol(coTF_binding_matrix)==0){
    for(yname in colnames(Y)){
        allout <- cbind(allout,matrix(rep("NA",3),nrow=1))
    }
}else{
    for(yname in colnames(Y)){
    single_linear_coeff <- as.matrix(t(apply(coTF_binding_matrix,2,single_round_linear,Y[,yname])))
    allout <- cbind(allout,single_linear_coeff)
    }
}
colnames(allout) <- paste(rep(colnames(Y),each=3),c('Estimate','Pval','adjRsq'),sep="_")
write.table(allout,file=paste0('elastic_mgaussian/',crid,'_outinfo.txt'),quote=F,sep="\t",row.names=T,col.names=T)









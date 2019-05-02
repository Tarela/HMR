#export PATH=/Users/sh8tv/Dropbox/CR/HMRpackage/initversion/refpackage/bwsummary:$PATH

./HMR -p /Users/sh8tv/Dropbox/CR/HMR_testdata/mESC_GSM1562337_CBX7.bed -s /Users/sh8tv/Dropbox/CR/HMR_testdata/mESC_GSM1399500_H3K27me3.bw,/Users/sh8tv/Dropbox/CR/HMR_testdata/mESC_GSM1954953_H3K9me3.bw -f /Users/sh8tv/Dropbox/CR/HMR_testdata/TFCRpeak_uniq/ -o test7 --overwrite
#Rscript tmpPackage/detectNonClassic.r test5 mESC_GSM1399500_H3K27me3,mESC_GSM1954953_H3K9me3 0.001 0.5 1se all /Users/sh8tv/Dropbox/CR/HMRpackage/initversion/test/tmpPackage/


ncHMR_detector -p /scratch/sh8tv/Project/tmp/HMRpackage/usedata/mESC_GSM1562337_CBX7.bed -s /scratch/sh8tv/Project/tmp/HMRpackage/usedata/mESC_GSM1399500_H3K27me3.bw,/scratch/sh8tv/Project/tmp/HMRpackage/usedata/mESC_GSM1954953_H3K9me3.bw -f /scratch/sh8tv/Project/tmp/HMRpackage/usedata/TFCRpeak_uniq/ -o test8 --overwrite



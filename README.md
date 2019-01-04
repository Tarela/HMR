# HMR package

HMR: Detecting non-canonical function of histone modification regulator

Usage: HMR -p HMRpeak.bed -s HMsig.bw -f TFovFolder -o outputname

Options:
  --version             show program's version number and exit
  -h, --help            Show this help message and exit.
  -p HMRPEAK, --HMRpeak=HMRPEAK
                        Peak file for genome-wide HMR binding sites (bed
                        format, at least 3 column)
  -s SIGNAL, --Signal=SIGNAL
                        Signal track for genome-wide HM substrate signal
                        (bigwig format), comma separated for multiple HM
                        substrate (each substrate accept ONLY 1 bw file)
  -f PEAKFOLDER, --peakFolder=PEAKFOLDER
                        Folder for peak files of potential co-factors,
                        ABSOLUTE path
  -o OUTNAME, --outname=OUTNAME
                        Name (prefix) of output results and output directory,
                        default is NA
  --Qvalue=QVALUE       [optional] Cutoff of Q-value, default is 1e-2 (0.01),
                        exclusive from P-value
  --Pvalue=PVALUE       [optional] Cutoff of P-value, default is 1e-3 (0.001),
                        exclusive from Q-value
  --Alpha=ALPHA         [optional] alpha parameters for elasticNet, choose
                        from 0~1, 1 for lasso and 0 for ridge, default is 0.5
  --LambdaChoice=LAMBDACHOICE
                        [optional] Solution to determine Lambda (choose from
                        1se and min, default is 1se)
  --TopNcofactors=TOPNCOFACTORS
                        [optional] TopN predicted co-factors with highest
                        association with non-classic function is reported
                        (choose any number or all(default) to report topN
                        predicted co-factors that pass the thresholds)
  --overwrite           [optional] force overwrite, this cmd will rm existing
                        result if set !!


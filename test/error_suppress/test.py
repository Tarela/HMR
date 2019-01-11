import subprocess
def sp(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac
def sperr(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stderr=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac

bwfile = "/Users/sh8tv/Dropbox/CR/HMR_testdata/mESC_GSM1399500_H3K27me3.bw"

print sp("/Users/sh8tv/Dropbox/CR/HMRpackage/initversion/refpackage/bwsummary/bigWigSummary_mac %s chrJH584304.1 48487 48767 1"%(bwfile))
print sp("/Users/sh8tv/Dropbox/CR/HMRpackage/initversion/refpackage/bwsummary/bigWigSummary_mac %s chr9 35305126 35305591 1"%(bwfile))


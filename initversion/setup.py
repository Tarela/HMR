#!/usr/bin/env python
"""Description
Setup script for HMR  -- Detecting non-canonical function of histone modification regulator
Copyright (c) 2019 Shengen Hu <tarelahu@gmail.com>
This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).
"""
import os
import sys
import subprocess
from distutils.core import setup, Extension


def sp(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stdout=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac

def compile_bedtools():
    curdir = os.getcwd()
    os.chdir('refpackage/bedtools')
    sp('make 1>/dev/null 2>&1 ')
    sp('chmod 755 *')
    os.chdir(curdir)
    
def check_bedtools():
    checkhandle = sp('which bedtools')
    if checkhandle[0].strip() == "":
        return 0
    else:
        return 1
def check_R():
    checkhandle = sp('which Rscript')
    if checkhandle[0].strip() == "":
        return 0
    else:
        return 1       
def main(): 
    if sys.version_info[0] != 2 or sys.version_info[1] < 7:
	    print >> sys.stderr, "ERROR: HMR requires Python 2.7"
	    sys.exit()
    has_R = check_R()
    if has_R == 0:
	    print >> sys.stderr, "ERROR: HMR requires R & Rscript under default PATH"
	    sys.exit()
        
    has_bedtools = check_bedtools()
    print 'Intalling HMR, may take serval minutes'
    if has_bedtools == 0:
        compile_bedtools()
        setup(name="HMRpipe",
              version="1.1.3",
              description="HMR: Detecting non-canonical function of histone modification regulator",
              author='Shengen Hu',
              author_email='Tarelahu@gmail.com',
              url='https://github.com/Tarela/HMR.git',
              package_dir={'HMRpipe' : 'lib'},
              packages=['HMRpipe'],
              package_data={'HMRpipe': [#'Config/Drseq_template.conf',
                                      #'Rscript/analysis.r',
                                      #'Rscript/individual_qc.r',
                                      #'Rscript/readsbulkQC.r',
                                      'Rscript/detectNonCanonical.r'
                                         ]},
              scripts=['bin/HMR','refpackage/bedtools/bin/bedtools','refpackage/bwsummary/bigWigSummary_linux','refpackage/bwsummary/bigWigSummary_mac'],
                        
              classifiers=[
            'Development Status :: version1.0 finish',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: pipeline',
            ],
              requires=['glmnet_py'],
          )
        print 'bedtools is not detected under default PATH, bedtools is also installed'
        print 'Installation of HMR is DONE'
    
    else:
        setup(name="HMRpipe",
              version="1.0",
              description="HMR: Detecting non-canonical function of histone modification regulator",
              author='Shengen Hu',
              author_email='Tarelahu@gmail.com',
              url='https://github.com/Tarela/HMR.git',
              package_dir={'HMRpipe' : 'lib'},
              packages=['HMRpipe'],
              package_data={'HMRpipe': [#'Config/Drseq_template.conf',
                                      #'Rscript/analysis.r',
                                      #'Rscript/individual_qc.r',
                                      #'Rscript/readsbulkQC.r',
                                      'Rscript/detectNonCanonical.r'
                                         ]},
              scripts=['bin/HMR','refpackage/bwsummary/bigWigSummary_linux','refpackage/bwsummary/bigWigSummary_mac'],
                        
              classifiers=[
            'Development Status :: version1.0 finish',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: pipeline',
            ],
              requires=['glmnet_py'],
          )
        print 'Installation of HMR is DONE'


if __name__ == '__main__':
    main()


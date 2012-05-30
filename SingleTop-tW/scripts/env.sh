source $VO_CMS_SW_DIR/cmsset_default.sh
export CVS_RSH=/usr/bin/ssh
PS1="\u@\h:\w ] "

alias python='~/Python-2.6.5/bin/python2.6'


#export ROOTSYS="/user/cmssoft/root_5.26.00e_iihe_default_dcap/root"
#export ROOTSYS="/user/cmssoft/root_5.28.00b_iihe_default_dcap/root"
#export ROOTSYS="/localgrid/cmstools/ROOT/5.26.00b_bak/slc4_ia32_gcc34/root"
#source /localgrid/cmstools/ROOT/5.26.00b_bak/slc4_ia32_gcc34/root/bin/thisroot.sh
export ROOTSYS="/user/cmssoft/root_5.30.02/root"
#export ROOTSYS=/user/cmssoft/root_old
#export ROOTSYS=/user/cmssoft/root_5.27.06_iihe



export PATH=$PATH:$ROOTSYS/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/lib

export CVS_RSH=/usr/bin/ssh
export CVSROOT=:ext:rebeca@cmscvs.cern.ch:/cvs_server/repositories/CMSSW

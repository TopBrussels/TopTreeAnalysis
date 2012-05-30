/// Made faster Freya-style
#include <iostream>
#include "TSystem.h"


void shapechain(int nsel, int mode, int syst = 0, bool up = true, bool silent = true){  

  gSystem->CompileMacro("shapelooper.C","k");

 
  // samples used
  double x_sec = 0.;
  char plotName[300];
  sprintf(plotName,"test");
  
  if (nsel == 0)                	{sprintf(plotName,"tt");}
  else if (nsel == 1)   		{sprintf(plotName,"twdr");}
  else if (nsel == -1)   		{sprintf(plotName,"twds");}
  else if (nsel == 2)   		{sprintf(plotName,"zjetsall");}
  else if (nsel == 3)   		{sprintf(plotName,"di");}
  else if (nsel == 4)			{sprintf(plotName, "st");}
  else if (nsel == 5)   		{sprintf(plotName,"wjets");}
  else if (nsel == 6)   		{sprintf(plotName,"qcd_mu");}
  else if (nsel == 7)                	{sprintf(plotName,"others");}
  else if (nsel == 8)                	{sprintf(plotName,"tt_scaleup");}
  else if (nsel == 9)                	{sprintf(plotName,"tt_scaledown");}
  else if (nsel == 10)                	{sprintf(plotName,"tw_sdo");}
  else if (nsel == 11)                	{sprintf(plotName,"tw_sup");}
  else if (nsel == 12)                	{sprintf(plotName,"tt_matchingup");}
  else if (nsel == 13)                	{sprintf(plotName,"tt_matchingdown");}

  else if (nsel == 555)                	{sprintf(plotName,"mc");}
  
  else if (nsel == 666)                	{sprintf(plotName,"data");}
  else if (nsel == 6661)                {sprintf(plotName,"data1");}
  else if (nsel == 6662)                {sprintf(plotName,"data2");}
  
  if (mode != 0 &&  mode !=1 && mode !=2) mode = 0;
  if (!silent){
    cout << "[Info:]" ;
    if      (mode == 0) 	cout << " emu channel, " ;
    else if (mode == 1) 	cout << " mumu channel, " ;
    else if (mode == 2) 	cout << " ee channel, " ;
  }
  

  bool RunA = false;
  bool RunB = false;
  char myRootFile[300];
  sprintf(myRootFile,"outputs/out_%d_%s.root", mode, plotName);
  if (RunA) sprintf(myRootFile,"outputs/out_runA_%d_%s.root", mode, plotName);
  if (RunB) sprintf(myRootFile,"outputs/out_runB_%d_%s.root", mode, plotName);
  if (syst == 1 && up) sprintf(myRootFile,"outputs/PUsysUp_%d_%s.root", mode, plotName);
  if (syst == 1 && !up) sprintf(myRootFile,"outputs/PUsysDown_%d_%s.root", mode, plotName);
  if (syst == 2 && up) sprintf(myRootFile,"outputs/JESsysUp_%d_%s.root", mode, plotName);
  if (syst == 2 && !up) sprintf(myRootFile,"outputs/JESsysDown_%d_%s.root", mode, plotName);
  if (syst == 3 && up) sprintf(myRootFile,"outputs/JERsysUp_%d_%s.root", mode, plotName);
  if (syst == 3 && !up) sprintf(myRootFile,"outputs/JERsysDown_%d_%s.root", mode, plotName);
  if (syst == 4 && up) sprintf(myRootFile,"outputs/METsysUp_%d_%s.root", mode, plotName);
  if (syst == 4 && !up) sprintf(myRootFile,"outputs/METsysDown_%d_%s.root", mode, plotName);
  
  
  TChain *myCh = new TChain("myTree","myTree");
  myCh->Add(myRootFile);

  shapelooper an(myCh);
  std::cout << plotName << " (" << myRootFile << ",  " << myCh->GetEntries() << " entries )" << std::endl;
  an.myLoop(nsel,mode,syst,up,silent);
  
}

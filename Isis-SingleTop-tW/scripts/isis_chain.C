/// isis 
#include <iostream>
#include "TSystem.h"



void isis_chain(int nsel = 0, int mode = 0, bool silent = false){  

  gSystem->CompileMacro("isis_looper.C","k");

 
  // samples used
  double x_sec = 0.;
  char plotName[300];
  
  
  if (nsel == 0)                	{sprintf(plotName,"tt");}
  else if (nsel == 1)   		{sprintf(plotName,"twdr");}
  else if (nsel == -1)   		{sprintf(plotName,"twds");}
  else if (nsel == 2)   		{sprintf(plotName,"zjetsall");}
  else if (nsel == 3)   		{sprintf(plotName,"di");}
  else if (nsel == 4)			{sprintf(plotName, "st");}
  else if (nsel == 5)   		{sprintf(plotName,"wjets");}
 
  else if (nsel == 7)                	{sprintf(plotName,"others");}

  else if (nsel == 555)                	{sprintf(plotName,"mc");}
  
  else if (nsel == 666)                	{sprintf(plotName,"data");}
  
  
  
  // JER  
  else if (nsel == -10)                   {sprintf(plotName,"tt");}  
  else if (nsel ==  10)                   {sprintf(plotName,"tt");}
  else if (nsel == -20)                   {sprintf(plotName,"twdr");}
  else if (nsel ==  20)                   {sprintf(plotName,"twdr");}
  else if (nsel == -30)                   {sprintf(plotName,"others");}
  else if (nsel ==  30)                   {sprintf(plotName,"others");}
  
  //JES
  else if (nsel == -11)                   {sprintf(plotName,"tt");}  
  else if (nsel ==  11)                   {sprintf(plotName,"tt");}
  else if (nsel == -21)                   {sprintf(plotName,"twdr");}
  else if (nsel ==  21)                   {sprintf(plotName,"twdr");}
  else if (nsel == -31)                   {sprintf(plotName,"others");}
  else if (nsel ==  31)                   {sprintf(plotName,"others");}
  //PU
  else if (nsel == -12)                   {sprintf(plotName,"tt");}  
  else if (nsel ==  12)                   {sprintf(plotName,"tt");}
  else if (nsel == -22)                   {sprintf(plotName,"twdr");}
  else if (nsel ==  22)                   {sprintf(plotName,"twdr");}
  else if (nsel == -32)                   {sprintf(plotName,"others");}
  else if (nsel ==  32)                   {sprintf(plotName,"others");}
  
  //SF
  else if (nsel == -13)                   {sprintf(plotName,"tt");}  
  else if (nsel ==  13)                   {sprintf(plotName,"tt");}
 else if (nsel == -23)                   {sprintf(plotName,"twdr");}
  else if (nsel ==  23)                   {sprintf(plotName,"twdr");}
  else if (nsel == -33)                   {sprintf(plotName,"others");}
  else if (nsel ==  33)                   {sprintf(plotName,"others");}
  
    //MET
  else if (nsel == -14)                   {sprintf(plotName,"tt");}  
  else if (nsel ==  14)                   {sprintf(plotName,"tt");}
  else if (nsel == -24)                   {sprintf(plotName,"twdr");}
  else if (nsel ==  24)                   {sprintf(plotName,"twdr");}
  else if (nsel == -34)                   {sprintf(plotName,"others");}
  else if (nsel ==  34)                   {sprintf(plotName,"others");}
  
      //topmass
  else if (nsel == -15)                   {sprintf(plotName,"tt");}  
  else if (nsel ==  15)                   {sprintf(plotName,"tt");}
  else if (nsel == -25)                   {sprintf(plotName,"twdr");}
  else if (nsel ==  25)                   {sprintf(plotName,"twdr");}
 
  
  
      //Q2
  else if (nsel == -16)                   {sprintf(plotName,"tt");}  
  else if (nsel ==  16)                   {sprintf(plotName,"tt");}
  else if (nsel == -26)                   {sprintf(plotName,"twdr");}
  else if (nsel ==  26)                   {sprintf(plotName,"twdr");}

  
  
      //LES
  else if (nsel == -17)                   {sprintf(plotName,"tt");}  
  else if (nsel ==  17)                   {sprintf(plotName,"tt");}
  else if (nsel == -27)                   {sprintf(plotName,"twdr");}
  else if (nsel ==  27)                   {sprintf(plotName,"twdr");}
  else if (nsel == -37)                   {sprintf(plotName,"others");}
  else if (nsel ==  37)                   {sprintf(plotName,"others");}
  
      //Matching
  else if (nsel == -18)                   {sprintf(plotName,"tt");}  
  else if (nsel ==  18)                   {sprintf(plotName,"tt");}
 
  // ZSF 
  else if (nsel == -19)                   {sprintf(plotName,"tt");}  
  else if (nsel ==  19)                   {sprintf(plotName,"tt");}
  else if (nsel == -29)                   {sprintf(plotName,"twdr");}
  else if (nsel ==  29)                   {sprintf(plotName,"twdr");}
  else if (nsel == -39)                   {sprintf(plotName,"others");}
  else if (nsel ==  39)                   {sprintf(plotName,"others");}
  
  if (mode != 0 &&  mode !=1 && mode !=2) mode = 0;
  if (!silent){
    cout << "[Info:]" ;
    if      (mode == 0) 	cout << " emu channel, " ;
    else if (mode == 1) 	cout << " mumu channel, " ;
    else if (mode == 2) 	cout << " ee channel, " ;
  }
  
  char myRootFile[300];
   if(nsel == -19 || nsel == -29 || nsel == -39){
  sprintf(myRootFile,"outputs/out_ZSFsysDown_%d_%s.root", mode, plotName);
  }else if(nsel == 19 || nsel == 29 || nsel == 39){
   sprintf(myRootFile,"outputs/out_ZSFsysUp_%d_%s.root", mode, plotName);
  }else if(nsel == -10 || nsel == -20 || nsel == -30){
  sprintf(myRootFile,"outputs/JERsysDown_%d_%s.root", mode, plotName);
  }else if(nsel == 10 || nsel == 20 || nsel == 30){
   sprintf(myRootFile,"outputs/JERsysUp_%d_%s.root", mode, plotName);
  }else if(nsel == -11 || nsel == -21 || nsel == -31){
   sprintf(myRootFile,"outputs/JESsysDown_%d_%s.root", mode, plotName);
  }else if(nsel == 11 || nsel == 21 || nsel == 31){
   sprintf(myRootFile,"outputs/JESsysUp_%d_%s.root", mode, plotName);
  }else if(nsel == -12 || nsel == -22 || nsel == -32){
   sprintf(myRootFile,"outputs/PUsysDown_%d_%s.root", mode, plotName);
  }else if(nsel == 12 || nsel == 22 || nsel == 32){
   sprintf(myRootFile,"outputs/PUsysUp_%d_%s.root", mode, plotName);
  }else if(nsel == -13 || nsel == -23 || nsel == -33 ){
   sprintf(myRootFile,"outputs/SFsysDown_%d_%s.root", mode, plotName);
  }else if(nsel == 13 || nsel == 23 || nsel == 33){
   sprintf(myRootFile,"outputs/SFsysUp_%d_%s.root", mode, plotName);
  }else if(nsel == -14 || nsel == -24 || nsel == -34){
   sprintf(myRootFile,"outputs/METsysDown_%d_%s.root", mode, plotName);
  }else if(nsel == 14 || nsel == 24 || nsel == 34){
   sprintf(myRootFile,"outputs/METsysUp_%d_%s.root", mode, plotName);
  }else if(nsel == -15 || nsel == -25 ){
   sprintf(myRootFile,"outputs/TopMassDown_%d_%s.root", mode, plotName);
  }else if(nsel == 15 || nsel == 25 ){
   sprintf(myRootFile,"outputs/TopMassUp_%d_%s.root", mode, plotName);
  }else if(nsel == -16 || nsel == -26 ){
   sprintf(myRootFile,"outputs/Q2Down_%d_%s.root", mode, plotName);
  }else if(nsel == 16 || nsel == 26) {
   sprintf(myRootFile,"outputs/Q2Up_%d_%s.root", mode, plotName);
  }else if(nsel == -17 || nsel == -27 || nsel == -37){
   sprintf(myRootFile,"outputs/eleSFsysDown_%d_%s.root", mode, plotName);
  }else if(nsel == 17 || nsel == 27 || nsel == 37){
   sprintf(myRootFile,"outputs/eleSFsysUp_%d_%s.root", mode, plotName);
  }else if(nsel == -18 ){
   sprintf(myRootFile,"outputs/matchingDown_%d_%s.root", mode, plotName);
  }else if(nsel == 18 ){
   sprintf(myRootFile,"outputs/matchingUp_%d_%s.root", mode, plotName);
  }  else{
   sprintf(myRootFile,"outputs/out_%d_%s.root", mode, plotName);
  }
  
  
  
  TChain *myCh = new TChain("myTree","myTree");
  myCh->Add(myRootFile);

  isis_looper an(myCh);
  std::cout << plotName << " (" << myRootFile << ",  " << myCh->GetEntries() << " entries )" << std::endl;
  an.myLoop(nsel,mode,silent);
  
}

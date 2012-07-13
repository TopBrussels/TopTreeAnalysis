#include "../interface/BTagName.h"

BTagName::BTagName(){    
  
}

BTagName::~BTagName(){

}

std::string BTagName::NameGiving(int TCHEbTagLoop, int NumberTCHEbTags,  int TCHPbTagLoop, int NumberTCHPbTags, int SSVHEbTagLoop, int NumberSSVHEbTags, int SSVHPbTagLoop, int NumberSSVHPbTags,int CSVbTagLoop, int NumberCSVbTags,int JPbTagLoop, int NumberJPbTags,int JBPbTagLoop, int NumberJBPbTags){

  int bTagSum = TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1+JPbTagLoop-1+JBPbTagLoop-1;
  
  if(TCHEbTagLoop == 1)
    bTagFileOutput[bTagSum] = "No bTag constraint applied";
  else if(TCHEbTagLoop == 2)
    bTagFileOutput[bTagSum] =" Loose TCHE working point bTag (1.7) for hadronic b-jet, no constraint on leptonic b-jet " ;
  else if(TCHEbTagLoop == 3)
    bTagFileOutput[bTagSum] =" Loose TCHE working point bTag (1.7) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if (TCHEbTagLoop == 4)
    bTagFileOutput[bTagSum] =" Loose TCHE working point bTag (1.7) for leptonic b-jet or hadronic b-jet ";
  else if (TCHEbTagLoop == 5)
    bTagFileOutput[bTagSum] =" Loose TCHE working point bTag (1.7) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == 6)
    bTagFileOutput[bTagSum] = " Medium TCHE working point bTag (3.3) for hadronic b-jet, no constraint on leptonic b-jet ";
  else if (TCHEbTagLoop == 7)
    bTagFileOutput[bTagSum] = " Medium TCHE working point bTag (3.3) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if (TCHEbTagLoop == 8)
    bTagFileOutput[bTagSum] = " Medium TCHE working point bTag (3.3) for leptonic b-jet or hadronic b-jet ";  
  else if (TCHEbTagLoop == 9)
    bTagFileOutput[bTagSum] = " Medium TCHE working point bTag (3.3) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == 10)
    bTagFileOutput[bTagSum] = " Tight TCHE working point bTag (10.2) for hadronic b-jet, no constraint on leptonic b-jet ( Tight TCHE is not supported )";	
  else if (TCHEbTagLoop == 11)
    bTagFileOutput[bTagSum] = " Tight TCHE working point bTag (10.2) for leptonic b-jet, no constraint on hadronic b-jet ( Tight TCHE is not supported )";		  
  else if (TCHEbTagLoop == 12)
    bTagFileOutput[bTagSum] = " Tight TCHE working point bTag (10.2) for leptonic b-jet or hadronic b-jet ( Tight TCHE is not supported )";  
  else if (TCHEbTagLoop == 13)
    bTagFileOutput[bTagSum] = " Tight TCHE working point bTag (10.2) for leptonic b-jet and hadronic b-jet ( Tight TCHE is not supported )";		
  
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 2)
    bTagFileOutput[bTagSum] = " Loose TCHP working point bTag (1.19) for hadronic b-jet, no constraint on leptonic b-jet (Loose TCHP is not supported)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 3)
    bTagFileOutput[bTagSum] = " Loose TCHP working point bTag (1.19) for leptonic b-jet, no constraint on hadronic b-jet (Loose TCHP is not supported)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 4)
    bTagFileOutput[bTagSum] = " Loose TCHP working point bTag (1.19) for leptonic b-jet or hadronic b-jet (Loose TCHP is not supported)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 5)
    bTagFileOutput[bTagSum] = " Loose TCHP working point bTag (1.19) for leptonic b-jet and hadronic b-jet (Loose TCHP is not supported)";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 6)
    bTagFileOutput[bTagSum] = " Medium TCHP working point bTag (1.93) for hadronic b-jet, no constraint on leptonic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 7)
    bTagFileOutput[bTagSum] = " Medium TCHP working point bTag (1.93) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 8)
    bTagFileOutput[bTagSum] = " Medium TCHP working point bTag (1.93) for leptonic b-jet or hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 9)
    bTagFileOutput[bTagSum] = " Medium TCHP working point bTag (1.93) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 10)
    bTagFileOutput[bTagSum] = " Tight TCHP working point bTag (3.41) for hadronic b-jet, no constraint on leptonic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 11)
    bTagFileOutput[bTagSum] = " Tight TCHP working point bTag (3.41) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 12)
    bTagFileOutput[bTagSum] = " Tight TCHP working point bTag (3.41) for leptonic b-jet or hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 13)
    bTagFileOutput[bTagSum] = " Tight TCHP working point bTag (3.41) for leptonic b-jet and hadronic b-jet ";
  
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop ==2)
    bTagFileOutput[bTagSum] = " Loose SSVHE working point bTag (0) for hadronic b-jet, no constraint on leptonic b-jet (SSVHE == 0 no POG wp)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 3)
    bTagFileOutput[bTagSum] = " Loose SSVHE working point bTag (0) for leptonic b-jet, no constraint on hadronic b-jet(SSVHE == 0 no POG wp) ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 4)
    bTagFileOutput[bTagSum] = " Loose SSVHE working point bTag (0) for leptonic b-jet or hadronic b-jet (SSVHE == 0 no POG wp)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 5)
    bTagFileOutput[bTagSum] = " Loose SSVHE working point bTag (0) for leptonic b-jet and hadronic b-jet (SSVHE == 0 no POG wp)";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 6)
    bTagFileOutput[bTagSum] = " Medium SSVHE working point bTag (1.74) for hadronic b-jet, no constraint on leptonic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 7)
    bTagFileOutput[bTagSum] = " Medium SSVHE working point bTag (1.74) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 8)
    bTagFileOutput[bTagSum] = " Medium SSVHE working point bTag (1.74) for leptonic b-jet or hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 9)
    bTagFileOutput[bTagSum] = " Medium SSVHE working point bTag (1.74) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 10)
    bTagFileOutput[bTagSum] = " Tight SSVHE working point bTag (3.05) for hadronic b-jet, no constraint on leptonic b-jet ( Tight SSVHE is not supported )";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 11)
    bTagFileOutput[bTagSum] = " Tight SSVHE working point bTag (3.05) for leptonic b-jet, no constraint on hadronic b-jet ( Tight SSVHE is not supported )";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 12)
    bTagFileOutput[bTagSum] = " Tight SSVHE working point bTag (3.05) for leptonic b-jet or hadronic b-jet ( Tight SSVHE is not supported )";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 13)
    bTagFileOutput[bTagSum] = " Tight SSVHE working point bTag (3.05) for leptonic b-jet and hadronic b-jet ( Tight SSVHE is not supported )";
  
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==2)
    bTagFileOutput[bTagSum] = " Loose SSVHP working point bTag (0) for hadronic b-jet, no constraint on leptonic b-jet (SSVHE == 0 no POG wp)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 3)
    bTagFileOutput[bTagSum] = " Loose SSVHP working point bTag (0) for leptonic b-jet, no constraint on hadronic b-jet(SSVHE == 0 no POG wp) ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 4)
    bTagFileOutput[bTagSum] = " Loose SSVHP working point bTag (0) for leptonic b-jet or hadronic b-jet (SSVHE == 0 no POG wp)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 5)
    bTagFileOutput[bTagSum] = " Loose SSVHP working point bTag (0) for leptonic b-jet and hadronic b-jet (SSVHE == 0 no POG wp)";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 6)
    bTagFileOutput[bTagSum] = " Medium SSVHP working point bTag 1) for hadronic b-jet, no constraint on leptonic b-jet (SSVHE == 1 no POG wp)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 7)
    bTagFileOutput[bTagSum] = " Medium SSVHP working point bTag (1) for leptonic b-jet, no constraint on hadronic b-jet (SSVHE == 1 no POG wp)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 8)
    bTagFileOutput[bTagSum] = " Medium SSVHP working point bTag (1) for leptonic b-jet or hadronic b-jet(SSVHE == 1 no POG wp) ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 9)
    bTagFileOutput[bTagSum] = " Medium SSVHP working point bTag (1) for leptonic b-jet and hadronic b-jet (SSVHE == 1 no POG wp)";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 10)
    bTagFileOutput[bTagSum] = " Tight SSVHP working point bTag (2) for hadronic b-jet, no constraint on leptonic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 11)
    bTagFileOutput[bTagSum] = " Tight SSVHP working point bTag (2) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 12)
    bTagFileOutput[bTagSum] = " Tight SSVHP working point bTag (2) for leptonic b-jet or hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 13)
    bTagFileOutput[bTagSum] = " Tight SSVHP working point bTag (2) for leptonic b-jet and hadronic b-jet ";

    //Combined secondary vertex loop:
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 2)
    bTagFileOutput[bTagSum] = " Loose CSV working point bTag (0.244) for hadronic b-jet, no constraint on leptonic b-jet ";
    else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 3)
    bTagFileOutput[bTagSum] = " Loose CSV working point bTag (0.244) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 4)
    bTagFileOutput[bTagSum] = " Loose CSV working point bTag (0.244) for leptonic b-jet or hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 5)
    bTagFileOutput[bTagSum] = " Loose CSV working point bTag (0.244) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 6)
    bTagFileOutput[bTagSum] = " Medium CSV working point bTag (0.679) for hadronic b-jet, no constraint on leptonic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 7)
    bTagFileOutput[bTagSum] = " Medium CSV working point bTag (0.679) for leptonic b-jet, no constraint on hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 8)
    bTagFileOutput[bTagSum] = " Medium CSV working point bTag (0.679) for leptonic b-jet or hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 9)
    bTagFileOutput[bTagSum] = " Medium CSV working point bTag (0.679) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 10)
    bTagFileOutput[bTagSum] = " Tight CSV working point bTag (0.898) for hadronic b-jet, no constraint on leptonic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 11)
    bTagFileOutput[bTagSum] = " Tight CSV working point bTag (0.898) for leptonic b-jet, no constraint on hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 12)
    bTagFileOutput[bTagSum] = " Tight CSV working point bTag (0.898) for leptonic b-jet or hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 13)
    bTagFileOutput[bTagSum] = " Tight CSV working point bTag (0.898) for leptonic b-jet and hadronic b-jet  ";

  //Jet probability loop:
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 2)
    bTagFileOutput[bTagSum] = " Loose JP working point bTag (0.275) for hadronic b-jet, no constraint on leptonic b-jet ";
    else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 3)
    bTagFileOutput[bTagSum] = " Loose JP working point bTag (0.275) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 4)
    bTagFileOutput[bTagSum] = " Loose JP working point bTag (0.275) for leptonic b-jet or hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 5)
    bTagFileOutput[bTagSum] = " Loose JP working point bTag (0.275) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 6)
    bTagFileOutput[bTagSum] = " Medium JP working point bTag (0.545) for hadronic b-jet, no constraint on leptonic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 7)
    bTagFileOutput[bTagSum] = " Medium JP working point bTag (0.545) for leptonic b-jet, no constraint on hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 8)
    bTagFileOutput[bTagSum] = " Medium JP working point bTag (0.545) for leptonic b-jet or hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 9)
    bTagFileOutput[bTagSum] = " Medium JP working point bTag (0.545) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 10)
    bTagFileOutput[bTagSum] = " Tight JP working point bTag (0.790) for hadronic b-jet, no constraint on leptonic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 11)
    bTagFileOutput[bTagSum] = " Tight JP working point bTag (0.790) for leptonic b-jet, no constraint on hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 12)
    bTagFileOutput[bTagSum] = " Tight JP working point bTag (0.790) for leptonic b-jet or hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 13)
    bTagFileOutput[bTagSum] = " Tight JP working point bTag (0.790) for leptonic b-jet and hadronic b-jet  ";

  //Jet B probability loop:
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 2)
    bTagFileOutput[bTagSum] = " Loose JBP working point bTag (1.33) for hadronic b-jet, no constraint on leptonic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 3)
    bTagFileOutput[bTagSum] = " Loose JBP working point bTag (1.33) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 4)
    bTagFileOutput[bTagSum] = " Loose JBP working point bTag (1.33) for leptonic b-jet or hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 5)
    bTagFileOutput[bTagSum] = " Loose JBP working point bTag (1.33) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 6)
    bTagFileOutput[bTagSum] = " Medium JBP working point bTag (2.55) for hadronic b-jet, no constraint on leptonic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 7)
    bTagFileOutput[bTagSum] = " Medium JBP working point bTag (2.55) for leptonic b-jet, no constraint on hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 8)
    bTagFileOutput[bTagSum] = " Medium JBP working point bTag (2.55) for leptonic b-jet or hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 9)
    bTagFileOutput[bTagSum] = " Medium JBP working point bTag (2.55) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 10)
    bTagFileOutput[bTagSum] = " Tight JBP working point bTag (3.74) for hadronic b-jet, no constraint on leptonic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 11)
    bTagFileOutput[bTagSum] = " Tight JBP working point bTag (3.74) for leptonic b-jet, no constraint on hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 12)
    bTagFileOutput[bTagSum] = " Tight JBP working point bTag (3.74) for leptonic b-jet or hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 13)
    bTagFileOutput[bTagSum] = " Tight JBP working point bTag (3.74) for leptonic b-jet and hadronic b-jet  ";

  return bTagFileOutput[bTagSum];
  
}

std::string BTagName::NameGivingPres(int TCHEbTagLoop, int NumberTCHEbTags,  int TCHPbTagLoop, int NumberTCHPbTags, int SSVHEbTagLoop, int NumberSSVHEbTags, int SSVHPbTagLoop, int NumberSSVHPbTags,int CSVbTagLoop, int NumberCSVbTags,int JPbTagLoop, int NumberJPbTags,int JBPbTagLoop, int NumberJBPbTags){
  
  int bTagSum = TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1+JPbTagLoop-1+JBPbTagLoop-1;

  //Track counting high efficiency loop:
  if(TCHEbTagLoop == 1)
    PresOutput[bTagSum] = "No bTag";
  else if(TCHEbTagLoop == 2)
    PresOutput[bTagSum] =" L TCHE hadr" ;
  else if(TCHEbTagLoop == 3)
    PresOutput[bTagSum] =" L TCHE lept";
  else if (TCHEbTagLoop == 4)
    PresOutput[bTagSum] =" 1 L TCHE ";
  else if (TCHEbTagLoop == 5)
    PresOutput[bTagSum] =" 2 L TCHE ";
  else if(TCHEbTagLoop == 6)
    PresOutput[bTagSum] = " M TCHE hadr";
  else if (TCHEbTagLoop == 7)
    PresOutput[bTagSum] = " M TCHE lept";
  else if (TCHEbTagLoop == 8)
    PresOutput[bTagSum] = " 1 M TCHE ";  
  else if (TCHEbTagLoop == 9)
    PresOutput[bTagSum] = " 2 M TCHE ";
  else if(TCHEbTagLoop == 10)
    PresOutput[bTagSum] = " T TCHE hadr";	
  else if (TCHEbTagLoop == 11)
    PresOutput[bTagSum] = " T TCHE lept";		  
  else if (TCHEbTagLoop == 12)
    PresOutput[bTagSum] = " 1 T TCHE ";  
  else if (TCHEbTagLoop == 13)
    PresOutput[bTagSum] = " 2 T TCHE ";		
  
  //Track counting high purity loop:  
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 2)
    PresOutput[bTagSum] = " L TCHP hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 3)
    PresOutput[bTagSum] = " L TCHP lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 4)
    PresOutput[bTagSum] = " 1 L TCHP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 5)
    PresOutput[bTagSum]= " 2 L TCHP ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 6)
    PresOutput[bTagSum] = " M TCHP hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 7)
    PresOutput[bTagSum] = " M TCHP lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 8)
    PresOutput[bTagSum] = " 1 M TCHP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 9)
    PresOutput[bTagSum] = " 2 M TCHP";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 10)
    PresOutput[bTagSum] = " T TCHP hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 11)
    PresOutput[bTagSum] = " T TCHP lept ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 12)
    PresOutput[bTagSum] = " 1 T TCHP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 13)
    PresOutput[bTagSum] = " 2 T TCHP ";
  
  //Simple secondary vertex high efficiency loop:
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop ==2)
    PresOutput[bTagSum] = " L SSVHE hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 3)
    PresOutput[bTagSum] = " L SSVHE lept ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 4)
    PresOutput[bTagSum] = " 1 L SSVHE";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 5)
    PresOutput[bTagSum] = " 2 L SSVHE";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 6)
    PresOutput[bTagSum] = " M SSVHE hadr ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 7)
    PresOutput[bTagSum] = " M SSVHE lept ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 8)
    PresOutput[bTagSum] = " 1 M SSVHE";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 9)
    PresOutput[bTagSum] = " 2 M SSVHE";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 10)
    PresOutput[bTagSum] = " T SSVHE hadr ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 11)
    PresOutput[bTagSum] = " T SSVHE lept ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 12)
    PresOutput[bTagSum] = " 1 T SSVHE ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 13)
    PresOutput[bTagSum] = " 2 T SSVHE";
  
  //Simple Secondary vertex high purity loop:
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==2)
    PresOutput[bTagSum] = " L SSVHP hadr ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 3)
    PresOutput[bTagSum] = " L SSVHP lept ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 4)
    PresOutput[bTagSum] = " 1 L SSVHP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 5)
    PresOutput[bTagSum] = " 2 L SSVHP";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 6)
    PresOutput[bTagSum] = " M SSVHP hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 7)
    PresOutput[bTagSum] = " M SSVHP lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 8)
    PresOutput[bTagSum] = " 1 M SSVHP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 9)
    PresOutput[bTagSum] = " 2 M SSVHP";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 10)
    PresOutput[bTagSum] = " T SSVHP hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 11)
    PresOutput[bTagSum] = " T SSVHP lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 12)
    PresOutput[bTagSum] = " 1 T SSVHP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 13)
    PresOutput[bTagSum] = " 2 T SSVHP";

  //Combined Secondary vertex values:
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 2)
    PresOutput[bTagSum] = " L CSV hadr ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 3)
    PresOutput[bTagSum] = " L CSV lept ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 4)
    PresOutput[bTagSum] = " 1 L CSV ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 5)
    PresOutput[bTagSum] = " 2 L CSV";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 6)
    PresOutput[bTagSum] = " M CSV hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 7)
    PresOutput[bTagSum] = " M CSV lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 8)
    PresOutput[bTagSum] = " 1 M CSV";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 9)
    PresOutput[bTagSum] = " 2 M CSV";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 10)
    PresOutput[bTagSum] = " T CSV hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 11)
    PresOutput[bTagSum] = " T CSV lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 12)
    PresOutput[bTagSum] = " 1 T CSV";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 13)
    PresOutput[bTagSum] = " 2 T CSV";

  //Jet probabibility values:
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 2)
    PresOutput[bTagSum] = " L JP hadr ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 3)
    PresOutput[bTagSum] = " L JP lept ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 4)
    PresOutput[bTagSum] = " 1 L JP ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 5)
    PresOutput[bTagSum] = " 2 L JP";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1  && CSVbTagLoop == NumberCSVbTags+1 & JPbTagLoop == 6)
    PresOutput[bTagSum] = " M JP hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 7)
    PresOutput[bTagSum] = " M JP lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 8)
    PresOutput[bTagSum] = " 1 M JP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 9)
    PresOutput[bTagSum] = " 2 M JP";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 10)
    PresOutput[bTagSum] = " T JP hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 11)
    PresOutput[bTagSum] = " T JP lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 12)
    PresOutput[bTagSum] = " 1 T JP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == 13)
    PresOutput[bTagSum] = " 2 T JP";

  //Jet B probabibility values:
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 2)
    PresOutput[bTagSum] = " L JBP hadr ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 3)
    PresOutput[bTagSum] = " L JBP lept ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 4)
    PresOutput[bTagSum] = " 1 L JBP ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 5)
    PresOutput[bTagSum] = " 2 L JBP";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1  && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 & JBPbTagLoop == 6)
    PresOutput[bTagSum] = " M JBP hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 7)
    PresOutput[bTagSum] = " M JBP lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 8)
    PresOutput[bTagSum] = " 1 M JBP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 9)
    PresOutput[bTagSum] = " 2 M JBP";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 10)
    PresOutput[bTagSum] = " T JBP hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 11)
    PresOutput[bTagSum] = " T JBP lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 12)
    PresOutput[bTagSum] = " 1 T JBP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == NumberCSVbTags+1 && JPbTagLoop == NumberJPbTags+1 && JBPbTagLoop == 13)
    PresOutput[bTagSum] = " 2 T JBP";

  return PresOutput[bTagSum];
  
}
  

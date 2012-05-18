#include "../interface/BTagName.h"

BTagName::BTagName(){    
  
}

BTagName::~BTagName(){

}

std::string BTagName::NameGiving(int TCHEbTagLoop, int NumberTCHEbTags,  int TCHPbTagLoop, int NumberTCHPbTags, int SSVHEbTagLoop, int NumberSSVHEbTags, int SSVHPbTagLoop, int NumberSSVHPbTags,int CSVbTagLoop, int NumberCSVbTags){
  
  if(TCHEbTagLoop == 1)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = "No bTag constraint applied";
  else if(TCHEbTagLoop == 2)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] =" Loose TCHE working point bTag (1.7) for hadronic b-jet, no constraint on leptonic b-jet " ;
  else if(TCHEbTagLoop == 3)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] =" Loose TCHE working point bTag (1.7) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if (TCHEbTagLoop == 4)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] =" Loose TCHE working point bTag (1.7) for leptonic b-jet or hadronic b-jet ";
  else if (TCHEbTagLoop == 5)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] =" Loose TCHE working point bTag (1.7) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == 6)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium TCHE working point bTag (3.3) for hadronic b-jet, no constraint on leptonic b-jet ";
  else if (TCHEbTagLoop == 7)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium TCHE working point bTag (3.3) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if (TCHEbTagLoop == 8)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium TCHE working point bTag (3.3) for leptonic b-jet or hadronic b-jet ";  
  else if (TCHEbTagLoop == 9)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium TCHE working point bTag (3.3) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == 10)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight TCHE working point bTag (10.2) for hadronic b-jet, no constraint on leptonic b-jet ( Tight TCHE is not supported )";	
  else if (TCHEbTagLoop == 11)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight TCHE working point bTag (10.2) for leptonic b-jet, no constraint on hadronic b-jet ( Tight TCHE is not supported )";		  
  else if (TCHEbTagLoop == 12)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight TCHE working point bTag (10.2) for leptonic b-jet or hadronic b-jet ( Tight TCHE is not supported )";  
  else if (TCHEbTagLoop == 13)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight TCHE working point bTag (10.2) for leptonic b-jet and hadronic b-jet ( Tight TCHE is not supported )";		
  
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 2)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose TCHP working point bTag (1.19) for hadronic b-jet, no constraint on leptonic b-jet (Loose TCHP is not supported)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 3)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose TCHP working point bTag (1.19) for leptonic b-jet, no constraint on hadronic b-jet (Loose TCHP is not supported)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 4)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose TCHP working point bTag (1.19) for leptonic b-jet or hadronic b-jet (Loose TCHP is not supported)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 5)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose TCHP working point bTag (1.19) for leptonic b-jet and hadronic b-jet (Loose TCHP is not supported)";
  
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 6)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium TCHP working point bTag (1.93) for hadronic b-jet, no constraint on leptonic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 7)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium TCHP working point bTag (1.93) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 8)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium TCHP working point bTag (1.93) for leptonic b-jet or hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 9)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium TCHP working point bTag (1.93) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 10)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight TCHP working point bTag (3.41) for hadronic b-jet, no constraint on leptonic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 11)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight TCHP working point bTag (3.41) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 12)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight TCHP working point bTag (3.41) for leptonic b-jet or hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 13)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight TCHP working point bTag (3.41) for leptonic b-jet and hadronic b-jet ";
  
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop ==2)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose SSVHE working point bTag (0) for hadronic b-jet, no constraint on leptonic b-jet (SSVHE == 0 no POG wp)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 3)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose SSVHE working point bTag (0) for leptonic b-jet, no constraint on hadronic b-jet(SSVHE == 0 no POG wp) ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 4)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose SSVHE working point bTag (0) for leptonic b-jet or hadronic b-jet (SSVHE == 0 no POG wp)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 5)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose SSVHE working point bTag (0) for leptonic b-jet and hadronic b-jet (SSVHE == 0 no POG wp)";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 6)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium SSVHE working point bTag (1.74) for hadronic b-jet, no constraint on leptonic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 7)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium SSVHE working point bTag (1.74) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 8)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium SSVHE working point bTag (1.74) for leptonic b-jet or hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 9)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium SSVHE working point bTag (1.74) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 10)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight SSVHE working point bTag (3.05) for hadronic b-jet, no constraint on leptonic b-jet ( Tight SSVHE is not supported )";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 11)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight SSVHE working point bTag (3.05) for leptonic b-jet, no constraint on hadronic b-jet ( Tight SSVHE is not supported )";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 12)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight SSVHE working point bTag (3.05) for leptonic b-jet or hadronic b-jet ( Tight SSVHE is not supported )";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 13)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight SSVHE working point bTag (3.05) for leptonic b-jet and hadronic b-jet ( Tight SSVHE is not supported )";
  
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==2)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose SSVHP working point bTag (0) for hadronic b-jet, no constraint on leptonic b-jet (SSVHE == 0 no POG wp)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 3)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose SSVHP working point bTag (0) for leptonic b-jet, no constraint on hadronic b-jet(SSVHE == 0 no POG wp) ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 4)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose SSVHP working point bTag (0) for leptonic b-jet or hadronic b-jet (SSVHE == 0 no POG wp)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 5)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose SSVHP working point bTag (0) for leptonic b-jet and hadronic b-jet (SSVHE == 0 no POG wp)";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 6)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium SSVHP working point bTag 1) for hadronic b-jet, no constraint on leptonic b-jet (SSVHE == 1 no POG wp)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 7)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium SSVHP working point bTag (1) for leptonic b-jet, no constraint on hadronic b-jet (SSVHE == 1 no POG wp)";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 8)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium SSVHP working point bTag (1) for leptonic b-jet or hadronic b-jet(SSVHE == 1 no POG wp) ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 9)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium SSVHP working point bTag (1) for leptonic b-jet and hadronic b-jet (SSVHE == 1 no POG wp)";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 10)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight SSVHP working point bTag (2) for hadronic b-jet, no constraint on leptonic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 11)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight SSVHP working point bTag (2) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 12)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight SSVHP working point bTag (2) for leptonic b-jet or hadronic b-jet ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 13)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight SSVHP working point bTag (2) for leptonic b-jet and hadronic b-jet ";

  //Combined secondary vertex loop:
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 2)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose CSV working point bTag (0.244) for hadronic b-jet, no constraint on leptonic b-jet ";
    else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 3)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose CSV working point bTag (0.244) for leptonic b-jet, no constraint on hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 4)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose CSV working point bTag (0.244) for leptonic b-jet or hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 5)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Loose CSV working point bTag (0.244) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 6)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium CSV working point bTag (0.679) for hadronic b-jet, no constraint on leptonic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 7)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium CSV working point bTag (0.679) for leptonic b-jet, no constraint on hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 8)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium CSV working point bTag (0.679) for leptonic b-jet or hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 9)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Medium CSV working point bTag (0.679) for leptonic b-jet and hadronic b-jet ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 10)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight CSV working point bTag (0.898) for hadronic b-jet, no constraint on leptonic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 11)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight CSV working point bTag (0.898) for leptonic b-jet, no constraint on hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 12)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight CSV working point bTag (0.898) for leptonic b-jet or hadronic b-jet  ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 13)
    bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " Tight CSV working point bTag (0.898) for leptonic b-jet and hadronic b-jet  ";

  return bTagFileOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1];
  
}

std::string BTagName::NameGivingPres(int TCHEbTagLoop, int NumberTCHEbTags,  int TCHPbTagLoop, int NumberTCHPbTags, int SSVHEbTagLoop, int NumberSSVHEbTags, int SSVHPbTagLoop, int NumberSSVHPbTags,int CSVbTagLoop, int NumberCSVbTags){
  
  //Track counting high efficiency loop:
  if(TCHEbTagLoop == 1)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = "No bTag";
  else if(TCHEbTagLoop == 2)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] =" L TCHE hadr" ;
  else if(TCHEbTagLoop == 3)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] =" L TCHE lept";
  else if (TCHEbTagLoop == 4)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] =" 1 L TCHE ";
  else if (TCHEbTagLoop == 5)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] =" 2 L TCHE ";
  else if(TCHEbTagLoop == 6)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " M TCHE hadr";
  else if (TCHEbTagLoop == 7)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " M TCHE lept";
  else if (TCHEbTagLoop == 8)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 1 M TCHE ";  
  else if (TCHEbTagLoop == 9)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 2 M TCHE ";
  else if(TCHEbTagLoop == 10)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " T TCHE hadr";	
  else if (TCHEbTagLoop == 11)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " T TCHE lept";		  
  else if (TCHEbTagLoop == 12)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 1 T TCHE ";  
  else if (TCHEbTagLoop == 13)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 2 T TCHE ";		
  
  //Track counting high purity loop:  
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 2)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " L TCHP hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 3)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " L TCHP lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 4)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 1 L TCHP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 5)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1]= " 2 L TCHP ";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 6)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " M TCHP hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 7)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " M TCHP lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 8)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 1 M TCHP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 9)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 2 M TCHP";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 10)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " T TCHP hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 11)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " T TCHP lept ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 12)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 1 T TCHP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == 13)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 2 T TCHP ";
  
  //Simple secondary vertex high efficiency loop:
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop ==2)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " L SSVHE hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 3)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " L SSVHE lept ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 4)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 1 L SSVHE";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 5)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 2 L SSVHE";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 6)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " M SSVHE hadr ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 7)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " M SSVHE lept ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 8)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 1 M SSVHE";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 9)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 2 M SSVHE";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 10)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " T SSVHE hadr ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 11)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " T SSVHE lept ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 12)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 1 T SSVHE ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == 13)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 2 T SSVHE";
  
  //Simple Secondary vertex high purity loop:
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==2)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " L SSVHP hadr ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 3)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " L SSVHP lept ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 4)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 1 L SSVHP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 5)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 2 L SSVHP";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 6)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " M SSVHP hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 7)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " M SSVHP lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 8)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 1 M SSVHP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 9)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 2 M SSVHP";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 10)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " T SSVHP hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 11)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " T SSVHP lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 12)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 1 T SSVHP";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == 13)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 2 T SSVHP";

  //Combined Secondary vertex values:
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop ==NumberSSVHPbTags+1 && CSVbTagLoop == 2)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " L CSV hadr ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 3)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " L CSV lept ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 4)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 1 L CSV ";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 5)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 2 L CSV";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 6)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " M CSV hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 7)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " M CSV lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 8)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 1 M CSV";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 9)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 2 M CSV";
  else if(TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 10)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " T CSV hadr";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 11)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " T CSV lept";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 12)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 1 T CSV";
  else if (TCHEbTagLoop == NumberTCHEbTags+1 && TCHPbTagLoop == NumberTCHPbTags+1 && SSVHEbTagLoop == NumberSSVHEbTags+1 && SSVHPbTagLoop == NumberSSVHPbTags+1 && CSVbTagLoop == 13)
    PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1] = " 2 T CSV";


  return PresOutput[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1];
  
}
  

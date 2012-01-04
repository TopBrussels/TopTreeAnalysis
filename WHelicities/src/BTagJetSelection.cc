#include "../interface/BTagJetSelection.h"

using namespace std;

BTagJetSelection::BTagJetSelection(){
 
}

BTagJetSelection::~BTagJetSelection(){

}

int BTagJetSelection::HighestProbSelection(int bTagLoop, int ConsideredBTagger, vector<float> KinFitProb, vector<float> MlbProb, vector<float> btagTCHE, vector<float> btagTCHP, vector<float> btagSSVHE, vector<float> btagSSVHP){

  float ProbBTag[4][4][4];
  int BLeptIndex[12];
  int BHadrIndex[12];
  
   //Fist index represents used btagger (0 = tche, 1 = tchp, 2 = ssvhe & 3 = ssvhp)
  //Second index represents working point (0 = no btag, 1 = loose, 2 = medium & 3 = tight)
  //Third index represents jet numbering
  for(int ii=0; ii<4; ii++){
    for(int jj=0;jj<4;jj++){
      for(int kk=0;kk<4;kk++){
	ProbBTag[ii][jj][kk]=1.;
      }
    }
  }
  
  BLeptIndex[0]=3;
  BLeptIndex[1]=2;
  BLeptIndex[2]=3;
  BLeptIndex[3]=1;
  BLeptIndex[4]=2;
  BLeptIndex[5]=1;
  BLeptIndex[6]=3;
  BLeptIndex[7]=0;
  BLeptIndex[8]=2;
  BLeptIndex[9]=0;
  BLeptIndex[10]=1;
  BLeptIndex[11]=0;
  BHadrIndex[0]=2;
  BHadrIndex[1]=3;
  BHadrIndex[2]=1;
  BHadrIndex[3]=3;
  BHadrIndex[4]=1;
  BHadrIndex[5]=2;
  BHadrIndex[6]=0;
  BHadrIndex[7]=3;
  BHadrIndex[8]=0;
  BHadrIndex[9]=2;
  BHadrIndex[10]=0;
  BHadrIndex[11]=1;
  
   //Values that not pass the btag constraint get probability == 0
  for(int ii=0; ii<4; ii++){
    if(btagTCHE[ii]<1.7){ProbBTag[0][1][ii]=0.;}
    if(btagTCHE[ii]<3.3){ProbBTag[0][2][ii]=0.;}
    if(btagTCHE[ii]<10.2){ProbBTag[0][3][ii]=0.;}
    if(btagTCHP[ii]<1.19){ProbBTag[1][1][ii]=0.;}
    if(btagTCHP[ii]<1.93){ProbBTag[1][2][ii]=0.;}
    if(btagTCHP[ii]<3.41){ProbBTag[1][3][ii]=0.;}
    if(btagSSVHE[ii]<0.){ProbBTag[2][1][ii]=0.;}
    if(btagSSVHE[ii]<1.74){ProbBTag[2][2][ii]=0.;}
    if(btagSSVHE[ii]<3.05){ProbBTag[2][3][ii]=0.;}
    if(btagSSVHP[ii]<0.){ProbBTag[3][1][ii]=0.;}
    if(btagSSVHP[ii]<1.){ProbBTag[3][2][ii]=0.;}
    if(btagSSVHP[ii]<2.){ProbBTag[3][3][ii]=0.;}
  }
  
  int MaximumProbability=0;
  int WorkingPointNumber=0;
  int SelectedJetCombination = 999;
  float TotalProbability[12];
  if(bTagLoop >1 && bTagLoop <= 5){WorkingPointNumber=1;}
  else if (bTagLoop >=6 && bTagLoop <= 9){WorkingPointNumber=2;}
  else if(bTagLoop >=10 && bTagLoop <=13){WorkingPointNumber=3;}
  
  if(bTagLoop == 1){  //No btag applied
    for(int ii=0;ii<12;ii++){
      TotalProbability[ii] = KinFitProb[ii]*MlbProb[BLeptIndex[ii]]*ProbBTag[ConsideredBTagger][WorkingPointNumber][0];
      if(MaximumProbability<TotalProbability[ii]){
	MaximumProbability=TotalProbability[ii];
	SelectedJetCombination=ii;
      }
    }
  }
  else if(bTagLoop == 2 || bTagLoop == 6 || bTagLoop == 10){  //1 btag on hadr jet
    for(int ii=0;ii<12;ii++){
      TotalProbability[ii] = KinFitProb[ii]*MlbProb[BLeptIndex[ii]]*ProbBTag[ConsideredBTagger][WorkingPointNumber][BHadrIndex[ii]];
      if(MaximumProbability<TotalProbability[ii]){
	MaximumProbability=TotalProbability[ii];
	SelectedJetCombination=ii;
      }
    }
  }
  else if(bTagLoop == 3 || bTagLoop == 7 || bTagLoop == 11){  //1 btag on lept jet
    for(int ii=0;ii<12;ii++){
      TotalProbability[ii] = KinFitProb[ii]*MlbProb[BLeptIndex[ii]]*ProbBTag[ConsideredBTagger][WorkingPointNumber][BLeptIndex[ii]];
      if(MaximumProbability<TotalProbability[ii]){
	MaximumProbability=TotalProbability[ii];
	SelectedJetCombination=ii;
      }
    }
  }
  else if(bTagLoop == 4 || bTagLoop == 8 || bTagLoop == 12){  //1 btag 
    for(int ii=0;ii<12;ii++){
      float bTagProb=0.;
      if(ProbBTag[ConsideredBTagger][WorkingPointNumber][BHadrIndex[ii]] == 1 || ProbBTag[ConsideredBTagger][WorkingPointNumber][BLeptIndex[ii]] ==1)
	bTagProb = 1;
      TotalProbability[ii] = KinFitProb[ii]*MlbProb[BLeptIndex[ii]]*bTagProb;
      if(MaximumProbability<TotalProbability[ii]){
	MaximumProbability=TotalProbability[ii];
	SelectedJetCombination=ii;
      }
    }
  }
  else if(bTagLoop == 5 || bTagLoop == 9 || bTagLoop == 13){  //2 btags
    for(int ii=0;ii<12;ii++){
      TotalProbability[ii] = KinFitProb[ii]*MlbProb[BLeptIndex[ii]]*ProbBTag[ConsideredBTagger][WorkingPointNumber][BHadrIndex[ii]]*ProbBTag[ConsideredBTagger][WorkingPointNumber][BLeptIndex[ii]];
      if(MaximumProbability<TotalProbability[ii]){
	MaximumProbability=TotalProbability[ii];
	SelectedJetCombination=ii;
      }
    }
  }
  
  return SelectedJetCombination;
}

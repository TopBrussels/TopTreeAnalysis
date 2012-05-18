
#include "../interface/BTagCosThetaCalculation.h"

BTagCosThetaCalculation::BTagCosThetaCalculation(){

}

BTagCosThetaCalculation::~BTagCosThetaCalculation(){

}

float BTagCosThetaCalculation::Calculation(TLorentzVector lepton, TLorentzVector Neutrino, TLorentzVector leptonicBJet){

  float CosTheta = 999;	  
  //----------------------------------------------
  //  Calculating cos theta value
  //----------------------------------------------
  
  TLorentzVector WLeptonic = (Neutrino+lepton);
  TLorentzVector TopLeptonic = (Neutrino+lepton+leptonicBJet);
  
  //Reboost the particles to rest frames
  TLorentzVector MuonWZMF = lepton; // In W Zero Mass Frame (WZMF)
  TLorentzVector WLeptonicTZMF = WLeptonic;  // In Top Zero Mass Frame (TZMF)	  
  
  MuonWZMF.Boost(-WLeptonic.BoostVector());
  WLeptonicTZMF.Boost(-TopLeptonic.BoostVector());
  
  //Calculating cos:	      
  CosTheta = ((WLeptonicTZMF.Vect()).Dot(MuonWZMF.Vect()))/(((WLeptonicTZMF.Vect()).Mag())*((MuonWZMF.Vect()).Mag()));       	
  
  return CosTheta;	
}        

float BTagCosThetaCalculation::CalcOrigKins(int BLeptonicIndex, int BHadronicIndex,  TLorentzVector muon, vector<TLorentzVector> selectedJets, float MassW, float MassTop){
  
  float CosThetaOrigKins = 999; 
  //---------------------------------------------------------------
  //     Calculating MET_Pz() (equation ax² + bx + c = 0 ):     
  //---------------------------------------------------------------
  //
  // MET_Px() and MET_Py() is known since there is no Px() and Py() component before the collision.
  // The Pz() component before the collision is not known for proton-proton collisions.
  //
  float NeutrinoPx = -(muon +selectedJets[0]+ selectedJets[1] + selectedJets[2] + selectedJets[3] ).Px();
  float NeutrinoPy = -(muon +selectedJets[0]+ selectedJets[1] + selectedJets[2] + selectedJets[3] ).Py();	
  float NeutrinoPz=999;
  float NeutrinoE=999;
  Neutrino.SetPxPyPzE(NeutrinoPx,NeutrinoPy,NeutrinoPz,NeutrinoE);  //Initialize in case no neutrino can be reconstructed!!	  
      
  //Calculating solutions for quadratic equation  
  float aCoefficient = 4*pow(muon.E(),2)-4*pow(muon.Pz(),2);
  float bCoefficient = 4*(muon.Pz())*(pow(muon.M(),2)-pow(MassW,2)-2*(muon.Px())*NeutrinoPx-2*(muon.Py())*NeutrinoPy);
  float cCoefficient = -pow(muon.M(),4)-pow(MassW,4)-4*pow(muon.Px(),2)*pow(NeutrinoPx,2)-4*pow(muon.Py(),2)*pow(NeutrinoPy,2)+4*(pow(muon.M(),2)-pow(MassW,2))*((muon.Px())*NeutrinoPx+(muon.Py())*NeutrinoPy)-8*(muon.Px())*NeutrinoPx*(muon.Py())*NeutrinoPy+4*pow(muon.E(),2)*pow(NeutrinoPx,2)+4*pow(muon.E(),2)*pow(NeutrinoPy,2)+2*(pow(MassW,2))*(pow(muon.M(),2));
      
  float DCoefficient = pow(bCoefficient,2)-4*aCoefficient*cCoefficient;           
  if(DCoefficient>=0){	  //Still need to find solution for D<0 !!!	      
    float NeutrinoPzOne = ((-bCoefficient + sqrt(DCoefficient))/(aCoefficient*2));
    float NeutrinoPzTwo = ((-bCoefficient - sqrt(DCoefficient))/(aCoefficient*2));	  
	
    float NeutrinoEOne = sqrt(NeutrinoPx*NeutrinoPx+NeutrinoPy*NeutrinoPy+NeutrinoPzOne*NeutrinoPzOne);
    float NeutrinoETwo = sqrt(NeutrinoPx*NeutrinoPx+NeutrinoPy*NeutrinoPy+NeutrinoPzTwo*NeutrinoPzTwo);
    TLorentzVector NeutrinoOne;
    NeutrinoOne.SetPxPyPzE(NeutrinoPx,NeutrinoPy,NeutrinoPzOne,NeutrinoEOne);  //Has to be initialized like this !
    TLorentzVector NeutrinoTwo;
    NeutrinoTwo.SetPxPyPzE(NeutrinoPx,NeutrinoPy,NeutrinoPzTwo,NeutrinoETwo);
	
    float TopMassDiffOne = (MassTop -((NeutrinoOne+muon+selectedJets[BLeptonicIndex]).M()) );	  
    float TopMassDiffTwo = (MassTop -((NeutrinoTwo+muon+selectedJets[BLeptonicIndex]).M()) );	   
	
    //-----   Selected neutrino has the smallest deviation from the top mass   -----
    //
    if(fabs(TopMassDiffOne)<fabs(TopMassDiffTwo)){
      NeutrinoPz=NeutrinoPzOne;
      NeutrinoE = NeutrinoEOne;
    }
    else{
      NeutrinoPz=NeutrinoPzTwo;
      NeutrinoE = NeutrinoETwo;
    }	  	  
    Neutrino.SetPxPyPzE(NeutrinoPx,NeutrinoPy,NeutrinoPz,NeutrinoE);

    //----------------------------------------------
    //  Calculating cos theta value
    //----------------------------------------------
	
    TLorentzVector WLeptonic = (Neutrino+muon);
    TLorentzVector TopLeptonic = (Neutrino+muon+selectedJets[BLeptonicIndex]);

    //histo1D["TopMassLeptReco"]->Fill(TopLeptonic.M());
    //histo1D["WMassLeptReco"]->Fill(WLeptonic.M());
	
    //Reboost the particles to rest frames 
    //
    TLorentzVector MuonWZMF = muon; // In W Zero Mass Frame (WZMF)
    TLorentzVector WLeptonicTZMF = WLeptonic;  // In Top Zero Mass Frame (TZMF)	  
	
    MuonWZMF.Boost(-WLeptonic.BoostVector());
    WLeptonicTZMF.Boost(-TopLeptonic.BoostVector());
	
    //Calculating cos:	      
    CosThetaOrigKins = ((WLeptonicTZMF.Vect()).Dot(MuonWZMF.Vect()))/(((WLeptonicTZMF.Vect()).Mag())*((MuonWZMF.Vect()).Mag()));       	

  }  //end of D>=0 loop   
	
  return CosThetaOrigKins;
	
}        




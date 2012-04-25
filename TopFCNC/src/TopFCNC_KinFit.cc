#include "../interface/TopFCNC_KinFit.h"

TopFCNC_KinFit::TopFCNC_KinFit(Dataset* dataset, ResolutionFit *resFitLeptons, ResolutionFit *resFitLightJets, ResolutionFit *resFitBJets){

  dataset_ = dataset;
  
  resFitLeptons_   = resFitLeptons;
  resFitLightJets_ = resFitLightJets;
  resFitBJets_     = resFitBJets;

  Prob_ = -1.;
  Chi2_ =  0.;
  Ndof_ =  0;

  verbose_ = false;

  if(dataset_)
  {
    // Make some plots
/*
    string mWUnCorrtitle = "mWUnCorr_"+dataset_->Name();
    string mTopUnCorrtitle = "mTopUnCorr_"+dataset_->Name();
    string mWCorrtitle = "mWCorr_"+dataset_->Name();
    string mTopCorrtitle = "mTopCorr_"+dataset_->Name();
    histo1D_[mWUnCorrtitle] = new TH1F(mWUnCorrtitle.c_str(),mWUnCorrtitle.c_str(),50,0,150);
    histo1D_[mTopUnCorrtitle] = new TH1F(mTopUnCorrtitle.c_str(),mTopUnCorrtitle.c_str(),50,100,250);
    histo1D_[mWCorrtitle] = new TH1F(mWCorrtitle.c_str(),mWCorrtitle.c_str(),50,0,150);
    histo1D_[mTopCorrtitle] = new TH1F(mTopCorrtitle.c_str(),mTopCorrtitle.c_str(),50,100,250);
*/
  }
}

TopFCNC_KinFit::~TopFCNC_KinFit(){
}

void TopFCNC_KinFit::FitEvent(TopFCNC_Evt *topFCNC_Evt, float wMass, float zMass, float topMass)
{

  Prob_ = -1.;
  Chi2_ =  0.;
  Ndof_ =  0;

  vector<TRootJet> jets  = topFCNC_Evt->selectedJets();
  if(verbose_)
	  cout<<"- Number of selected jets : "<<jets.size()<<endl;

  UInt_t *numbers = new UInt_t[jets.size()];
  for(UInt_t i=0;i<jets.size();i++) numbers[i]=i;

  TLorentzVector lepton1FromZ = topFCNC_Evt->lepton1FromZ();
  TLorentzVector lepton2FromZ = topFCNC_Evt->lepton2FromZ();

  if(topFCNC_Evt->isDiLeptonic()){
    // Topology to reconstruct : tt-> bW + qZ -> bqq + qll

    UInt_t NofJets = 4 ;

    UInt_t *comb = new UInt_t[NofJets];
    for(UInt_t i=0;i <NofJets;i++) comb[i]=i;
	  UInt_t *Permutation = new UInt_t[NofJets];
    for(UInt_t i=0;i<NofJets;i++) Permutation[i]=i;

    Double_t Prob_tmp = 0.;
    TMatrixD Ml1(3,3), Ml2(3,3); // for leptons
    TMatrixD Mb(3,3),  Mq(3,3), Mj1(3,3), Mj2(3,3); // for jets

    do{
      if(verbose_)
		    cout<<"-- Jet combination considered : "<<comb[0]<<"/"<<comb[1]<<"/"<<comb[2]<<"/"<<comb[3]<<endl;

		  Prob_tmp = 0.;

	    do{
		    if(verbose_)
		  	  cout<<"--- Permutations : "<<comb[0]<<"/"<<comb[1]<<"/"<<comb[2]<<"/"<<comb[3]<<endl;
  
        Ml1.Zero();
        Ml2.Zero();
        Mb.Zero();
        Mq.Zero();
        Mj1.Zero();
        Mj2.Zero();
  
        TLorentzVector bJet      = jets[comb[0]];
        TLorentzVector qJet      = jets[comb[1]];
        TLorentzVector lightJet1 = jets[comb[2]];
        TLorentzVector lightJet2 = jets[comb[3]];
  
        Ml1(0,0) = pow(resFitLeptons_->EtResolution(&lepton1FromZ), 2);
        Ml1(1,1) = pow(resFitLeptons_->ThetaResolution(&lepton1FromZ), 2);
        Ml1(2,2) = pow(resFitLeptons_->PhiResolution(&lepton1FromZ), 2);
    
        Ml2(0,0) = pow(resFitLeptons_->EtResolution(&lepton2FromZ), 2);
        Ml2(1,1) = pow(resFitLeptons_->ThetaResolution(&lepton2FromZ), 2);
        Ml2(2,2) = pow(resFitLeptons_->PhiResolution(&lepton2FromZ), 2);

        Mb(0,0)  = pow(resFitBJets_->EtResolution(&bJet), 2);
        Mb(1,1)  = pow(resFitBJets_->ThetaResolution(&bJet), 2);
        Mb(2,2)  = pow(resFitBJets_->PhiResolution(&bJet), 2);
 
        Mq(0,0)  = pow(resFitLightJets_->EtResolution(&qJet), 2);
        Mq(1,1)  = pow(resFitLightJets_->ThetaResolution(&qJet), 2);
        Mq(2,2)  = pow(resFitLightJets_->PhiResolution(&qJet), 2);

        Mj1(0,0) = pow(resFitLightJets_->EtResolution(&lightJet1), 2);
        Mj1(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet1), 2);
        Mj1(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet1), 2);
    
        Mj2(0,0) = pow(resFitLightJets_->EtResolution(&lightJet2), 2);
        Mj2(1,1) = pow(resFitLightJets_->ThetaResolution(&lightJet2), 2);
        Mj2(2,2) = pow(resFitLightJets_->PhiResolution(&lightJet2), 2);
		    if(verbose_)
		  	  cout<<"---- Correlation matrices instantiated "<<endl;

        TKinFitter *kinfit_    = new TKinFitter("topFCNC_Fit", "topFCNC_Fit");
        kinfit_->setVerbosity(0);

        TFitParticleEtThetaPhiEMomFix *fitLep1   = new TFitParticleEtThetaPhiEMomFix("fitLep1", "fitLep1", &lepton1FromZ, &Ml1);
        TFitParticleEtThetaPhiEMomFix *fitLep2   = new TFitParticleEtThetaPhiEMomFix("fitLep2", "fitLep2", &lepton2FromZ, &Ml1);

        TFitParticleEtThetaPhiEMomFix *fitB      = new TFitParticleEtThetaPhiEMomFix("bJet", "bJet", &bJet, &Mb);
        TFitParticleEtThetaPhiEMomFix *fitQ      = new TFitParticleEtThetaPhiEMomFix("qJet", "qJet", &qJet, &Mq);
        TFitParticleEtThetaPhiEMomFix *fitLight1 = new TFitParticleEtThetaPhiEMomFix("lightJet1", "lightJet1", &lightJet1, &Mj1);
        TFitParticleEtThetaPhiEMomFix *fitLight2 = new TFitParticleEtThetaPhiEMomFix("lightJet2", "lightJet2", &lightJet2, &Mj2);

        kinfit_->addMeasParticles(fitLep1,fitLep2,fitB,fitQ,fitLight1,fitLight2);

        TFitConstraintM *consHadW     = new TFitConstraintM("HadWMass",      "MassConstraint", 0, 0,   wMass );
        TFitConstraintM *consLepZ     = new TFitConstraintM("LepZMass",      "MassConstraint", 0, 0,   zMass );
        TFitConstraintM *consHadTop   = new TFitConstraintM("Had_TopMass",   "MassConstraint", 0, 0, topMass );
        TFitConstraintM *consFcncTop  = new TFitConstraintM("FCNC_TopMass",  "MassConstraint", 0, 0, topMass );
        TFitConstraintM *consEqualTop = new TFitConstraintM("EqualTopMasses","EqualTopMasses", 0, 0,       0.);

        consHadW->addParticles1(fitLight1,fitLight2);
        consLepZ->addParticles1(fitLep1,fitLep2);

        consHadTop->addParticles1(fitB,fitLight1,fitLight2);

        consFcncTop->addParticles1(fitQ,fitLep1,fitLep2);

        consEqualTop->addParticles1(fitB,fitLight1,fitLight2);
        consEqualTop->addParticles2(fitQ,fitLep1,fitLep2);
                    
        kinfit_->addConstraint(consHadW);
        kinfit_->addConstraint(consLepZ);
        kinfit_->addConstraint(consHadTop);
        kinfit_->addConstraint(consFcncTop);
        kinfit_->addConstraint(consEqualTop);
        kinfit_->setMaxNbIter(30);
        kinfit_->setMaxDeltaS(5e-5);
        kinfit_->setMaxF(1e-4);
      
        //do the fit!
        kinfit_->fit();
        if(kinfit_->getStatus() == 0) // if the fitter converged
          Prob_tmp = TMath::Prob(kinfit_->getS(), kinfit_->getNDF());
        else
          Prob_tmp = -1.;

        if(Prob_tmp>Prob_){
          for(UInt_t i=0;i<NofJets;i++) Permutation[i] = comb[i];
          Prob_ = Prob_tmp;
          Chi2_ = kinfit_->getS();
          Ndof_ = kinfit_->getNDF();
        }

        delete kinfit_;

        delete fitLep1;
        delete fitLep2;

        delete fitB;
        delete fitQ;
        delete fitLight1;
        delete fitLight2;

        delete consHadW;
        delete consLepZ;
        delete consHadTop;
        delete consFcncTop;
        delete consEqualTop;

			}
		  while(next_permutation(comb,comb+NofJets));
   	}
    while(next_combination(numbers,numbers+jets.size(),comb,comb+NofJets));

    topFCNC_Evt->SetB(TRootJet(jets[Permutation[0]]));
    topFCNC_Evt->SetQ(TRootJet(jets[Permutation[1]]));
    topFCNC_Evt->SetQuarkFromW(TRootJet(jets[Permutation[2]]));
    topFCNC_Evt->SetQuarkBarFromW(TRootJet(jets[Permutation[3]]));
    
    delete Permutation;
    delete comb;
  }
  else{
    // Topology to reconstruct : tt-> bW + qZ -> blv + qll
    Prob_ = 0;
    Chi2_ = 0;
    Ndof_ = 0;
  }
/*
  if(writePNG)
  {
    stringstream s1; s1 << event->runId();
    stringstream s2; s2 << event->lumiBlockId();
    stringstream s3; s3 << event->eventId();
    stringstream s4; s4 << jetCombi;
    string monsterName = "Monster_Data_" + s1.str() + "_" + s2.str() + "_" + s3.str() + "_" + s4.str();
    
    mkdir("PlotsJES",0777);
    string path = "PlotsJES/Monsters/";
    mkdir(path.c_str(),0777);
    
    if(measureTopMassDiff_)
    {
      TH1F* histo1D = (TH1F*) histo->ProjectionY(monsterName.c_str());
      TCanvas* tempCanvas = TCanvasCreator(histo1D, monsterName);
//      tempCanvas->SaveAs( (path + monsterName + ".root").c_str() );
      tempCanvas->SaveAs( (path + monsterName + ".png").c_str() );
      delete histo1D;
    }
    else
    {
      TCanvas* tempCanvas = TCanvasCreator(histo, monsterName);
//      tempCanvas->SaveAs( (path + monsterName + ".root").c_str() );
      tempCanvas->SaveAs( (path + monsterName + ".png").c_str() );
    }
  }
*/
}

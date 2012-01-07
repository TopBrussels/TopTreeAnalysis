#include "../interface/LightJetCombiner.h"

bool sortMVAValuesAgain (MVAValues f1, MVAValues f2)
{
  return f1.MVAValue > f2.MVAValue;
}

LightJetCombiner::LightJetCombiner(bool trainMVA, float Luminosity, const vector<Dataset*>& datasets, bool measureTopMass, string MVAMethod, string postfix)
{
  trainMVA_ = trainMVA;
	postfix_ = postfix;
  trainer_ = 0;
  computer_ = 0;
  
  measureTopMass_ = measureTopMass;
  EndJobWasRun_ = false;
  
  Luminosity_ = Luminosity;
  datasets_ = datasets;
  
  string MVAPrefix = "LightJetCombTrainer";
  string MVAOut = "MVAOutput_TTbarJES.root";
	
  FullKinFit* kinFit = new FullKinFit(0, 0, 0, measureTopMass_, false, true);
  dummyMonster_ = kinFit->DummyMonster();
  delete kinFit;
  
  if(trainMVA_)
  {
   	cout << "The MVA stuff will be trained!" << endl;
   	trainer_ = new MVATrainer(MVAMethod,MVAPrefix,MVAOut);
//    trainer_->addMethod("LikelihoodD");
//    trainer_->addMethod("LikelihoodPCA");
//    trainer_->addMethod("MLP");
//    trainer_->addMethod("MLPD");

    trainer_->bookInputVar("btag");
//    trainer_->bookInputVar("ThPtOverSumPt");
//    trainer_->bookInputVar("AngleThWh");
    trainer_->bookInputVar("AngleBlMu");
//    trainer_->bookInputVar("AngleThBl");
//    trainer_->bookInputVar("AngleThBh");
//    trainer_->bookInputVar("AngleThMu");
  }
  else
  {
    cout << "The MVA stuff will be computed!" << endl;
    
    // First prepare some plots
    MSPlot_["MVA_input_btag"] = new MultiSamplePlot(datasets, "MVA_input_btag", 100, -0.0001, 1, "btag");
    MSPlot_["MVA_input_ThPtOverSumPt"] = new MultiSamplePlot(datasets, "MVA_input_ThPtOverSumPt", 100, -0.0001, 1, "ThPtOverSumPt");
    MSPlot_["MVA_input_AngleThWh"] = new MultiSamplePlot(datasets, "MVA_input_AngleThWh", 80, -0.0001, 3.2, "AngleThWh");
    MSPlot_["MVA_input_AngleBlMu"] = new MultiSamplePlot(datasets, "MVA_input_AngleBlMu", 80, -0.0001, 3.2, "AngleBlMu");
    MSPlot_["MVA_input_AngleThBl"] = new MultiSamplePlot(datasets, "MVA_input_AngleThBl", 80, -0.0001, 3.2, "AngleThBl");
    MSPlot_["MVA_input_AngleThBh"] = new MultiSamplePlot(datasets, "MVA_input_AngleThBh", 80, -0.0001, 3.2, "AngleThBh");
    MSPlot_["MVA_input_AngleThMu"] = new MultiSamplePlot(datasets, "MVA_input_AngleThMu", 80, -0.0001, 3.2, "AngleThMu");
    
    MSPlot_["maxMVA_btag"] = new MultiSamplePlot(datasets, "maxMVA_btag", 100, -0.0001, 1, "btag");
    MSPlot_["maxMVA_ThPtOverSumPt"] = new MultiSamplePlot(datasets, "maxMVA_ThPtOverSumPt", 100, -0.0001, 1, "ThPtOverSumPt");
    MSPlot_["maxMVA_AngleThWh"] = new MultiSamplePlot(datasets, "maxMVA_AngleThWh", 80, -0.0001, 3.2, "AngleThWh");
    MSPlot_["maxMVA_AngleBlMu"] = new MultiSamplePlot(datasets, "maxMVA_AngleBlMu", 80, -0.0001, 3.2, "AngleBlMu");
    MSPlot_["maxMVA_AngleThBl"] = new MultiSamplePlot(datasets, "maxMVA_AngleThBl", 80, -0.0001, 3.2, "AngleThBl");
    MSPlot_["maxMVA_AngleThBh"] = new MultiSamplePlot(datasets, "maxMVA_AngleThBh", 80, -0.0001, 3.2, "AngleThBh");
    MSPlot_["maxMVA_AngleThMu"] = new MultiSamplePlot(datasets, "maxMVA_AngleThMu", 80, -0.0001, 3.2, "AngleThMu");
    MSPlot_["maxMVA_HadrWMass"] = new MultiSamplePlot(datasets, "maxMVA_HadrWMass", 100, -0.0001, 200, "HadrWmass");
    MSPlot_["maxMVA_HadrTopMass"] = new MultiSamplePlot(datasets, "maxMVA_HadrTopMass", 100, -0.0001, 1000, "Hadronic Top Mass");
    
	  vector<string> MVAvars;
	  MVAvars.push_back("btag");
//	  MVAvars.push_back("ThPtOverSumPt");
//    MVAvars.push_back("AngleThWh");
    MVAvars.push_back("AngleBlMu");
//    MVAvars.push_back("AngleThBl");
//    MVAvars.push_back("AngleThBh");
//	  MVAvars.push_back("AngleThMu");
	  
	  computer_ = new MVAComputer(MVAMethod,MVAOut,MVAPrefix,MVAvars,postfix_);
	//    computer_->addMethod("LikelihoodD");
	//    computer_->addMethod("LikelihoodPCA");
	//    computer_->addMethod("MLP");
	//    computer_->addMethod("MLPD");
  }
  
  // Plots for the correlation
  histo2D_["hadronicWMassVSbtag"] = new TH2F("hadronicWMassVSbtag","Hadronic W Mass VS btag product",100,0,150,100,-0.0001,1);
  histo2D_["hadronicWMassVSAngleThWh"] = new TH2F("hadronicWMassVSAngleThWh","Hadronic W Mass VS Angle Th Wh",100,0,150,100,-0.0001,3.2);
  histo2D_["hadronicWMassVSAngleBlMu"] = new TH2F("hadronicWMassVSAngleBlMu","Hadronic W Mass VS Angle Bl Mu",100,0,150,100,-0.0001,3.2);
  histo2D_["hadronicWMassVSAngleThBl"] = new TH2F("hadronicWMassVSAngleThBl","Hadronic W Mass VS Angle Th Bl",100,0,150,100,-0.0001,3.2);
  histo2D_["hadronicWMassVSAngleThBh"] = new TH2F("hadronicWMassVSAngleThBh","Hadronic W Mass VS Angle Th Bh",100,0,150,100,-0.0001,3.2);
  histo2D_["hadronicWMassVSAngleThMu"] = new TH2F("hadronicWMassVSAngleThMu","Hadronic W Mass VS Angle Th Mu",100,0,150,100,-0.0001,3.2);
  histo2D_["hadronicWMassVSThPtOverSumPt"] = new TH2F("hadronicWMassVSThPtOverSumPt","Hadronic W Mass VS ThPtOverSumPt",100,0,150,100,-0.0001,1);
      
  histo2D_["hadronicTopMassVSbtag"] = new TH2F("hadronicTopMassVSbtag","Hadronic Top Mass VS btag product",100,100,250,100,-0.0001,1);
  histo2D_["hadronicTopMassVSAngleThWh"] = new TH2F("hadronicTopMassVSAngleThWh","Hadronic Top Mass VS Angle Th Wh",100,100,250,100,-0.0001,3.2);
  histo2D_["hadronicTopMassVSAngleBlMu"] = new TH2F("hadronicTopMassVSAngleBlMu","Hadronic Top Mass VS Angle Bl Mu",100,100,250,100,-0.0001,3.2);
  histo2D_["hadronicTopMassVSAngleThBl"] = new TH2F("hadronicTopMassVSAngleThBl","Hadronic Top Mass VS Angle Th Bl",100,100,250,100,-0.0001,3.2);
  histo2D_["hadronicTopMassVSAngleThBh"] = new TH2F("hadronicTopMassVSAngleThBh","Hadronic Top Mass VS Angle Th Bh",100,100,250,100,-0.0001,3.2);
  histo2D_["hadronicTopMassVSAngleThMu"] = new TH2F("hadronicTopMassVSAngleThMu","Hadronic Top Mass VS Angle Th Mu",100,100,250,100,-0.0001,3.2);
  histo2D_["hadronicTopMassVSThPtOverSumPt"] = new TH2F("hadronicTopMassVSThPtOverSumPt","Hadronic Top Mass VS ThPtOverSumPt",100,100,250,100,-0.0001,1);
}

LightJetCombiner::~LightJetCombiner()
{
  if(dummyMonster_) delete dummyMonster_;
  if(trainer_) delete trainer_;
  if(computer_) delete computer_;
}

void LightJetCombiner::ProcessEvent(Dataset* dataSet, Monster* monster, float scaleFactor)
{
  vectorMVA_.clear();
  string dataSetName = dataSet->Name();
  
  ///////////////////////
  // Train/Compute MVA //
  ///////////////////////
  
  float previousMaxMVA = -9999;
  
  float maxMVA_btag = -9999;
  float maxMVA_AngleThWh = -9999;
  float maxMVA_AngleBlMu = -9999;
  float maxMVA_AngleThBl = -9999;
  float maxMVA_AngleThBh = -9999;
  float maxMVA_AngleThMu = -9999;
  float maxMVA_ThPtOverSumPt = -9999;
//  vector<float> bTag = monster->bTagTCHE();
  vector<float> bTag = monster->bTagSSVHE();
  vector<TLorentzVector> selectedJets = monster->selectedJets();
  TLorentzVector Mu = monster->lepton();
  
  // index convention -> i,j: jets from Hadronic W  k: Hadronic b and l: leptonic b
  for (int i=0; i<4; i++)
  {
    for (int j=0; j<4; j++)
    {
 	    for (int k=0; k<4; k++)
 	    {
 	      for (int l=0; l<4; l++)
 	      {
          if (i < j && i != j && i != k && i != l && j != k && j != l && k != l)
          {
//            D = log(1 + | L3D |/σL3D )

            //calculate the vars
            // btag
            float btag_i = 0;
            if(bTag[i] > -0.090) btag_i = exp(bTag[i]) - 1;
            float btag_j = 0;
            if(bTag[j] > -0.090) btag_j = exp(bTag[j]) - 1;
            float btag_k = 0;
            if(bTag[k] > -0.090) btag_k = exp(bTag[k]) - 1;
            float btag_l = 0;
            if(bTag[l] > -0.090) btag_l = exp(bTag[l]) - 1;
            
//            float btag = btag_k + btag_l;
//            btag = btag / ( btag_i + btag_j + btag_k + btag_l );
            
//            cout << btag << endl;

            float btag = pow(btag_k,2) + pow(btag_l,2);
            btag = btag / ( pow(btag_i,2) + pow(btag_j,2) + pow(btag_k,2) + pow(btag_l,2) );
            
            // build the lorentz-vectors
            TLorentzVector Wh = selectedJets[i]+selectedJets[j];
            TLorentzVector Bh = selectedJets[k];
            TLorentzVector Th = Wh+Bh;
            TLorentzVector Bl = selectedJets[l];
            
            float AngleThWh = Th.Angle(Wh.Vect());
            float AngleBlMu = Bl.Angle(Mu.Vect());
            float AngleThBl = Th.Angle(Bl.Vect());
            float AngleThBh = Th.Angle(Bh.Vect());
            float AngleThMu = Th.Angle(Mu.Vect());
              
            // pt(Th)/SumPt 3j - possible 3 jet combinations are 012 013 023 and 123
            float sumPt = (selectedJets[0]+selectedJets[1]+selectedJets[2]).Pt();
            sumPt += (selectedJets[0]+selectedJets[1]+selectedJets[3]).Pt();
            sumPt += (selectedJets[0]+selectedJets[2]+selectedJets[3]).Pt();
            sumPt += (selectedJets[1]+selectedJets[2]+selectedJets[3]).Pt();
              
            float ThPtOverSumPt = Th.Pt()/sumPt;
            
            if( dataSetName.find("TTbarJets_SemiMu") == 0 )
            {
              if ( (i == monster->hadrLJet1() || i == monster->hadrLJet2()) && (j == monster->hadrLJet1() || j == monster->hadrLJet2()) && k == monster->hadrBJet() )
              {
  	            // Do some things only for the good jet combination on the hadronic side
                if( l == monster->leptBJet() ) // fully matched event
                {
                  histo2D_["hadronicWMassVSbtag"]->Fill(Wh.M(), btag);
                  histo2D_["hadronicWMassVSAngleThWh"]->Fill(Wh.M(), AngleThWh);
                  histo2D_["hadronicWMassVSAngleBlMu"]->Fill(Wh.M(), AngleBlMu);
                  histo2D_["hadronicWMassVSAngleThBl"]->Fill(Wh.M(), AngleThBl);
                  histo2D_["hadronicWMassVSAngleThBh"]->Fill(Wh.M(), AngleThBh);
                  histo2D_["hadronicWMassVSAngleThMu"]->Fill(Wh.M(), AngleThMu);
                  histo2D_["hadronicWMassVSThPtOverSumPt"]->Fill(Wh.M(), ThPtOverSumPt);
                      
                  histo2D_["hadronicTopMassVSbtag"]->Fill(Th.M(), btag);
                  histo2D_["hadronicTopMassVSAngleThWh"]->Fill(Th.M(), AngleThWh);
                  histo2D_["hadronicTopMassVSAngleBlMu"]->Fill(Th.M(), AngleBlMu);
                  histo2D_["hadronicTopMassVSAngleThBl"]->Fill(Th.M(), AngleThBl);
                  histo2D_["hadronicTopMassVSAngleThBh"]->Fill(Th.M(), AngleThBh);
                  histo2D_["hadronicTopMassVSAngleThMu"]->Fill(Th.M(), AngleThMu);
                  histo2D_["hadronicTopMassVSThPtOverSumPt"]->Fill(Th.M(), ThPtOverSumPt);
         		    }
         		  }
         		  
              if( trainMVA_ && monster->all4JetsMCMatched() ) // Only train for events where you are sure what is S and B
              {
                // additional cuts on jet combinations to use for training
                int oldCombiIndex = -1;
                for(unsigned int jCombi=0; jCombi<12; jCombi++)
                {
         	        unsigned int *oldMVAres = monster->mvaResult(jCombi);
         	        if( ( oldMVAres[0] == (unsigned int)i || oldMVAres[0] == (unsigned int)j ) && 
         	          ( oldMVAres[1] == (unsigned int)i || oldMVAres[1] == (unsigned int)j ) && 
         	          oldMVAres[2] == (unsigned int)k && oldMVAres[3] == (unsigned int)l )
                    oldCombiIndex = jCombi;
                }
                if(oldCombiIndex < 0 || oldCombiIndex > 11) cout << "LightJetCombiner:  oldCombiIndex = " << oldCombiIndex << endl;
          
                if( maxProb(monster, oldCombiIndex) > 0.98 && probNoCorr(monster, dummyMonster_, measureTopMass_, oldCombiIndex) > 0.02 )
                {
                  // fill the MVA trainer variables and check what should be the good combination
                  if ( (i == monster->hadrLJet1() || i == monster->hadrLJet2()) && (j == monster->hadrLJet1() || j == monster->hadrLJet2()) && k == monster->hadrBJet() && l == monster->leptBJet() )
                  {
    	              trainer_->Fill("S","btag",btag);
//    	              trainer_->Fill("S","ThPtOverSumPt",ThPtOverSumPt);
//     	              trainer_->Fill("S","AngleThWh",AngleThWh);
     	              trainer_->Fill("S","AngleBlMu",AngleBlMu);
//  	                trainer_->Fill("S","AngleThBl",AngleThBl);
//     	              trainer_->Fill("S","AngleThBh",AngleThBh);
//    	              trainer_->Fill("S","AngleThMu",AngleThMu);
    					    }
 	  	            else
 	  	            {
 	  	              trainer_->Fill("B","btag",btag);
//   		              trainer_->Fill("B","ThPtOverSumPt",ThPtOverSumPt);
//   		              trainer_->Fill("B","AngleThWh",AngleThWh);
 	  	              trainer_->Fill("B","AngleBlMu",AngleBlMu);
// 	  	              trainer_->Fill("B","AngleThBl",AngleThBl);
//   		              trainer_->Fill("B","AngleThBh",AngleThBh);
// 	  	              trainer_->Fill("B","AngleThMu",AngleThMu);
 	  	            }
 	  	          }
              }
            }
            
            if( !trainMVA_ )
            {
              // compute MVA stuff
              MSPlot_["MVA_input_btag"]->Fill(btag, dataSet, true, Luminosity_*scaleFactor);
              MSPlot_["MVA_input_ThPtOverSumPt"]->Fill(ThPtOverSumPt, dataSet, true, Luminosity_*scaleFactor);
              MSPlot_["MVA_input_AngleThWh"]->Fill(AngleThWh, dataSet, true, Luminosity_*scaleFactor);
              MSPlot_["MVA_input_AngleBlMu"]->Fill(AngleBlMu, dataSet, true, Luminosity_*scaleFactor);
              MSPlot_["MVA_input_AngleThBl"]->Fill(AngleThBl, dataSet, true, Luminosity_*scaleFactor);
              MSPlot_["MVA_input_AngleThBh"]->Fill(AngleThBh, dataSet, true, Luminosity_*scaleFactor);
              MSPlot_["MVA_input_AngleThMu"]->Fill(AngleThMu, dataSet, true, Luminosity_*scaleFactor);
              
              computer_->FillVar("btag",btag);
//              computer_->FillVar("ThPtOverSumPt",ThPtOverSumPt);
//              computer_->FillVar("AngleThWh",AngleThWh);
              computer_->FillVar("AngleBlMu",AngleBlMu);
//              computer_->FillVar("AngleThBl",AngleThBl);
//              computer_->FillVar("AngleThBh",AngleThBh);
//              computer_->FillVar("AngleThMu",AngleThMu);
              
              // get map of  MVAValues  
              std::map<std::string,Float_t> MVAVals = computer_->GetMVAValues();
              
              // vectorMVA_
              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
              {
                MVAValues tempMVA;
                tempMVA.MVAAlgo = it->first;
                tempMVA.MVAValue = it->second;
                tempMVA.WJet1 = i;
                tempMVA.WJet2 = j;
                tempMVA.HadrBJet = k;
                tempMVA.LeptBJet = l;
                
                vectorMVA_.push_back(tempMVA);
                
                if(tempMVA.MVAValue > previousMaxMVA)
                {
                  previousMaxMVA = tempMVA.MVAValue;
  
                  maxMVA_btag = btag;
                  maxMVA_AngleThWh = AngleThWh;
                  maxMVA_AngleBlMu = AngleBlMu;
                  maxMVA_AngleThBl = AngleThBl;
                  maxMVA_AngleThBh = AngleThBh;
                  maxMVA_AngleThMu = AngleThMu;
                  maxMVA_ThPtOverSumPt = ThPtOverSumPt;
                }
              }
            }
          }
        }
      }
    }
  } // end of jetcomb loop
  
  if( !trainMVA_ )
  {
    vector<string> MVAnames = computer_->GetAllMethods();
    for(unsigned int i=0; i<MVAnames.size(); i++)
    {
      pair<float, vector<unsigned int> > MVAResult = getMVAValue(MVAnames[i], 1); // 1 means the highest MVA value
      if(MVAResult.first != previousMaxMVA)
        cout << "MVAResult.first != previousMaxMVA ???" << endl;
      
      //only the jets from the hadronic top should be matched
      bool hadronictopJetsMatched_MVAdef = false;
      if ( (( (int)MVAResult.second[0] == monster->hadrLJet1() && (int)MVAResult.second[1] == monster->hadrLJet2()) 
        || ( (int)MVAResult.second[0] == monster->hadrLJet2() && (int)MVAResult.second[1] == monster->hadrLJet1())) 
        && (int)MVAResult.second[2] == monster->hadrBJet() )
        hadronictopJetsMatched_MVAdef = true;
      
      float mW_maxMVA = (selectedJets[MVAResult.second[0]] + selectedJets[MVAResult.second[1]]).M();
      float mTop_maxMVA = (selectedJets[MVAResult.second[0]] + selectedJets[MVAResult.second[1]] + selectedJets[MVAResult.second[2]]).M();
      MSPlot_["maxMVA_btag"]->Fill(maxMVA_btag, dataSet, true, Luminosity_*scaleFactor);
      MSPlot_["maxMVA_ThPtOverSumPt"]->Fill(maxMVA_ThPtOverSumPt, dataSet, true, Luminosity_*scaleFactor);
      MSPlot_["maxMVA_AngleThWh"]->Fill(maxMVA_AngleThWh, dataSet, true, Luminosity_*scaleFactor);
      MSPlot_["maxMVA_AngleBlMu"]->Fill(maxMVA_AngleBlMu, dataSet, true, Luminosity_*scaleFactor);
      MSPlot_["maxMVA_AngleThBl"]->Fill(maxMVA_AngleThBl, dataSet, true, Luminosity_*scaleFactor);
      MSPlot_["maxMVA_AngleThBh"]->Fill(maxMVA_AngleThBh, dataSet, true, Luminosity_*scaleFactor);
      MSPlot_["maxMVA_AngleThMu"]->Fill(maxMVA_AngleThMu, dataSet, true, Luminosity_*scaleFactor);
      MSPlot_["maxMVA_HadrWMass"]->Fill(mW_maxMVA, dataSet, true, Luminosity_*scaleFactor);
      MSPlot_["maxMVA_HadrTopMass"]->Fill(mTop_maxMVA, dataSet, true, Luminosity_*scaleFactor);
      
      string title = "MaxMVA_"+MVAnames[i]+"_allcomb";
      string titleSemiMuBG = "MaxMVA_"+MVAnames[i]+"_allSemiMuBGcomb"; //alSemiMulBGcomb contains badcomb events and also events where the mcParticle matching was not successful
      string axisTitle = "MVA "+MVAnames[i]+" response";
      
      // check if plots exist, else create them
      if (MSPlot_.find(title) == MSPlot_.end())
        MSPlot_[title] = new MultiSamplePlot(datasets_, title.c_str(), 20, 0, 1, axisTitle.c_str());
      if(histo1D_.find(titleSemiMuBG) == histo1D_.end())
        histo1D_[titleSemiMuBG] = new TH1F(titleSemiMuBG.c_str(),axisTitle.c_str(),50,0,1);
      
      MSPlot_[title]->Fill(MVAResult.first, dataSet, true, Luminosity_*scaleFactor);
      if( hadronictopJetsMatched_MVAdef == false && dataSetName.find("TTbarJets_SemiMu") == 0 )
        histo1D_[titleSemiMuBG]->Fill(MVAResult.first);
      
      if( dataSetName.find("TTbarJets_SemiMu") == 0 )
      {
        std::string titleGood = "MaxMVA_"+MVAnames[i]+"_goodcomb";
   	    std::string titleBad = "MaxMVA_"+MVAnames[i]+"_badcomb";
   	    std::string axisTitle = "MVA "+MVAnames[i]+" response";
        std::string titleGoodmWVSmtop = "MaxMVA_mWVSmtop_"+MVAnames[i]+"_goodcomb";
   	    std::string titleBadmWVSmtop = "MaxMVA_mWVSmtop_"+MVAnames[i]+"_badcomb";
      	
 	      // check if goodcomb badcomb histo exists, else create it
   	    if (histo1D_.find(titleGood) == histo1D_.end())
   	    {
   	      histo1D_[titleGood] = new TH1F(titleGood.c_str(),axisTitle.c_str(),50,0,1);
   	      histo1D_[titleBad] = new TH1F(titleBad.c_str(),axisTitle.c_str(),50,0,1);
   	      histo2D_[titleGoodmWVSmtop] = new TH2F(titleGoodmWVSmtop.c_str(),axisTitle.c_str(),100,0,200,200,0,400);
   	      histo2D_[titleBadmWVSmtop] = new TH2F(titleBadmWVSmtop.c_str(),axisTitle.c_str(),100,0,200,200,0,400);
   	    }
        
   	    //cout << MVAnames[i] << " " << MVAResult->first << " " << hadronictopJetsMatched_MVAdef << endl;
   	    if ( hadronictopJetsMatched_MVAdef )
   	    {
		      histo1D_[titleGood]->Fill( MVAResult.first );
			    histo2D_[titleGoodmWVSmtop]->Fill(mW_maxMVA,mTop_maxMVA);
   	    }
		    else if ( monster->allHadronicJetsMCMatched() )
		    {
   	      histo1D_[titleBad]->Fill( MVAResult.first );
			    histo2D_[titleBadmWVSmtop]->Fill(mW_maxMVA,mTop_maxMVA);
      	}
	    }
	    
	    // put it into the monster!
	    if( MVAnames[i] == "Likelihood" )
	    {
//	      cout << " --> Putting the stuff into the monster!" << endl;
//	      cout << "  Old kinfit example: " << monster->fitResult(0, 21, 21) <<  endl;
	      float weight = monster->eventWeight();
	      vector<TH2F> oldMonsterHistos, newMonsterHistos;
	      for(unsigned int iCombi=0; iCombi<12; iCombi++)
//          oldMonsterHistos.push_back( TH2F() );
	        oldMonsterHistos.push_back( *histoMonster(monster, dummyMonster_, iCombi) );
	      
//	      cout << "  Old kinfit example again: " << oldMonsterHistos[0].GetBinContent(22,22) <<  endl;
	      
   	    vector< float > mvaValsVector;
   	    vector< vector<unsigned int> > mvaResultsVector;
   	    for(unsigned int iCombi=0; iCombi<12; iCombi++)
   	    {
   	      pair<float, vector<unsigned int> > tmpMvaVals = getMVAValue("Likelihood", iCombi+1);
   	      mvaResultsVector.push_back(tmpMvaVals.second);
   	      mvaValsVector.push_back(tmpMvaVals.first);
   	      
   	      // check which old combination this was
   	      int nCombisFound=0;
   	      for(unsigned int jCombi=0; jCombi<12; jCombi++)
   	      {
   	        unsigned int *oldMVAres = monster->mvaResult(jCombi);
   	        if( ( oldMVAres[0] == tmpMvaVals.second[0] || oldMVAres[0] == tmpMvaVals.second[1] ) &&
   	          ( oldMVAres[1] == tmpMvaVals.second[0] || oldMVAres[1] == tmpMvaVals.second[1] ) &&
   	          oldMVAres[2] == tmpMvaVals.second[2] && oldMVAres[3] == tmpMvaVals.second[3] )
   	        {
//   	          cout << " old " << jCombi << " " << monster->mvaVal(jCombi) << " : " << oldMVAres[0] << " " << oldMVAres[1] << " " << oldMVAres[2] << " " << oldMVAres[3] << " | new "
//   	            << iCombi << " " << tmpMvaVals.first << " : " << tmpMvaVals.second[0] << " " << tmpMvaVals.second[1] << " " << tmpMvaVals.second[2] << " " << tmpMvaVals.second[3] << endl;
   	          nCombisFound++;
   	          newMonsterHistos.push_back( oldMonsterHistos[jCombi] );
//   	          if(iCombi == 0)
//   	            cout << "  Old kinfit example again2: " << oldMonsterHistos[0].GetBinContent(22,22) << "  Changing to new value: " << oldMonsterHistos[jCombi].GetBinContent(22,22) << endl;
   	        }
   	      }
   	      if(nCombisFound != 1)
   	      {
   	        cout << "LightJetCombiner:  nCombisFound = " << nCombisFound << endl;
   	        cout << monster->eventID() << " " << monster->runID() << " " << monster->lumiBlockID() << endl;
   	      }
        }
//        cout << "size: " << mvaResultsVector.size() << endl;
        monster->setMvaVals(mvaValsVector);
        monster->setMvaResults(mvaResultsVector);
        monster->setFitResults(newMonsterHistos, measureTopMass_);
        monster->setEventWeight(weight);
//        cout << "  New kinfit example: " << monster->fitResult(0, 21, 21) <<  endl;
	    }
    }
  }
}

pair<float, vector<unsigned int> > LightJetCombiner::getMVAValue(string MVAMethod, int rank)
{
  if(vectorMVA_.size() == 0)
    cout << "LightJetCombiner:  vectorMVA_.size() == 0   Was this event processed by the LightJetCombiner??" << endl;
  if(rank < 1 && rank > 12)
    cout << "LightJetCombiner: Wrong rank!!   rank = " << rank << endl;
    
  vector<MVAValues> MVAvalues;
  
  for(unsigned int i=0; i<vectorMVA_.size(); i++)
    if( vectorMVA_[i].MVAAlgo.find(MVAMethod) == 0 )
      MVAvalues.push_back(vectorMVA_[i]);
  
  sort( MVAvalues.begin(), MVAvalues.end(), sortMVAValuesAgain );
  
  pair<float, vector<unsigned int> > out;
  vector<unsigned int> tempVec;
  tempVec.push_back(MVAvalues[rank-1].WJet1);
  tempVec.push_back(MVAvalues[rank-1].WJet2);
  tempVec.push_back(MVAvalues[rank-1].HadrBJet);
  tempVec.push_back(MVAvalues[rank-1].LeptBJet);
  out.second = tempVec;
  out.first = MVAvalues[rank-1].MVAValue;
  return out;
}

void LightJetCombiner::Write(TFile* fout, bool savePNG, string pathPNG)
{
  if(trainMVA_)
  {
    cout << " - Training the MVA!" << endl;
    trainer_->TrainMVA("Block","",0,0,"",0,0,postfix_);
  }
  
  cout << " - Writing out the LightJetCombiner stuff ..." << endl;
  
//  string pathPNGExpCorr = pathPNG + "../ExpCorr/";
//	if( ! trainMVA_ ) expCorrIncl_->Write(fout, true, pathPNGExpCorr);

  fout->cd();
  string dirname = "1D_histograms_JetCombiner";
	if(fout->Get(dirname.c_str())==0)
		fout->mkdir(dirname.c_str());
	
	fout->cd(dirname.c_str());
  // Write MVA distributions
  if( !trainMVA_ )
  {
    map<string,Float_t> MVAVals = computer_->GetMVAValues();
    vector< pair<TGraph*, string> > effpurplots;
    vector< pair<TGraph*, string> > effpurplotsAllSemiMuBG;
    
    for (map<string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
    {
      vector< pair<TH1F*,string> > histos;

      string titleGood = "MaxMVA_"+it->first+"_goodcomb";
      string titleBad = "MaxMVA_"+it->first+"_badcomb";
      string titleSemiMuBG = "MaxMVA_"+it->first+"_allSemiMuBGcomb";
      
      if (histo1D_.find(titleGood) != histo1D_.end())
      {
	      string label = "Maximal MVA value using method: "+it->first;
	
	      /////////////
	      // S/(S+B) //
	      /////////////
	      TH1F* newB = (TH1F*)histo1D_[titleBad]->Clone();
	      newB->Add(histo1D_[titleGood]);
	
	      TGraphAsymmErrors* SoverSB = new TGraphAsymmErrors(histo1D_[titleGood],newB);
	
	      SoverSB->SetTitle((it->first+" response: S/(S+B)").c_str());
	      SoverSB->GetXaxis()->SetTitle((it->first+" response").c_str());
	      SoverSB->GetYaxis()->SetTitle("S/(S+B)");
	
	      TCanvas* tempCanvas2 = TCanvasCreator(SoverSB, "c1");
	      tempCanvas2->SaveAs( (pathPNG+"SoverSB_"+it->first+".png").c_str() );
	
	      tempCanvas2->Write();
	
	      ///////////////////////////////////
	      // make efficiency vs purity plot
        ///////////////////////////////////
	      double x[100000], yBadComb[100000], yAllSemiMuBGComb[100000];
	      int nVals = 0;

	      // calculate efficiency and purity
	      for (double cut=0; cut < 1; cut+=0.01)
	      {
	        Int_t cutBin = histo1D_[titleGood]->GetXaxis()->FindBin(cut);

	        Double_t nSignalAbove = histo1D_[titleGood]->Integral(cutBin,histo1D_[titleGood]->GetNbinsX());
	        Double_t nSignalBelow = histo1D_[titleGood]->Integral(0,cutBin);
	        Double_t nBadCombAbove = histo1D_[titleBad]->Integral(cutBin,histo1D_[titleBad]->GetNbinsX());
	        Double_t nAllSemiMuBGCombAbove = histo1D_[titleSemiMuBG]->Integral(cutBin,histo1D_[titleSemiMuBG]->GetNbinsX());
	        
	        x[nVals] = (Double_t) nSignalAbove / (Double_t) (nSignalAbove + nSignalBelow); //efficiency
	        yBadComb[nVals] = (Double_t) nSignalAbove / (Double_t) (nSignalAbove + nBadCombAbove); //purity
	        yAllSemiMuBGComb[nVals] = (Double_t) nSignalAbove / (Double_t) (nSignalAbove + nAllSemiMuBGCombAbove); //purity
	        
	        nVals++;
	        //cout <<"Cut " << cut << " -> " <<  cutBin << " " << nSignalBelow << " " << nSignalAbove << " " << efficiency << " " << purity << endl;
        }
	      TGraph* plot = new TGraph(nVals,x,yBadComb);
	      effpurplots.push_back(pair<TGraph*,string>(plot,it->first));
	      TGraph* plotallSemiMuBG = new TGraph(nVals,x,yAllSemiMuBGComb);
	      effpurplotsAllSemiMuBG.push_back(pair<TGraph*,string>(plotallSemiMuBG,it->first));
	      
	      cout << "   - MVA method: " << it->first << " ||  Eff for SemiMu-Matched: " << yBadComb[0] << " ||  Eff for All SemiMu: " << yAllSemiMuBGComb[0] << endl;
	      cout << "                 nGoodCombi: " << histo1D_[titleGood]->GetEffectiveEntries() << " || nBadCombi: " << histo1D_[titleBad]->GetEffectiveEntries() << endl;
	
	      //////////////////////////////////////
	      // good overlayed with bad jet combi
	      //////////////////////////////////////
	      histo1D_[titleGood]->GetXaxis()->SetTitle((it->first+" Response").c_str());
	      histo1D_[titleBad]->GetXaxis()->SetTitle((it->first+" Response").c_str());
	
	      histo1D_[titleGood]->GetYaxis()->SetTitle("# events");
	      histo1D_[titleBad]->GetYaxis()->SetTitle("# events");
	
	      histos.push_back(std::pair<TH1F*,std::string>(histo1D_[titleGood],"Good Jet Combinations"));
	      histos.push_back(std::pair<TH1F*,std::string>(histo1D_[titleBad],"Bad Jet Combinations"));
	
	      TCanvas* tempCanvas = TCanvasCreator(histos, "LF",label,false);
	      tempCanvas->SaveAs( (pathPNG+"MaxMVA_"+it->first+".png").c_str() );
	
	      tempCanvas->Write();
	
	      delete tempCanvas;
	      delete tempCanvas2;
	
	      histos.clear();
      }
    }
    
    ///////////////////////////////
    // WRITE EFF PUR plots canvas
    ///////////////////////////////
    string effpurtitle= "Eff_Pur_MVA";

    if (effpurplots.size() > 0)
    {
      effpurplots[0].first->GetXaxis()->SetTitle("Efficiency");
      effpurplots[0].first->GetYaxis()->SetTitle("Purity");
      effpurplots[0].first->SetTitle("Efficiency VS purity MVA discriminator cut, only successfull mc-matched hadronic tops");
      effpurplots[0].first->SetLineWidth(3);
    }
    TCanvas* tempCanvas3 = TCanvasCreator(effpurplots, effpurtitle.c_str());
    tempCanvas3->SaveAs( (pathPNG+effpurtitle+".png").c_str() );
    tempCanvas3->Write();
    
    effpurtitle= "Eff_Pur_MVA_AllSemiMuBG";
    if (effpurplotsAllSemiMuBG.size() > 0)
    {
      effpurplotsAllSemiMuBG[0].first->GetXaxis()->SetTitle("Efficiency");
      effpurplotsAllSemiMuBG[0].first->GetYaxis()->SetTitle("Purity");
      effpurplotsAllSemiMuBG[0].first->SetTitle("Efficiency VS purity MVA discriminator cut, all selected TTbar semi-mu events");
      effpurplotsAllSemiMuBG[0].first->SetLineWidth(3);
    }
    TCanvas* tempCanvas5 = TCanvasCreator(effpurplotsAllSemiMuBG, effpurtitle.c_str());
    tempCanvas5->SaveAs( (pathPNG+effpurtitle+".png").c_str() );
    tempCanvas5->Write();
  } // done with MVA plotting

  mkdir( (pathPNG+"MSPlot/").c_str(),0777);
  
  // Write out the MSPlot
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot_.begin(); it != MSPlot_.end(); it++)
  {
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(false, name, true, true, true, true, true);
    temp->Write(fout, name, true, pathPNG+"MSPlot/");
  }
  
  for(std::map<std::string,TH1F*>::const_iterator it = histo1D_.begin(); it != histo1D_.end(); it++)
	{
		TH1F *temp = it->second;
		int N = temp->GetNbinsX();
  	temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
  	temp->SetBinContent(N+1,0);
		temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}

  // 2D
  string dirname2D = "2D_histograms_graphs_JetCombiner";
	if(fout->Get(dirname.c_str())==0)
		fout->mkdir(dirname.c_str());
	fout->cd(dirname.c_str());
	
	for(std::map<std::string,TH2F*>::const_iterator it = histo2D_.begin(); it != histo2D_.end(); it++)
	{
		TH2F *temp = it->second;
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}
	
	//Write TGraphAsymmErrors
	for(map<string,TGraphAsymmErrors*>::const_iterator it = graphAsymmErr_.begin(); it != graphAsymmErr_.end(); it++)
	{
	  TGraphAsymmErrors *temp = it->second;
	  temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}
}

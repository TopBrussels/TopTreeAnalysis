#include "TopTreeAnalysis/JESMeasurement/IdeogramTools/interface/MassFitInfoFiller.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

void fillMassFitInfoHeader(Dataset* dataSet, ideogram::MassFitInfoHeader& header, bool isSemiMu, bool containsPDFsyst)
{
  header.channel_ = isSemiMu;
  string dataSetName = dataSet->Name();
  if( dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 || dataSetName.find("DAtaDriven") != string::npos )
    header.isMC_ = false;
  else header.isMC_ = true;
  header.datasetName_ = dataSetName;
  header.datasetTitle_ = dataSet->Title();
  header.datasetCrossSection_ = dataSet->Xsection();
  header.datasetIntegratedLuminosity_ = dataSet->EquivalentLumi();
  header.pdf_ = containsPDFsyst;
}

void fillMassFitInfoEvent(LightMonster* monster, const ideogram::MassFitInfoHeader& header, ideogram::MassFitInfoEvent& event, float PUweight, float PUup, float PUdown)
{
  event.runNumber = monster->runID();
  event.luminositySection = monster->lumiBlockID();
  event.eventNumber = monster->eventID();
  vector<TLorentzVector> selJets = monster->selectedJets();
  vector<float> bTag = monster->bTagSSVHE();
  event.nJet = selJets.size();
  if(header.channel_ == false)
  {
    event.nElectron = 1;
    event.nMuon = 0;
  }
  else
  {
    event.nElectron = 0;
    event.nMuon = 1;
  }
  event.nVertex = monster->nPV();
  event.nPU = monster->nPU();
  event.leptonCharge = monster->leptonCharge();
  TLorentzVector lepton = monster->lepton();
  event.leptonPt = lepton.Pt();
  event.leptonEta = lepton.Eta();
  event.leptonPhi = lepton.Phi();
  TLorentzVector met =  monster->met();
  event.METx = met.Px();
  event.METy = met.Py();
  for (size_t j = 0 ; j != ideogram::MassFitInfoEvent::kNJET ; ++j)
  {
    event.jetPt[j] = selJets[j].Pt();
    event.jetEta[j] = selJets[j].Eta();
    event.jetPhi[j] = selJets[j].Phi();
    event.jetBTagDiscriminant[j] = bTag[j];
  }
  event.eventWeight = monster->eventWeight();
  event.lumiWeight = PUweight;
  event.PUdown = PUdown;
  event.PUup = PUup;
}

void fillMassFitResult(LightMonster* monster, ideogram::MassFitInfoEvent& event, float chi2cut)
{
  size_t nFit(0);

  // Calculate total weight, for normalization
  double sumChi2Weight(0.0);
  double sumBTagWeight(0.0);
  double sumChi2BTagWeight(0.0);

  for(unsigned int iCombi=0; iCombi<12; iCombi++)
  {
    unsigned int* combi = monster->mvaResult(iCombi);
    if( monster->chi2MTopFit(iCombi) < 10 )
    {
      // fill the stuff!
      event.fittedTopMass[nFit] = monster->mTopFit(iCombi);
      event.fittedTopMassSigma[nFit] = monster->sigmaMTopFit(iCombi);
      event.fitChi2[nFit] = monster->chi2MTopFit(iCombi);
      for (size_t j = 0 ; j!= 4 ; ++j) event.fitJetType[nFit*(ideogram::MassFitInfoEvent::kNJET)+j] = combi[j];
      event.fitChi2Weight[nFit] = exp(-0.5*event.fitChi2[nFit]);
      event.fitBTagWeight[nFit] = 1.0; // No weight for the moment
      event.fitWeight[nFit] = event.fitChi2Weight[nFit]; // No account is taken for BTagWeight for now
      event.fitConverge[nFit] = true;
      sumChi2Weight += event.fitChi2Weight[nFit];
      sumBTagWeight += event.fitBTagWeight[nFit];
      sumChi2BTagWeight += event.fitWeight[nFit];
      nFit++;
    }
  }

  event.nFit = nFit;
  event.nFitConverge = nFit; // Implement ad-hoc calculation of fit convergence criteria here
  
  // normalize the weight
  for (size_t f = 0 ; f != nFit ; ++f)
  {
      event.fitChi2Weight[f] = event.fitChi2Weight[f]/sumChi2Weight;
      event.fitBTagWeight[f] = event.fitBTagWeight[f]/sumBTagWeight;
      event.fitWeight[f] = event.fitWeight[f]/sumChi2BTagWeight;
  }

  // fill the sorted chi-squared index
  std::vector< std::pair<size_t, double> > Chi2Index;
  for (size_t f = 0 ; f != nFit ; ++f) Chi2Index.push_back(std::pair<size_t,double>(f,event.fitChi2[f]));
  std::sort(Chi2Index.begin(),Chi2Index.end(),ideogram::IndexedQuantityAbsLessThan<double>);
  for (size_t f = 0 ; f != nFit ; ++f) event.fitSortedChi2Index[f] = Chi2Index[f].first;
}

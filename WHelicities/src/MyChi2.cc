 #include "../interface/MyChi2.h" 

using namespace Annik;

void MyChi2::SetNTuple(TNtuple *genttbarhisto[CosThetaNbins])
{
  genttbarhisto_ = genttbarhisto;
}

double MyChi2::DoEval(const double* x) const {
    double f0Fit = x[0];
    double flFit = x[1];
    double frFit = 0;
    double nn=1.;
    if (ndimen_==3) { 
      frFit = 1-x[0]-x[1]; 
      nn = x[2];
    }
    else if (ndimen_==2){
      flFit = 1.-x[0];
      nn = x[1];
    }

    double ChiSquaredAllBins =0.;
    double chi2 = 0.;
    
    int nncells = signalh_->GetSize()-2;  // -1 and not -2, to get the last bin

    // initialization
    double nmcaux[nncells];
    for (int ibin=0; ibin< nncells; ibin++){
      nmcaux[ibin]=0.;
    }
  
    // loop over ntuples containing gen-level costh info relative to each costhRec bin
    Int_t nGENentries[nncells];
    for (int ihgen=0; ihgen<nncells; ihgen++){ 
 
      // Here I am assuming that the input array only contains entries that 
      // have costheta* reconstructed
      // If no costhetarec for this igen, one should skip this entry in the sum
      // or send it to an underflow/overflow bin
      
      Float_t costhgen, evtweight;
      nGENentries[ihgen] = genttbarhisto_[ihgen]->GetEntries();
      
      genttbarhisto_[ihgen]->SetBranchAddress("costhgen",&costhgen);
      genttbarhisto_[ihgen]->SetBranchAddress("evtweight",&evtweight);
    
      if (nGENentries[ihgen] ==0) {
	nmcaux[ihgen] =0;
	//unrwgt[ihgen]=0; //Double which has the function of summing up the number of generated events after normalization, but before reweighting in the fitter.
	// Double is necessary because cos theta * only needs to be reweighted for generated Muon leptonic decays
	// --> Not needed in my code since we can split in ttbarSemiMu ??
      }
      else{
        for (Int_t igen=0; igen<nGENentries[ihgen]; ++igen) {
          genttbarhisto_[ihgen]->GetEntry(igen);
          double xx = costhgen;
          double thisevweight = evtweight;

          double SM = 0.3325*(1-xx)*(1-xx)*3/8+0.6671*(1-xx*xx)*3/4 + 0.0004*(1+xx)*(1+xx)*3/8 ;
          double newmodel = flFit*(1-xx)*(1-xx)*3/8+f0Fit*(1-xx*xx)*3/4 + frFit*(1+xx)*(1+xx)*3/8 ;
          
          nmcaux[ihgen] +=thisevweight*newmodel/SM;
        }
      }	  
    }
    for (int ibin=0; ibin< nncells; ibin++){
      double ndata_i = datah_->GetBinContent(ibin+1);       
      double nbkg_i = bkgh_->GetBinContent(ibin+1);     
      double  nmc_i =  nbkg_i +  nn*nmcaux[ibin] ;   // nn--> free normalization found by the Fit
      
      if (nmc_i>0){
	chi2 += ( nmc_i - ndata_i * TMath::Log(nmc_i) ); 
	ChiSquaredAllBins=ChiSquaredAllBins+(((nmc_i-ndata_i)/(sqrt(ndata_i)))*((nmc_i-ndata_i)/(sqrt(ndata_i))));//nmc_i = all samples ndata_ = semiMu	
      }      
    }    
    return ChiSquaredAllBins;
}

#ifndef VJetEstimation_h
#define VJetEstimation_h

#include <iostream>
#include <iomanip>
#include <vector>

#include "Math/VectorUtil.h"
#include "Math/Polynomial.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TList.h"
#include "TObject.h"


#include "TopTreeAnalysis/Content/interface/MCObsExpectation.h"
#include "TopTreeProducer/interface/TRootJet.h"

  // RooFit librairies
#include "RooArgSet.h"
#include "RooAddition.h"
#include "RooCategory.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooHist.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
  //#include "RooSimultaneous.h"

#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
  //#include "RooNLLVar.h"
  //#include "RooProfileLL.h"
#include "RooPlot.h"

#include "RooFitResult.h"

using namespace std;
using namespace TopTree;
using namespace RooFit ;

class VJetEstimation : public TObject
{
  
public:
  VJetEstimation();
  /** Constructor assigning NULL to pointers. */
	VJetEstimation(UInt_t NofBtagWorkingPoInt_t, Float_t* BtagWorkingPoInt_t, UInt_t NofJets, UInt_t NofJetBins, UInt_t NofDatasets, vector<Int_t> iDTTLike, vector<Int_t> iDVLike, vector<Int_t> iDVbLike);
	/** Constructor */
	VJetEstimation(UInt_t NofBtagWorkingPoInt_t, Float_t* BtagWorkingPoInt_t, UInt_t NofJets, UInt_t NofJetBins, Double_t** EffXbq, UInt_t NofDatasets, vector<Int_t> iDTTLike, vector<Int_t> iDVLike, vector<Int_t> iDVbLike);
	/** Constructor */
	VJetEstimation(const VJetEstimation& vjet);
 	/** Copy constructor (not yet implemented) */
	~VJetEstimation();
	/** Destructor */
	
    //	void FillInputs(Double_t**** n, Double_t*** MultiJets_Estimated_Nj);
 	/** Method used to fill the histograms containing the number of events with 0,1,2 and 3 b-jets for all datasets. Needs also the estimated number of multijet events */
    //	void FillInputs(Double_t**** n);
 	/** Method used to fill the histograms containing the number of events with 0,1,2 and 3 b-jets for all datasets.*/
	void Fill(vector<TRootJet*> &SelectedJets, UInt_t idx, Int_t btagAlgo, Double_t weight=1.);
 	/** Method to be used for each event to fill the histograms containing the number of events with 0,1,2 and 3 b-jets for the dataset idx. */
	void FillInputs(vector <vector< vector< vector <Double_t> > > > Inputs);
	/** Method to be used to fill the histograms containing the number of events with 0,1,2 and 3 b-jets for all datasets. */
	void FillInputs(vector <vector< vector< vector <Double_t> > > > Inputs, vector <vector< vector<Double_t> > > Multijet_est_nj);
	/** Method to be used to fill the histograms containing the number of events with 0,1,2 and 3 b-jets for all datasets. */
  
	void ReScaleInputs(UInt_t idx, Float_t Ntot, Double_t factor);
	/** Rescale the input number of events from dataset idx. Method to be used after the event loop.*/
	Float_t WilsonScoreIntervalLow(Float_t Non, Float_t Ntot);
	/** Upper limit of the Wilson score Int_terval for binomial parameter (being the selection efficiency here) (ArXiv:hep-ph/0905.3831v2)*/
	Float_t WilsonScoreIntervalHigh(Float_t Non, Float_t Ntot);
	/** Lower limit of the Wilson score Int_terval for binomial parameter (being the selection efficiency here) (ArXiv:hep-ph/0905.3831v2)*/
	Float_t WilsonScoreIntervalMean(Float_t Non, Float_t Ntot);
	/** Mean between the upper and lower limit of the Wilson score Int_terval for binomial parameter (being the selection efficiency here) (ArXiv:hep-ph/0905.3831v2)*/
	void SumOverAllInputs();
	/** Sum the weighted contribution of all the datasets. Method to be used after having looped over all events of all datasets.*/
	void ComputeEffbqFromMC(UInt_t idx);
	/** Compute ebq efficiencies from MC. Method to be used after having looped over all events of all datasets.*/
	void ComputeEffFromMC();
	/** Compute the b/mis-tagging efficiencies from MC. Method to be used after having looped over all events of all datasets.*/
  
    //void ComputeCondProb(bool Verbose = false);
    // /** Compute the conditional probabilities */
	
	void FillSummaryHistos();
	/** Fill summary histograms. */
	void BckgdSubstraction(vector<MCObsExpectation*> &hists, vector<string> &name, Float_t lumi);
	/** Substract background b-jet ditribution provided as input */
    //void Write(TFile* file, string label = string(""), bool verbose = false);
 	/** Method to book histograms */
	
    /////////////////////////////////
    //Main Methods to be executed ...
    /////////////////////////////////
	void SetInitialValues_Nttlike(vector<Double_t> init);
	/** Set initial values of Nttlike for the minizer. */
	void SetInitialValues_Nttlike(UInt_t njets,Double_t init);
	/** Set initial values of Nttlike for the minizer. */
  
	void SetInitialValues_Nvlike(vector<Double_t> init);
	/** Set initial values of Nvlike for the minizer. */
	void SetInitialValues_Nvlike(UInt_t njets,Double_t init);
	/** Set initial values of Nvlike for the minizer. */
  
	void SetInitialValues_Eb(vector< vector<Double_t> > init);
	/** Set initial values of Eb for the minizer. */
	void SetInitialValues_Eb(UInt_t njets,vector<Double_t> init);
	/** Set initial values of Eb for the minizer. */
	void SetInitialValues_Eb(UInt_t njets,UInt_t btagIdx, Double_t init);
	/** Set initial values of Eb for the minizer. */
  
	void SetInitialValues_Eudsc(vector< vector<Double_t> > init);
	/** Set initial values of Eudsc for the minizer. */
	void SetInitialValues_Eudsc(UInt_t njets,vector<Double_t> init);
	/** Set initial values of Eudsc for the minizer. */
	void SetInitialValues_Eudsc(UInt_t njets,UInt_t btagIdx, Double_t init);
	/** Set initial values of Eudsc for the minizer. */
  
	void SetInitialValues_Euds(vector< vector<Double_t> > init);
	/** Set initial values of Euds for the minizer. */
	void SetInitialValues_Euds(UInt_t njets,vector<Double_t> init);
	/** Set initial values of Euds for the minizer. */
	void SetInitialValues_Euds(UInt_t njets,UInt_t btagIdx, Double_t init);
	/** Set initial values of Euds for the minizer. */
  
	void UnBinnedMaximumLikelihoodEst(UInt_t btag_wp_idx, UInt_t njets, vector<Int_t> &FixedVarIdx, bool doMinos, bool doToyMC, bool Verbose);
	/** Maximum likelihood estimation with RooFit. */
  
	void UnBinnedMaximumLikelihoodEst(UInt_t njets, vector<Int_t> &FixedVarIdx, bool doMinos, bool doToyMC, bool Verbose){
		for(UInt_t i=0;i<NbOfBtagWorkingPoint_;i++) UnBinnedMaximumLikelihoodEst(i, njets, FixedVarIdx, doMinos, doToyMC, Verbose);
	};
	/** Maximum likelihood estimation with RooFit (solve for all b-tagging working points). */
  
	void UnBinnedMaximumLikelihoodEst(vector<Int_t> &FixedVarIdx, bool doMinos, bool doToyMC, bool Verbose){
		for(UInt_t i=0;i<NbOfJetsBins_;i++) UnBinnedMaximumLikelihoodEst(i, FixedVarIdx, doMinos, doToyMC, Verbose);
	};
	/** Maximum likelihood estimation with RooFit (solve for all jet multiplicities). */
	
	void UnBinnedMaximumJointWPLikelihoodEst(UInt_t njets, vector<Int_t> &FixedVarIdx, bool doMinos, bool doToyMC, bool Verbose);
	/** Maximum likelihood estimation with RooFit, combining all the b-tagging working points. */
  
	void UnBinnedMaximumJointWPLikelihoodEst(vector<Int_t> &FixedVarIdx, bool doMinos, bool doToyMC, bool Verbose){
		for(UInt_t i=0;i<NbOfJetsBins_;i++) UnBinnedMaximumJointWPLikelihoodEst(i, FixedVarIdx, doMinos, doToyMC, Verbose);
	};
	/** Maximum likelihood estimation with RooFit, combining all the b-tagging working points (solve for all jet multiplicities). */
  
    //void UnBinnedMaximumNjetsLikelihoodEst(Int_t btag_wp_idx, vector<Int_t> &FixedVarIdx, bool doMinos, bool Verbose);
    // /** Maximum likelihood estimation with RooFit, combining all the jet multiplicities. */
	
    //void UnBinnedMaximumJoInt_tLikelihoodEst(bool Verbose);
    // /** Maximum likelihood estimation with RooFit, combining both all the jet multiplicities and all b-tagging working points. */
	
    /////////////////////
    // PrInt_t  methods
    ////////////////////
	void PrintInputs();
	/** Method to prInt_t the number of events passed as imputs per dataset*/
	void PrintInputs(UInt_t njets);
	/** Method to prInt_t the number of events passed as imputs per dataset, per jet multiplicity*/
	void PrintResults();
	/** Method to prInt_t the output of the estimation method*/
	void PrintResults(UInt_t njets);
	/** Method to prInt_t the output of the estimation method, per jet multiplicity*/
	void PrintResults_LatexFormat(ofstream &ofile);
	/** Method to prInt_t a latex-formatted output of the estimation method*/
	
    /////////////////////
    // Cross-check methods
    /////////////////////
	bool CheckEstimation(Double_t threshold);
	/** Returns true if the value of the minimizer function at the minimum is below the threshold. */
	bool CheckEstimation(Double_t threshold, UInt_t njets);
	/** Returns true if the value of the minimizer function at the minimum is below the threshold. */
	void CheckEstimationLinearity(vector<Int_t> &FixedVarIdx, bool doMinos, Int_t nbOfRescaleFact, Double_t **RescaleFact, Int_t Idx);
	/** Method to check if the response of the estimation scales linearly with the inputs. */
  
    //////////////////////////////////////
    // Access to the class members
    //////////////////////////////////////
  
	UInt_t GetNbOfBtagWorkingPoInt_t()    const {return NbOfBtagWorkingPoint_;};
	/** returns the number of b-tagging working points used by the estimation method */
	UInt_t GetNbOfDatasets()            const {return NbOfDatasets_;};
	/** returns the number of datasets used by the estimation method */
	UInt_t GetNbOfBJetsBins()           const {return NbOfBJetsBins_;};
	/** returns the number of b-jet bins used by the estimation method */
	UInt_t GetNbOfJetsBins()            const {return NbOfJetsBins_;};
	/** returns the number of jet bins used by the estimation method */
	UInt_t GetNjets(UInt_t jetidx)const {return (jetidx<NbOfJetsBins_ ? Njets_[jetidx] : -999);}
  
	vector<Int_t> GetiDatasetsTTLike() const {return iDatasetsTTLike_;};
	vector<Int_t> GetiDatasetsVLike()  const {return iDatasetsVLike_;};
	vector<Int_t> GetiDatasetsVbLike() const {return iDatasetsVbLike_;};
  
	vector< vector< vector< vector< Double_t > > > > GetN()                                              const{ return N_;};
	vector< vector< vector< Double_t > > >           GetN(UInt_t wp)                                     const{ return N_[wp];};
	vector< vector< Double_t > >                     GetN(UInt_t wp,UInt_t njets)                        const{ return N_[wp][njets];};
	vector< Double_t >                               GetN(UInt_t wp,UInt_t njets,Int_t nbjets)             const{ return N_[wp][njets][nbjets];};
	Double_t                                         GetN(UInt_t wp,UInt_t njets,Int_t nbjets,Int_t dataset) const{ return N_[wp][njets][nbjets][dataset];};
  
	vector< vector< vector< Double_t > > > GetNbjets()                                       const{ return Nbjets_;};
	vector< vector< Double_t > >           GetNbjets(UInt_t wp)                        	 const{ return Nbjets_[wp];};
	vector< Double_t >                     GetNbjets(UInt_t wp, UInt_t njets)             	 const{ return Nbjets_[wp][njets];};
	Double_t                               GetNbjets(UInt_t wp, UInt_t njets, UInt_t nbjets) const{ return Nbjets_[wp][njets][nbjets];};
	
	Double_t   GetMinLL(UInt_t njets) const{ return minValue_[njets];};
  /*
   Double_t*  GetMultiJet_Est(UInt_t wp) const;
   Double_t*  GetMultiJet_Est(UInt_t wp, UInt_t njets) const;
   Double_t*  GetMultiJet_Est(UInt_t wp, UInt_t njets, UInt_t nbjets) const;
   */
	Double_t   GetPredEb(Int_t wp) const;
	Double_t   GetPredEb(Int_t wp, Int_t njets) const;
	Double_t   GetPredEbErr(Int_t wp) const;
	Double_t   GetPredEbErr(Int_t wp, Int_t njets) const;
	void     SetPredEb(Int_t wp, Int_t njets, Double_t value) ;
  
	Double_t   GetPredEudsc(Int_t wp) const;
	Double_t   GetPredEudsc(Int_t wp, Int_t njets) const;
	Double_t   GetPredEudscErr(Int_t wp) const;
	Double_t   GetPredEudscErr(Int_t wp, Int_t njets) const;
	void     SetPredEudsc(Int_t wp, Int_t njets, Double_t value) ;
  
	Double_t   GetPredEuds(Int_t wp) const;
	Double_t   GetPredEuds(Int_t wp, Int_t njets) const;
	Double_t   GetPredEudsErr(Int_t wp) const;
	Double_t   GetPredEudsErr(Int_t wp, Int_t njets) const;
	void     SetPredEuds(Int_t wp, Int_t njets, Double_t value) ;
  
	Double_t   GetPredNv(Int_t wp) const;
	Double_t   GetPredNv(Int_t wp, Int_t njets) const;
	Double_t   GetPredNv(Int_t wp, Int_t njets, Int_t nbjets) const;
	Double_t   GetPredNvErr(Int_t wp) const;
	Double_t   GetPredNvErr(Int_t wp, Int_t njets) const;
	Double_t   GetPredNvErr(Int_t wp, Int_t njets, Int_t nbjets) const;
  
	Double_t   GetPredNvb(Int_t wp) const;
	Double_t   GetPredNvb(Int_t wp, Int_t njets) const;
	Double_t   GetPredNvb(Int_t wp, Int_t njets, Int_t nbjets) const;
	Double_t   GetPredNvbErr(Int_t wp) const;
	Double_t   GetPredNvbErr(Int_t wp, Int_t njets) const;
	Double_t   GetPredNvbErr(Int_t wp, Int_t njets, Int_t nbjets) const;
  
	Double_t   GetPredNtt(Int_t wp) const;
	Double_t   GetPredNtt(Int_t wp, Int_t njets) const;
	Double_t   GetPredNtt(Int_t wp, Int_t njets, Int_t nbjets) const;
	Double_t   GetPredNttErr(Int_t wp) const;
	Double_t   GetPredNttErr(Int_t wp, Int_t njets) const;
	Double_t   GetPredNttErr(Int_t wp, Int_t njets, Int_t nbjets) const;
  
	Double_t   GetPredNtotal(Int_t wp) const;
	Double_t   GetPredNtotal(Int_t wp, Int_t njets) const;
	Double_t   GetPredNtotal(Int_t wp, Int_t njets, Int_t nbjets) const;
    //	Double_t   GetPredNtotalErr(Int_t wp) const;
    //	Double_t   GetPredNtotalErr(Int_t wp, Int_t njets) const;
    //	Double_t   GetPredNtotalErr(Int_t wp, Int_t njets, Int_t nbjets) const;
  
	Double_t   GetPredN(Int_t idx) const;
	Double_t   GetPredN(Int_t idx, Int_t njets) const;
  
    //////////////////////////////////////
    // Access to the parameter estimated values
    //////////////////////////////////////
  
	Double_t GetEstEb(Int_t wp) const;
	Double_t GetEstEb(Int_t wp, Int_t njets) const;
	Double_t GetEstEbErr(Int_t wp) const;
	Double_t GetEstEbErr(Int_t wp, Int_t njets) const;
	Double_t GetEstEbErrHigh(Int_t wp, Int_t njets) const;
	Double_t GetEstEbErrLow(Int_t wp, Int_t njets) const;
  
	Double_t GetEstEudsc(Int_t wp) const;
	Double_t GetEstEudsc(Int_t wp, Int_t njets) const;
	Double_t GetEstEudscErr(Int_t wp) const;
	Double_t GetEstEudscErr(Int_t wp, Int_t njets) const;
	Double_t GetEstEudscErrHigh(Int_t wp, Int_t njets) const;
	Double_t GetEstEudscErrLow(Int_t wp, Int_t njets) const;
  
	Double_t GetEstEuds(Int_t wp) const;
	Double_t GetEstEuds(Int_t wp, Int_t njets) const;
	Double_t GetEstEudsErr(Int_t wp) const;
	Double_t GetEstEudsErr(Int_t wp, Int_t njets) const;
	Double_t GetEstEudsErrHigh(Int_t wp, Int_t njets) const;
	Double_t GetEstEudsErrLow(Int_t wp, Int_t njets) const;
  
	Double_t GetEstNv(Int_t wp) const;
	Double_t GetEstNv(Int_t wp, Int_t njets) const;
	Double_t GetEstNv(Int_t wp, Int_t njets, Int_t nbjets) const;
	Double_t GetEstNvErr(Int_t wp) const;
	Double_t GetEstNvErr(Int_t wp, Int_t njets) const;
	Double_t GetEstNvErrUp(  Int_t wp, Int_t njets) const;
	Double_t GetEstNvErrDown(Int_t wp, Int_t njets) const;
	Double_t GetEstNvErr(Int_t wp, Int_t njets, Int_t nbjets) const;
  
	Double_t GetEstNvb(Int_t wp) const;
	Double_t GetEstNvb(Int_t wp, Int_t njets) const;
	Double_t GetEstNvb(Int_t wp, Int_t njets, Int_t nbjets) const;
	Double_t GetEstNvbErr(Int_t wp) const;
	Double_t GetEstNvbErr(Int_t wp, Int_t njets) const;
	Double_t GetEstNvbErr(Int_t wp, Int_t njets, Int_t nbjets) const;
  
	Double_t GetEstNtt(Int_t wp) const;
	Double_t GetEstNtt(Int_t wp, Int_t njets) const;
	Double_t GetEstNtt(Int_t wp, Int_t njets, Int_t nbjets) const;
	Double_t GetEstNttErr(Int_t wp) const;
	Double_t GetEstNttErr(Int_t wp, Int_t njets) const;
	Double_t GetEstNttErrUp(  Int_t wp, Int_t njets) const;
	Double_t GetEstNttErrDown(Int_t wp, Int_t njets) const;
	Double_t GetEstNttErr(Int_t wp, Int_t njets, Int_t nbjets) const;
  
	Double_t GetEstNtotal(Int_t wp) const;
	Double_t GetEstNtotal(Int_t wp, Int_t njets) const;
	Double_t GetEstNtotal(Int_t wp, Int_t njets, Int_t nbjets) const;
  
	Double_t GetTTEff2bq(Int_t njets) const{return e2bq_[njets];}
	void   SetTTEff2bq(Int_t njets, Double_t value) const{e2bq_[njets] = value;};
  
	Double_t GetTTEff1bq(Int_t njets) const{ return e1bq_[njets];}
	void   SetTTEff1bq(Int_t njets, Double_t value) const{e1bq_[njets] = value;};
  
	Double_t GetTTEff0bq(Int_t njets) const{ return e0bq_[njets];}
	void   SetTTEff0bq(Int_t njets, Double_t value) const{e0bq_[njets] = value;};
  
	Double_t GetBtagWorkingPoint(Int_t idx)             const{return BtagWorkingPoint_[idx];};
	void   SetBtagWorkingPoint(Int_t idx, Double_t wp)  const{BtagWorkingPoint_[idx] = wp;};
  
	void   SetNbOfPE(Int_t NbOfPE) {NbOfPE_ = NbOfPE;};
    // private functions
private:
	
    //Important:
    //in the following, njets denotes the index of the Njets_ array :
  
    ///////////////////////////////////////////////////////////////////	
    // Methods to calculates the estimators ...  (from equations)
    ///////////////////////////////////////////////////////////////////
	Double_t Nbjets   (Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t njets, Int_t nbjets) const;
	Double_t Ntt_bjets(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t njets, Int_t nbjets) const;
	Double_t Nvb_bjets(Double_t Nvb, Double_t eb, Double_t euds,  Int_t njets, Int_t nbjets) const;
	Double_t Nv_bjets( Double_t Nv,  Double_t euds, Int_t njets, Int_t nbjets) const;
	Double_t Ntt_err_bjets(Double_t Ntt, Double_t Ntt_err, Double_t eb, Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t njets, Int_t nbjets) const;
	Double_t Nv_err_bjets( Double_t Nv,  Double_t Nv_err,  Double_t euds, Double_t euds_err, Int_t njets, Int_t nbjets) const;
  
	Double_t N0bjet   (Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t n) const;
	Double_t Ntt_0bjet(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) const;
	Double_t Nvb_0bjet(Double_t Nvb, Double_t eb, Double_t euds,  Int_t n) const;
	Double_t Nv_0bjet (Double_t Nv,  Double_t euds,  Int_t n) const;
	Double_t Ntt_err_0bjet(Double_t Ntt, Double_t Ntt_err, Double_t eb,   Double_t eb_err,   Double_t eudsc, Double_t eudsc_err, Int_t n) const;
	Double_t Nv_err_0bjet (Double_t Nv,  Double_t Nv_err,  Double_t euds, Double_t euds_err, Int_t n) const;
	
	Double_t N1bjet   (Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t n) const;
	Double_t Ntt_1bjet(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) const;
	Double_t Nvb_1bjet(Double_t Nvb, Double_t eb, Double_t euds,  Int_t n) const;
	Double_t Nv_1bjet (Double_t Nv,  Double_t euds,  Int_t n) const;
	Double_t Ntt_err_1bjet(Double_t Ntt, Double_t Ntt_err, Double_t eb,   Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) const;
	Double_t Nv_err_1bjet (Double_t Nv,  Double_t Nv_err,  Double_t euds, Double_t euds_err, Int_t n) const;
	
	Double_t N2bjets   (Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t n) const;
	Double_t Ntt_2bjets(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) const;
	Double_t Nvb_2bjets(Double_t Nvb, Double_t eb, Double_t euds,  Int_t n) const;
	Double_t Nv_2bjets (Double_t Nv,  Double_t euds,  Int_t n) const;
	Double_t Ntt_err_2bjets(Double_t Ntt, Double_t Ntt_err, Double_t eb,   Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) const;
	Double_t Nv_err_2bjets (Double_t Nv,  Double_t Nv_err,  Double_t euds, Double_t euds_err, Int_t n) const;
	
	Double_t N3bjets   (Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t n) const;
	Double_t Ntt_3bjets(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) const;
	Double_t Nvb_3bjets(Double_t Nvb, Double_t eb, Double_t euds,  Int_t n) const;
	Double_t Nv_3bjets (Double_t Nv,  Double_t euds,  Int_t n) const;
	Double_t Ntt_err_3bjets(Double_t Ntt, Double_t Ntt_err, Double_t eb,   Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) const;
	Double_t Nv_err_3bjets (Double_t Nv,  Double_t Nv_err,  Double_t euds, Double_t euds_err, Int_t n) const;
  /*	
   Double_t eudsc_fromN3bjets(Double_t N3bjets, Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc_old, Double_t euds, Int_t n) const;
   Double_t    eb_fromN2bjets(Double_t N2bjets, Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb_old, Double_t eudsc, Double_t euds, Int_t n) const;
   //Double_t Ntt_fromN2bjets(Double_t N2bjets, Double_t Nv, Double_t eb, Double_t eudsc, Int_t n) const;
   Double_t   Ntt_fromN1bjet( Double_t N1bjet,  Double_t Nv,  Double_t Nvb, Double_t eb,  Double_t eudsc, Double_t euds,  Int_t n) const;
   Double_t   Nvb_fromN1bjet( Double_t N1bjet,  Double_t Ntt, Double_t Nv,  Double_t eb,  Double_t eudsc, Double_t euds,  Int_t n) const;
   Double_t  euds_fromN1bjet( Double_t N1bjet,  Double_t Ntt, Double_t Nv,  Double_t Nvb, Double_t eb,    Double_t eudsc, Double_t euds_old, Int_t n) const;
   Double_t  euds_fromN0bjet( Double_t N0bjet,  Double_t Ntt, Double_t Nv,  Double_t Nvb, Double_t eb,    Double_t eudsc, Double_t euds_old, Int_t n) const;
   Double_t    Nv_fromN0bjet( Double_t N0bjet,  Double_t Ntt, Double_t Nvb, Double_t eb,  Double_t eudsc, Double_t euds,  Int_t n) const;
   */	
private:
  Bool_t MCdata_;
	/** Boolean indicating if running on real data or Monte Carlo simulations. */
	UInt_t  NbOfDatasets_;
	/** Number of datasets (in case of Monte Carlo simulations). */
	vector<Int_t> iDatasetsTTLike_;
	/** List of datasets indeces for what is considered as tt-like events (in case of Monte Carlo simulations). */
	vector<Int_t> iDatasetsVLike_;
	/** List of datasets indeces for what is considered as V-like events (in case of Monte Carlo simulations). */
	vector<Int_t> iDatasetsVbLike_;
	/** List of datasets indeces for what is considered as Vb-like events (in case of Monte Carlo simulations). */
  
	UInt_t  NbOfBJetsBins_;
	/** Number of b-jet multiplicity considered ( == 4 for 0,1,2 and >=3 jets). */
	UInt_t  NbOfJetsBins_;
	/** Number of jet multiplicity considered ( == 3 for 4,5 and >=6 jets). */
	UInt_t *Njets_; //[NbOfJetsBins_]
	/** Array containing the jet multiplicities. */
	UInt_t  NbOfBtagWorkingPoint_;
	/** Number of different b-tagging working point used by the equation solver. */
	Float_t   *BtagWorkingPoint_; //[NbOfBtagWorkingPoint_]
	/** Values of the b-tagging working points used by the equation solver. */
	vector< vector< vector< vector< Double_t > > > > N_;
	/** initial number of events (per btagging working point, per jet multiplicity and per dataset). */
	vector< vector< vector< vector< Double_t > > > > N_err_;
	/** statistical errors on the initial number of events (per btagging working point, per jet multiplicity and per dataset). */
	vector< vector< vector< vector< Double_t > > > > N_err_hi_;
	/** statistical upper errors on the initial number of events (per btagging working point, per jet multiplicity and per dataset). */
	vector< vector< vector< vector< Double_t > > > > N_err_lo_;
	/** statistical lower errors on the initial number of events (per btagging working point, per jet multiplicity and per dataset). */
	vector< vector< vector< Double_t > > > Nbjets_;
	/** initial number of events (per btagging working point and per jet multiplicity). */
	vector< vector< vector< Double_t > > > MultiJet_Est_;
  /** Estimated number of multi-jet events (per btagging working point and per jet multiplicity) */
	
  vector<vector<Double_t> >  eb_mc_;
	/** B-tagging efficiency for tt-like events, calculated from MC. Depends on the btagging working point (per jet multiplicity). */
  vector<vector<Double_t> >  eb_err_mc_;
	/** associated error */
	vector<vector<Double_t> >  eudsc_mc_;
	/** Mis-tagging efficiency for tt-like events, calculated from MC. Depends on the btagging working point (per jet multiplicity). */
	vector<vector<Double_t> >  eudsc_err_mc_;
	/** associated error */
	vector<vector<Double_t> >  euds_mc_;
	/** Mis-tagging efficiency for V-like events, calculated from MC. Depends on the btagging working point (per jet multiplicity). */
	vector<vector<Double_t> >  euds_err_mc_;
	/** associated error */
  
  vector<vector<Double_t> >  eb_;
	/** B-tagging estimation parameter for tt-like events. Depends on the btagging working point (per jet multiplicity). */
  vector<vector<Double_t> >  eb_err_;
	/** Error on the B-tagging estimation parameter for tt-like events. Depends on the btagging working point (per jet multiplicity). */
	vector<vector<Double_t> >  eb_err_up_;
	vector<vector<Double_t> >  eb_err_down_;
	vector<vector<Double_t> >  eudsc_;
	/** Mis-tagging estimation parameter for tt-like events. Depends on the btagging working point (per jet multiplicity). */
	vector<vector<Double_t> >  eudsc_err_;
	/** Error on the Mis-tagging estimation parameter for tt-like events. Depends on the btagging working point (per jet multiplicity). */
	vector<vector<Double_t> >  eudsc_err_up_;
	vector<vector<Double_t> >  eudsc_err_down_;
	vector<vector<Double_t> >  euds_;
	/** Mis-tagging estimation parameter for V-like events. Depends on the btagging working point (per jet multiplicity). */
	vector<vector<Double_t> >  euds_err_;
	/** Error on the Mis-tagging estimation parameter for V-like events. Depends on the btagging working point (per jet multiplicity). */
	vector<vector<Double_t> >  euds_err_up_;
	vector<vector<Double_t> >  euds_err_down_;
  
  /*** Global estimation parameters depending on the jet multiplicity (per b-jet multiplicity) */
  Double_t *Ntt_ ; //[NbOfJetsBins_]
	/** Estimation parameter for the number of selected tt-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
  Double_t *Ntt_err_ ; //[NbOfJetsBins_]
	/** Error on the estimation parameter for the number of selected tt-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
  Double_t *Ntt_err_up_ ; //[NbOfJetsBins_]
	/** Upper limit error on the estimation parameter for the number of selected tt-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
  Double_t *Ntt_err_down_ ; //[NbOfJetsBins_]
	/** Lower limit error on the estimation parameter for the number of selected tt-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
	Double_t *Nv_; //[NbOfJetsBins_]
	/** Estimation parameter for the number of selected V-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
	Double_t *Nv_err_; //[NbOfJetsBins_]
	/** Error on the estimation parameter for the number of selected V-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
  Double_t *Nv_err_up_ ; //[NbOfJetsBins_]
	/** Upper limit error on the estimation parameter for the number of selected V-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
  Double_t *Nv_err_down_ ; //[NbOfJetsBins_]
	/** Lower limit error on the estimation parameter for the number of selected V-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
	Double_t *Nvb_; //[NbOfJetsBins_]
	/** Estimation parameter for the number of selected Vb-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
	Double_t *Nvb_err_; //[NbOfJetsBins_]
	/** Estimation parameter for the number of selected Vb-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
  
    // Parameter defining the probability for a b-quark from ttbar final state to lead to a selected b-jet
	Double_t *ebq_; //[NbOfJetsBins_]
	Double_t *e0bq_; //[NbOfJetsBins_]
	/** Parameter defining the probability to select 0 b-quark from ttbar decay in the final state */
	Double_t *e1bq_; //[NbOfJetsBins_]
	/** Parameter defining the probability to select 1 b-quark from ttbar decay in the final state */
	Double_t *e2bq_; //[NbOfJetsBins_]
	/** Parameter defining the probability to select 2 b-quarks from ttbar decay in the final state */
	
	Int_t NbOfPE_;
	/** Number of Pseudo-experiments to be run by RooFit when switching doToyMC boolean to true) */
  
	vector<Double_t> init_Nttlike_;
	/** Initial values of Nttlike for the MLE */
	vector<Double_t> init_Nvlike_;
	/** Initial values of Nvlike for the MLE */
	vector< vector<Double_t> > init_Eb_;
	/** Initial values of eb for the MLE */
	vector< vector<Double_t> > init_Eudsc_;
	/** Initial values of eudsc for the MLE */
	vector< vector<Double_t> > init_Euds_;
	/** Initial values of euds for the MLE */
	
	Double_t *minValue_; //[NbOfJetsBins_]
	/** Threshold on the estimation function value */
  
  /* Dans CheckEstimationLinearity(...) */
  TGraphErrors* RescaledTTLikeEstimation; //
	/** Histograms for the stability check : */
  TGraphErrors* RescaledVLikeEstimation; //
	/** Histograms for the stability check : */
  TGraphErrors* RescaledVbLikeEstimation; //
	/** Histograms for the stability check ; nulle part */
	
  /* Dans CheckEstimationLinearity(...) */
  TCanvas* tCanva_RescaledTTLikeEstimation; //
	/** Canva for the stability check histogram */
  TCanvas* tCanva_RescaledVLikeEstimation; //
	/** Canva for the stability check histogram */
  TCanvas* tCanva_RescaledVbLikeEstimation; //
	/** Canva for the stability check histogram */
  
    //TH2F***hFlavorHistory_;
    ///** Flavor history */
	/*  TH1F***hNbOfBGenJets_; //NbOfDatasets_]*/
	vector< vector <TH1F> > hNbOfBGenJets_;
	/** Histograms for b-genjets */
	vector< TH3F > hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_;   /*NbOfDatasets_]*/
	/** Histograms for b-jets (b-tagging) */
	vector< vector <TH3F> > hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;   /*NbOfDatasets_]NbOfBtagWorkingPoint_]*/
	/** Histograms for b-jets, tagged as b-jets (b-tagging) */
	vector< TH3F > hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_;   /*NbOfDatasets_]*/
	/** Histograms for c-jets (mis-tagging) */
	vector< vector <TH3F> >  hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;   /*NbOfDatasets_]NbOfBtagWorkingPoint_]*/
	/** Histograms for c-jets, tagged as b-jets (mis-tagging) */
	vector< TH3F > hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_;   /*NbOfDatasets_]*/
	/** Histograms for uds-jets (mis-tagging) */
	vector< vector <TH3F> > hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;   /*NbOfDatasets_]NbOfBtagWorkingPoint_]*/
	/** Histograms for uds-jets, tagged as b-jets (mis-tagging) */
	vector< TH3F > hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_;   /*NbOfDatasets_]*/
	/** Histograms for non b-jets (mis-tagging) */
	vector< vector <TH3F> > hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;   /*NbOfDatasets_]NbOfBtagWorkingPoint_]*/
	/** Histograms for non b-jets, tagged as b-jets (mis-tagging) */
  
  /*
   Dans ComputeEffFromMC()
   */
  vector< vector< TEfficiency* > > bTagEff_vs_Njets_; //
	/** TEfficiency for the calculated b-tagging efficiency as a function of the jet multiplicity per dataset and per b-tagging working points */
  vector< vector< TEfficiency* > > cTagEff_vs_Njets_; //
	/** TEfficiency for the calculated c-tagging efficiency as a function of the jet multiplicity per dataset and per b-tagging working points */
  vector< vector< TEfficiency* > > udsTagEff_vs_Njets_; //
	/** TEfficiency for the calculated uds-tagging efficiency as a function of the jet multiplicity per dataset and per b-tagging working points */
  vector< vector< TEfficiency* > > misTagEff_vs_Njets_; //
	/** TEfficiency for the calculated mis-tagging efficiency as a function of the jet multiplicity per dataset and per b-tagging working points */
  
  vector< TGraphAsymmErrors* > bTagEff_vs_Njets_TTlike_; //
	/** TGraphAsymmErrors for the calculated b-tagging efficiency as a function of the jet multiplicity per b-tagging working points */
  vector< TGraphAsymmErrors* > misTagEff_vs_Njets_TTlike_; //
	/** TGraphAsymmErrors for the calculated mis-tagging efficiency as a function of the jet multiplicity per b-tagging working points */
  vector< TGraphAsymmErrors* > misTagEff_vs_Njets_Vlike_; //
	/** TGraphAsymmErrors for the calculated mis-tagging efficiency as a function of the jet multiplicity per b-tagging working points */
	
	vector< vector < vector <TH1F> > > hNbjets_mc_; 
	/**Monte-Carlo normalized b-jets distribution per b-tagging working points, per jet multiplicity (for ttlike, Vblike and Vlike events) */
	vector< vector < vector <TH1F> > > hNbjets_pdf_mc_; 
	/**Monte-Carlo pdf for the b-jets distribution per b-tagging working points, per jet multiplicity (for ttlike, Vblike and Vlike events)*/
	vector< vector < vector <TH1F> > > hNbjets_pdf_est_; 
	/**Estimated pdf for the b-jets distribution per b-tagging working points, per jet multiplicity (for ttlike, Vblike and Vlike events)*/
	
    // Estimation Summary histograms
	vector< vector < vector <TH1F> > > hNjetsEstSummary;
	vector< vector < vector <TH1F> > > hNbjetsEstSummary;
	vector< vector < vector <TH1F> > > hNjetsMCSummary;
	vector< vector < vector <TH1F> > > hNbjetsMCSummary;
  /* Dans Constructeur */
  TLegend* MyLeg; 
  
  /*
   Dans Nulle part
   */
    //	vector< vector < THStack* > > hsNjets_MC; //
    //	vector< vector < THStack* > > hsNjets_Est; //
    //	vector< vector < THStack* > > hsNbjets_MC; //
    //	vector< vector < THStack* > > hsNbjets_Est; //
	
  /* Dans constructeur */
  vector< vector < TCanvas* > > tCanva_Njets_Summary; //1
  vector< vector < TCanvas* > > tCanva_Nbjets_Summary; //1
  
	ClassDef (VJetEstimation,1);
};

#endif

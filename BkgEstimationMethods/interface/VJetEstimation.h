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

#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"

#include "TopTreeAnalysis/Content/interface/MCObsExpectation.h"
#include "TopTreeProducer/interface/TRootJet.h"

// RooFit librairies
#include "RooArgSet.h"
#include "RooAddition.h"
#include "RooCategory.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooHist.h"
//#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
//#include "RooSimultaneous.h"

#include "RooMCStudy.h"
#include "RooMinuit.h"
//#include "RooMinimizer.h"
//#include "RooNLLVar.h"
//#include "RooProfileLL.h"
#include "RooPlot.h"

#include "RooFitResult.h"

using namespace std;
using namespace TopTree;
using namespace RooFit ;

class VJetEstimation {

  public:
  	VJetEstimation(){};
	/** Empty constructor */
	VJetEstimation(int NofBtagWorkingPoint, float* BtagWorkingPoint, int NofJets, int NofJetBins, int NofDatasets, vector<int> iDTTLike, vector<int> iDVLike, vector<int> iDVbLike);
	/** Constructor */
	VJetEstimation(int NofBtagWorkingPoint, float* BtagWorkingPoint, int NofJets, int NofJetBins, double** EffXbq, int NofDatasets, vector<int> iDTTLike, vector<int> iDVLike, vector<int> iDVbLike);
	/** Constructor */
	~VJetEstimation();
	/** Destructor */

	VJetEstimation(const VJetEstimation& vjet);
 	/** Copy constructor (not yet implemented) */
	
	void FillInputs(double**** n, double*** MultiJets_Estimated_Nj);
 	/** Method used to fill the histograms containing the number of events with 0,1,2 and 3 b-jets for all datasets. Needs also the estimated number of multijet events */
	void FillInputs(double**** n);
 	/** Method used to fill the histograms containing the number of events with 0,1,2 and 3 b-jets for all datasets.*/
	void FillInputs(vector<TRootJet*> &SelectedJets, int idx, int btagAlgo);
 	/** Method to used for each event to fill the histograms containing the number of events with 0,1,2 and 3 b-jets for the dataset idx. */
	void FillInputs(vector<TRootJet> &SelectedJets, unsigned int idx, int btagAlgo);
 	/** Method to used for each event to fill the histograms containing the number of events with 0,1,2 and 3 b-jets for the dataset idx. */
	void ReScaleInputs(int idx, float Ntot, double factor);
	/** Rescale the input number of events from dataset idx. Method to be used after the event loop.*/
	float WilsonScoreIntervalLow(float Non, float Ntot);
	/** Upper limit of the Wilson score interval for binomial parameter (being the selection efficiency here) (ArXiv:hep-ph/0905.3831v2)*/
	float WilsonScoreIntervalHigh(float Non, float Ntot);
	/** Lower limit of the Wilson score interval for binomial parameter (being the selection efficiency here) (ArXiv:hep-ph/0905.3831v2)*/
	float WilsonScoreIntervalMean(float Non, float Ntot);
	/** Mean between the upper and lower limit of the Wilson score interval for binomial parameter (being the selection efficiency here) (ArXiv:hep-ph/0905.3831v2)*/
	void SumOverAllInputs();
	/** Sum the weighted contribution of all the datasets. Method to be used after having looped over all events of all datasets.*/
	void ComputeEffbqFromMC(int idx);
	/** Compute ebq efficiencies from MC. Method to be used after having looped over all events of all datasets.*/
	void ComputeEffFromMC();
	/** Compute the b/mis-tagging efficiencies from MC. Method to be used after having looped over all events of all datasets.*/
        void ComputeCondProb(bool Verbose = false);
	/** Compute the conditional probabilities */
	
	void FillSummaryHistos();
	/** Fill summary histograms. */
	void BckgdSubstraction(vector<MCObsExpectation*> &hists, vector<string> &name, float lumi);
	/** Substract background b-jet ditribution provided as input */
	void Write(TFile* file, string label = string(""), bool verbose = false);
 	/** Method to book histograms */
	
	/////////////////////////////////
	//Main Methods to be executed ...
	/////////////////////////////////
	void SetInitialValues(double **init) { init_ = init;};
	/** Set initial values for the minizer. */
	void BinnedMaximumLikelihoodEst(int btag_wp_idx, int njets, string &method, string &option, vector<int> &FixedVarIdx, bool doScan, bool Verbose);
	/** Equation solver based on a likelihood function minimization. */
	void BinnedMaximumLikelihoodEst(int njets, string &method, string &option, vector<int> &FixedVarIdx, bool doScan, bool Verbose){for(int	i=0;i<NbOfBtagWorkingPoint_;i++) BinnedMaximumLikelihoodEst(i, njets, method, option, FixedVarIdx, doScan, Verbose);};
	/** Equation solver based on a likelihood function minimization (solve for all b-tagging working points). */
	void BinnedMaximumLikelihoodEst(string &method, string &option, vector<int> &FixedVarIdx, bool doScan, bool Verbose){for(int i=0;i<NbOfJetsBins_;i++) BinnedMaximumLikelihoodEst(i, method, option, FixedVarIdx, doScan, Verbose);};
	/** Equation solver based on a likelihood function minimization (solve for all jet multiplicity). */

	void BinnedMaximumJointWPLikelihoodEst(int njets, string &method, string &option, vector<int> &FixedVarIdx, bool doScan, bool Verbose);
	/** Equation solver based on a joint likelihood function minimization, combining all the b-tagging working points. */
	void BinnedMaximumJointWPLikelihoodEst(string &method, string &option, vector<int> &FixedVarIdx, bool doScan, bool Verbose){for(int i=0;i<NbOfJetsBins_;i++) BinnedMaximumJointWPLikelihoodEst(i, method, option, FixedVarIdx, doScan, Verbose);};
	/** Equation solver based on a joint likelihood function minimization, combining all the b-tagging working points (solve for all jet multiplicity). */

	void BinnedMaximumJointLikelihoodEst(string &method, string &option, bool doScan, bool Verbose);
	/** Equation solver based on a joint likelihood function minimization, combining all the jet multiplicities. */


	void UnBinnedMaximumLikelihoodEst(int btag_wp_idx, int njets, vector<int> &FixedVarIdx, bool doMinos, bool doToyMC, bool Verbose);
	/** Maximum likelihood estimation with RooFit. */
	void UnBinnedMaximumLikelihoodEst(int njets, vector<int> &FixedVarIdx, bool doMinos, bool doToyMC, bool Verbose){for(int i=0;i<NbOfBtagWorkingPoint_;i++) UnBinnedMaximumLikelihoodEst(i, njets, FixedVarIdx, doMinos, doToyMC, Verbose);};
	/** Maximum likelihood estimation with RooFit (solve for all b-tagging working points). */
	void UnBinnedMaximumLikelihoodEst(vector<int> &FixedVarIdx, bool doMinos, bool doToyMC, bool Verbose){for(int i=0;i<NbOfJetsBins_;i++) UnBinnedMaximumLikelihoodEst(i, FixedVarIdx, doMinos, doToyMC, Verbose);};
	/** Maximum likelihood estimation with RooFit (solve for all jet multiplicities). */
	
	void UnBinnedMaximumJointWPLikelihoodEst(int njets, vector<int> &FixedVarIdx, bool doMinos, bool Verbose);
	/** Maximum likelihood estimation with RooFit, combining all the b-tagging working points. */
	void UnBinnedMaximumJointWPLikelihoodEst(vector<int> &FixedVarIdx, bool doMinos, bool Verbose){for(int i=0;i<NbOfJetsBins_;i++) UnBinnedMaximumJointWPLikelihoodEst(i, FixedVarIdx, doMinos, Verbose);};
	/** Maximum likelihood estimation with RooFit, combining all the b-tagging working points. */

	void UnBinnedMaximumNjetsLikelihoodEst(int btag_wp_idx, vector<int> &FixedVarIdx, bool doMinos, bool Verbose);
	/** Maximum likelihood estimation with RooFit, combining all the jet multiplicities. */
	
	void UnBinnedMaximumJointLikelihoodEst(bool Verbose);
	/** Maximum likelihood estimation with RooFit, combining both all the jet multiplicities and all b-tagging working points. */
	
	/////////////////////
	// Print  methods
	////////////////////
	void PrintInputs();
	/** Method to print the number of events passed as imputs per dataset*/
	void PrintInputs(int njets);
	/** Method to print the number of events passed as imputs per dataset, per jet multiplicity*/
	void PrintResults();
	/** Method to print the output of the estimation method*/
	void PrintResults(int njets);
	/** Method to print the output of the estimation method, per jet multiplicity*/
	void PrintResults_LatexFormat(ofstream &ofile);
	/** Method to print a latex-formatted output of the estimation method*/
	
	/////////////////////
	// Cross-check methods
	/////////////////////
	bool CheckEstimation(double threshold);
	/** Returns true if the value of the minimizer function at the minimum is below the threshold. */
	bool CheckEstimation(double threshold, int njets);
	/** Returns true if the value of the minimizer function at the minimum is below the threshold. */
	void CheckEstimationLinearity(int nbOfRescaleFact, double **RescaleFact, int Idx);
	/** Method to check if the response of the estimation scales linearly with the inputs. */

	//////////////////////////////////////
	// Access to the class members
	//////////////////////////////////////

	int GetNbOfBtagWorkingPoint()    const {return NbOfBtagWorkingPoint_;};
	/** returns the number of b-tagging working points used by the estimation method */
	int GetNbOfDatasets()            const {return NbOfDatasets_;};
	/** returns the number of datasets used by the estimation method */
	int GetNbOfBJetsBins()           const {return NbOfBJetsBins_;};
	/** returns the number of b-jet bins used by the estimation method */
	int GetNbOfJetsBins()            const {return NbOfJetsBins_;};
	/** returns the number of jet bins used by the estimation method */
	int GetNjets(int jetidx)         const {return (jetidx<NbOfJetsBins_ ? Njets_[jetidx] : -999);}
	double ***** GetCondProb()       const {return condProb_;};
	/** returns the conditional probabilities (needed by the VJetEstPseudoExp class*/

	vector<int> GetiDatasetsTTLike() const {return iDatasetsTTLike_;};
	vector<int> GetiDatasetsVLike()  const {return iDatasetsVLike_;};
	vector<int> GetiDatasetsVbLike() const {return iDatasetsVbLike_;};

	double****GetN()                                           const{ return N_;};
	double*** GetN(int wp)                                     const{ return N_[wp];};
	double**  GetN(int wp, int njets)                          const{ return N_[wp][njets];};
	double*   GetN(int wp, int njets, int nbjets)              const{ return N_[wp][njets][nbjets];};
	double    GetN(int wp, int njets, int nbjets, int dataset) const{ return N_[wp][njets][nbjets][dataset];};

	double*** GetNbjets()                              const{ return Nbjets_;};
	double**  GetNbjets(int wp)                        const{ return Nbjets_[wp];};
	double*   GetNbjets(int wp, int njets)             const{ return Nbjets_[wp][njets];};
	double    GetNbjets(int wp, int njets, int nbjets) const{ return Nbjets_[wp][njets][nbjets];};
	
	double   GetMinLL(int njets) const{ return minValue_[njets];};

	double*  GetMultiJet_Est(int wp) const;
	double*  GetMultiJet_Est(int wp, int njets) const;
	double*  GetMultiJet_Est(int wp, int njets, int nbjets) const;

	double   GetPredEb(int wp) const;
	double   GetPredEb(int wp, int njets) const;
	double   GetPredEbErr(int wp) const;
	double   GetPredEbErr(int wp, int njets) const;
	void     SetPredEb(int wp, int njets, double value) const;

	double   GetPredEudsc(int wp) const;
	double   GetPredEudsc(int wp, int njets) const;
	double   GetPredEudscErr(int wp) const;
	double   GetPredEudscErr(int wp, int njets) const;
	void     SetPredEudsc(int wp, int njets, double value) const;

	double   GetPredEuds(int wp) const;
	double   GetPredEuds(int wp, int njets) const;
	double   GetPredEudsErr(int wp) const;
	double   GetPredEudsErr(int wp, int njets) const;
	void     SetPredEuds(int wp, int njets, double value) const;

	double   GetPredNv(int wp) const;
	double   GetPredNv(int wp, int njets) const;
	double   GetPredNv(int wp, int njets, int nbjets) const;
	double   GetPredNvErr(int wp) const;
	double   GetPredNvErr(int wp, int njets) const;
	double   GetPredNvErr(int wp, int njets, int nbjets) const;

	double   GetPredNvb(int wp) const;
	double   GetPredNvb(int wp, int njets) const;
	double   GetPredNvb(int wp, int njets, int nbjets) const;
	double   GetPredNvbErr(int wp) const;
	double   GetPredNvbErr(int wp, int njets) const;
	double   GetPredNvbErr(int wp, int njets, int nbjets) const;

	double   GetPredNtt(int wp) const;
	double   GetPredNtt(int wp, int njets) const;
	double   GetPredNtt(int wp, int njets, int nbjets) const;
	double   GetPredNttErr(int wp) const;
	double   GetPredNttErr(int wp, int njets) const;
	double   GetPredNttErr(int wp, int njets, int nbjets) const;

	double   GetPredNtotal(int wp) const;
	double   GetPredNtotal(int wp, int njets) const;
	double   GetPredNtotal(int wp, int njets, int nbjets) const;
	double   GetPredNtotalErr(int wp) const;
	double   GetPredNtotalErr(int wp, int njets) const;
	double   GetPredNtotalErr(int wp, int njets, int nbjets) const;

	double   GetPredN(int idx) const;
	double   GetPredN(int idx, int njets) const;

	//////////////////////////////////////
	// Access to the parameter estimated values
	//////////////////////////////////////

	double GetEstEb(int wp) const;
	double GetEstEb(int wp, int njets) const;
	double GetEstEbErr(int wp) const;
	double GetEstEbErr(int wp, int njets) const;

	double GetEstEudsc(int wp) const;
	double GetEstEudsc(int wp, int njets) const;
	double GetEstEudscErr(int wp) const;
	double GetEstEudscErr(int wp, int njets) const;

	double GetEstEuds(int wp) const;
	double GetEstEuds(int wp, int njets) const;
	double GetEstEudsErr(int wp) const;
	double GetEstEudsErr(int wp, int njets) const;

	double GetEstNv(int wp) const;
	double GetEstNv(int wp, int njets) const;
	double GetEstNv(int wp, int njets, int nbjets) const;
	double GetEstNvErr(int wp) const;
	double GetEstNvErr(int wp, int njets) const;
	double GetEstNvErrUp(  int wp, int njets) const;
	double GetEstNvErrDown(int wp, int njets) const;
	double GetEstNvErr(int wp, int njets, int nbjets) const;

	double GetEstNvb(int wp) const;
	double GetEstNvb(int wp, int njets) const;
	double GetEstNvb(int wp, int njets, int nbjets) const;
	double GetEstNvbErr(int wp) const;
	double GetEstNvbErr(int wp, int njets) const;
	double GetEstNvbErr(int wp, int njets, int nbjets) const;

	double GetEstNtt(int wp) const;
	double GetEstNtt(int wp, int njets) const;
	double GetEstNtt(int wp, int njets, int nbjets) const;
	double GetEstNttErr(int wp) const;
	double GetEstNttErr(int wp, int njets) const;
	double GetEstNttErrUp(  int wp, int njets) const;
	double GetEstNttErrDown(int wp, int njets) const;
	double GetEstNttErr(int wp, int njets, int nbjets) const;

	double GetEstNtotal(int wp) const;
	double GetEstNtotal(int wp, int njets) const;
	double GetEstNtotal(int wp, int njets, int nbjets) const;

	double GetTTEff2bq(int njets) const{return e2bq_[njets];}
	void   SetTTEff2bq(int njets, double value) const{e2bq_[njets] = value;};

	double GetTTEff1bq(int njets) const{ return e1bq_[njets];}
	void   SetTTEff1bq(int njets, double value) const{e1bq_[njets] = value;};

	double GetTTEff0bq(int njets) const{ return e0bq_[njets];}
	void   SetTTEff0bq(int njets, double value) const{e0bq_[njets] = value;};

	double GetBtagWorkingPoint(int idx)             const{return BtagWorkingPoint_[idx];};
	void   SetBtagWorkingPoint(int idx, double wp)  const{BtagWorkingPoint_[idx] = wp;};

	///////////////////////////////////////////////////////////////////
	// Getting and setting statistical and systematic errors
	///////////////////////////////////////////////////////////////////

	// Statistical and systematic erros (retrieved from WJetsPseudoExp)

	void SetNvVar (int wp, int njets, int nbjets, double error) {NvVar_[wp][njets][nbjets] = error;}
	void SetNvbVar(int wp, int njets, int nbjets, double error) {NvbVar_[wp][njets][nbjets] = error;}
	void SetNttVar(int wp, int njets, int nbjets, double error) {NttVar_[wp][njets][nbjets] = error;}

	void SetNvBias (int wp, int njets, int nbjets, double error)  {NvBias_[wp][njets][nbjets]  = error;}
	void SetNvbBias(int wp, int njets, int nbjets, double error)  {NvbBias_[wp][njets][nbjets]  = error;}
	void SetNttBias(int wp, int njets, int nbjets, double error)  {NttBias_[wp][njets][nbjets]  = error;}

	double GetNvVar   (int wp) const;
	double GetNvBias  (int wp) const;
	double GetNvMSE   (int wp) const;
	double GetNvbVar  (int wp) const;
	double GetNvbBias (int wp) const;
	double GetNvbMSE  (int wp) const;
	double GetNttVar  (int wp) const;
	double GetNttBias (int wp) const;
	double GetNttMSE  (int wp) const;

	double GetNvVar   (int wp, int njets) const {return NvVar_[wp][njets][NbOfBJetsBins_];}
	double GetNvBias  (int wp, int njets) const {return NvBias_[wp][njets][NbOfBJetsBins_];}
	double GetNvMSE   (int wp, int njets) const {return sqrt(pow(GetNvBias(wp,njets),2)+pow(GetNvVar(wp,njets),2));}
	double GetNvbVar  (int wp, int njets) const {return NvbVar_[wp][njets][NbOfBJetsBins_];}
	double GetNvbBias (int wp, int njets) const {return NvbBias_[wp][njets][NbOfBJetsBins_];}
	double GetNvbMSE  (int wp, int njets) const {return sqrt(pow(GetNvbBias(wp,njets),2)+pow(GetNvbVar(wp,njets),2));}
	double GetNttVar  (int wp, int njets) const {return NttVar_[wp][njets][NbOfBJetsBins_];}
	double GetNttBias (int wp, int njets) const {return NttBias_[wp][njets][NbOfBJetsBins_];}
	double GetNttMSE  (int wp, int njets) const {return sqrt(pow(GetNttBias(wp,njets),2)+pow(GetNttVar(wp,njets),2));}

	double GetNvVar   (int wp, int njets, int nbjets) const {return (njets==-1 ?  NvVar_[wp][NbOfJetsBins_][nbjets] :  NvVar_[wp][njets][nbjets]);}
	double GetNvBias  (int wp, int njets, int nbjets) const {return (njets==-1 ? NvBias_[wp][NbOfJetsBins_][nbjets] : NvBias_[wp][njets][nbjets]);}
	double GetNvMSE   (int wp, int njets, int nbjets) const {return sqrt(pow(GetNvBias(wp,njets,nbjets),2)+pow(GetNvVar(wp,njets,nbjets),2));}
	double GetNvbVar  (int wp, int njets, int nbjets) const {return (njets==-1 ?  NvbVar_[wp][NbOfJetsBins_][nbjets] :  NvbVar_[wp][njets][nbjets]);}
	double GetNvbBias (int wp, int njets, int nbjets) const {return (njets==-1 ? NvbBias_[wp][NbOfJetsBins_][nbjets] : NvbBias_[wp][njets][nbjets]);}
	double GetNvbMSE  (int wp, int njets, int nbjets) const {return sqrt(pow(GetNvbBias(wp,njets,nbjets),2)+pow(GetNvbVar(wp,njets,nbjets),2));}
	double GetNttVar  (int wp, int njets, int nbjets) const {return (njets==-1 ?  NttVar_[wp][NbOfJetsBins_][nbjets] :  NttVar_[wp][njets][nbjets]);}
	double GetNttBias (int wp, int njets, int nbjets) const {return (njets==-1 ? NttBias_[wp][NbOfJetsBins_][nbjets] : NttBias_[wp][njets][nbjets]);}
	double GetNttMSE  (int wp, int njets, int nbjets) const {return sqrt(pow(GetNttBias(wp,njets,nbjets),2)+pow(GetNttVar(wp,njets,nbjets),2));}

   friend class VJetEstPseudoExp;

   // private functions
   private:
	
	//Important:
	//in the following, njets denotes the index of the Njets_ array :

	//Calculate a chisq value between the estimated and predicted values of Nbjets
	double LikelihoodFunct(          const double *xx) const ;
	double JointNjetsLikelihoodFunct(const double *xx) const ;
	double JointWPLikelihoodFunct(   const double *xx) const ;
	double JointLikelihoodFunct(     const double *xx) const ;
	double MinimizerFunct(           const double *xx) const ;

	///////////////////////////////////////////////////////////////////	
	// Methods to calculates the estimators ...  (from equations)
	///////////////////////////////////////////////////////////////////
	double Nbjets   (double &Ntt, double &Nv, double &Nvb, double &eb, double &eudsc,  double &euds, int &njets, int &nbjets) const;
	double Ntt_bjets(double &Ntt, double &eb, double &eudsc, int &njets, int &nbjets) const;
	double Nvb_bjets(double &Nvb, double &eb, double &euds,  int &njets, int &nbjets) const;
	double Nv_bjets( double &Nv,  double &euds, int &njets, int &nbjets) const;
	double Ntt_err_bjets(double &Ntt, double &Ntt_err, double &eb, double &eb_err, double &eudsc, double &eudsc_err, int &njets, int &nbjets) const;
	double Nv_err_bjets( double &Nv,  double &Nv_err,  double &euds, double &euds_err, int &njets, int &nbjets) const;

	double N0bjet   (double &Ntt, double &Nv, double &Nvb, double &eb, double &eudsc,  double &euds, int &n) const;
	double Ntt_0bjet(double &Ntt, double &eb, double &eudsc, int &n) const;
	double Nvb_0bjet(double &Nvb, double &eb, double &euds,  int &n) const;
	double Nv_0bjet (double &Nv,  double &euds,  int &n) const;
	double Ntt_err_0bjet(double &Ntt, double &Ntt_err, double &eb,   double &eb_err,   double &eudsc, double &eudsc_err, int &n) const;
	double Nv_err_0bjet (double &Nv,  double &Nv_err,  double &euds, double &euds_err, int &n) const;
	
	double N1bjet   (double &Ntt, double &Nv, double &Nvb, double &eb, double &eudsc,  double &euds, int &n) const;
	double Ntt_1bjet(double &Ntt, double &eb, double &eudsc, int &n) const;
	double Nvb_1bjet(double &Nvb, double &eb, double &euds,  int &n) const;
	double Nv_1bjet (double &Nv,  double &euds,  int &n) const;
	double Ntt_err_1bjet(double &Ntt, double &Ntt_err, double &eb,   double &eb_err, double &eudsc, double &eudsc_err, int &n) const;
	double Nv_err_1bjet (double &Nv,  double &Nv_err,  double &euds, double &euds_err, int &n) const;
	
	double N2bjets   (double &Ntt, double &Nv, double &Nvb, double &eb, double &eudsc,  double &euds, int &n) const;
	double Ntt_2bjets(double &Ntt, double &eb, double &eudsc, int &n) const;
	double Nvb_2bjets(double &Nvb, double &eb, double &euds,  int &n) const;
	double Nv_2bjets (double &Nv,  double &euds,  int &n) const;
	double Ntt_err_2bjets(double &Ntt, double &Ntt_err, double &eb,   double &eb_err, double &eudsc, double &eudsc_err, int &n) const;
	double Nv_err_2bjets (double &Nv,  double &Nv_err,  double &euds, double &euds_err, int &n) const;
	
	double N3bjets   (double &Ntt, double &Nv, double &Nvb, double &eb, double &eudsc,  double &euds, int &n) const;
	double Ntt_3bjets(double &Ntt, double &eb, double &eudsc, int &n) const;
	double Nvb_3bjets(double &Nvb, double &eb, double &euds,  int &n) const;
	double Nv_3bjets (double &Nv,  double &euds,  int &n) const;
	double Ntt_err_3bjets(double &Ntt, double &Ntt_err, double &eb,   double &eb_err, double &eudsc, double &eudsc_err, int &n) const;
	double Nv_err_3bjets (double &Nv,  double &Nv_err,  double &euds, double &euds_err, int &n) const;
	
	double eudsc_fromN3bjets(double &N3bjets, double &Ntt, double &Nv, double &Nvb, double &eb, double &eudsc_old, double &euds, int &n) const;
	double    eb_fromN2bjets(double &N2bjets, double &Ntt, double &Nv, double &Nvb, double &eb_old, double &eudsc, double &euds, int &n) const;
	//double Ntt_fromN2bjets(double &N2bjets, double &Nv, double &eb, double &eudsc, int &n) const;
	double   Ntt_fromN1bjet( double &N1bjet,  double &Nv,  double &Nvb, double &eb,  double &eudsc, double &euds,  int &n) const;
	double   Nvb_fromN1bjet( double &N1bjet,  double &Ntt, double &Nv,  double &eb,  double &eudsc, double &euds,  int &n) const;
	double  euds_fromN1bjet( double &N1bjet,  double &Ntt, double &Nv,  double &Nvb, double &eb,    double &eudsc, double &euds_old, int &n) const;
	double  euds_fromN0bjet( double &N0bjet,  double &Ntt, double &Nv,  double &Nvb, double &eb,    double &eudsc, double &euds_old, int &n) const;
	double    Nv_fromN0bjet( double &N0bjet,  double &Ntt, double &Nvb, double &eb,  double &eudsc, double &euds,  int &n) const;
	
   private:
   	bool MCdata_;
	/** Boolean indicating if running on real data or Monte Carlo simulations. */
	int  NbOfDatasets_;
	/** Number of datasets (in case of Monte Carlo simulations). */
	vector<int> iDatasetsTTLike_;
	/** List of datasets indeces for what is considered as tt-like events (in case of Monte Carlo simulations). */
	vector<int> iDatasetsVLike_;
	/** List of datasets indeces for what is considered as V-like events (in case of Monte Carlo simulations). */
	vector<int> iDatasetsVbLike_;
	/** List of datasets indeces for what is considered as Vb-like events (in case of Monte Carlo simulations). */

	int  NbOfBJetsBins_;
	/** Number of b-jet multiplicity considered ( == 4 for 0,1,2 and >=3 jets). */
	int  NbOfJetsBins_;
	/** Number of jet multiplicity considered ( == 3 for 4,5 and >=6 jets). */
	int* Njets_;
	/** Array containing the jet multiplicities. */
	int  NbOfBtagWorkingPoint_;
	/** Number of different b-tagging working point used by the equation solver. */
	float*   BtagWorkingPoint_;
	/** Values of the b-tagging working points used by the equation solver. */
	double*** Nbjets_; 
	/** initial number of events (per btagging working point and per jet multiplicity). */
	double**** N_;
	/** initial number of events (per btagging working point, per jet multiplicity and per dataset). */
	double**** N_err_;
	/** statistical errors on the initial number of events (per btagging working point, per jet multiplicity and per dataset). */
	double*** MultiJet_Est_;
        /** Estimated number of multi-jet events (per btagging working point and per jet multiplicity) */
	
        double** eb_mc_;
	/** B-tagging efficiency for tt-like events, calculated from MC. Depends on the btagging working point (per jet multiplicity). */
        double** eb_err_mc_;
	/** associated error */
	double** eudsc_mc_;
	/** Mis-tagging efficiency for tt-like events, calculated from MC. Depends on the btagging working point (per jet multiplicity). */
	double** eudsc_err_mc_;
	/** associated error */
	double** euds_mc_;
	/** Mis-tagging efficiency for V-like events, calculated from MC. Depends on the btagging working point (per jet multiplicity). */
	double** euds_err_mc_;
	/** associated error */

        double** eb_;
	/** B-tagging estimation parameter for tt-like events. Depends on the btagging working point (per jet multiplicity). */
        double** eb_err_;
	/** Error on the B-tagging estimation parameter for tt-like events. Depends on the btagging working point (per jet multiplicity). */
	double** eudsc_;
	/** Mis-tagging estimation parameter for tt-like events. Depends on the btagging working point (per jet multiplicity). */
	double** eudsc_err_;
	/** Error on the Mis-tagging estimation parameter for tt-like events. Depends on the btagging working point (per jet multiplicity). */
	double** euds_;
	/** Mis-tagging estimation parameter for V-like events. Depends on the btagging working point (per jet multiplicity). */
	double** euds_err_;
	/** Error on the Mis-tagging estimation parameter for V-like events. Depends on the btagging working point (per jet multiplicity). */

	// Global estimation parameters depending on the jet multiplicity (per b-jet multiplicity)
        double* Ntt_ ;
	/** Estimation parameter for the number of selected tt-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
        double* Ntt_err_ ;
	/** Error on the estimation parameter for the number of selected tt-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
        double* Ntt_err_up_ ;
	/** Upper limit error on the estimation parameter for the number of selected tt-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
        double* Ntt_err_down_ ;
	/** Lower limit error on the estimation parameter for the number of selected tt-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
	double* Nv_;
	/** Estimation parameter for the number of selected V-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
	double* Nv_err_;
	/** Error on the estimation parameter for the number of selected V-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
        double* Nv_err_up_ ;
	/** Upper limit error on the estimation parameter for the number of selected V-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
        double* Nv_err_down_ ;
	/** Lower limit error on the estimation parameter for the number of selected V-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
	double* Nvb_;
	/** Estimation parameter for the number of selected Vb-like events. Depends on the jet multiplicity (per b-jet multiplicity). */
	double* Nvb_err_;
	/** Estimation parameter for the number of selected Vb-like events. Depends on the jet multiplicity (per b-jet multiplicity). */

	// Parameter defining the probability for a b-quark from ttbar final state to lead to a selected b-jet
	double* ebq_;
	double* e0bq_;
	/** Parameter defining the probability to select 0 b-quark from ttbar decay in the final state */
	double* e1bq_;
	/** Parameter defining the probability to select 1 b-quark from ttbar decay in the final state */
	double* e2bq_;
	/** Parameter defining the probability to select 2 b-quarks from ttbar decay in the final state */
	
	double***** condProb_;
	/** Conditional probability to have X b-jets in the final state with a certain b-tagging working point, given that there were
	Y b-jets with a looser b-tagging working point. */
	
	ROOT::Math::IMultiGenFunction * fFunc;
	/** Pointer to function to minimize when using the nonlinear solvers. */

	double **init_;
	/** Initial values for the non-linear equation solver */
	
	double *minValue_;
	/** Threshold on the estimation function value */
	
	//Statistical and systematic errors (to be retrieved from VJetsEstPseudoExp) for each b-tagging working point
	double ***NvVar_;
	/** Variance of Nv estimator (row : nb of jets (from 4 to >6), column : nb of b-jets (from 0 to 3, last column=inclusive)) */
	double ***NvbVar_;
	/** Variance of Nvb estimator (row : nb of jets (from 4 to >6), column : nb of b-jets (from 0 to 3, last column=inclusive)) */
	double ***NttVar_;
	/** Variance of Ntt estimator (row : nb of jets (from 4 to >6), column : nb of b-jets (from 0 to 3, last column=inclusive)) */
	double ***NvBias_;
	/** Bias of Nv estimator (row : nb of jets (from 4 to >6), column : nb of b-jets (from 0 to 3, last column=inclusive)) */
	double ***NvbBias_;
	/** Bias of Nvb estimator (row : nb of jets (from 4 to >6), column : nb of b-jets (from 0 to 3, last column=inclusive)) */
	double ***NttBias_;
	/** Bias of Ntt estimator (row : nb of jets (from 4 to >6), column : nb of b-jets (from 0 to 3, last column=inclusive)) */

	// Histograms for the ChiSq-based equation solver
	TGraph***  hScanMin_Eb;
	TGraph***  hScanMin_Eudsc;
	TGraph***  hScanMin_Euds;
	TGraph**   hScanMin_Ntt;
	TGraph**   hScanMin_Nv;
	TGraph**   hScanMin_Nvb;

	// Histograms for the stability check
	TGraphErrors *RescaledTTLikeEstimation;
	TGraphErrors *RescaledVLikeEstimation;
	TGraphErrors *RescaledVbLikeEstimation;
	
	TCanvas *tCanva_RescaledTTLikeEstimation;
	TCanvas *tCanva_RescaledVLikeEstimation;
	TCanvas *tCanva_RescaledVbLikeEstimation;

	TH2F***hFlavorHistory_;
	/** Flavor history */
	TH1F***hNbOfBGenJets_;
	/** Histograms for b-genjets */
	TH3F** hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_;
	/** Histograms for b-jets (b-tagging) */
	TH3F*** hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;
	/** Histograms for b-jets, tagged as b-jets (b-tagging) */
	TH3F** hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_;
	/** Histograms for c-jets (mis-tagging) */
	TH3F*** hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;
	/** Histograms for c-jets, tagged as b-jets (mis-tagging) */
	TH3F** hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_;
	/** Histograms for uds-jets (mis-tagging) */
	TH3F*** hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;
	/** Histograms for uds-jets, tagged as b-jets (mis-tagging) */
	TH3F** hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_;
	/** Histograms for non b-jets (mis-tagging) */
	TH3F*** hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;
	/** Histograms for non b-jets, tagged as b-jets (mis-tagging) */
	TGraphAsymmErrors*** bTagEff_vs_Njets_;
	/** TGraph for the calculated b-tagging efficiency as a function of the jet multiplicity */
	TGraphAsymmErrors*** cTagEff_vs_Njets_;
	/** TGraph for the calculated c-tagging efficiency as a function of the jet multiplicity */
	TGraphAsymmErrors*** udsTagEff_vs_Njets_;
	/** TGraph for the calculated uds-tagging efficiency as a function of the jet multiplicity */
	TGraphAsymmErrors*** misTagEff_vs_Njets_;
	/** TGraph for the calculated mis-tagging efficiency as a function of the jet multiplicity */
	
	TH1F**** hNbjets_mc_;
	/**Monte-Carlo normalized b-jets distribution*/
	TH1F**** hNbjets_pdf_mc_;
	/**Monte-Carlo pdf for the b-jets distribution*/
	TH1F**** hNbjets_pdf_est_;
	/**Estimated pdf for the b-jets distribution*/
	
	// Estimation Summary histograms
	TH1F**** hNjetsEstSummary;
	TH1F**** hNbjetsEstSummary;
	TH1F**** hNjetsMCSummary;
	TH1F**** hNbjetsMCSummary;
	TLegend* MyLeg;

	THStack*** hsNjets_MC;
	THStack*** hsNjets_Est;
	THStack*** hsNbjets_MC;
	THStack*** hsNbjets_Est;
	
	TCanvas*** tCanva_Njets_Summary;
	TCanvas*** tCanva_Nbjets_Summary;	

};

#endif

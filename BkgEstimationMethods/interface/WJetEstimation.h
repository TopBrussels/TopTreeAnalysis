#ifndef WJetEstimation_h
#define WJetEstimation_h

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include "Math/VectorUtil.h"
#include "Math/Polynomial.h"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"

using namespace std;

        /**
        //	Aim: Estimation of the W+jets bkg (#events) with a ME method using b-tagging
	//	This method use 4 estimators in a system of 4 equations for the events where
	//	there is 4, 5, >=6 jets.\n
	//      This estimators are: Eff(b) - Eff(udsc) - NofEvents(W) - NofEvents(ttbar)\n
	//	
	//	It's an iterative process. 20 iterations should be enough in principle.
	//	The results are summed up.
	//	
	//	Example of usage:
	//	- call the constructor (if the number of datasets is not specified, the code won't make use or display any MC related info)
	//	- Fill the vectors and arrays used in the FillInputs
	//	- FillInputs(...)
	//	- Estimation()
	//	- PrintResults() can be called later or not called at all
	//	- Draw()
	//	- Write()
	//
	//	NB: the FillInputs methods take as arguments the number of QCD events estimated in the ABCD method per bin of b-tagging
	//	Missing pieces:
	//	- Implement a Fill method which will run over events and fill automatically
	//	the arrays and vector which are given as arguments to the FillInputs method
	*/

class WJetEstimation {

  public:
	//WJetEstimation();
	//WJetEstimation(int);
	WJetEstimation(int NofDatasets, vector<int> iDTTLike, vector<int> iDWLike);
	WJetEstimation(int Iterations, int NofDatasets, vector<int> iDTTLike, vector<int> iDWLike);
	WJetEstimation(int Iterations);
	~WJetEstimation();

 	//copy constructor
	WJetEstimation(const WJetEstimation& wjet);
	
	//MultiJets_Estimated_Nj: nof events of QCD estimation from data (ABCD) for the bins 0-1-2-3 b-jets !!
	// This method has to be called before all the others to expect results ... ;-)
	//N0 to N3 are the nof events with X bjets per jet multiplicity and datatset  (ex: N0[5jets][ttjets])
	//MultiJets_Estimated_Nj is the nof QCD events estimated per jet multiplicity and per bjet multiplicity (ex: MultiJets_Estimated_Nj[5jets][2bjets])
	//iDatasetsTTLike contains the iterator for TTjets-like events in the previously described array (idem for WJets-like events)
	//void FillInputs(int NbOfDatasets, double** N0, double** N1, double** N2, double** N3, double** MultiJets_Estimated_Nj, vector<int> iDatasetsTTLike, vector<int> iDatasetsWLike);
	void FillInputs(double** n0, double** n1, double** n2, double** n3, double** MultiJets_Estimated_Nj);
	void FillInputs(double** n0, double** n1, double** n2, double** n3);
	
	void Write(TFile* file, string label = string(""));
	
	/////////////////////
	// Print  methods
	////////////////////
	void PrintInputs();
	void PrintResults();
	void PrintResults_Latex();
	void Print4Estimators();
	
	void CheckConsistency();
	bool CheckEstimation(int njets);
	bool CheckEstimation(bool);
	void CheckEstimationLinearity(const int &, int, double **, int, int);
	void CheckEstimationStability(double, int, double, int, double, int, double, int, bool);

	//Main Method to be executed ...
	void Estimation(bool Print, bool DoFillSummary);
	//void Estimation(int NofIterations, bool Print, bool DoFillSummary);

   	//Rescale to an other lumi ...
	void ReScale(float factor);

	//////////////////////////////////////
	// Access to the 4 estimators
	// to be used after
	//////////////////////////////////////
	double GetEffudscEstimated() const;
	double GetEffbEstimated() const;
	double GetNWEventsEstimated() const;
	double GetNttEventsEstimated() const;
	double GetNtotalEventsEstimated() const;

        //nbjets refer to the number of bjets: 0,1,2,3
	double GetNWEventsEstimatedNbjets(int nbjets) const;
	double GetNttEventsRatioEstimatedNbjets(int nbjets) const;
	double GetNttEventsEstimatedNbjets(int nbjets) const;
	double GetNWEventsRatioEstimatedNbjets(int nbjets) const;
	
	//Give the number of events - integral - per #jets - per #b-jets
	double NofEvents() const { double a = 0; for(int i=0;i<NofJetsBins;i++) a+=NofEvents(i); return a;}
	double NofEvents(int njets) const { return (njets<NofJetsBins? N0bjets_[njets]+N1bjets_[njets]+N2bjets_[njets]+N3bjets_[njets]: -9999.);}
	double NofEvents(int njets, int nbjets) const;
	double NofEventsNbjets(int nbjets) const;

	///////////////////////////////////////////////////////////////////
	// Additional getters
	///////////////////////////////////////////////////////////////////
	int GetNbOfDatasets() const { return NbOfDatasets;}
	vector<int> GetiDatasetsTTLike()const {return iDatasetsTTLike;}
	vector<int> GetiDatasetsWLike() const {return iDatasetsWLike;}
	double*  GetMultiJet_Est_N0() const {return MultiJet_Est_N0_;}
	double*  GetMultiJet_Est_N1() const {return MultiJet_Est_N1_;}
	double*  GetMultiJet_Est_N2() const {return MultiJet_Est_N2_;}
	double*  GetMultiJet_Est_N3() const {return MultiJet_Est_N3_;}
	double** GetN0() const {return N0_;}
	double** GetN1() const {return N1_;}
	double** GetN2() const {return N2_;}
	double** GetN3() const {return N3_;}	
	double   GetEffudscEstimated(int njets) const;
	double   GetEffbEstimated(int njets) const;
	double   GetNWEventsEstimated(int njets) const;
	double   GetNttEventsEstimated(int njets) const;
	double   GetNtotalEventsEstimated(int njets) const;
	double   GetNTTlikeMC() const;
	double   GetNTTlikeMC(int njets) const;
	double   GetNWlikeMC() const;
	double   GetNWlikeMC(int njets) const;

	// Statistical and systematic erros (retrieved from WJetsPseudoExp) 
	void SetStatisticalWLikeError( int njets, int nbjets, double error) {statisticalWLikeError_[njets][nbjets] = error;}
	void SetStatisticalTTLikeError(int njets, int nbjets, double error) {statisticalTTLikeError_[njets][nbjets] = error;}

	void SetSystematicWLikeError( int njets, int nbjets, double error)  {systematicWLikeError_[njets][nbjets]  = error;}
	void SetSystematicTTLikeError(int njets, int nbjets, double error)  {systematicTTLikeError_[njets][nbjets]  = error;}

	double GetStatisticalWLikeError()  const;
	double GetSystematicWLikeError()   const;
	double GetStatisticalTTLikeError() const;
	double GetSystematicTTLikeError()  const;

	double GetStatisticalWLikeError(int njets)  const {return statisticalWLikeError_[njets][4];}
	double GetSystematicWLikeError( int njets)  const {return systematicWLikeError_[njets][4];}
	double GetStatisticalTTLikeError(int njets) const {return statisticalTTLikeError_[njets][4];}
	double GetSystematicTTLikeError( int njets) const {return systematicTTLikeError_[njets][4];}

	double GetStatisticalWLikeError(int njets, int nbjets)  const {return statisticalWLikeError_[njets][nbjets];}
	double GetSystematicWLikeError( int njets, int nbjets)  const {return systematicWLikeError_[njets][nbjets];}
	double GetStatisticalTTLikeError(int njets, int nbjets) const {return statisticalTTLikeError_[njets][nbjets];}
	double GetSystematicTTLikeError( int njets, int nbjets) const {return systematicTTLikeError_[njets][nbjets];}

   friend class WJetEstPseudoExp;

               	// private functions
   private:
	
	void Initialisation(bool &Print);
	void Iteration(int &iterator, bool &Print);
	//Fill Summary histograms :
	void FillSummaryHistos();
	//Important:
	//in the following methods, njets means the iterator on the NofJets contained in the Njets_ array and not this number itself :
	// 0-> 4jets 1-> 5jets 2->6jets or more
	void PrintInputs(int njets);
	void PrintResults(int njets);
	void Print4Estimators(int njets);
	void CheckConsistency(int njets);
	double EstimationChiSq();
	double EstimationChiSq(int njets);
	//methods called by Estimation()
	void Initialisation(bool &Print, int njets);
	void Iteration(int &It, bool &Print, int njets);
	//access to estimator
	double GetNWEventsEstimatedNbjets(int nbjets, int njets) const;
	double GetNttEventsEstimatedNbjets(int nbjets, int njets) const;

	///////////////////////////////////////////////////////////////////	
	// Methods to calculates the estimators ...  (from equations)
	///////////////////////////////////////////////////////////////////
	double N0bjet   (double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const;
	double Ntt_0bjet(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const;
	double Nw_0bjet (double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const;
	
	double N1bjet   (double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const;
	double Ntt_1bjet(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const;
	double Nw_1bjet (double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const;
	
	double N2bjets   (double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const;
	double Ntt_2bjets(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const;
	double Nw_2bjets (double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const;
	
	double N3bjets   (double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const;
	double Ntt_3bjets(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const;
	double Nw_3bjets (double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const;
	
	double eb_fromN2bjets(double &N2bjets, double &Ntt, double &Nw, double &eb_old, double &eudsc, int &n) const;
	vector<double> eb_cand_fromN3bjets(double &N3bjets, double &Ntt, double &Nw, double &eudsc, int &n) const;
	vector<double> eb_cand_fromN2bjets(double &N2bjets, double &Ntt, double &Nw, double &eudsc, int &n) const;
	vector<double> eb_cand_fromN1bjet (double &N1bjet,  double &Ntt, double &Nw, double &eudsc, int &n) const;
	vector<double> eb_cand_fromN0bjet (double &N0bjet,  double &Ntt, double &Nw, double &eudsc, int &n) const;
	
	double eudsc_fromN3bjets(double &N3bjets, double &Ntt, double &Nw, double &eb, double &eudsc_old, int &n) const;
	vector<double> eudsc_cand_fromN3bjets(double &N3bjets, double &Ntt, double &Nw, double &eb, int &n) const;
	vector<double> eudsc_cand_fromN2bjets(double &N2bjets, double &Ntt, double &Nw, double &eb, int &n) const;
	vector<double> eudsc_cand_fromN1bjet (double &N1bjet,  double &Ntt, double &Nw, double &eb, int &n) const;
	vector<double> eudsc_cand_fromN0bjet (double &N0bjet,  double &Ntt, double &Nw, double &eb, int &n) const;
	
	double Ntt_fromN2bjets(double &N2bjets, double &Nw, double &eb, double &eudsc, int &n) const;
	double Ntt_fromN1bjet( double &N1bjet,  double &Nw, double &eb, double &eudsc, int &n) const;
	double Nw_fromNtotal(  double &Ntotal,  double &Ntt) const;
	double Nw_fromN0bjet(  double &N0bjet,  double &Ntt, double &eb, double &eudsc, int &n) const;
	
	double GetGlobalEffudscEstimate(double&, double&, double&, double&, double&, double&, double&, double&, int&, bool&) const;
	double GetGlobalEffbEstimate(   double&, double&, double&, double&, double&, double&, double&, double&, int&, bool&) const;

   private:
   	bool MCdata_;
	int NbOfDatasets;
	int It_;
	TH1F*** histos;
	TH1F*** hNjetsEstSummary;
	TH1F*** hNbjetsEstSummary;
	TH1F*** hNjetsMCSummary;
	TH1F*** hNbjetsMCSummary;
	TH2D**  EstimationChiSq_Eb_vs_Eudsc;
	TH2D**  EstimationChiSq_Ntt_vs_Nw;
	TLegend *MyLeg;

	THStack** hsNjets_MC;
	THStack** hsNjets_Est;
	THStack** hsNbjets_MC;
	THStack** hsNbjets_Est;
	
	TCanvas **tCanva_Njets_Summary;
	TCanvas **tCanva_Nbjets_Summary;

	TGraphErrors *RescaledTTlikeEstimation;
	TGraphErrors *RescaledWlikeEstimation;
	
	TCanvas *tCanva_RescaledTTlikeEstimation;
	TCanvas *tCanva_RescaledWlikeEstimation;
	
	vector<int> iDatasetsTTLike;
	vector<int> iDatasetsWLike;

	int  NofJetsBins;
	//all the following members are array of jet multiplicities
	int* Njets_;
	//initial number of b-jets 
	double* N0bjets_; 
	double* N1bjets_;
	double* N2bjets_;
	double* N3bjets_; 
	//initial number of b-jets per datasets (NbOfDatasets)
	// Multi-jets/tt+jets/W+jets/Z+jets/st-s/st-t/st-tW
	double** N0_;
	double** N1_;
	double** N2_;
	double** N3_;
        //Number of events estimated for multi-jets QCD process	
	double* MultiJet_Est_N0_;
	double* MultiJet_Est_N1_;
	double* MultiJet_Est_N2_;
	double* MultiJet_Est_N3_;
	
	//The 4 estimators
	double* eudsc_;
        double* eb_;
        double* Ntt_ ;
	double* Nw_;

	//Statistical and systematic errors (to be retrieved from WJetsEstPseudoExp)
	double **statisticalWLikeError_; // row : nb of jets (from 4 to >6), column : nb of b-jets (from 0 to 3, last column=inclusive)
	double **statisticalTTLikeError_;// row : nb of jets (from 4 to >6), column : nb of b-jets (from 0 to 3, last column=inclusive)
	double **systematicWLikeError_;  // row : nb of jets (from 4 to >6), column : nb of b-jets (from 0 to 3, last column=inclusive)
	double **systematicTTLikeError_; // row : nb of jets (from 4 to >6), column : nb of b-jets (from 0 to 3, last column=inclusive)

};

#endif

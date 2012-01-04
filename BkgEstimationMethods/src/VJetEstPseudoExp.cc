#include "../interface/VJetEstPseudoExp.h"

VJetEstPseudoExp::VJetEstPseudoExp(int NbOfPseudoExp, VJetEstimation *vjEst){
	cout<<"Object from the class VJetEstPseudoExp being instantiated"<<endl;
	NbOfPseudoExp_    = NbOfPseudoExp;
	vjEst_            = vjEst;
	NbOfDataSets_     = vjEst_->GetNbOfDatasets();
	NbOfJetsBins_     = vjEst_->GetNbOfJetsBins();
	NbOfBJetsBins_    = vjEst_->GetNbOfBJetsBins();
	NbOfBtagWorkingPoint_ = vjEst_->GetNbOfBtagWorkingPoint();

	NbOfSuccessfulPE_ = 0;
	NbOfRealPE_= 0;

	rand_ = new ROOT::Math::Random<ROOT::Math::GSLRngMT>();
	condProb_ = 0;
	
	/// Initial numbers of events for each b-jet multiplicity
	/// -- per b-tagging working point, jet multiplicity and per dataset
	initNbjets_ = new double***[NbOfBtagWorkingPoint_];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		initNbjets_[i] = new double**[NbOfJetsBins_];//vjEst_->GetN0(0);
		for(int j=0; j<NbOfJetsBins_; j++){
			initNbjets_[i][j] = new double*[NbOfBJetsBins_];
			for(int k=0; k<NbOfBJetsBins_; k++){
			initNbjets_[i][j][k] = new double[NbOfDataSets_];
				for(int l=0; l<NbOfDataSets_; l++){
					initNbjets_[i][j][k][l] = 0;//(double)vjEst_->GetN(i)[j][k][l];
				}
			}
		}
	}

	/// Randomized numbers of events with 0,1,2 and 3 b-jets
	/// -- per b-tagging working points, per jet multiplicity and per dataset
	randN0bjets_ = new double***[NbOfBtagWorkingPoint_];
	randN1bjets_ = new double***[NbOfBtagWorkingPoint_];
	randN2bjets_ = new double***[NbOfBtagWorkingPoint_];
	randN3bjets_ = new double***[NbOfBtagWorkingPoint_];
	for(int j=0;j<NbOfBtagWorkingPoint_;j++){
		randN0bjets_[j] = new double**[NbOfJetsBins_];
		randN1bjets_[j] = new double**[NbOfJetsBins_];
		randN2bjets_[j] = new double**[NbOfJetsBins_];
		randN3bjets_[j] = new double**[NbOfJetsBins_];
		for(int k=0;k<NbOfJetsBins_;k++){
			randN0bjets_[j][k] = new double*[NbOfDataSets_];
			randN1bjets_[j][k] = new double*[NbOfDataSets_];
			randN2bjets_[j][k] = new double*[NbOfDataSets_];
			randN3bjets_[j][k] = new double*[NbOfDataSets_];
			for(int l=0;l<NbOfDataSets_;l++){
				randN0bjets_[j][k][l] = new double[NbOfPseudoExp];
				randN1bjets_[j][k][l] = new double[NbOfPseudoExp];
				randN2bjets_[j][k][l] = new double[NbOfPseudoExp];
				randN3bjets_[j][k][l] = new double[NbOfPseudoExp];
				for(int i=0; i<NbOfPseudoExp; i++){
					randN0bjets_[j][k][l][i] = 0;
					randN1bjets_[j][k][l][i] = 0;
					randN2bjets_[j][k][l][i] = 0;
					randN3bjets_[j][k][l][i] = 0;
				}
			}
		}
	}

	/// Histograms for the estimation parameters 
	/// -- per b-tagging working point and jet multiplicity
	hEudsc_ = new TH1F**[NbOfBtagWorkingPoint_];
	hEuds_  = new TH1F**[NbOfBtagWorkingPoint_];
	hEb_    = new TH1F**[NbOfBtagWorkingPoint_];

	hPullEudsc_ = new TH1F**[NbOfBtagWorkingPoint_];
	hPullEuds_  = new TH1F**[NbOfBtagWorkingPoint_];
	hPullEb_    = new TH1F**[NbOfBtagWorkingPoint_];

	hDiffEudsc_ = new TH1F**[NbOfBtagWorkingPoint_];
	hDiffEuds_  = new TH1F**[NbOfBtagWorkingPoint_];
	hDiffEb_    = new TH1F**[NbOfBtagWorkingPoint_];

	hNtt_       = new TH1F***[NbOfBtagWorkingPoint_];
	hNvb_       = new TH1F***[NbOfBtagWorkingPoint_];
	hNv_        = new TH1F***[NbOfBtagWorkingPoint_];

	hPullNtt_            = new TH1F***[NbOfBtagWorkingPoint_];
	hPullNvb_            = new TH1F***[NbOfBtagWorkingPoint_];
	hPullNv_             = new TH1F***[NbOfBtagWorkingPoint_];

	hDiffNtt_            = new TH1F***[NbOfBtagWorkingPoint_];
	hDiffNvb_            = new TH1F***[NbOfBtagWorkingPoint_];
	hDiffNv_             = new TH1F***[NbOfBtagWorkingPoint_];

	hDiffNtt_vs_MinLL_   = new TH2F***[NbOfBtagWorkingPoint_];
	hDiffNvb_vs_MinLL_   = new TH2F***[NbOfBtagWorkingPoint_];
	hDiffNv_vs_MinLL_    = new TH2F***[NbOfBtagWorkingPoint_];
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		hEudsc_[i] = new TH1F*[NbOfJetsBins_];
		hEuds_[i]  = new TH1F*[NbOfJetsBins_];
		hEb_[i]    = new TH1F*[NbOfJetsBins_];
		hPullEudsc_[i] = new TH1F*[NbOfJetsBins_];
		hPullEuds_[i]  = new TH1F*[NbOfJetsBins_];
		hPullEb_[i]    = new TH1F*[NbOfJetsBins_];
		hDiffEudsc_[i] = new TH1F*[NbOfJetsBins_];
		hDiffEuds_[i]  = new TH1F*[NbOfJetsBins_];
		hDiffEb_[i]    = new TH1F*[NbOfJetsBins_];
		/// -- per jet multiplicity
		hNtt_[i]     = new TH1F**[NbOfJetsBins_+1]; /// the last histogram is (jet) inclusive 
		hNvb_[i]     = new TH1F**[NbOfJetsBins_+1];
		hNv_[i]      = new TH1F**[NbOfJetsBins_+1];
		hPullNtt_[i] = new TH1F**[NbOfJetsBins_+1]; /// the last histogram is (jet) inclusive 
		hPullNvb_[i] = new TH1F**[NbOfJetsBins_+1];
		hPullNv_[i]  = new TH1F**[NbOfJetsBins_+1];
		hDiffNtt_[i] = new TH1F**[NbOfJetsBins_+1]; /// the last histogram is (jet) inclusive 
		hDiffNvb_[i] = new TH1F**[NbOfJetsBins_+1];
		hDiffNv_[i]  = new TH1F**[NbOfJetsBins_+1];
		hDiffNtt_vs_MinLL_[i] = new TH2F**[NbOfJetsBins_+1]; /// the last histogram is (jet) inclusive 
		hDiffNvb_vs_MinLL_[i] = new TH2F**[NbOfJetsBins_+1];
		hDiffNv_vs_MinLL_[i]  = new TH2F**[NbOfJetsBins_+1];
		for(int j=0;j<=NbOfJetsBins_;j++){
			hEb_[i][j]    = 0;
			hEudsc_[i][j] = 0;
			hEuds_[i][j]  = 0;
			hPullEb_[i][j]    = 0;
			hPullEudsc_[i][j] = 0;
			hPullEuds_[i][j]  = 0;
			hDiffEb_[i][j]    = 0;
			hDiffEudsc_[i][j] = 0;
			hDiffEuds_[i][j]  = 0;
			hNtt_[i][j] = new TH1F*[NbOfBJetsBins_+1]; /// the last histogram is (b-jet) inclusive 
			hNvb_[i][j] = new TH1F*[NbOfBJetsBins_+1];
			hNv_[i][j]  = new TH1F*[NbOfBJetsBins_+1];
			hPullNtt_[i][j] = new TH1F*[NbOfBJetsBins_+1]; /// the last histogram is (b-jet) inclusive 
			hPullNvb_[i][j] = new TH1F*[NbOfBJetsBins_+1];
			hPullNv_[i][j]  = new TH1F*[NbOfBJetsBins_+1];
			hDiffNtt_[i][j] = new TH1F*[NbOfBJetsBins_+1]; /// the last histogram is (b-jet) inclusive 
			hDiffNvb_[i][j] = new TH1F*[NbOfBJetsBins_+1];
			hDiffNv_[i][j]  = new TH1F*[NbOfBJetsBins_+1];
			hDiffNtt_vs_MinLL_[i][j] = new TH2F*[NbOfBJetsBins_+1]; /// the last histogram is (b-jet) inclusive 
			hDiffNvb_vs_MinLL_[i][j] = new TH2F*[NbOfBJetsBins_+1];
			hDiffNv_vs_MinLL_[i][j]  = new TH2F*[NbOfBJetsBins_+1];
			for(int k=0;k<=NbOfBJetsBins_;k++){
				hNtt_[i][j][k] = 0;
				hNvb_[i][j][k] = 0;
				hNv_[i][j][k]  = 0;
				hPullNtt_[i][j][k] = 0;
				hPullNvb_[i][j][k] = 0;
				hPullNv_[i][j][k]  = 0;
				hDiffNtt_[i][j][k] = 0;
				hDiffNvb_[i][j][k] = 0;
				hDiffNv_[i][j][k]  = 0;
				hDiffNtt_vs_MinLL_[i][j][k] = 0;
				hDiffNvb_vs_MinLL_[i][j][k] = 0;
				hDiffNv_vs_MinLL_[i][j][k]  = 0;
			}
		}
	}

	/* Histograms for the randomized numbers of events with 0,1,2 and 3 b-jets
	-- per b-tagging working points, per jet multiplicity and per dataset */
	hRandN0bjets_ = new TH1F***[NbOfBtagWorkingPoint_];
	hRandN1bjets_ = new TH1F***[NbOfBtagWorkingPoint_];
	hRandN2bjets_ = new TH1F***[NbOfBtagWorkingPoint_];
	hRandN3bjets_ = new TH1F***[NbOfBtagWorkingPoint_];

	/*Histograms for the randomized numbers of events with 0,1,2 and 3 b-jets vs the min value of the LogL function */
	hRandNbjets_vs_MinLL_ = new TH3F***[NbOfBtagWorkingPoint_];
	
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		hRandN0bjets_[i] = new TH1F**[NbOfJetsBins_];
		hRandN1bjets_[i] = new TH1F**[NbOfJetsBins_];
		hRandN2bjets_[i] = new TH1F**[NbOfJetsBins_];
		hRandN3bjets_[i] = new TH1F**[NbOfJetsBins_];
		hRandNbjets_vs_MinLL_[i] = new TH3F**[NbOfJetsBins_];
		for(int j=0;j<NbOfJetsBins_;j++)
		{
			hRandN0bjets_[i][j] = new TH1F*[NbOfDataSets_];
			hRandN1bjets_[i][j] = new TH1F*[NbOfDataSets_];
			hRandN2bjets_[i][j] = new TH1F*[NbOfDataSets_];
			hRandN3bjets_[i][j] = new TH1F*[NbOfDataSets_];
			hRandNbjets_vs_MinLL_[i][j] = new TH3F*[NbOfDataSets_];
				for(int k=0;k<NbOfDataSets_;k++){
					hRandN0bjets_[i][j][k] = 0;
					hRandN1bjets_[i][j][k] = 0;
					hRandN2bjets_[i][j][k] = 0;
					hRandN3bjets_[i][j][k] = 0;
					hRandNbjets_vs_MinLL_[i][j][k] = 0;
				}
		}
	}
	cout<<"Object from the class VJetEstPseudoExp correctly instantiated"<<endl;
}

/**________________________________________________________________________________________________________________*/
VJetEstPseudoExp::~VJetEstPseudoExp(){
	delete rand_;

	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		for(int j=0;j<NbOfJetsBins_;j++){
			for(int k=0;k<NbOfBJetsBins_; k++) delete [] initNbjets_[i][j][k];
			for(int k=0;k<NbOfDataSets_;k++){
				delete [] randN0bjets_[i][j][k];
				delete [] randN1bjets_[i][j][k];
				delete [] randN2bjets_[i][j][k];
				delete [] randN3bjets_[i][j][k];
			}
			delete [] initNbjets_[i][j];
			delete [] randN0bjets_[i][j];
			delete [] randN1bjets_[i][j];
			delete [] randN2bjets_[i][j];
			delete [] randN3bjets_[i][j];
		}
		delete [] initNbjets_[i];
		delete [] randN0bjets_[i];
		delete [] randN1bjets_[i];
		delete [] randN2bjets_[i];
		delete [] randN3bjets_[i];
	}
	delete [] initNbjets_;
	delete [] randN0bjets_;
	delete [] randN1bjets_;
	delete [] randN2bjets_;
	delete [] randN3bjets_;

	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		for(int j=0;j<NbOfJetsBins_;j++){
			delete hEudsc_[i][j];
			delete hEuds_[i][j];
			delete hEb_[i][j];
			delete hPullEudsc_[i][j];
			delete hPullEuds_[i][j];
			delete hPullEb_[i][j];
			delete hDiffEudsc_[i][j];
			delete hDiffEuds_[i][j];
			delete hDiffEb_[i][j];
		}
		delete [] hEudsc_[i];
		delete [] hEuds_[i];
		delete [] hEb_[i];
		delete [] hPullEudsc_[i];
		delete [] hPullEuds_[i];
		delete [] hPullEb_[i];
		delete [] hDiffEudsc_[i];
		delete [] hDiffEuds_[i];
		delete [] hDiffEb_[i];
		for(int j=0;j<=NbOfJetsBins_;j++){
			for(int k=0;k<=NbOfBJetsBins_;k++){
				delete hNtt_[i][j][k];
				delete hNvb_[i][j][k];
				delete hNv_[i][j][k] ;
				delete hPullNtt_[i][j][k];
				delete hPullNvb_[i][j][k];
				delete hPullNv_[i][j][k] ;
				delete hDiffNtt_[i][j][k];
				delete hDiffNvb_[i][j][k];
				delete hDiffNv_[i][j][k] ;
				delete hDiffNtt_vs_MinLL_[i][j][k];
				delete hDiffNvb_vs_MinLL_[i][j][k];
				delete hDiffNv_vs_MinLL_[i][j][k] ;
			}
			delete [] hNtt_[i][j];
			delete [] hNvb_[i][j];
			delete [] hNv_[i][j];
			delete [] hPullNtt_[i][j];
			delete [] hPullNvb_[i][j];
			delete [] hPullNv_[i][j];
			delete [] hDiffNtt_[i][j];
			delete [] hDiffNvb_[i][j];
			delete [] hDiffNv_[i][j];
			delete [] hDiffNtt_vs_MinLL_[i][j];
			delete [] hDiffNvb_vs_MinLL_[i][j];
			delete [] hDiffNv_vs_MinLL_[i][j];
		}
		delete [] hNtt_[i];
		delete [] hNvb_[i];
		delete [] hNv_[i];
		delete [] hPullNtt_[i];
		delete [] hPullNvb_[i];
		delete [] hPullNv_[i];
		delete [] hDiffNtt_[i];
		delete [] hDiffNvb_[i];
		delete [] hDiffNv_[i];
		delete [] hDiffNtt_vs_MinLL_[i];
		delete [] hDiffNvb_vs_MinLL_[i];
		delete [] hDiffNv_vs_MinLL_[i];
	}
	delete [] hEudsc_;
	delete [] hEuds_;
	delete [] hEb_;
	delete [] hPullEudsc_;
	delete [] hPullEuds_;
	delete [] hPullEb_;
	delete [] hDiffEudsc_;
	delete [] hDiffEuds_;
	delete [] hDiffEb_;

	delete [] hNtt_;
	delete [] hNvb_;
	delete [] hNv_ ;
	delete [] hDiffNtt_;
	delete [] hDiffNvb_;
	delete [] hDiffNv_ ;
	delete [] hPullNtt_;
	delete [] hPullNvb_;
	delete [] hPullNv_ ;
	delete [] hDiffNtt_vs_MinLL_;
	delete [] hDiffNvb_vs_MinLL_;
	delete [] hDiffNv_vs_MinLL_;

	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		for(int j=0;j<NbOfJetsBins_;j++)
		{
			for(int k=0;k<NbOfDataSets_;k++){
				delete hRandN0bjets_[i][j][k];
				delete hRandN1bjets_[i][j][k];
				delete hRandN2bjets_[i][j][k];
				delete hRandN3bjets_[i][j][k];
				delete hRandNbjets_vs_MinLL_[i][j][k];
			}
			delete [] hRandN0bjets_[i][j];
			delete [] hRandN1bjets_[i][j];
			delete [] hRandN2bjets_[i][j];
			delete [] hRandN3bjets_[i][j];
			delete [] hRandNbjets_vs_MinLL_[i][j];
		}
		delete [] hRandN0bjets_[i];
		delete [] hRandN1bjets_[i];
		delete [] hRandN2bjets_[i];
		delete [] hRandN3bjets_[i];
		delete [] hRandNbjets_vs_MinLL_[i];
	}
	delete [] hRandN0bjets_;
	delete [] hRandN1bjets_;
	delete [] hRandN2bjets_;
	delete [] hRandN3bjets_;
	delete [] hRandNbjets_vs_MinLL_;
}

/**________________________________________________________________________________________________________________*/
void VJetEstPseudoExp::PrintResults(){

	cout<<" --------------------------------------------------"<<endl;
	cout<<" ------------ W/Z background estimation -----------"<<endl;
	cout<<" -- Required number of pseudo-experiments : "<<NbOfPseudoExp_<<endl;
	cout<<" -- Number of pseudo-experiments that were actually thrown : "<<NbOfRealPE_<<endl;
	cout<<" -- Number of successful P.E. : "<<NbOfSuccessfulPE_<<endl;
	cout<<" -- Statistical errors calculated from the pseudo-experience :"<<endl;
	
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		cout<<"For the b-tagging working point (nr "<<i<<")"<<endl;
		for(int j = 0; j<NbOfJetsBins_; j++){
		cout<<" --- estimates with "<<vjEst_->GetNjets(j)<<" jets :"<<endl;
		cout<<" ---- V-like estimate :"<<endl;
		cout<<" ----- stat. error (Variance) : "<<GetNvVar(i,j,false,0,0)<<endl;
		cout<<" ----- syst. error (Mean)     : "<<GetNvBias(i,j,false,0,0)<<endl;
		cout<<" ---- TT-like estimate :"<<endl;
		cout<<" ----- stat. error (Variance) : "<<GetNttVar(i,j,false,0,0)<<endl;
		cout<<" ----- syst. error (Mean)     : "<<GetNttBias(i,j,false,0,0)<<endl;
		cout<<" ----------------------"<<endl;
		}

		cout<<" --- inclusive estimates :"<<endl;
		cout<<" ---- V-like estimate :"<<endl;
		cout<<" ----- stat. error (Variance) : "<<GetNvVar(i,false,0,0)<<endl;
		cout<<" ----- syst. error (Mean)     : "<<GetNvBias(i,false,0,0)<<endl;
		cout<<" ---- TT-like estimate :"<<endl;
		cout<<" ----- stat. error (Variance) : "<<GetNttVar(i,false,0,0)<<endl;
		cout<<" ----- syst. error (Mean)     : "<<GetNttBias(i,false,0,0)<<endl;
		cout<<" --------------------------------------------------"<<endl;
	}
}

/**________________________________________________________________________________________________________________*/
void VJetEstPseudoExp::PrintResults_latex(ofstream &ofile){

	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		ofile<<"\\begin{table}"<<endl;
		ofile<<"	\\centering"<<endl;
		ofile<<"	 \\begin{tabular}{|l||c|c|c|c|}"<<endl;
		ofile<<"	\\hline"<<endl;
		ofile<<"\\multicolumn{5}{|l}{V-like estimate : }\\\\"<<endl;
		ofile<<"&";for(int j=0;j<NbOfJetsBins_;j++){ofile<<"$"<<vjEst_->GetNjets(j)<<"$ jets & ";}ofile<<"Inclusive \\\\"<<endl;
		//ofile<<" & $4$ jets & $5$ jets & $6$ jets & inclusive \\\\"<<endl;
		ofile<<"Statistical error &$";for(int j=0;j<vjEst_->GetNbOfBtagWorkingPoint();j++){ofile<<GetNvVar( i,j,false,0,0)<<"&";}ofile<<GetNvVar( i,false,0,0)<<"\\\\"<<endl;
		ofile<<"Systematic error  &$";for(int j=0;j<vjEst_->GetNbOfBtagWorkingPoint();j++){ofile<<GetNvBias(i,j,false,0,0)<<"&";}ofile<<GetNvBias(i,false,0,0)<<"\\\\"<<endl;
		ofile<<"	\\hline"<<endl;
		ofile<<"\\multicolumn{5}{|l}{$\\ttbar$-like estimate : }\\\\"<<endl;
		ofile<<"&";for(int j=0;j<NbOfJetsBins_;j++){ofile<<"$"<<vjEst_->GetNjets(j)<<"$ jets & ";}ofile<<"Inclusive \\\\"<<endl;
		//ofile<<" & $4$ jets & $5$ jets & $6$ jets & inclusive \\\\"<<endl;
		ofile<<"Statistical error &$";for(int j=0;j<vjEst_->GetNbOfBtagWorkingPoint();j++){ofile<<GetNttVar( i,j,false,0,0)<<"&";}ofile<<GetNttVar( i,false,0,0)<<"\\\\"<<endl;
		ofile<<"Systematic error  &$";for(int j=0;j<vjEst_->GetNbOfBtagWorkingPoint();j++){ofile<<GetNttBias(i,j,false,0,0)<<"&";}ofile<<GetNttBias(i,false,0,0)<<"\\\\"<<endl;
		ofile<<"	\\hline"<<endl;
		ofile<<"  \\end{tabular}"<<endl;
		ofile<<"  \\caption{}"<<endl;
		ofile<<"  \\label{tab:}"<<endl;
		ofile<<"\\end{table}"<<endl;	
	}
}

/**________________________________________________________________________________________________________________*/
void VJetEstPseudoExp::RollTheDice(string &method, string &option, vector<int> SetFixedVarIdx, bool doMinos, bool UseUnBinMLE, bool UseMJLE, bool Print){
	//bool   Verbose = false;
	condProb_ = vjEst_->GetCondProb();
	if(condProb_ == 0) {cout<<"Can't perform the randomization. Please use the method ComputeCondProb of the VJetEstimation class..."<<endl; return;}
	
	NbOfSuccessfulPE_ = 0;
	RooMsgService::instance().Print();

	double ***PseudoExp_EstEb    = new double**[NbOfBtagWorkingPoint_];
	double ***PseudoExp_EstEudsc = new double**[NbOfBtagWorkingPoint_];
	double ***PseudoExp_EstEuds  = new double**[NbOfBtagWorkingPoint_];

	double ***PseudoExp_PullEstEb    = new double**[NbOfBtagWorkingPoint_];
	double ***PseudoExp_PullEstEudsc = new double**[NbOfBtagWorkingPoint_];
	double ***PseudoExp_PullEstEuds  = new double**[NbOfBtagWorkingPoint_];

	double ***PseudoExp_DiffEstEb    = new double**[NbOfBtagWorkingPoint_];
	double ***PseudoExp_DiffEstEudsc = new double**[NbOfBtagWorkingPoint_];
	double ***PseudoExp_DiffEstEuds  = new double**[NbOfBtagWorkingPoint_];

	double ****PseudoExp_EstNvlike  = new double***[NbOfBtagWorkingPoint_];
	double ****PseudoExp_EstNvblike = new double***[NbOfBtagWorkingPoint_];
	double ****PseudoExp_EstNttlike = new double***[NbOfBtagWorkingPoint_];

	double ****PseudoExp_PullEstNvlike  = new double***[NbOfBtagWorkingPoint_];
	double ****PseudoExp_PullEstNvblike = new double***[NbOfBtagWorkingPoint_];
	double ****PseudoExp_PullEstNttlike = new double***[NbOfBtagWorkingPoint_];

	double ****PseudoExp_DiffEstNvlike  = new double***[NbOfBtagWorkingPoint_];
	double ****PseudoExp_DiffEstNvblike = new double***[NbOfBtagWorkingPoint_];
	double ****PseudoExp_DiffEstNttlike = new double***[NbOfBtagWorkingPoint_];

	double ****tempNbjets = new double***[NbOfBtagWorkingPoint_];
	double   **minLL      = new double*[NbOfJetsBins_];
	for(int i=0;i<NbOfJetsBins_;i++) minLL[i] = new double[NbOfPseudoExp_];

	for(int i=0;i<NbOfBtagWorkingPoint_;i++){ 
		PseudoExp_EstNvlike[i]  = new double**[NbOfJetsBins_+1];
		PseudoExp_EstNvblike[i] = new double**[NbOfJetsBins_+1];
		PseudoExp_EstNttlike[i] = new double**[NbOfJetsBins_+1];
		PseudoExp_PullEstNvlike[i]  = new double**[NbOfJetsBins_+1];
		PseudoExp_PullEstNvblike[i] = new double**[NbOfJetsBins_+1];
		PseudoExp_PullEstNttlike[i] = new double**[NbOfJetsBins_+1];
		PseudoExp_DiffEstNvlike[i]  = new double**[NbOfJetsBins_+1];
		PseudoExp_DiffEstNvblike[i] = new double**[NbOfJetsBins_+1];
		PseudoExp_DiffEstNttlike[i] = new double**[NbOfJetsBins_+1];
		for(int j=0;j<NbOfJetsBins_+1;j++){
			PseudoExp_EstNvlike[i][j]   = new double*[NbOfBJetsBins_+1];
			PseudoExp_EstNvblike[i][j]  = new double*[NbOfBJetsBins_+1];
			PseudoExp_EstNttlike[i][j]  = new double*[NbOfBJetsBins_+1];
			PseudoExp_PullEstNvlike[i][j]   = new double*[NbOfBJetsBins_+1];
			PseudoExp_PullEstNvblike[i][j]  = new double*[NbOfBJetsBins_+1];
			PseudoExp_PullEstNttlike[i][j]  = new double*[NbOfBJetsBins_+1];
			PseudoExp_DiffEstNvlike[i][j]   = new double*[NbOfBJetsBins_+1];
			PseudoExp_DiffEstNvblike[i][j]  = new double*[NbOfBJetsBins_+1];
			PseudoExp_DiffEstNttlike[i][j]  = new double*[NbOfBJetsBins_+1];
			for(int k=0;k<NbOfBJetsBins_+1;k++){
				PseudoExp_EstNvlike[i][j][k]   = new double[NbOfPseudoExp_];
				PseudoExp_EstNvblike[i][j][k]  = new double[NbOfPseudoExp_];
				PseudoExp_EstNttlike[i][j][k]  = new double[NbOfPseudoExp_];
				PseudoExp_PullEstNvlike[i][j][k]   = new double[NbOfPseudoExp_];
				PseudoExp_PullEstNvblike[i][j][k]  = new double[NbOfPseudoExp_];
				PseudoExp_PullEstNttlike[i][j][k]  = new double[NbOfPseudoExp_];
				PseudoExp_DiffEstNvlike[i][j][k]   = new double[NbOfPseudoExp_];
				PseudoExp_DiffEstNvblike[i][j][k]  = new double[NbOfPseudoExp_];
				PseudoExp_DiffEstNttlike[i][j][k]  = new double[NbOfPseudoExp_];
			}
		}
		PseudoExp_EstEb[i]    = new double*[NbOfJetsBins_];
		PseudoExp_EstEudsc[i] = new double*[NbOfJetsBins_];
		PseudoExp_EstEuds[i]  = new double*[NbOfJetsBins_];
		PseudoExp_PullEstEb[i]    = new double*[NbOfJetsBins_];
		PseudoExp_PullEstEudsc[i] = new double*[NbOfJetsBins_];
		PseudoExp_PullEstEuds[i]  = new double*[NbOfJetsBins_];
		PseudoExp_DiffEstEb[i]    = new double*[NbOfJetsBins_];
		PseudoExp_DiffEstEudsc[i] = new double*[NbOfJetsBins_];
		PseudoExp_DiffEstEuds[i]  = new double*[NbOfJetsBins_];

		tempNbjets[i] = new double**[NbOfJetsBins_];
		for(int j=0;j<NbOfJetsBins_;j++){
			PseudoExp_EstEb[i][j]    = new double[NbOfPseudoExp_];
			PseudoExp_EstEudsc[i][j] = new double[NbOfPseudoExp_];
			PseudoExp_EstEuds[i][j]  = new double[NbOfPseudoExp_];
			PseudoExp_PullEstEb[i][j]    = new double[NbOfPseudoExp_];
			PseudoExp_PullEstEudsc[i][j] = new double[NbOfPseudoExp_];
			PseudoExp_PullEstEuds[i][j]  = new double[NbOfPseudoExp_];
			PseudoExp_DiffEstEb[i][j]    = new double[NbOfPseudoExp_];
			PseudoExp_DiffEstEudsc[i][j] = new double[NbOfPseudoExp_];
			PseudoExp_DiffEstEuds[i][j]  = new double[NbOfPseudoExp_];

			tempNbjets[i][j] = new double*[NbOfBJetsBins_];

			for(int k=0;k<NbOfBJetsBins_;k++) tempNbjets[i][j][k] = new double[NbOfDataSets_];
		}
	}
	
	double minfact = 0.9;
	double maxfact = 1.1;
	double sum = 0;
	double InitEstNv[NbOfJetsBins_];
	double InitEstNtt[NbOfJetsBins_];
	for(int i=0;i<NbOfJetsBins_;i++){
		InitEstNv[i]  = vjEst_->GetPredNv( 0,i);
		InitEstNtt[i] = vjEst_->GetPredNtt(0,i);
	}

	//vector<unsigned int> outcomes[NbOfBtagWorkingPoint_-1][NbOfJetsBins_][NbOfDataSets_];
	vector<unsigned int> outcomes[NbOfBJetsBins_];
	for(int m=0;m<NbOfBJetsBins_;m++){outcomes[m].resize(4);};
	
	vector<double> prob[NbOfBtagWorkingPoint_-1][NbOfJetsBins_][NbOfDataSets_][NbOfBJetsBins_];
	for(int j=NbOfBtagWorkingPoint_-2;j>=0;j--){
		for(int k=0;k<NbOfJetsBins_;k++){
			for(int l=0;l<NbOfDataSets_;l++){
				for(int m=0;m<NbOfBJetsBins_;m++){
					prob[j][k][l][m].resize(NbOfBJetsBins_);
					//outcomes[j][k][l][m].resize(4);
					for(int n=0;n<NbOfBJetsBins_;n++) prob[j][k][l][m][n] = condProb_[j][k][l][m][n];
					if(Print) cout<<"Sum of the prob = "<<accumulate(prob[j][k][l][m].begin(),prob[j][k][l][m].begin()+NbOfBJetsBins_,0.0)<<endl;
					sum = accumulate(prob[j][k][l][m].begin(),prob[j][k][l][m].begin()+NbOfBJetsBins_,0.0);
					if( sum <0.9999 ){
						cout<<"Sum of the prob for(j,k,l,m) = ("<<j<<""<<k<<""<<l<<""<<m<<") is not equal to 1 ("<<sum<<")"<<endl;
						return;
					}
				}
			}
		}
	}

	if(Print) std::cout<<"// RollTheDice method initialized..."<<std::endl;

	std::cout<<"// Let's start the RollTheDice method..."<<std::endl;
	for(int i=0; i<NbOfPseudoExp_; i++){
		NbOfRealPE_++;
		if(i%500==0) std::cout<<"RollTheDice method : running the "<<i+1<<"th pseudo-exp"<<std::endl;
		/** As the b-jet multiplicity distributions, when working with more than one b-tagging working points, are correlated,
		special care is needed when randomizing	these distributions.
		The distribution corresponding to the tighest working point is randomized according to a Poisson distribution
		and this randomization is propagated to the others via multinomial distribution */

		for(int k=0;k<NbOfJetsBins_;k++){
	 		for(int l=0;l<NbOfDataSets_;l++){
				//cout<<"initNbjets_["<<NbOfBtagWorkingPoint_-1<<"]["<<k<<"]["<<0<<"]["<<l<<"] = "<<initNbjets_[NbOfBtagWorkingPoint_-1][k][0][l]<<endl;
				randN0bjets_[NbOfBtagWorkingPoint_-1][k][l][i] = rand_->Poisson(initNbjets_[NbOfBtagWorkingPoint_-1][k][0][l]);
				tempNbjets[NbOfBtagWorkingPoint_-1][k][0][l] = randN0bjets_[NbOfBtagWorkingPoint_-1][k][l][i];
				//cout<<"initNbjets_["<<NbOfBtagWorkingPoint_-1<<"]["<<k<<"]["<<1<<"]["<<l<<"] = "<<initNbjets_[NbOfBtagWorkingPoint_-1][k][1][l]<<endl;
				randN1bjets_[NbOfBtagWorkingPoint_-1][k][l][i] = rand_->Poisson(initNbjets_[NbOfBtagWorkingPoint_-1][k][1][l]);
				tempNbjets[NbOfBtagWorkingPoint_-1][k][1][l] = randN1bjets_[NbOfBtagWorkingPoint_-1][k][l][i];
				//cout<<"initNbjets_["<<NbOfBtagWorkingPoint_-1<<"]["<<k<<"]["<<2<<"]["<<l<<"] = "<<initNbjets_[NbOfBtagWorkingPoint_-1][k][2][l]<<endl;
				randN2bjets_[NbOfBtagWorkingPoint_-1][k][l][i] = rand_->Poisson(initNbjets_[NbOfBtagWorkingPoint_-1][k][2][l]);
				tempNbjets[NbOfBtagWorkingPoint_-1][k][2][l] = randN2bjets_[NbOfBtagWorkingPoint_-1][k][l][i];
				//cout<<"initNbjets_["<<NbOfBtagWorkingPoint_-1<<"]["<<k<<"]["<<3<<"]["<<l<<"] = "<<initNbjets_[NbOfBtagWorkingPoint_-1][k][3][l]<<endl;
				randN3bjets_[NbOfBtagWorkingPoint_-1][k][l][i] = rand_->Poisson(initNbjets_[NbOfBtagWorkingPoint_-1][k][3][l]);
				tempNbjets[NbOfBtagWorkingPoint_-1][k][3][l] = randN3bjets_[NbOfBtagWorkingPoint_-1][k][l][i];
			}
		}

		for(int j=NbOfBtagWorkingPoint_-2;j>=0;j--){
			for(int k=0;k<NbOfJetsBins_;k++){
	 			for(int l=0;l<NbOfDataSets_;l++){
					outcomes[0] = rand_->Multinomial(randN0bjets_[j+1][k][l][i],prob[j][k][l][0]);
					outcomes[1] = rand_->Multinomial(randN1bjets_[j+1][k][l][i],prob[j][k][l][1]);
					outcomes[2] = rand_->Multinomial(randN2bjets_[j+1][k][l][i],prob[j][k][l][2]);
					outcomes[3] = rand_->Multinomial(randN3bjets_[j+1][k][l][i],prob[j][k][l][3]);
					
					randN0bjets_[j][k][l][i] = outcomes[0][0]+outcomes[1][0]+outcomes[2][0]+outcomes[3][0];
					tempNbjets[j][k][0][l] = randN0bjets_[j][k][l][i];
					randN1bjets_[j][k][l][i] = outcomes[0][1]+outcomes[1][1]+outcomes[2][1]+outcomes[3][1];
					tempNbjets[j][k][1][l] = randN1bjets_[j][k][l][i];
					randN2bjets_[j][k][l][i] = outcomes[0][2]+outcomes[1][2]+outcomes[2][2]+outcomes[3][2];
					tempNbjets[j][k][2][l] = randN2bjets_[j][k][l][i];
					randN3bjets_[j][k][l][i] = outcomes[0][3]+outcomes[1][3]+outcomes[2][3]+outcomes[3][3];
					tempNbjets[j][k][3][l] = randN3bjets_[j][k][l][i];
				}
			}
			
		}

		vjEst_->FillInputs(tempNbjets);
		
		if(!UseMJLE){
			if(!UseUnBinMLE) vjEst_->  BinnedMaximumLikelihoodEst(method, option, SetFixedVarIdx, false, Print);
			else{
				RooMsgService::instance().setStreamStatus(1,kFALSE);
				vjEst_->UnBinnedMaximumLikelihoodEst(SetFixedVarIdx, doMinos, false, Print);
			}
		}
		else{
			if(!UseUnBinMLE) vjEst_->  BinnedMaximumJointWPLikelihoodEst(method, option, SetFixedVarIdx, false, Print); 
			else{
				RooMsgService::instance().setStreamStatus(1,kFALSE);
				vjEst_->UnBinnedMaximumJointWPLikelihoodEst(SetFixedVarIdx, doMinos, Print);
			}
		}
		for(int j=0;j<NbOfJetsBins_;j++) minLL[j][i] = vjEst_->GetMinLL(j);
		NbOfSuccessfulPE_++;

		/// Fill arrays with values of :
		for(int j=0;j<NbOfBtagWorkingPoint_;j++){
			for(int k=0;k<NbOfJetsBins_;k++){
				/// -- eb :
				PseudoExp_EstEb[j][k][i]        = vjEst_->GetEstEb(j,k);
				//PseudoExp_PullEstEb[j][k][i]    = vjEst_->GetEstEb(j,k)    - vjEst_->GetPredEb(j,k);
				PseudoExp_DiffEstEb[j][k][i]    = vjEst_->GetEstEb(j,k)    - vjEst_->GetPredEb(j,k);
				/// -- eudsc :
				PseudoExp_EstEudsc[j][k][i]     = vjEst_->GetEstEudsc(j,k);
				PseudoExp_DiffEstEudsc[j][k][i] = vjEst_->GetEstEudsc(j,k) - vjEst_->GetPredEudsc(j,k);
				/// -- euds :
				PseudoExp_EstEuds[j][k][i]      = vjEst_->GetEstEuds(j,k);
				PseudoExp_DiffEstEuds[j][k][i]  = vjEst_->GetEstEuds(j,k)  - vjEst_->GetPredEuds(j,k);

				for(int l=0;l<NbOfBJetsBins_;l++){
					if(k==0) PseudoExp_EstNvlike[j][NbOfJetsBins_][l][i]      = vjEst_->GetEstNv( j,-1,l);
					if(k==0) PseudoExp_EstNvblike[j][NbOfJetsBins_][l][i]     = vjEst_->GetEstNvb(j,-1,l);
					if(k==0) PseudoExp_EstNttlike[j][NbOfJetsBins_][l][i]     = vjEst_->GetEstNtt(j,-1,l);
					if(k==0) PseudoExp_DiffEstNvlike[j][NbOfJetsBins_][l][i]  = vjEst_->GetEstNv( j,-1,l) - vjEst_->GetPredNv( j,-1,l);
					if(k==0) PseudoExp_DiffEstNvblike[j][NbOfJetsBins_][l][i] = vjEst_->GetEstNvb(j,-1,l) - vjEst_->GetPredNvb(j,-1,l);
					if(k==0) PseudoExp_DiffEstNttlike[j][NbOfJetsBins_][l][i] = vjEst_->GetEstNtt(j,-1,l) - vjEst_->GetPredNtt(j,-1,l);
					PseudoExp_EstNvlike[j][k][l][i]      = vjEst_->GetEstNv( j,k,l);
					PseudoExp_EstNvblike[j][k][l][i]     = vjEst_->GetEstNvb(j,k,l);
					PseudoExp_EstNttlike[j][k][l][i]     = vjEst_->GetEstNtt(j,k,l);
					PseudoExp_DiffEstNvlike[j][k][l][i]  = vjEst_->GetEstNv( j,k,l) - vjEst_->GetPredNv( j,k,l);
					PseudoExp_DiffEstNvblike[j][k][l][i] = vjEst_->GetEstNvb(j,k,l) - vjEst_->GetPredNvb(j,k,l);
					PseudoExp_DiffEstNttlike[j][k][l][i] = vjEst_->GetEstNtt(j,k,l) - vjEst_->GetPredNtt(j,k,l);
				}
				PseudoExp_EstNvlike[j][k][NbOfBJetsBins_][i]      = vjEst_->GetEstNv( j,k);
				PseudoExp_EstNvblike[j][k][NbOfBJetsBins_][i]     = vjEst_->GetEstNvb(j,k);
				PseudoExp_EstNttlike[j][k][NbOfBJetsBins_][i]     = vjEst_->GetEstNtt(j,k);
				if(doMinos){
					if((vjEst_->GetEstNv( j,k) - InitEstNv[k] > 0 ? vjEst_->GetEstNvErrUp(j,k) : fabs(vjEst_->GetEstNvErrDown(j,k)))!= 0){
						PseudoExp_PullEstNvlike[j][k][NbOfBJetsBins_][i]  = (vjEst_->GetEstNv( j,k) - InitEstNv[k] > 0 ? (vjEst_->GetEstNv(j,k)-InitEstNv[k])/vjEst_->GetEstNvErrUp(j,k) : (vjEst_->GetEstNv(j,k)-InitEstNv[k])/fabs(vjEst_->GetEstNvErrDown(j,k)));
					}else   PseudoExp_PullEstNvlike[j][k][NbOfBJetsBins_][i]  = (vjEst_->GetEstNv( j,k) - InitEstNv[k])/ vjEst_->GetEstNvErr( j,k);
					PseudoExp_PullEstNvblike[j][k][NbOfBJetsBins_][i] = (vjEst_->GetEstNvb(j,k) - 0);
					if((vjEst_->GetEstNtt( j,k) - InitEstNtt[k] > 0 ? vjEst_->GetEstNttErrUp(j,k) : fabs(vjEst_->GetEstNttErrDown(j,k)))!= 0){
						PseudoExp_PullEstNttlike[j][k][NbOfBJetsBins_][i]  = (vjEst_->GetEstNtt( j,k) - InitEstNtt[k] > 0 ? (vjEst_->GetEstNtt(j,k)-InitEstNtt[k])/vjEst_->GetEstNttErrUp(j,k) : (vjEst_->GetEstNtt(j,k)-InitEstNtt[k])/fabs(vjEst_->GetEstNttErrDown(j,k)));
					}else   PseudoExp_PullEstNttlike[j][k][NbOfBJetsBins_][i]  = (vjEst_->GetEstNtt( j,k) - InitEstNtt[k])/ vjEst_->GetEstNttErr( j,k);
				}else{
					PseudoExp_PullEstNvlike[j][k][NbOfBJetsBins_][i]  = (vjEst_->GetEstNv( j,k) - InitEstNv[k])/ vjEst_->GetEstNvErr( j,k);
					PseudoExp_PullEstNvblike[j][k][NbOfBJetsBins_][i] = (vjEst_->GetEstNvb(j,k) - 0);
					PseudoExp_PullEstNttlike[j][k][NbOfBJetsBins_][i] = (vjEst_->GetEstNtt(j,k) - InitEstNtt[k])/vjEst_->GetEstNttErr(j,k);
				}
				PseudoExp_DiffEstNvlike[j][k][NbOfBJetsBins_][i]  = vjEst_->GetEstNv( j,k) - vjEst_->GetPredNv( j,k);
				PseudoExp_DiffEstNvblike[j][k][NbOfBJetsBins_][i] = vjEst_->GetEstNvb(j,k) - vjEst_->GetPredNvb(j,k);
				PseudoExp_DiffEstNttlike[j][k][NbOfBJetsBins_][i] = vjEst_->GetEstNtt(j,k) - vjEst_->GetPredNtt(j,k);
			}
			PseudoExp_EstNvlike[j][NbOfJetsBins_][NbOfBJetsBins_][i]      = vjEst_->GetEstNv( j);
			PseudoExp_EstNvblike[j][NbOfJetsBins_][NbOfBJetsBins_][i]     = vjEst_->GetEstNvb(j);
			PseudoExp_EstNttlike[j][NbOfJetsBins_][NbOfBJetsBins_][i]     = vjEst_->GetEstNtt(j);
			PseudoExp_DiffEstNvlike[j][NbOfJetsBins_][NbOfBJetsBins_][i]  = vjEst_->GetEstNv( j) - vjEst_->GetPredNv( j);
			PseudoExp_DiffEstNvblike[j][NbOfJetsBins_][NbOfBJetsBins_][i] = vjEst_->GetEstNvb(j) - vjEst_->GetPredNvb(j);
			PseudoExp_DiffEstNttlike[j][NbOfJetsBins_][NbOfBJetsBins_][i] = vjEst_->GetEstNtt(j) - vjEst_->GetPredNtt(j);
		}
	}

	char name[100];
	double Xmax[4];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		for(int j=0; j<NbOfJetsBins_; j++){
			sprintf(name,"Eb_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hEb_[i][j]    = new TH1F(name,"",250,*min_element(PseudoExp_EstEb[i][j],   PseudoExp_EstEb[i][j]+NbOfPseudoExp_),   *max_element(PseudoExp_EstEb[i][j],   PseudoExp_EstEb[i][j]+NbOfPseudoExp_));
			hEb_[i][j]    ->GetXaxis()->SetTitle("#epsilon_{b}^{Est.}");
			sprintf(name,"Eudsc_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hEudsc_[i][j] = new TH1F(name,"",250,*min_element(PseudoExp_EstEudsc[i][j],PseudoExp_EstEudsc[i][j]+NbOfPseudoExp_),*max_element(PseudoExp_EstEudsc[i][j],PseudoExp_EstEudsc[i][j]+NbOfPseudoExp_));
			hEudsc_[i][j] ->GetXaxis()->SetTitle("#epsilon_{udsc}^{Est.}");
			sprintf(name,"Euds_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hEuds_[i][j] = new TH1F(name,"",250,*min_element( PseudoExp_EstEuds[i][j], PseudoExp_EstEuds[i][j]+NbOfPseudoExp_), *max_element(PseudoExp_EstEuds[i][j], PseudoExp_EstEuds[i][j]+NbOfPseudoExp_));
			hEuds_[i][j] ->GetXaxis()->SetTitle("#epsilon_{uds}^{Est.}");

			sprintf(name,"PullEb_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hPullEb_[i][j]    = new TH1F(name,"",250,*min_element(PseudoExp_PullEstEb[i][j],   PseudoExp_PullEstEb[i][j]+NbOfPseudoExp_),   *max_element(PseudoExp_PullEstEb[i][j],   PseudoExp_PullEstEb[i][j]+NbOfPseudoExp_));
			hPullEb_[i][j]    ->GetXaxis()->SetTitle("#epsilon_{b}^{Est.} : pull");
			sprintf(name,"PullEudsc_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hPullEudsc_[i][j] = new TH1F(name,"",250,*min_element(PseudoExp_PullEstEudsc[i][j],PseudoExp_PullEstEudsc[i][j]+NbOfPseudoExp_),*max_element(PseudoExp_PullEstEudsc[i][j],PseudoExp_PullEstEudsc[i][j]+NbOfPseudoExp_));
			hPullEudsc_[i][j] ->GetXaxis()->SetTitle("#epsilon_{udsc}^{Est.} : pull");
			sprintf(name,"PullEuds_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hPullEuds_[i][j]  = new TH1F(name,"",250,*min_element( PseudoExp_PullEstEuds[i][j], PseudoExp_PullEstEuds[i][j]+NbOfPseudoExp_), *max_element(PseudoExp_PullEstEuds[i][j], PseudoExp_PullEstEuds[i][j]+NbOfPseudoExp_));
			hPullEuds_[i][j]  ->GetXaxis()->SetTitle("#epsilon_{uds}^{Est.} : pull");

			sprintf(name,"DiffEb_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hDiffEb_[i][j]    = new TH1F(name,"",250,*min_element(PseudoExp_DiffEstEb[i][j],   PseudoExp_DiffEstEb[i][j]+NbOfPseudoExp_),   *max_element(PseudoExp_DiffEstEb[i][j],   PseudoExp_DiffEstEb[i][j]+NbOfPseudoExp_));
			hDiffEb_[i][j]    ->GetXaxis()->SetTitle("#Delta(#epsilon_{b}^{Est.} - #epsilon_{b}^{Exp.})");
			sprintf(name,"DiffEudsc_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hDiffEudsc_[i][j] = new TH1F(name,"",250,*min_element(PseudoExp_DiffEstEudsc[i][j],PseudoExp_DiffEstEudsc[i][j]+NbOfPseudoExp_),*max_element(PseudoExp_DiffEstEudsc[i][j],PseudoExp_DiffEstEudsc[i][j]+NbOfPseudoExp_));
			hDiffEudsc_[i][j] ->GetXaxis()->SetTitle("#Delta(#epsilon_{udsc}^{Est.} - #epsilon_{udsc}^{Exp.})");
			sprintf(name,"DiffEuds_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hDiffEuds_[i][j] = new TH1F(name,"",250,*min_element( PseudoExp_DiffEstEuds[i][j], PseudoExp_DiffEstEuds[i][j]+NbOfPseudoExp_), *max_element(PseudoExp_DiffEstEuds[i][j], PseudoExp_DiffEstEuds[i][j]+NbOfPseudoExp_));
			hDiffEuds_[i][j] ->GetXaxis()->SetTitle("#Delta(#epsilon_{uds}^{Est.} - #epsilon_{uds}^{Exp.})");

			for(int k=0;k<NbOfDataSets_;k++){
				sprintf(name,"randN0bjets_%d_wp_%d_jets_%d_dataset",i,vjEst_->GetNjets(j),k);
				hRandN0bjets_[i][j][k] = new TH1F(name,"",250,*min_element(randN0bjets_[i][j][k],randN0bjets_[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(randN0bjets_[i][j][k],randN0bjets_[i][j][k]+NbOfPseudoExp_)*maxfact);
				Xmax[0] = hRandN0bjets_[i][j][k]->GetXaxis()->GetXmax();
				sprintf(name,"randN1bjets_%d_wp_%d_jets_%d_dataset",i,vjEst_->GetNjets(j),k);
				hRandN1bjets_[i][j][k] = new TH1F(name,"",250,*min_element(randN1bjets_[i][j][k],randN1bjets_[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(randN1bjets_[i][j][k],randN1bjets_[i][j][k]+NbOfPseudoExp_)*maxfact);
				Xmax[1] = hRandN1bjets_[i][j][k]->GetXaxis()->GetXmax();
				sprintf(name,"randN2bjets_%d_wp_%d_jets_%d_dataset",i,vjEst_->GetNjets(j),k);
				hRandN2bjets_[i][j][k] = new TH1F(name,"",250,*min_element(randN2bjets_[i][j][k],randN2bjets_[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(randN2bjets_[i][j][k],randN2bjets_[i][j][k]+NbOfPseudoExp_)*maxfact);
				Xmax[2] = hRandN2bjets_[i][j][k]->GetXaxis()->GetXmax();
				sprintf(name,"randN3bjets_%d_wp_%d_jets_%d_dataset",i,vjEst_->GetNjets(j),k);
				hRandN3bjets_[i][j][k] = new TH1F(name,"",250,*min_element(randN3bjets_[i][j][k],randN3bjets_[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(randN3bjets_[i][j][k],randN3bjets_[i][j][k]+NbOfPseudoExp_)*maxfact);
				Xmax[3] = hRandN3bjets_[i][j][k]->GetXaxis()->GetXmax();
				sprintf(name,"randNbjets_vs_minLL_%d_wp_%d_jets_%d_dataset",i,vjEst_->GetNjets(j),k);
				hRandNbjets_vs_MinLL_[i][j][k] = new TH3F(name,"",4,0,4,(int)*max_element(Xmax,Xmax+4),0,*max_element(Xmax,Xmax+4),250,*min_element(minLL[j],minLL[j]+NbOfPseudoExp_),*max_element(minLL[j],minLL[j]+NbOfPseudoExp_));
			}

			for(int k=0; k<NbOfBJetsBins_; k++){
				sprintf(name,"Nv_%d_wp_%d_jets_%d_bjets_PseudoExp",i,vjEst_->GetNjets(j),k);
				hNv_[i][j][k] = new TH1F(name,"",250,*min_element(PseudoExp_EstNvlike[i][j][k], PseudoExp_EstNvlike[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_EstNvlike[i][j][k], PseudoExp_EstNvlike[i][j][k]+NbOfPseudoExp_)*maxfact);
				hNv_[i][j][k] ->SetStats();
				hNv_[i][j][k] ->GetXaxis()->SetTitle("N^{est}_{V-like}");
				sprintf(name,"Nvb_%d_wp_%d_jets_%d_bjets_PseudoExp",i,vjEst_->GetNjets(j),k);
				hNvb_[i][j][k] = new TH1F(name,"",250,*min_element(PseudoExp_EstNvblike[i][j][k], PseudoExp_EstNvblike[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_EstNvblike[i][j][k], PseudoExp_EstNvblike[i][j][k]+NbOfPseudoExp_)*maxfact);
				hNvb_[i][j][k] ->SetStats();
				hNvb_[i][j][k] ->GetXaxis()->SetTitle("N^{est}_{Vb-like}");
				sprintf(name,"Ntt_%d_wp_%d_jets_%d_bjets_PseudoExp",i,vjEst_->GetNjets(j),k);
				hNtt_[i][j][k] = new TH1F(name,"",250,*min_element(PseudoExp_EstNttlike[i][j][k], PseudoExp_EstNttlike[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_EstNttlike[i][j][k], PseudoExp_EstNttlike[i][j][k]+NbOfPseudoExp_)*maxfact);
				hNtt_[i][j][k] ->SetStats();
				hNtt_[i][j][k] ->GetXaxis()->SetTitle("N^{est}_{#bar{t}t-like}");

				sprintf(name,"PullNv_%d_wp_%d_jets_%d_bjets_PseudoExp",i,vjEst_->GetNjets(j),k);
				hPullNv_[i][j][k]  = new TH1F(name,"",250,*min_element(PseudoExp_PullEstNvlike[i][j][k], PseudoExp_PullEstNvlike[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_PullEstNvlike[i][j][k], PseudoExp_PullEstNvlike[i][j][k]+NbOfPseudoExp_)*maxfact);
				hPullNv_[i][j][k]  ->SetStats();
				hPullNv_[i][j][k]  ->GetXaxis()->SetTitle("N^{est}_{V-like} : pull");
				sprintf(name,"PullNvb_%d_wp_%d_jets_%d_bjets_PseudoExp",i,vjEst_->GetNjets(j),k);
				hPullNvb_[i][j][k] = new TH1F(name,"",250,*min_element(PseudoExp_PullEstNvblike[i][j][k], PseudoExp_PullEstNvblike[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_PullEstNvblike[i][j][k], PseudoExp_PullEstNvblike[i][j][k]+NbOfPseudoExp_)*maxfact);
				hPullNvb_[i][j][k] ->SetStats();
				hPullNvb_[i][j][k] ->GetXaxis()->SetTitle("N^{est}_{Vb-like} : pull");
				sprintf(name,"PullNtt_%d_wp_%d_jets_%d_bjets_PseudoExp",i,vjEst_->GetNjets(j),k);
				hPullNtt_[i][j][k] = new TH1F(name,"",250,*min_element(PseudoExp_PullEstNttlike[i][j][k], PseudoExp_PullEstNttlike[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_PullEstNttlike[i][j][k], PseudoExp_PullEstNttlike[i][j][k]+NbOfPseudoExp_)*maxfact);
				hPullNtt_[i][j][k] ->SetStats();
				hPullNtt_[i][j][k] ->GetXaxis()->SetTitle("N^{est}_{#bar{t}t-like} : pull");

				sprintf(name,"DiffNv_%d_wp_%d_jets_%d_bjets_PseudoExp",i,vjEst_->GetNjets(j),k);
				hDiffNv_[i][j][k] = new TH1F(name,"",250,*min_element(PseudoExp_DiffEstNvlike[i][j][k], PseudoExp_DiffEstNvlike[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_DiffEstNvlike[i][j][k], PseudoExp_DiffEstNvlike[i][j][k]+NbOfPseudoExp_)*maxfact);
				hDiffNv_[i][j][k] ->SetStats();
				hDiffNv_[i][j][k] ->GetXaxis()->SetTitle("#Delta(N^{est}_{V-like}-N^{exp}_{V-like})");
				sprintf(name,"DiffNvb_%d_wp_%d_jets_%d_bjets_PseudoExp",i,vjEst_->GetNjets(j),k);
				hDiffNvb_[i][j][k] = new TH1F(name,"",250,*min_element(PseudoExp_DiffEstNvblike[i][j][k], PseudoExp_DiffEstNvblike[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_DiffEstNvblike[i][j][k], PseudoExp_DiffEstNvblike[i][j][k]+NbOfPseudoExp_)*maxfact);
				hDiffNvb_[i][j][k] ->SetStats();
				hDiffNvb_[i][j][k] ->GetXaxis()->SetTitle("#Delta(N^{est}_{Vb-like}-N^{exp}_{Vb-like})");
				sprintf(name,"DiffNtt_%d_wp_%d_jets_%d_bjets_PseudoExp",i,vjEst_->GetNjets(j),k);
				hDiffNtt_[i][j][k] = new TH1F(name,"",250,*min_element(PseudoExp_DiffEstNttlike[i][j][k], PseudoExp_DiffEstNttlike[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_DiffEstNttlike[i][j][k], PseudoExp_DiffEstNttlike[i][j][k]+NbOfPseudoExp_)*maxfact);
				hDiffNtt_[i][j][k] ->SetStats();
				hDiffNtt_[i][j][k] ->GetXaxis()->SetTitle("#Delta(N^{est}_{#bar{t}t-like}-N^{exp}_{#bar{t}t-like})");

				sprintf(name,"DiffNv_vs_MinLL_%d_wp_%d_jets_%d_bjets_PseudoExp",i,vjEst_->GetNjets(j),k);
				hDiffNv_vs_MinLL_[i][j][k] = new TH2F(name,"",250,*min_element(PseudoExp_DiffEstNvlike[i][j][k],PseudoExp_DiffEstNvlike[i][j][k]+NbOfPseudoExp_)*maxfact,*max_element(PseudoExp_DiffEstNvlike[i][j][k],PseudoExp_DiffEstNvlike[i][j][k]+NbOfPseudoExp_)*maxfact,250,*min_element(minLL[j],minLL[j]+NbOfPseudoExp_),*max_element(minLL[j],minLL[j]+NbOfPseudoExp_));
				hDiffNv_vs_MinLL_[i][j][k] ->SetStats();
				hDiffNv_vs_MinLL_[i][j][k] ->GetXaxis()->SetTitle("#Delta(N^{est}_{V-like}-N^{exp}_{V-like})");
				hDiffNv_vs_MinLL_[i][j][k] ->GetYaxis()->SetTitle("min(-log(L))");
				sprintf(name,"DiffNvb_vs_MinLL_%d_wp_%d_jets_%d_bjets_PseudoExp",i,vjEst_->GetNjets(j),k);
				hDiffNvb_vs_MinLL_[i][j][k] = new TH2F(name,"",250,*min_element(PseudoExp_DiffEstNvblike[i][j][k],PseudoExp_DiffEstNvblike[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_DiffEstNvblike[i][j][k],PseudoExp_DiffEstNvblike[i][j][k]+NbOfPseudoExp_)*maxfact,250,*min_element(minLL[j],minLL[j]+NbOfPseudoExp_),*max_element(minLL[j],minLL[j]+NbOfPseudoExp_));
				hDiffNvb_vs_MinLL_[i][j][k] ->SetStats();
				hDiffNvb_vs_MinLL_[i][j][k] ->GetXaxis()->SetTitle("#Delta(N^{est}_{Vb-like}-N^{exp}_{Vb-like})");
				hDiffNvb_vs_MinLL_[i][j][k] ->GetYaxis()->SetTitle("min(-log(L))");
				sprintf(name,"DiffNtt_vs_MinLL_%d_wp_%d_jets_%d_bjets_PseudoExp",i,vjEst_->GetNjets(j),k);
				hDiffNtt_vs_MinLL_[i][j][k] = new TH2F(name,"",250,*min_element(PseudoExp_DiffEstNttlike[i][j][k],PseudoExp_DiffEstNttlike[i][j][k]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_DiffEstNttlike[i][j][k],PseudoExp_DiffEstNttlike[i][j][k]+NbOfPseudoExp_)*maxfact,250,*min_element(minLL[j],minLL[j]+NbOfPseudoExp_),*max_element(minLL[j],minLL[j]+NbOfPseudoExp_));
				hDiffNtt_vs_MinLL_[i][j][k] ->SetStats();
				hDiffNtt_vs_MinLL_[i][j][k] ->GetXaxis()->SetTitle("#Delta(N^{est}_{#bar{t}t-like}-N^{exp}_{#bar{t}t-like})");
				hDiffNtt_vs_MinLL_[i][j][k] ->GetYaxis()->SetTitle("min(-log(L))");
			}
		}
	
		for(int j=0; j<NbOfJetsBins_; j++){
			sprintf(name,"Nv_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hNv_[i][j][NbOfBJetsBins_] = new TH1F(name,"",250,*min_element(PseudoExp_EstNvlike[i][j][NbOfBJetsBins_], PseudoExp_EstNvlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_EstNvlike[i][j][NbOfBJetsBins_], PseudoExp_EstNvlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*maxfact);
			hNv_[i][j][NbOfBJetsBins_] ->SetStats();
			hNv_[i][j][NbOfBJetsBins_] ->GetXaxis()->SetTitle("N^{est}_{V-like}");
			sprintf(name,"Nvb_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hNvb_[i][j][NbOfBJetsBins_] = new TH1F(name,"",250,*min_element(PseudoExp_EstNvblike[i][j][NbOfBJetsBins_], PseudoExp_EstNvblike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_EstNvblike[i][j][NbOfBJetsBins_], PseudoExp_EstNvblike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*maxfact);
			hNvb_[i][j][NbOfBJetsBins_] ->SetStats();
			hNvb_[i][j][NbOfBJetsBins_] ->GetXaxis()->SetTitle("N^{est}_{Vb-like}");
			sprintf(name,"Ntt_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hNtt_[i][j][NbOfBJetsBins_] = new TH1F(name,"",250,*min_element(PseudoExp_EstNttlike[i][j][NbOfBJetsBins_], PseudoExp_EstNttlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_EstNttlike[i][j][NbOfBJetsBins_], PseudoExp_EstNttlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*maxfact);
			hNtt_[i][j][NbOfBJetsBins_] ->SetStats();
			hNtt_[i][j][NbOfBJetsBins_] ->GetXaxis()->SetTitle("N^{est}_{#bar{t}t-like}");

			sprintf(name,"PullNv_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hPullNv_[i][j][NbOfBJetsBins_]  = new TH1F(name,"",250,*min_element(PseudoExp_PullEstNvlike[i][j][NbOfBJetsBins_], PseudoExp_PullEstNvlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_PullEstNvlike[i][j][NbOfBJetsBins_], PseudoExp_PullEstNvlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*maxfact);
			hPullNv_[i][j][NbOfBJetsBins_]  ->SetStats();
			hPullNv_[i][j][NbOfBJetsBins_]  ->GetXaxis()->SetTitle("N^{est}_{V-like} : pull");
			sprintf(name,"PullNvb_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hPullNvb_[i][j][NbOfBJetsBins_] = new TH1F(name,"",250,*min_element(PseudoExp_PullEstNvblike[i][j][NbOfBJetsBins_], PseudoExp_PullEstNvblike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_PullEstNvblike[i][j][NbOfBJetsBins_], PseudoExp_PullEstNvblike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*maxfact);
			hPullNvb_[i][j][NbOfBJetsBins_] ->SetStats();
			hPullNvb_[i][j][NbOfBJetsBins_] ->GetXaxis()->SetTitle("N^{est}_{Vb-like} : pull");
			sprintf(name,"PullNtt_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hPullNtt_[i][j][NbOfBJetsBins_] = new TH1F(name,"",250,*min_element(PseudoExp_PullEstNttlike[i][j][NbOfBJetsBins_], PseudoExp_PullEstNttlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_PullEstNttlike[i][j][NbOfBJetsBins_], PseudoExp_PullEstNttlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*maxfact);
			hPullNtt_[i][j][NbOfBJetsBins_] ->SetStats();
			hPullNtt_[i][j][NbOfBJetsBins_] ->GetXaxis()->SetTitle("N^{est}_{#bar{t}t-like} : pull");

			sprintf(name,"DiffNv_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hDiffNv_[i][j][NbOfBJetsBins_] = new TH1F(name,"",250,*min_element(PseudoExp_DiffEstNvlike[i][j][NbOfBJetsBins_], PseudoExp_DiffEstNvlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_DiffEstNvlike[i][j][NbOfBJetsBins_], PseudoExp_DiffEstNvlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*maxfact);
			hDiffNv_[i][j][NbOfBJetsBins_] ->SetStats();
			hDiffNv_[i][j][NbOfBJetsBins_] ->GetXaxis()->SetTitle("#Delta(N^{est}_{V-like}-N^{exp}_{V-like})");
			sprintf(name,"DiffNvb_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hDiffNvb_[i][j][NbOfBJetsBins_] = new TH1F(name,"",250,*min_element(PseudoExp_DiffEstNvblike[i][j][NbOfBJetsBins_], PseudoExp_DiffEstNvblike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_DiffEstNvblike[i][j][NbOfBJetsBins_], PseudoExp_DiffEstNvblike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*maxfact);
			hDiffNvb_[i][j][NbOfBJetsBins_] ->SetStats();
			hDiffNvb_[i][j][NbOfBJetsBins_] ->GetXaxis()->SetTitle("#Delta(N^{est}_{Vb-like}-N^{exp}_{Vb-like})");
			sprintf(name,"DiffNtt_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hDiffNtt_[i][j][NbOfBJetsBins_] = new TH1F(name,"",250,*min_element(PseudoExp_DiffEstNttlike[i][j][NbOfBJetsBins_], PseudoExp_DiffEstNttlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_DiffEstNttlike[i][j][NbOfBJetsBins_], PseudoExp_DiffEstNttlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*maxfact);
			hDiffNtt_[i][j][NbOfBJetsBins_] ->SetStats();
			hDiffNtt_[i][j][NbOfBJetsBins_] ->GetXaxis()->SetTitle("#Delta(N^{est}_{#bar{t}t-like}-N^{exp}_{#bar{t}t-like})");

			sprintf(name,"DiffNv_vs_MinLL_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hDiffNv_vs_MinLL_[i][j][NbOfBJetsBins_] = new TH2F(name,"",250,*min_element(PseudoExp_DiffEstNvlike[i][j][NbOfBJetsBins_],PseudoExp_DiffEstNvlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_DiffEstNvlike[i][j][NbOfBJetsBins_],PseudoExp_DiffEstNvlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*maxfact,250,*min_element(minLL[j],minLL[j]+NbOfPseudoExp_),*max_element(minLL[j],minLL[j]+NbOfPseudoExp_));
			hDiffNv_vs_MinLL_[i][j][NbOfBJetsBins_] ->SetStats();
			hDiffNv_vs_MinLL_[i][j][NbOfBJetsBins_] ->GetXaxis()->SetTitle("#Delta(N^{est}_{V-like}-N^{exp}_{V-like})");
			hDiffNv_vs_MinLL_[i][j][NbOfBJetsBins_] ->GetYaxis()->SetTitle("min(-log(L))");
			sprintf(name,"DiffNvb_vs_MinLL_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hDiffNvb_vs_MinLL_[i][j][NbOfBJetsBins_] = new TH2F(name,"",250,*min_element(PseudoExp_DiffEstNvblike[i][j][NbOfBJetsBins_],PseudoExp_DiffEstNvblike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_DiffEstNvblike[i][j][NbOfBJetsBins_],PseudoExp_DiffEstNvblike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*maxfact,250,*min_element(minLL[j],minLL[j]+NbOfPseudoExp_),*max_element(minLL[j],minLL[j]+NbOfPseudoExp_));
			hDiffNvb_vs_MinLL_[i][j][NbOfBJetsBins_] ->SetStats();
			hDiffNvb_vs_MinLL_[i][j][NbOfBJetsBins_] ->GetXaxis()->SetTitle("#Delta(N^{est}_{Vb-like}-N^{exp}_{Vb-like})");
			hDiffNvb_vs_MinLL_[i][j][NbOfBJetsBins_] ->GetYaxis()->SetTitle("min(-log(L))");
			sprintf(name,"DiffNtt_vs_MinLL_%d_wp_%d_jets_PseudoExp",i,vjEst_->GetNjets(j));
			hDiffNtt_vs_MinLL_[i][j][NbOfBJetsBins_] = new TH2F(name,"",250,*min_element(PseudoExp_DiffEstNttlike[i][j][NbOfBJetsBins_],PseudoExp_DiffEstNttlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_DiffEstNttlike[i][j][NbOfBJetsBins_],PseudoExp_DiffEstNttlike[i][j][NbOfBJetsBins_]+NbOfPseudoExp_)*maxfact,250,*min_element(minLL[j],minLL[j]+NbOfPseudoExp_),*max_element(minLL[j],minLL[j]+NbOfPseudoExp_));
			hDiffNtt_vs_MinLL_[i][j][NbOfBJetsBins_] ->SetStats();
			hDiffNtt_vs_MinLL_[i][j][NbOfBJetsBins_] ->GetXaxis()->SetTitle("#Delta(N^{est}_{#bar{t}t-like}-N^{exp}_{#bar{t}t-like})");
			hDiffNtt_vs_MinLL_[i][j][NbOfBJetsBins_] ->GetYaxis()->SetTitle("min(-log(L))");
		}

		for(int j=0; j<NbOfBJetsBins_; j++){
			sprintf(name,"Nv_%d_wp_%d_b-jets_PseudoExp",i,j);
			hNv_[i][NbOfJetsBins_][j] = new TH1F(name,"",250,*min_element(PseudoExp_EstNvlike[i][NbOfJetsBins_][j], PseudoExp_EstNvlike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_EstNvlike[i][NbOfJetsBins_][j], PseudoExp_EstNvlike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*maxfact);
			hNv_[i][NbOfJetsBins_][j] ->SetStats();
			hNv_[i][NbOfJetsBins_][j] ->GetXaxis()->SetTitle("N^{est}_{V-like}");
			sprintf(name,"Nvb_%d_wp_%d_b-jets_PseudoExp",i,j);
			hNvb_[i][NbOfJetsBins_][j] = new TH1F(name,"",250,*min_element(PseudoExp_EstNvblike[i][NbOfJetsBins_][j], PseudoExp_EstNvblike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_EstNvblike[i][NbOfJetsBins_][j], PseudoExp_EstNvblike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*maxfact);
			hNvb_[i][NbOfJetsBins_][j] ->SetStats();
			hNvb_[i][NbOfJetsBins_][j] ->GetXaxis()->SetTitle("N^{est}_{Vb-like}");
			sprintf(name,"Ntt_%d_wp_%d_b-jets_PseudoExp",i,j);
			hNtt_[i][NbOfJetsBins_][j] = new TH1F(name,"",250,*min_element(PseudoExp_EstNttlike[i][NbOfJetsBins_][j], PseudoExp_EstNttlike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_EstNttlike[i][NbOfJetsBins_][j], PseudoExp_EstNttlike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*maxfact);
			hNtt_[i][NbOfJetsBins_][j] ->SetStats();
			hNtt_[i][NbOfJetsBins_][j] ->GetXaxis()->SetTitle("N^{est}_{#bar{t}t-like}");

			sprintf(name,"PullNv_%d_wp_%d_b-jets_PseudoExp",i,j);
			hPullNv_[i][NbOfJetsBins_][j] = new TH1F(name,"",250,*min_element(PseudoExp_PullEstNvlike[i][NbOfJetsBins_][j], PseudoExp_PullEstNvlike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_PullEstNvlike[i][NbOfJetsBins_][j], PseudoExp_PullEstNvlike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*maxfact);
			hPullNv_[i][NbOfJetsBins_][j] ->SetStats();
			hPullNv_[i][NbOfJetsBins_][j] ->GetXaxis()->SetTitle("N^{est}_{V-like} : pull");
			sprintf(name,"PullNvb_%d_wp_%d_b-jets_PseudoExp",i,j);
			hPullNvb_[i][NbOfJetsBins_][j] = new TH1F(name,"",250,*min_element(PseudoExp_PullEstNvblike[i][NbOfJetsBins_][j], PseudoExp_PullEstNvblike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_PullEstNvblike[i][NbOfJetsBins_][j], PseudoExp_PullEstNvblike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*maxfact);
			hPullNvb_[i][NbOfJetsBins_][j] ->SetStats();
			hPullNvb_[i][NbOfJetsBins_][j] ->GetXaxis()->SetTitle("N^{est}_{Vb-like} : pull");
			sprintf(name,"PullNtt_%d_wp_%d_b-jets_PseudoExp",i,j);
			hPullNtt_[i][NbOfJetsBins_][j] = new TH1F(name,"",250,*min_element(PseudoExp_PullEstNttlike[i][NbOfJetsBins_][j], PseudoExp_PullEstNttlike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_PullEstNttlike[i][NbOfJetsBins_][j], PseudoExp_PullEstNttlike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*maxfact);
			hPullNtt_[i][NbOfJetsBins_][j] ->SetStats();
			hPullNtt_[i][NbOfJetsBins_][j] ->GetXaxis()->SetTitle("N^{est}_{#bar{t}t-like} : pull");

			sprintf(name,"DiffNv_%d_wp_%d_b-jets_PseudoExp",i,j);
			hDiffNv_[i][NbOfJetsBins_][j] = new TH1F(name,"",250,*min_element(PseudoExp_DiffEstNvlike[i][NbOfJetsBins_][j], PseudoExp_DiffEstNvlike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_DiffEstNvlike[i][NbOfJetsBins_][j], PseudoExp_DiffEstNvlike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*maxfact);
			hDiffNv_[i][NbOfJetsBins_][j] ->SetStats();
			hDiffNv_[i][NbOfJetsBins_][j] ->GetXaxis()->SetTitle("#Delta(N^{est}_{V-like}-N^{exp}_{V-like})");
			sprintf(name,"DiffNvb_%d_wp_%d_b-jets_PseudoExp",i,j);
			hDiffNvb_[i][NbOfJetsBins_][j] = new TH1F(name,"",250,*min_element(PseudoExp_DiffEstNvblike[i][NbOfJetsBins_][j], PseudoExp_DiffEstNvblike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_DiffEstNvblike[i][NbOfJetsBins_][j], PseudoExp_DiffEstNvblike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*maxfact);
			hDiffNvb_[i][NbOfJetsBins_][j] ->SetStats();
			hDiffNvb_[i][NbOfJetsBins_][j] ->GetXaxis()->SetTitle("#Delta(N^{est}_{Vb-like}-N^{exp}_{Vb-like})");
			sprintf(name,"DiffNtt_%d_wp_%d_b-jets_PseudoExp",i,j);
			hDiffNtt_[i][NbOfJetsBins_][j] = new TH1F(name,"",250,*min_element(PseudoExp_DiffEstNttlike[i][NbOfJetsBins_][j], PseudoExp_DiffEstNttlike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*minfact,*max_element(PseudoExp_DiffEstNttlike[i][NbOfJetsBins_][j], PseudoExp_DiffEstNttlike[i][NbOfJetsBins_][j]+NbOfPseudoExp_)*maxfact);
			hDiffNtt_[i][NbOfJetsBins_][j] ->SetStats();
			hDiffNtt_[i][NbOfJetsBins_][j] ->GetXaxis()->SetTitle("#Delta(N^{est}_{#bar{t}t-like}-N^{exp}_{#bar{t}t-like})");
		}
	}

	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		for(int j=0; j<NbOfJetsBins_; j++){
			for(int k=0; k<NbOfPseudoExp_; k++){
				hEb_[i][j]   ->Fill(PseudoExp_EstEb[i][j][k]);
				hEudsc_[i][j]->Fill(PseudoExp_EstEudsc[i][j][k]);
				hEuds_[i][j] ->Fill(PseudoExp_EstEuds[i][j][k]);
				hPullEb_[i][j]   ->Fill(PseudoExp_PullEstEb[i][j][k]);
				hPullEudsc_[i][j]->Fill(PseudoExp_PullEstEudsc[i][j][k]);
				hPullEuds_[i][j] ->Fill(PseudoExp_PullEstEuds[i][j][k]);
				hDiffEb_[i][j]   ->Fill(PseudoExp_DiffEstEb[i][j][k]);
				hDiffEudsc_[i][j]->Fill(PseudoExp_DiffEstEudsc[i][j][k]);
				hDiffEuds_[i][j] ->Fill(PseudoExp_DiffEstEuds[i][j][k]);
			}
			for(int k=0; k<NbOfDataSets_; k++){
				for(int l=0; l<NbOfPseudoExp_; l++){
					hRandN0bjets_[i][j][k]->Fill(randN0bjets_[i][j][k][l]);
					hRandN1bjets_[i][j][k]->Fill(randN1bjets_[i][j][k][l]);
					hRandN2bjets_[i][j][k]->Fill(randN2bjets_[i][j][k][l]);
					hRandN3bjets_[i][j][k]->Fill(randN3bjets_[i][j][k][l]);
					hRandNbjets_vs_MinLL_[i][j][k]->Fill(0.,randN0bjets_[i][j][k][l],minLL[j][l]);
					hRandNbjets_vs_MinLL_[i][j][k]->Fill(1.,randN1bjets_[i][j][k][l],minLL[j][l]);
					hRandNbjets_vs_MinLL_[i][j][k]->Fill(2.,randN2bjets_[i][j][k][l],minLL[j][l]);
					hRandNbjets_vs_MinLL_[i][j][k]->Fill(3.,randN3bjets_[i][j][k][l],minLL[j][l]);
				}
			}
			for(int k=0; k<NbOfBJetsBins_; k++){
				for(int l=0; l<NbOfPseudoExp_; l++){
					hNv_[i][j][k] ->Fill(PseudoExp_EstNvlike[i][j][k][l]);
					hNvb_[i][j][k]->Fill(PseudoExp_EstNvblike[i][j][k][l]);
					hNtt_[i][j][k]->Fill(PseudoExp_EstNttlike[i][j][k][l]);
					hPullNv_[i][j][k] ->Fill(PseudoExp_PullEstNvlike[i][j][k][l]);
					hPullNvb_[i][j][k]->Fill(PseudoExp_PullEstNvblike[i][j][k][l]);
					hPullNtt_[i][j][k]->Fill(PseudoExp_PullEstNttlike[i][j][k][l]);
					hDiffNv_[i][j][k] ->Fill(PseudoExp_DiffEstNvlike[i][j][k][l]);
					hDiffNvb_[i][j][k]->Fill(PseudoExp_DiffEstNvblike[i][j][k][l]);
					hDiffNtt_[i][j][k]->Fill(PseudoExp_DiffEstNttlike[i][j][k][l]);
					hDiffNv_vs_MinLL_[i][j][k] ->Fill(PseudoExp_DiffEstNvlike[i][j][k][l],minLL[j][l]);
					hDiffNvb_vs_MinLL_[i][j][k]->Fill(PseudoExp_DiffEstNvblike[i][j][k][l],minLL[j][l]);
					hDiffNtt_vs_MinLL_[i][j][k]->Fill(PseudoExp_DiffEstNttlike[i][j][k][l],minLL[j][l]);
				}
			}
		}
		for(int j=0;j<NbOfJetsBins_;j++){
			for(int k=0; k<NbOfPseudoExp_; k++){
				hNv_[i][j][NbOfBJetsBins_] ->Fill(PseudoExp_EstNvlike[i][j][NbOfBJetsBins_][k]);
				hNvb_[i][j][NbOfBJetsBins_]->Fill(PseudoExp_EstNvblike[i][j][NbOfBJetsBins_][k]);
				hNtt_[i][j][NbOfBJetsBins_]->Fill(PseudoExp_EstNttlike[i][j][NbOfBJetsBins_][k]);
				hPullNv_[i][j][NbOfBJetsBins_] ->Fill(PseudoExp_PullEstNvlike[i][j][NbOfBJetsBins_][k]);
				hPullNvb_[i][j][NbOfBJetsBins_]->Fill(PseudoExp_PullEstNvblike[i][j][NbOfBJetsBins_][k]);
				hPullNtt_[i][j][NbOfBJetsBins_]->Fill(PseudoExp_PullEstNttlike[i][j][NbOfBJetsBins_][k]);
				hDiffNv_[i][j][NbOfBJetsBins_] ->Fill(PseudoExp_DiffEstNvlike[i][j][NbOfBJetsBins_][k]);
				hDiffNvb_[i][j][NbOfBJetsBins_]->Fill(PseudoExp_DiffEstNvblike[i][j][NbOfBJetsBins_][k]);
				hDiffNtt_[i][j][NbOfBJetsBins_]->Fill(PseudoExp_DiffEstNttlike[i][j][NbOfBJetsBins_][k]);
				hDiffNv_vs_MinLL_[i][j][NbOfBJetsBins_] ->Fill(PseudoExp_DiffEstNvlike[i][j][NbOfBJetsBins_][k],minLL[j][k]);
				hDiffNvb_vs_MinLL_[i][j][NbOfBJetsBins_]->Fill(PseudoExp_DiffEstNvblike[i][j][NbOfBJetsBins_][k],minLL[j][k]);
				hDiffNtt_vs_MinLL_[i][j][NbOfBJetsBins_]->Fill(PseudoExp_DiffEstNttlike[i][j][NbOfBJetsBins_][k],minLL[j][k]);
			}
		}
	
		for(int k=0;k<NbOfBJetsBins_;k++){
			for(int l=0; l<NbOfPseudoExp_; l++){
				hNv_[i][NbOfJetsBins_][k] ->Fill(PseudoExp_EstNvlike[i][NbOfJetsBins_][k][l]);
				hNvb_[i][NbOfJetsBins_][k]->Fill(PseudoExp_EstNvblike[i][NbOfJetsBins_][k][l]);
				hNtt_[i][NbOfJetsBins_][k]->Fill(PseudoExp_EstNttlike[i][NbOfJetsBins_][k][l]);
				hPullNv_[i][NbOfJetsBins_][k] ->Fill(PseudoExp_PullEstNvlike[i][NbOfJetsBins_][k][l]);
				hPullNvb_[i][NbOfJetsBins_][k]->Fill(PseudoExp_PullEstNvblike[i][NbOfJetsBins_][k][l]);
				hPullNtt_[i][NbOfJetsBins_][k]->Fill(PseudoExp_PullEstNttlike[i][NbOfJetsBins_][k][l]);
				hDiffNv_[i][NbOfJetsBins_][k] ->Fill(PseudoExp_DiffEstNvlike[i][NbOfJetsBins_][k][l]);
				hDiffNvb_[i][NbOfJetsBins_][k]->Fill(PseudoExp_DiffEstNvblike[i][NbOfJetsBins_][k][l]);
				hDiffNtt_[i][NbOfJetsBins_][k]->Fill(PseudoExp_DiffEstNttlike[i][NbOfJetsBins_][k][l]);
			}
		}
	}

	for(int i=0;i<NbOfJetsBins_;i++) delete [] minLL[i];
	delete [] minLL;
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		for(int j=0;j<NbOfJetsBins_;j++){
			for(int k=0;k<NbOfBJetsBins_;k++){
				delete [] tempNbjets[i][j][k];
			}
			for(int k=0;k<NbOfBJetsBins_+1;k++){
				delete [] PseudoExp_EstNvlike[i][j][k];
				delete [] PseudoExp_EstNvblike[i][j][k];
				delete [] PseudoExp_EstNttlike[i][j][k];
				delete [] PseudoExp_PullEstNvlike[i][j][k];
				delete [] PseudoExp_PullEstNvblike[i][j][k];
				delete [] PseudoExp_PullEstNttlike[i][j][k];
				delete [] PseudoExp_DiffEstNvlike[i][j][k];
				delete [] PseudoExp_DiffEstNvblike[i][j][k];
				delete [] PseudoExp_DiffEstNttlike[i][j][k];
			}
			delete [] PseudoExp_EstNvlike[i][j];
			delete [] PseudoExp_EstNvblike[i][j];
			delete [] PseudoExp_EstNttlike[i][j];
			delete [] PseudoExp_PullEstNvlike[i][j];
			delete [] PseudoExp_PullEstNvblike[i][j];
			delete [] PseudoExp_PullEstNttlike[i][j];
			delete [] PseudoExp_DiffEstNvlike[i][j];
			delete [] PseudoExp_DiffEstNvblike[i][j];
			delete [] PseudoExp_DiffEstNttlike[i][j];

			delete [] PseudoExp_EstEb[i][j];
			delete [] PseudoExp_EstEudsc[i][j];
			delete [] PseudoExp_EstEuds[i][j];
			delete [] PseudoExp_PullEstEb[i][j];
			delete [] PseudoExp_PullEstEudsc[i][j];
			delete [] PseudoExp_PullEstEuds[i][j];
			delete [] PseudoExp_DiffEstEb[i][j];
			delete [] PseudoExp_DiffEstEudsc[i][j];
			delete [] PseudoExp_DiffEstEuds[i][j];
			delete [] tempNbjets[i][j];
		}
		delete [] PseudoExp_EstNvlike[i];
		delete [] PseudoExp_EstNvblike[i];
		delete [] PseudoExp_EstNttlike[i];
		delete [] PseudoExp_PullEstNvlike[i];
		delete [] PseudoExp_PullEstNvblike[i];
		delete [] PseudoExp_PullEstNttlike[i];
		delete [] PseudoExp_DiffEstNvlike[i];
		delete [] PseudoExp_DiffEstNvblike[i];
		delete [] PseudoExp_DiffEstNttlike[i];

		delete [] PseudoExp_EstEb[i];
		delete [] PseudoExp_EstEudsc[i];
		delete [] PseudoExp_EstEuds[i];
		delete [] PseudoExp_PullEstEb[i];
		delete [] PseudoExp_PullEstEudsc[i];
		delete [] PseudoExp_PullEstEuds[i];
		delete [] PseudoExp_DiffEstEb[i];
		delete [] PseudoExp_DiffEstEudsc[i];
		delete [] PseudoExp_DiffEstEuds[i];
		delete [] tempNbjets[i];
	}
	delete [] PseudoExp_EstNvlike;
	delete [] PseudoExp_EstNvblike;
	delete [] PseudoExp_EstNttlike;
	delete [] PseudoExp_PullEstNvlike;
	delete [] PseudoExp_PullEstNvblike;
	delete [] PseudoExp_PullEstNttlike;
	delete [] PseudoExp_DiffEstNvlike;
	delete [] PseudoExp_DiffEstNvblike;
	delete [] PseudoExp_DiffEstNttlike;

	delete [] PseudoExp_EstEb;
	delete [] PseudoExp_EstEudsc;
	delete [] PseudoExp_EstEuds;
	delete [] PseudoExp_PullEstEb;
	delete [] PseudoExp_PullEstEudsc;
	delete [] PseudoExp_PullEstEuds;
	delete [] PseudoExp_DiffEstEb;
	delete [] PseudoExp_DiffEstEudsc;
	delete [] PseudoExp_DiffEstEuds;
	delete [] tempNbjets;
}

/**________________________________________________________________________________________________________________*/
void VJetEstPseudoExp::ResetInitialValues(){
	vjEst_->FillInputs(initNbjets_);
}

/**________________________________________________________________________________________________________________*/
void VJetEstPseudoExp::SetInitNbjets(double**** initNbjets){
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		initNbjets_[i] = new double**[NbOfJetsBins_];
		for(int j=0; j<NbOfJetsBins_; j++){
			initNbjets_[i][j] = new double*[NbOfBJetsBins_];
			for(int k=0; k<NbOfBJetsBins_; k++){
			initNbjets_[i][j][k] = new double[NbOfDataSets_];
				for(int l=0; l<NbOfDataSets_; l++){
					initNbjets_[i][j][k][l] = initNbjets[i][j][k][l];
				}
			}
		}
	}
}

/**________________________________________________________________________________________________________________*/
void VJetEstPseudoExp::GetPseudoExperiments(double**** randN0bjets, double**** randN1bjets, double**** randN2bjets, double**** randN3bjets){
	randN0bjets = randN0bjets_;
	randN1bjets = randN1bjets_;
	randN2bjets = randN2bjets_;
	randN3bjets = randN3bjets_;
}


/**________________________________________________________________________________________________________________*/
double VJetEstPseudoExp::GetNvBias(int wp, bool doFit, double xlow, double xhigh) const{
	if(!hNv_) return -1;
	double Bias = 0;
	for(int i=0;i<NbOfJetsBins_;i++) Bias += pow(GetNvBias(wp,i,doFit,xlow,xhigh),2);
	return sqrt(Bias);
}
/**________________________________________________________________________________________________________________*/
double VJetEstPseudoExp::GetNvBias(int wp, int ijets, bool doFit, double xlow, double xhigh) const{
	if(!hNv_[wp][ijets][NbOfBJetsBins_]) return -1;
	if(doFit){
		if(    hDiffNv_[wp][ijets][NbOfBJetsBins_]->GetFunction("gaus")==0) hDiffNv_[wp][ijets][NbOfBJetsBins_]->Fit("gaus","IREM","",xlow,xhigh);
		return hDiffNv_[wp][ijets][NbOfBJetsBins_]->GetFunction("gaus")->GetParameter(1);
	}
	else    return hDiffNv_[wp][ijets][NbOfBJetsBins_]->GetMean();
}
/**________________________________________________________________________________________________________________*/
double VJetEstPseudoExp::GetNvBias(int wp, int ijets, int ibjets, bool doFit, double xlow, double xhigh) const{
	if(!hNv_[wp][ijets][ibjets]) return -1;
	if(doFit){
		if(    hDiffNv_[wp][ijets][ibjets]->GetFunction("gaus")==0) hDiffNv_[wp][ijets][ibjets]->Fit("gaus","IREM","",xlow,xhigh);
		return hDiffNv_[wp][ijets][ibjets]->GetFunction("gaus")->GetParameter(1);
	}
	else    return hDiffNv_[wp][ijets][ibjets]->GetMean();
}

/**________________________________________________________________________________________________________________*/
double VJetEstPseudoExp::GetNttBias(int wp, bool doFit, double xlow, double xhigh) const{
	if(!hNtt_) return -1;
	double Bias = 0;
	for(int i=0;i<NbOfJetsBins_;i++) Bias += pow(GetNttBias(wp,i,doFit,xlow,xhigh),2);
	return sqrt(Bias);
}
/**________________________________________________________________________________________________________________*/
double VJetEstPseudoExp::GetNttBias(int wp, int ijets, bool doFit, double xlow, double xhigh) const{
	if(!hNtt_[wp][ijets][NbOfBJetsBins_]) return -1;
	if(doFit){
		if(    hDiffNtt_[wp][ijets][NbOfBJetsBins_]->GetFunction("gaus")==0) hDiffNtt_[wp][ijets][NbOfBJetsBins_]->Fit("gaus","IREM","",xlow,xhigh);
		return hDiffNtt_[wp][ijets][NbOfBJetsBins_]->GetFunction("gaus")->GetParameter(1);
	}
	else    return hDiffNtt_[wp][ijets][NbOfBJetsBins_]->GetMean();
}
/**________________________________________________________________________________________________________________*/
double VJetEstPseudoExp::GetNttBias(int wp, int ijets, int ibjets, bool doFit, double xlow, double xhigh) const{
	if(!hNtt_[wp][ijets][ibjets]) return -1;
	if(doFit){
		if(    hDiffNtt_[wp][ijets][ibjets]->GetFunction("gaus")==0) hDiffNtt_[wp][ijets][ibjets]->Fit("gaus","IREM","",xlow,xhigh);
		return hDiffNtt_[wp][ijets][ibjets]->GetFunction("gaus")->GetParameter(1);
	}
	else    return hDiffNtt_[wp][ijets][ibjets]->GetMean();
}

/**________________________________________________________________________________________________________________*/
double VJetEstPseudoExp::GetNvVar(int wp, bool doFit, double xlow, double xhigh) const{
	if(!hNv_) return -1;
	double Var = 0;
	for(int i=0;i<NbOfJetsBins_;i++) Var += pow(GetNvVar(wp,i,doFit,xlow,xhigh),2);
	return sqrt(Var);
}
/**________________________________________________________________________________________________________________*/
double VJetEstPseudoExp::GetNvVar(int wp, int ijets, bool doFit, double xlow, double xhigh) const{
	if(!hNv_[wp][ijets][NbOfBJetsBins_]) return -1;
	if(doFit){
		if(    hDiffNv_[wp][ijets][NbOfBJetsBins_]->GetFunction("gaus")==0) hDiffNv_[wp][ijets][NbOfBJetsBins_]->Fit("gaus","IREM","",xlow,xhigh);
		return hDiffNv_[wp][ijets][NbOfBJetsBins_]->GetFunction("gaus")->GetParameter(2);
	}
	else    return hDiffNv_[wp][ijets][NbOfBJetsBins_]->GetRMS();
}
/**________________________________________________________________________________________________________________*/
double VJetEstPseudoExp::GetNvVar(int wp, int ijets, int ibjets, bool doFit, double xlow, double xhigh) const{
	if(!hNv_[wp][ijets][ibjets]) return -1;
	if(doFit){
		if(    hDiffNv_[wp][ijets][ibjets]->GetFunction("gaus")==0) hDiffNv_[wp][ijets][ibjets]->Fit("gaus","IREM","",xlow,xhigh);
		return hDiffNv_[wp][ijets][ibjets]->GetFunction("gaus")->GetParameter(2);
	}
	else    return hDiffNv_[wp][ijets][ibjets]->GetRMS();
}
/**________________________________________________________________________________________________________________*/

double VJetEstPseudoExp::GetNttVar(int wp, bool doFit, double xlow, double xhigh) const{
	double Var = 0;
	for(int i=0;i<NbOfJetsBins_;i++) Var += pow(GetNttVar(wp,i,doFit,xlow,xhigh),2);
	return sqrt(Var);
}
/**________________________________________________________________________________________________________________*/
double VJetEstPseudoExp::GetNttVar(int wp, int ijets, bool doFit, double xlow, double xhigh) const{
	if(!hNtt_[wp][ijets][NbOfBJetsBins_]) return -1;
	if(doFit){
		if(    hDiffNtt_[wp][ijets][NbOfBJetsBins_]->GetFunction("gaus")==0) hDiffNtt_[wp][ijets][NbOfBJetsBins_]->Fit("gaus","IREM","",xlow,xhigh);
		return hDiffNtt_[wp][ijets][NbOfBJetsBins_]->GetFunction("gaus")->GetParameter(2);
	}
	else    return hDiffNtt_[wp][ijets][NbOfBJetsBins_]->GetRMS();
}
/**________________________________________________________________________________________________________________*/
double VJetEstPseudoExp::GetNttVar(int wp, int ijets, int ibjets, bool doFit, double xlow, double xhigh) const{
	if(!hNtt_[wp][ijets][ibjets]) return -1;
	if(doFit){
		if(    hDiffNtt_[wp][ijets][ibjets]->GetFunction("gaus")==0) hDiffNtt_[wp][ijets][ibjets]->Fit("gaus","IREM","",xlow,xhigh);
		return hDiffNtt_[wp][ijets][ibjets]->GetFunction("gaus")->GetParameter(2);
	}
	else    return hDiffNtt_[wp][ijets][ibjets]->GetRMS();
}




/**________________________________________________________________________________________________________________*/
TH1F* VJetEstPseudoExp::Get_hNv(int wp, int ijets, int ibjets){
	return hNv_[wp][ijets][ibjets];
}
/**________________________________________________________________________________________________________________*/
TH1F* VJetEstPseudoExp::Get_hNv(int wp, int ijets){
	return hNv_[wp][ijets][NbOfBJetsBins_];
}

/**________________________________________________________________________________________________________________*/
TH1F* VJetEstPseudoExp::Get_hNv(int wp){
	return hNv_[wp][NbOfJetsBins_][NbOfBJetsBins_];
}

/**________________________________________________________________________________________________________________*/
TH1F* VJetEstPseudoExp::Get_hDiffNv(int wp, int ijets, int ibjets){
	return hDiffNv_[wp][ijets][ibjets];
}
/**________________________________________________________________________________________________________________*/
TH1F* VJetEstPseudoExp::Get_hDiffNv(int wp, int ijets){
	return hDiffNv_[wp][ijets][NbOfBJetsBins_];
}

/**________________________________________________________________________________________________________________*/
TH1F* VJetEstPseudoExp::Get_hDiffNv(int wp){
	return hDiffNv_[wp][NbOfJetsBins_][NbOfBJetsBins_];
}

/**________________________________________________________________________________________________________________*/
TH1F* VJetEstPseudoExp::Get_hNtt(int wp, int ijets, int ibjets){
	return hNtt_[wp][ijets][ibjets];
}
/**________________________________________________________________________________________________________________*/
TH1F* VJetEstPseudoExp::Get_hNtt(int wp, int ijets){
	return hNtt_[wp][ijets][NbOfBJetsBins_];
}
/**________________________________________________________________________________________________________________*/
TH1F* VJetEstPseudoExp::Get_hNtt(int wp){
	return hDiffNtt_[wp][NbOfJetsBins_][NbOfBJetsBins_];
}

/**________________________________________________________________________________________________________________*/
TH1F* VJetEstPseudoExp::Get_hDiffNtt(int wp, int ijets, int ibjets){
	return hDiffNtt_[wp][ijets][ibjets];
}
/**________________________________________________________________________________________________________________*/
TH1F* VJetEstPseudoExp::Get_hDiffNtt(int wp, int ijets){
	return hDiffNtt_[wp][ijets][NbOfBJetsBins_];
}

/**________________________________________________________________________________________________________________*/
TH1F* VJetEstPseudoExp::Get_hDiffNtt(int wp){
	return hDiffNtt_[wp][NbOfJetsBins_][NbOfBJetsBins_];
}

/**________________________________________________________________________________________________________________*/
void VJetEstPseudoExp::Write(TFile* file, string label){
	char name[100];
	file->cd();
	string dirname = "VJetEstPseudoExp"+label;
	TDirectory *Mydir = file->mkdir(dirname.c_str());
	file->cd(dirname.c_str());
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		sprintf(name,"WorkingPoint_%d",i);
		Mydir->mkdir(name);
		Mydir->cd();
		for(int j=0;j<NbOfJetsBins_;j++){
			if(!hEb_[i][j]       == 0) hEb_[i][j]       ->Write();
			if(!hEudsc_[i][j]    == 0) hEudsc_[i][j]    ->Write();
			if(!hEuds_[i][j]     == 0) hEuds_[i][j]     ->Write();
			if(!hPullEb_[i][j]   == 0) hPullEb_[i][j]   ->Write();
			if(!hPullEudsc_[i][j]== 0) hPullEudsc_[i][j]->Write();
			if(!hPullEuds_[i][j] == 0) hPullEuds_[i][j] ->Write();
			if(!hDiffEb_[i][j]   == 0) hDiffEb_[i][j]   ->Write();
			if(!hDiffEudsc_[i][j]== 0) hDiffEudsc_[i][j]->Write();
			if(!hDiffEuds_[i][j] == 0) hDiffEuds_[i][j] ->Write();
			for(int k=0; k<NbOfDataSets_; k++){
				if(!hRandN0bjets_[i][j][k]==0)         hRandN0bjets_[i][j][k]->Write();
				if(!hRandN1bjets_[i][j][k]==0)         hRandN1bjets_[i][j][k]->Write();
				if(!hRandN2bjets_[i][j][k]==0)         hRandN2bjets_[i][j][k]->Write();
				if(!hRandN3bjets_[i][j][k]==0)         hRandN3bjets_[i][j][k]->Write();
				if(!hRandNbjets_vs_MinLL_[i][j][k]==0) hRandNbjets_vs_MinLL_[i][j][k]->Write();
			}
		}
		for(int j=0;j<=NbOfJetsBins_;j++){
			for(int k=0;k<=NbOfBJetsBins_;k++){
			if(j==NbOfJetsBins_ && k==NbOfBJetsBins_) continue;
			if(!hNv_[i][j][k]==0)  hNv_[i][j][k] ->Write();
			if(!hNtt_[i][j][k]==0) hNtt_[i][j][k]->Write();
			if(!hPullNv_[i][j][k]==0)  hPullNv_[i][j][k] ->Write();
			if(!hDiffNv_[i][j][k]==0)  hDiffNv_[i][j][k] ->Write();
			if(!hDiffNv_vs_MinLL_[i][j][k]==0) hDiffNv_vs_MinLL_[i][j][k]->Write();
			if(!hPullNtt_[i][j][k]==0) hPullNtt_[i][j][k]->Write();
			if(!hDiffNtt_[i][j][k]==0) hDiffNtt_[i][j][k]->Write();
			if(!hDiffNtt_vs_MinLL_[i][j][k]==0)hDiffNtt_vs_MinLL_[i][j][k]->Write();
			}
		}
	}
}

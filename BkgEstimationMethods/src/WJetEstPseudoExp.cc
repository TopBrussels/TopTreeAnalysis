#include "../interface/WJetEstPseudoExp.h"

WJetEstPseudoExp::WJetEstPseudoExp(int NbOfPseudoExp, WJetEstimation *wjEst){

	NbOfPseudoExp_    = NbOfPseudoExp;
	wjEst_            = wjEst;
	NbOfDataSets_     = wjEst->GetNbOfDatasets();

	NbOfSuccessfulPE_ = 0;
	rand              = new TRandom3(0);
	
	multijets_est_N0_ = wjEst_->GetMultiJet_Est_N0();
	multijets_est_N1_ = wjEst_->GetMultiJet_Est_N1();
	multijets_est_N2_ = wjEst_->GetMultiJet_Est_N2();
	multijets_est_N3_ = wjEst_->GetMultiJet_Est_N3();
	
	randN0bjets_ = new double**[NbOfPseudoExp];
	randN1bjets_ = new double**[NbOfPseudoExp];
	randN2bjets_ = new double**[NbOfPseudoExp];
	randN3bjets_ = new double**[NbOfPseudoExp];

	hEudsc  = new TH1F*[3];
	hEb     = new TH1F*[3];
	hNtt_Njets = new TH1F*[3];
	hNw_Njets  = new TH1F*[3];
	hNtt_Nbjets= new TH1F*[4];
	hNw_Nbjets = new TH1F*[4];
	hNtt_Incl  = 0;
	hNw_Incl   = 0;
	hNtotal = new TH1F*[3];

	hRandN0bjets = new TH1F**[3];
	hRandN1bjets = new TH1F**[3];
	hRandN2bjets = new TH1F**[3];
	hRandN3bjets = new TH1F**[3];
	for(int i=0;i<3;i++)
	{
		hRandN0bjets[i] = new TH1F*[2];
		hRandN1bjets[i] = new TH1F*[2];
		hRandN2bjets[i] = new TH1F*[2];
		hRandN3bjets[i] = new TH1F*[2];
	}
	
	for(int i=0; i<NbOfPseudoExp; i++){
		randN0bjets_[i] = new double*[3];
		randN1bjets_[i] = new double*[3];
		randN2bjets_[i] = new double*[3];
		randN3bjets_[i] = new double*[3];
		for(int j=0;j<3;j++){
			randN0bjets_[i][j] = new double[NbOfDataSets_];
			randN1bjets_[i][j] = new double[NbOfDataSets_];
			randN2bjets_[i][j] = new double[NbOfDataSets_];
			randN3bjets_[i][j] = new double[NbOfDataSets_];
			for(int k=0;k<NbOfDataSets_;k++){
				randN0bjets_[i][j][k] = 0;
				randN1bjets_[i][j][k] = 0;
				randN2bjets_[i][j][k] = 0;
				randN3bjets_[i][j][k] = 0;
			}			
		}
	}
}

WJetEstPseudoExp::~WJetEstPseudoExp(){	
	for(int i=0; i<NbOfPseudoExp_; i++){
		for(int j=0;j<3;j++){
			delete [] randN0bjets_[i][j];
			delete [] randN1bjets_[i][j];
			delete [] randN2bjets_[i][j];
			delete [] randN3bjets_[i][j];
		}
		delete [] randN0bjets_[i];
		delete [] randN1bjets_[i];
		delete [] randN2bjets_[i];
		delete [] randN3bjets_[i];
	}

	delete [] randN0bjets_;
	delete [] randN1bjets_;
	delete [] randN2bjets_;
	delete [] randN3bjets_;

	for(int i=0; i<3 ; i++){
		delete hEudsc[i];
		delete hEb[i];
		delete hNtt_Njets[i];
		delete hNw_Njets[i];
		delete hNtotal[i];

		for(int j=0; j<2; j++){
			delete hRandN0bjets[i][j];
			delete hRandN1bjets[i][j];
			delete hRandN2bjets[i][j];
			delete hRandN3bjets[i][j];
		}
	}

	for(int i=0; i<4 ; i++){
		delete hNtt_Nbjets[i];
		delete hNw_Nbjets[i];
	}
	delete hNtt_Incl;
	delete hNw_Incl;

	delete [] hEudsc;
	delete [] hEb;
	delete [] hNtt_Njets;
	delete [] hNw_Njets;
	delete [] hNtt_Nbjets;
	delete [] hNw_Nbjets;
	
	delete [] hRandN0bjets;
	delete [] hRandN1bjets;
	delete [] hRandN2bjets;
	delete [] hRandN3bjets;

	//delete wjEst_;
	//commented to avoid crash
}

void WJetEstPseudoExp::PrintResults(){
	cout<<" --------------------------------------------------"<<endl;
	cout<<" ------------ W/Z background estimation -----------"<<endl;
	cout<<" -- Number of pseudo-experiments : "<<NbOfPseudoExp_<<endl;
	cout<<" -- Number of successful P.E. : "<<NbOfSuccessfulPE_<<endl;
	cout<<" -- Statistical errors calculated from the pseudo-experience :"<<endl;
	for(int i = 0; i<3; i++){
		cout<<" --- estimates with "<<i+4<<" jets :"<<endl;
		cout<<" ---- W-like estimate :"<<endl;
		cout<<" ----- stat. error (RMS) : "<<GetNWError(i,false)<<endl;
		cout<<" ----- syst. error (Mean) : "<<GetNWBias(i,false)<<endl;
		cout<<" ---- TT-like estimate :"<<endl;
		cout<<" ----- stat. error (RMS) : "<<GetNTTError(i,false)<<endl;
		cout<<" ----- syst. error (Mean) : "<<GetNTTBias(i,false)<<endl;
		cout<<" ----------------------"<<endl;
	}
	cout<<" --- inclusive estimates :"<<endl;
	cout<<" ---- W-like estimate :"<<endl;
	cout<<" ----- stat. error (RMS) : "<< GetNWError(false)<<endl;
	cout<<" ----- syst. error (Mean) : "<<GetNWBias(false)<<endl;
	cout<<" ---- TT-like estimate :"<<endl;
	cout<<" ----- stat. error (RMS) : "<< GetNTTError(false)<<endl;
	cout<<" ----- syst. error (Mean) : "<<GetNTTBias(false)<<endl;
	cout<<" --------------------------------------------------"<<endl;
	cout<<"\\begin{table}"<<endl;
	cout<<"	\\centering"<<endl;
	cout<<"	 \\begin{tabular}{|l||c|c|c|c|}"<<endl;
	cout<<"	\\hline"<<endl;
	cout<<"\\multicolumn{5}{|l}{W-like estimate : \\\\}"<<endl;
	cout<<" & $4$ jets & $5$ jets & $6$ jets & inclusive \\\\"<<endl;
	cout<<"Statistical error &"<<GetNWError(0,false)<<"&"<<GetNWError(1,false)<<"&"<<GetNWError(2,false)<<"&"<<GetNWError(false)<<"\\\\"<<endl;
	cout<<"Systematic error  &"<<GetNWBias( 0,false)<<"&"<<GetNWBias( 1,false)<<"&"<<GetNWBias( 2,false)<<"&"<<GetNWBias( false)<<"\\\\"<<endl;
	cout<<"	\\hline"<<endl;
	cout<<"\\multicolumn{5}{|l}{$\\ttbar$-like estimate : \\\\}"<<endl;
	cout<<" & $4$ jets & $5$ jets & $6$ jets & inclusive \\\\"<<endl;
	cout<<"Statistical error &"<<GetNTTError(0,false)<<"&"<<GetNTTError(1,false)<<"&"<<GetNTTError(2,false)<<"&"<<GetNTTError(false)<<"\\\\"<<endl;
	cout<<"Systematic error  &"<<GetNTTBias( 0,false)<<"&"<<GetNTTBias( 1,false)<<"&"<<GetNTTBias( 2,false)<<"&"<<GetNTTBias( false)<<"\\\\"<<endl;
	cout<<"	\\hline"<<endl;
	cout<<"  \\end{tabular}"<<endl;
	cout<<"  \\caption{}"<<endl;
	cout<<"  \\label{tab:}"<<endl;
	cout<<"\\end{table}"<<endl;

}

void WJetEstPseudoExp::RollTheDice(bool Print){
	
	//double InitValues[4][3][NbOfDataSets_];
	double **InitValues[4];
	InitValues[0] = new double*[3];
	InitValues[1] = new double*[3];
	InitValues[2] = new double*[3];
	InitValues[3] = new double*[3];
	for(int i=0;i<3;i++){
		InitValues[0][i] = new double[NbOfDataSets_];
		InitValues[1][i] = new double[NbOfDataSets_];
		InitValues[2][i] = new double[NbOfDataSets_];
		InitValues[3][i] = new double[NbOfDataSets_];
		for(int j=0;j<NbOfDataSets_;j++){
			InitValues[0][i][j]= 0;
			InitValues[1][i][j]= 0;
			InitValues[2][i][j]= 0;
			InitValues[3][i][j]= 0;
		}
	}
	double PseudoExpResults[NbOfPseudoExp_][3][3];
	double PseudoExpResults_Wlike[NbOfPseudoExp_][3][4];
	double PseudoExpResults_WlikeIncl = 0;
	double PseudoExpResults_TTlike[NbOfPseudoExp_][3][4];
	double PseudoExpResults_TTlikeIncl = 0;
	double PseudoExpResults_Effb[NbOfPseudoExp_][3];
	double PseudoExpResults_Effudsc[NbOfPseudoExp_][3];

	double EstMinMax[3][3][2];
	double EstMinMax_Wlike[3][4][2];
	double EstMinMax_TTlike[3][4][2];
	double EffMinMax[3][2][2];

	double MCtruth[3][4][2];
	bool   Verbose = false;
	
	NbOfSuccessfulPE_ = 0;

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<NbOfDataSets_;j++)
		{
			InitValues[0][i][j] = wjEst_->GetN0()[i][j];
			InitValues[1][i][j] = wjEst_->GetN1()[i][j];
			InitValues[2][i][j] = wjEst_->GetN2()[i][j];
			InitValues[3][i][j] = wjEst_->GetN3()[i][j];
		}
		
		EstMinMax[i][0][0] = 99999; EstMinMax[i][0][1] = -99999;
		EstMinMax[i][1][0] = 99999; EstMinMax[i][1][1] = -99999;
		EstMinMax[i][2][0] = 99999; EstMinMax[i][2][1] = -99999;
		EffMinMax[i][0][0] = 99999; EffMinMax[i][0][1] = -99999;
		EffMinMax[i][1][0] = 99999; EffMinMax[i][1][1] = -99999;
		
		for(int j=0;j<4;j++)
		{
			EstMinMax_Wlike[i][j][0]  = 99999; EstMinMax_Wlike[i][j][1]  = -99999;
			EstMinMax_TTlike[i][j][0] = 99999; EstMinMax_TTlike[i][j][1] = -99999;
		}
	}
	std::cout<<"// Let's start the RollTheDice method..."<<std::endl;
	for(int i=0; i<NbOfPseudoExp_; i++)
	{
		if(Print && i%500==0) std::cout<<"RollTheDice method : running the "<<i+1<<"th pseudo-exp"<<std::endl;
		for(int j=0;j<3;j++)
		{
 			for(int k=0;k<NbOfDataSets_;k++)
			{
				randN0bjets_[i][j][k] = rand->PoissonD(InitValues[0][j][k]);
				randN1bjets_[i][j][k] = rand->PoissonD(InitValues[1][j][k]);
				randN2bjets_[i][j][k] = rand->PoissonD(InitValues[2][j][k]);
				randN3bjets_[i][j][k] = rand->PoissonD(InitValues[3][j][k]);
			}
		}
		
		wjEst_->FillInputs(randN0bjets_[i],randN1bjets_[i],randN2bjets_[i],randN3bjets_[i]);
		wjEst_->Estimation(Verbose,false);
		if(wjEst_->CheckEstimation(false)) NbOfSuccessfulPE_++;

		for(int j=0;j<3;j++)
		{
			for(int k=0;k<4;k++)
			{
				MCtruth[j][k][0] = 0;
				MCtruth[j][k][1] = 0;
			}
			for(unsigned int k=0;k<wjEst_->GetiDatasetsWLike().size();k++)
			{
				MCtruth[j][0][0] += randN0bjets_[i][j][wjEst_->GetiDatasetsWLike()[k]];
				MCtruth[j][1][0] += randN1bjets_[i][j][wjEst_->GetiDatasetsWLike()[k]];
				MCtruth[j][2][0] += randN2bjets_[i][j][wjEst_->GetiDatasetsWLike()[k]];
				MCtruth[j][3][0] += randN3bjets_[i][j][wjEst_->GetiDatasetsWLike()[k]];
			}
			for(unsigned int k=0;k<wjEst_->GetiDatasetsTTLike().size();k++)
			{
				MCtruth[j][0][1] += randN0bjets_[i][j][wjEst_->GetiDatasetsTTLike()[k]];
				MCtruth[j][1][1] += randN1bjets_[i][j][wjEst_->GetiDatasetsTTLike()[k]];
				MCtruth[j][2][1] += randN2bjets_[i][j][wjEst_->GetiDatasetsTTLike()[k]];
				MCtruth[j][3][1] += randN3bjets_[i][j][wjEst_->GetiDatasetsTTLike()[k]];
			}
		}

		for(int j=0;j<3;j++)
		{
			PseudoExpResults_Effudsc[i][j]    = wjEst_->GetEffudscEstimated(j);
			if(EffMinMax[j][0][0]>PseudoExpResults_Effudsc[i][j]) EffMinMax[j][0][0] = PseudoExpResults_Effudsc[i][j];
			if(EffMinMax[j][0][1]<PseudoExpResults_Effudsc[i][j]) EffMinMax[j][0][1] = PseudoExpResults_Effudsc[i][j];

			PseudoExpResults_Effb[i][j] = wjEst_->GetEffbEstimated(j);
			if(EffMinMax[j][1][0]>PseudoExpResults_Effb[i][j]) EffMinMax[j][1][0] = PseudoExpResults_Effb[i][j];
			if(EffMinMax[j][1][1]<PseudoExpResults_Effb[i][j]) EffMinMax[j][1][1] = PseudoExpResults_Effb[i][j];
			
			PseudoExpResults[i][j][0] = (MCtruth[j][0][0]+MCtruth[j][1][0]+MCtruth[j][2][0]+MCtruth[j][3][0])-wjEst_->GetNWEventsEstimated(j);
			if(EstMinMax[j][0][0]>PseudoExpResults[i][j][0]) EstMinMax[j][0][0] = PseudoExpResults[i][j][0];
			if(EstMinMax[j][0][1]<PseudoExpResults[i][j][0]) EstMinMax[j][0][1] = PseudoExpResults[i][j][0];

			PseudoExpResults[i][j][1] = (MCtruth[j][0][1]+MCtruth[j][1][1]+MCtruth[j][2][1]+MCtruth[j][3][1])-wjEst_->GetNttEventsEstimated(j);			
			if(EstMinMax[j][1][0]>PseudoExpResults[i][j][1]) EstMinMax[j][1][0] = PseudoExpResults[i][j][1];
			if(EstMinMax[j][1][1]<PseudoExpResults[i][j][1]) EstMinMax[j][1][1] = PseudoExpResults[i][j][1];

			PseudoExpResults[i][j][2] = wjEst_->NofEvents(j)-wjEst_->GetNtotalEventsEstimated(j);
			if(EstMinMax[j][2][0]>PseudoExpResults[i][j][2]) EstMinMax[j][2][0] = PseudoExpResults[i][j][2];
			if(EstMinMax[j][2][1]<PseudoExpResults[i][j][2]) EstMinMax[j][2][1] = PseudoExpResults[i][j][2];
			
			PseudoExpResults_Wlike[i][j][0]  = MCtruth[j][0][0]-wjEst_->GetNWEventsEstimatedNbjets(0,j);
			PseudoExpResults_Wlike[i][j][1]  = MCtruth[j][1][0]-wjEst_->GetNWEventsEstimatedNbjets(1,j);
			PseudoExpResults_Wlike[i][j][2]  = MCtruth[j][2][0]-wjEst_->GetNWEventsEstimatedNbjets(2,j);
			PseudoExpResults_Wlike[i][j][3]  = MCtruth[j][3][0]-wjEst_->GetNWEventsEstimatedNbjets(3,j);
			for(int k=0;k<4;k++)
			{			
				if(EstMinMax_Wlike[j][k][0]>PseudoExpResults_Wlike[i][j][k]) EstMinMax_Wlike[j][k][0] = PseudoExpResults_Wlike[i][j][k];
				if(EstMinMax_Wlike[j][k][1]<PseudoExpResults_Wlike[i][j][k]) EstMinMax_Wlike[j][k][1] = PseudoExpResults_Wlike[i][j][k];
			}
			
			PseudoExpResults_TTlike[i][j][0] = MCtruth[j][0][1]-wjEst_->GetNttEventsEstimatedNbjets(0,j);
			PseudoExpResults_TTlike[i][j][1] = MCtruth[j][1][1]-wjEst_->GetNttEventsEstimatedNbjets(1,j);
			PseudoExpResults_TTlike[i][j][2] = MCtruth[j][2][1]-wjEst_->GetNttEventsEstimatedNbjets(2,j);
			PseudoExpResults_TTlike[i][j][3] = MCtruth[j][3][1]-wjEst_->GetNttEventsEstimatedNbjets(3,j);
			for(int k=0;k<4;k++)
			{			
				if(EstMinMax_TTlike[j][k][0]>PseudoExpResults_TTlike[i][j][k]) EstMinMax_TTlike[j][k][0] = PseudoExpResults_TTlike[i][j][k];
				if(EstMinMax_TTlike[j][k][1]<PseudoExpResults_TTlike[i][j][k]) EstMinMax_TTlike[j][k][1] = PseudoExpResults_TTlike[i][j][k];
			}
		}
	}

	char name[100];
	for(int i=0; i<3; i++){	
		sprintf(name,"Eudsc_%d_jets_PseudoExp",i+4);
		hEudsc[i] = new TH1F(name,"",250,EffMinMax[i][0][0],EffMinMax[i][0][1]);
		sprintf(name,"Eb_%d_jets_PseudoExp",i+4);
		hEb[i]    = new TH1F(name,"",250,EffMinMax[i][1][0],EffMinMax[i][1][1]);
		sprintf(name,"Nw_%d_jets_PseudoExp",i+4);
		hNw_Njets[i] = new TH1F(name,"",250,EstMinMax[i][0][0],EstMinMax[i][0][1]);
		hNw_Njets[i]->SetStats();
		hNw_Njets[i]->GetXaxis()->SetTitle("#Delta(N^{exp}-N^{est})");
		sprintf(name,"Ntt_%d_jets_PseudoExp",i+4);
		hNtt_Njets[i] = new TH1F(name,"",250,EstMinMax[i][1][0],EstMinMax[i][1][1]);
		hNtt_Njets[i]->SetStats();
		hNtt_Njets[i]->GetXaxis()->SetTitle("#Delta(N^{exp}-N^{est})");
		sprintf(name,"Ntotal_%d_jets_PseudoExp",i+4);
		hNtotal[i] = new TH1F(name,"",250,EstMinMax[i][2][0],EstMinMax[i][2][1]);
		hNtotal[i]->SetStats();
		hNtotal[i]->GetXaxis()->SetTitle("#Delta(N^{exp}-N^{est})");
		
		sprintf(name,"Nw_N0bjets_%d_jets_PseudoExp",i+4);
		hRandN0bjets[i][0] = new TH1F(name,"",250,EstMinMax_Wlike[i][0][0],EstMinMax_Wlike[i][0][1]);
		sprintf(name,"Nw_N1bjets_%d_jets_PseudoExp",i+4);
		hRandN1bjets[i][0] = new TH1F(name,"",250,EstMinMax_Wlike[i][1][0],EstMinMax_Wlike[i][1][1]);
		sprintf(name,"Nw_N2bjets_%d_jets_PseudoExp",i+4);
		hRandN2bjets[i][0] = new TH1F(name,"",250,EstMinMax_Wlike[i][2][0],EstMinMax_Wlike[i][2][1]);
		sprintf(name,"Nw_N3bjets_%d_jets_PseudoExp",i+4);
		hRandN3bjets[i][0] = new TH1F(name,"",250,EstMinMax_Wlike[i][3][0],EstMinMax_Wlike[i][3][1]);

		sprintf(name,"Ntt_N0bjets_%d_jets_PseudoExp",i+4);
		hRandN0bjets[i][1] = new TH1F(name,"",250,EstMinMax_TTlike[i][0][0],EstMinMax_TTlike[i][0][1]);
		sprintf(name,"Ntt_N1bjets_%d_jets_PseudoExp",i+4);
		hRandN1bjets[i][1] = new TH1F(name,"",250,EstMinMax_TTlike[i][1][0],EstMinMax_TTlike[i][1][1]);
		sprintf(name,"Ntt_N2bjets_%d_jets_PseudoExp",i+4);
		hRandN2bjets[i][1] = new TH1F(name,"",250,EstMinMax_TTlike[i][2][0],EstMinMax_TTlike[i][2][1]);
		sprintf(name,"Ntt_N3bjets_%d_jets_PseudoExp",i+4);
		hRandN3bjets[i][1] = new TH1F(name,"",250,EstMinMax_TTlike[i][3][0],EstMinMax_TTlike[i][3][1]);
	}

	double EstMin_Wlike_Nbjets  = 99999; double EstMax_Wlike_Nbjets = -99999;
	double EstMin_TTlike_Nbjets = 99999; double EstMax_TTlike_Nbjets = -99999;

	double EstMin_Wlike_Incl  = 0; double EstMax_Wlike_Incl  = 0;
	double EstMin_TTlike_Incl = 0; double EstMax_TTlike_Incl = 0;

	for(int i=0; i<4; i++){
		EstMin_Wlike_Nbjets  = 99999; EstMax_Wlike_Nbjets  = -99999;
		EstMin_TTlike_Nbjets = 99999; EstMax_TTlike_Nbjets = -99999;
		for(int j=0;j<3;j++){
			if(EstMin_Wlike_Nbjets  < EstMinMax_Wlike[j][i][0])  EstMin_Wlike_Nbjets  = EstMinMax_Wlike[j][i][0];
			if(EstMin_TTlike_Nbjets < EstMinMax_TTlike[j][i][0]) EstMin_TTlike_Nbjets = EstMinMax_TTlike[j][i][0];
			if(EstMax_Wlike_Nbjets  < EstMinMax_Wlike[j][i][1])  EstMax_Wlike_Nbjets  = EstMinMax_Wlike[j][i][1];
			if(EstMax_TTlike_Nbjets < EstMinMax_TTlike[j][i][1]) EstMax_TTlike_Nbjets = EstMinMax_TTlike[j][i][1];
		}
		sprintf(name,"Nw_%d_bjets_PseudoExp",i);
		hNw_Nbjets[i]    = new TH1F(name,"",250,EstMin_Wlike_Nbjets,EstMax_Wlike_Nbjets);
		hNw_Nbjets[i]->SetStats();
		hNw_Nbjets[i]->GetXaxis()->SetTitle("#Delta(N^{exp}-N^{est})");
		sprintf(name,"Ntt_%d_bjets_PseudoExp",i);
		hNtt_Nbjets[i]   = new TH1F(name,"",250,EstMin_TTlike_Nbjets,EstMax_TTlike_Nbjets);
		hNtt_Nbjets[i]->SetStats();
		hNtt_Nbjets[i]->GetXaxis()->SetTitle("#Delta(N^{exp}-N^{est})");

		EstMin_Wlike_Incl  += EstMin_Wlike_Nbjets;
		EstMax_Wlike_Incl  += EstMax_Wlike_Nbjets;
		EstMin_TTlike_Incl += EstMin_TTlike_Nbjets;
		EstMax_TTlike_Incl += EstMax_TTlike_Nbjets;
	}

	sprintf(name,"Nw_Incl_PseudoExp");
	hNw_Incl  = new TH1F(name,"",250,EstMin_Wlike_Incl,EstMax_Wlike_Incl);
	hNw_Incl->SetStats();
	hNw_Incl->GetXaxis()->SetTitle("#Delta(N^{exp}-N^{est})");
	sprintf(name,"Ntt_Incl_PseudoExp");
	hNtt_Incl = new TH1F(name,"",250,EstMin_TTlike_Incl,EstMax_TTlike_Incl);
	hNtt_Incl->SetStats();
	hNtt_Incl->GetXaxis()->SetTitle("#Delta(N^{exp}-N^{est})");

	for(int i=0; i<NbOfPseudoExp_; i++){
		PseudoExpResults_WlikeIncl  = 0;
		PseudoExpResults_TTlikeIncl = 0;
		for(int j=0;j<3;j++){
			hEudsc[j]    ->Fill(PseudoExpResults_Effudsc[i][j]);
			hEb[j]       ->Fill(PseudoExpResults_Effb[i][j]);
			hNw_Njets[j] ->Fill(PseudoExpResults[i][j][0]);
			hNtt_Njets[j]->Fill(PseudoExpResults[i][j][1]);
			hNtotal[j]   ->Fill(PseudoExpResults[i][j][2]);
			
			hRandN0bjets[j][0]->Fill(PseudoExpResults_Wlike[i][j][0]);
			hRandN1bjets[j][0]->Fill(PseudoExpResults_Wlike[i][j][1]);
			hRandN2bjets[j][0]->Fill(PseudoExpResults_Wlike[i][j][2]);
			hRandN3bjets[j][0]->Fill(PseudoExpResults_Wlike[i][j][3]);

			hRandN0bjets[j][1]->Fill(PseudoExpResults_TTlike[i][j][0]);
			hRandN1bjets[j][1]->Fill(PseudoExpResults_TTlike[i][j][1]);
			hRandN2bjets[j][1]->Fill(PseudoExpResults_TTlike[i][j][2]);
			hRandN3bjets[j][1]->Fill(PseudoExpResults_TTlike[i][j][3]);
		}
		for(int j=0;j<4;j++){
			hNw_Nbjets[j] ->Fill(PseudoExpResults_Wlike[i][0][j] +PseudoExpResults_Wlike[i][1][j] +PseudoExpResults_Wlike[i][2][j]);
			hNtt_Nbjets[j]->Fill(PseudoExpResults_TTlike[i][0][j]+PseudoExpResults_TTlike[i][1][j]+PseudoExpResults_TTlike[i][2][j]);
			PseudoExpResults_WlikeIncl  += PseudoExpResults_Wlike[i][0][j] +PseudoExpResults_Wlike[i][1][j] +PseudoExpResults_Wlike[i][2][j];
			PseudoExpResults_TTlikeIncl += PseudoExpResults_TTlike[i][0][j]+PseudoExpResults_TTlike[i][1][j]+PseudoExpResults_TTlike[i][2][j];
		}
		hNw_Incl ->Fill(PseudoExpResults_WlikeIncl);
		hNtt_Incl->Fill(PseudoExpResults_TTlikeIncl);
	}
	wjEst_->FillInputs(InitValues[0],InitValues[1],InitValues[2],InitValues[3]);
}

void WJetEstPseudoExp::RollTheDice(double **multijets_est, bool Print){
	
	//double InitValues[4][3][NbOfDataSets_];
	double **InitValues[4];
	InitValues[0] = new double*[3];
	InitValues[1] = new double*[3];
	InitValues[2] = new double*[3];
	InitValues[3] = new double*[3];
	for(int i=0;i<3;i++){
		InitValues[0][i] = new double[NbOfDataSets_];
		InitValues[1][i] = new double[NbOfDataSets_];
		InitValues[2][i] = new double[NbOfDataSets_];
		InitValues[3][i] = new double[NbOfDataSets_];
		for(int j=0;j<NbOfDataSets_;j++){
			InitValues[0][i][j]= 0;
			InitValues[1][i][j]= 0;
			InitValues[2][i][j]= 0;
			InitValues[3][i][j]= 0;
		}
	}
	double PseudoExpResults[NbOfPseudoExp_][3][3];
	double PseudoExpResults_Wlike[NbOfPseudoExp_][3][4];
	double PseudoExpResults_WlikeIncl = 0;
	double PseudoExpResults_TTlike[NbOfPseudoExp_][3][4];
	double PseudoExpResults_TTlikeIncl = 0;
	double PseudoExpResults_Effb[NbOfPseudoExp_][3];
	double PseudoExpResults_Effudsc[NbOfPseudoExp_][3];

	double EstMinMax[3][3][2];
	double EstMinMax_Wlike[3][4][2];
	double EstMinMax_TTlike[3][4][2];
	double EffMinMax[3][2][2];

	double MCtruth[3][4][2];
	bool   Verbose = false;
	
	NbOfSuccessfulPE_ = 0;

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<NbOfDataSets_;j++)
		{
			InitValues[0][i][j] = wjEst_->GetN0()[i][j];
			InitValues[1][i][j] = wjEst_->GetN1()[i][j];
			InitValues[2][i][j] = wjEst_->GetN2()[i][j];
			InitValues[3][i][j] = wjEst_->GetN3()[i][j];
		}
		
		EstMinMax[i][0][0] = 99999; EstMinMax[i][0][1] = -99999;
		EstMinMax[i][1][0] = 99999; EstMinMax[i][1][1] = -99999;
		EstMinMax[i][2][0] = 99999; EstMinMax[i][2][1] = -99999;
		EffMinMax[i][0][0] = 99999; EffMinMax[i][0][1] = -99999;
		EffMinMax[i][1][0] = 99999; EffMinMax[i][1][1] = -99999;
		
		for(int j=0;j<4;j++)
		{
			EstMinMax_Wlike[i][j][0]  = 99999; EstMinMax_Wlike[i][j][1]  = -99999;
			EstMinMax_TTlike[i][j][0] = 99999; EstMinMax_TTlike[i][j][1] = -99999;
		}
	}
	std::cout<<"// Let's start the RollTheDice method..."<<std::endl;
	for(int i=0; i<NbOfPseudoExp_; i++)
	{
		if(Print && i%500==0) std::cout<<"RollTheDice method : running the "<<i+1<<"th pseudo-exp"<<std::endl;
		for(int j=0;j<3;j++)
		{
 			for(int k=0;k<NbOfDataSets_;k++)
			{
				randN0bjets_[i][j][k] = rand->PoissonD(InitValues[0][j][k]);
				randN1bjets_[i][j][k] = rand->PoissonD(InitValues[1][j][k]);
				randN2bjets_[i][j][k] = rand->PoissonD(InitValues[2][j][k]);
				randN3bjets_[i][j][k] = rand->PoissonD(InitValues[3][j][k]);
			}
		}
		
		wjEst_->FillInputs(randN0bjets_[i],randN1bjets_[i],randN2bjets_[i],randN3bjets_[i],multijets_est);
		wjEst_->Estimation(Verbose,false);
		if(wjEst_->CheckEstimation(false)) NbOfSuccessfulPE_++;

		for(int j=0;j<3;j++)
		{
			for(int k=0;k<4;k++)
			{
				MCtruth[j][k][0] = 0;
				MCtruth[j][k][1] = 0;
			}
			for(unsigned int k=0;k<wjEst_->GetiDatasetsWLike().size();k++)
			{
				MCtruth[j][0][0] += randN0bjets_[i][j][wjEst_->GetiDatasetsWLike()[k]];
				MCtruth[j][1][0] += randN1bjets_[i][j][wjEst_->GetiDatasetsWLike()[k]];
				MCtruth[j][2][0] += randN2bjets_[i][j][wjEst_->GetiDatasetsWLike()[k]];
				MCtruth[j][3][0] += randN3bjets_[i][j][wjEst_->GetiDatasetsWLike()[k]];
			}
			for(unsigned int k=0;k<wjEst_->GetiDatasetsTTLike().size();k++)
			{
				MCtruth[j][0][1] += randN0bjets_[i][j][wjEst_->GetiDatasetsTTLike()[k]];
				MCtruth[j][1][1] += randN1bjets_[i][j][wjEst_->GetiDatasetsTTLike()[k]];
				MCtruth[j][2][1] += randN2bjets_[i][j][wjEst_->GetiDatasetsTTLike()[k]];
				MCtruth[j][3][1] += randN3bjets_[i][j][wjEst_->GetiDatasetsTTLike()[k]];
			}
		}

		for(int j=0;j<3;j++)
		{
			PseudoExpResults_Effudsc[i][j]    = wjEst_->GetEffudscEstimated(j);
			if(EffMinMax[j][0][0]>PseudoExpResults_Effudsc[i][j]) EffMinMax[j][0][0] = PseudoExpResults_Effudsc[i][j];
			if(EffMinMax[j][0][1]<PseudoExpResults_Effudsc[i][j]) EffMinMax[j][0][1] = PseudoExpResults_Effudsc[i][j];

			PseudoExpResults_Effb[i][j] = wjEst_->GetEffbEstimated(j);
			if(EffMinMax[j][1][0]>PseudoExpResults_Effb[i][j]) EffMinMax[j][1][0] = PseudoExpResults_Effb[i][j];
			if(EffMinMax[j][1][1]<PseudoExpResults_Effb[i][j]) EffMinMax[j][1][1] = PseudoExpResults_Effb[i][j];
			
			PseudoExpResults[i][j][0] = (MCtruth[j][0][0]+MCtruth[j][1][0]+MCtruth[j][2][0]+MCtruth[j][3][0])-wjEst_->GetNWEventsEstimated(j);
			if(EstMinMax[j][0][0]>PseudoExpResults[i][j][0]) EstMinMax[j][0][0] = PseudoExpResults[i][j][0];
			if(EstMinMax[j][0][1]<PseudoExpResults[i][j][0]) EstMinMax[j][0][1] = PseudoExpResults[i][j][0];

			PseudoExpResults[i][j][1] = (MCtruth[j][0][1]+MCtruth[j][1][1]+MCtruth[j][2][1]+MCtruth[j][3][1])-wjEst_->GetNttEventsEstimated(j);			
			if(EstMinMax[j][1][0]>PseudoExpResults[i][j][1]) EstMinMax[j][1][0] = PseudoExpResults[i][j][1];
			if(EstMinMax[j][1][1]<PseudoExpResults[i][j][1]) EstMinMax[j][1][1] = PseudoExpResults[i][j][1];

			PseudoExpResults[i][j][2] = wjEst_->NofEvents(j)-wjEst_->GetNtotalEventsEstimated(j);
			if(EstMinMax[j][2][0]>PseudoExpResults[i][j][2]) EstMinMax[j][2][0] = PseudoExpResults[i][j][2];
			if(EstMinMax[j][2][1]<PseudoExpResults[i][j][2]) EstMinMax[j][2][1] = PseudoExpResults[i][j][2];
			
			PseudoExpResults_Wlike[i][j][0]  = MCtruth[j][0][0]-wjEst_->GetNWEventsEstimatedNbjets(0,j);
			PseudoExpResults_Wlike[i][j][1]  = MCtruth[j][1][0]-wjEst_->GetNWEventsEstimatedNbjets(1,j);
			PseudoExpResults_Wlike[i][j][2]  = MCtruth[j][2][0]-wjEst_->GetNWEventsEstimatedNbjets(2,j);
			PseudoExpResults_Wlike[i][j][3]  = MCtruth[j][3][0]-wjEst_->GetNWEventsEstimatedNbjets(3,j);
			for(int k=0;k<4;k++)
			{			
				if(EstMinMax_Wlike[j][k][0]>PseudoExpResults_Wlike[i][j][k]) EstMinMax_Wlike[j][k][0] = PseudoExpResults_Wlike[i][j][k];
				if(EstMinMax_Wlike[j][k][1]<PseudoExpResults_Wlike[i][j][k]) EstMinMax_Wlike[j][k][1] = PseudoExpResults_Wlike[i][j][k];
			}
			
			PseudoExpResults_TTlike[i][j][0] = MCtruth[j][0][1]-wjEst_->GetNttEventsEstimatedNbjets(0,j);
			PseudoExpResults_TTlike[i][j][1] = MCtruth[j][1][1]-wjEst_->GetNttEventsEstimatedNbjets(1,j);
			PseudoExpResults_TTlike[i][j][2] = MCtruth[j][2][1]-wjEst_->GetNttEventsEstimatedNbjets(2,j);
			PseudoExpResults_TTlike[i][j][3] = MCtruth[j][3][1]-wjEst_->GetNttEventsEstimatedNbjets(3,j);
			for(int k=0;k<4;k++)
			{			
				if(EstMinMax_TTlike[j][k][0]>PseudoExpResults_TTlike[i][j][k]) EstMinMax_TTlike[j][k][0] = PseudoExpResults_TTlike[i][j][k];
				if(EstMinMax_TTlike[j][k][1]<PseudoExpResults_TTlike[i][j][k]) EstMinMax_TTlike[j][k][1] = PseudoExpResults_TTlike[i][j][k];
			}
		}
	}

	char name[100];
	for(int i=0; i<3; i++){	
		sprintf(name,"Eudsc_%d_jets_PseudoExp",i+4);
		hEudsc[i] = new TH1F(name,"",250,EffMinMax[i][0][0],EffMinMax[i][0][1]);
		sprintf(name,"Eb_%d_jets_PseudoExp",i+4);
		hEb[i]    = new TH1F(name,"",250,EffMinMax[i][1][0],EffMinMax[i][1][1]);
		sprintf(name,"Nw_%d_jets_PseudoExp",i+4);
		hNw_Njets[i] = new TH1F(name,"",250,EstMinMax[i][0][0],EstMinMax[i][0][1]);
		hNw_Njets[i]->SetStats();
		hNw_Njets[i]->GetXaxis()->SetTitle("#Delta(N^{exp}-N^{est})");
		sprintf(name,"Ntt_%d_jets_PseudoExp",i+4);
		hNtt_Njets[i] = new TH1F(name,"",250,EstMinMax[i][1][0],EstMinMax[i][1][1]);
		hNtt_Njets[i]->SetStats();
		hNtt_Njets[i]->GetXaxis()->SetTitle("#Delta(N^{exp}-N^{est})");
		sprintf(name,"Ntotal_%d_jets_PseudoExp",i+4);
		hNtotal[i] = new TH1F(name,"",250,EstMinMax[i][2][0],EstMinMax[i][2][1]);
		hNtotal[i]->SetStats();
		hNtotal[i]->GetXaxis()->SetTitle("#Delta(N^{exp}-N^{est})");
		
		sprintf(name,"Nw_N0bjets_%d_jets_PseudoExp",i+4);
		hRandN0bjets[i][0] = new TH1F(name,"",250,EstMinMax_Wlike[i][0][0],EstMinMax_Wlike[i][0][1]);
		sprintf(name,"Nw_N1bjets_%d_jets_PseudoExp",i+4);
		hRandN1bjets[i][0] = new TH1F(name,"",250,EstMinMax_Wlike[i][1][0],EstMinMax_Wlike[i][1][1]);
		sprintf(name,"Nw_N2bjets_%d_jets_PseudoExp",i+4);
		hRandN2bjets[i][0] = new TH1F(name,"",250,EstMinMax_Wlike[i][2][0],EstMinMax_Wlike[i][2][1]);
		sprintf(name,"Nw_N3bjets_%d_jets_PseudoExp",i+4);
		hRandN3bjets[i][0] = new TH1F(name,"",250,EstMinMax_Wlike[i][3][0],EstMinMax_Wlike[i][3][1]);

		sprintf(name,"Ntt_N0bjets_%d_jets_PseudoExp",i+4);
		hRandN0bjets[i][1] = new TH1F(name,"",250,EstMinMax_TTlike[i][0][0],EstMinMax_TTlike[i][0][1]);
		sprintf(name,"Ntt_N1bjets_%d_jets_PseudoExp",i+4);
		hRandN1bjets[i][1] = new TH1F(name,"",250,EstMinMax_TTlike[i][1][0],EstMinMax_TTlike[i][1][1]);
		sprintf(name,"Ntt_N2bjets_%d_jets_PseudoExp",i+4);
		hRandN2bjets[i][1] = new TH1F(name,"",250,EstMinMax_TTlike[i][2][0],EstMinMax_TTlike[i][2][1]);
		sprintf(name,"Ntt_N3bjets_%d_jets_PseudoExp",i+4);
		hRandN3bjets[i][1] = new TH1F(name,"",250,EstMinMax_TTlike[i][3][0],EstMinMax_TTlike[i][3][1]);
	}

	double EstMin_Wlike_Nbjets  = 99999; double EstMax_Wlike_Nbjets = -99999;
	double EstMin_TTlike_Nbjets = 99999; double EstMax_TTlike_Nbjets = -99999;

	double EstMin_Wlike_Incl  = 0; double EstMax_Wlike_Incl  = 0;
	double EstMin_TTlike_Incl = 0; double EstMax_TTlike_Incl = 0;

	for(int i=0; i<4; i++){
		EstMin_Wlike_Nbjets  = 99999; EstMax_Wlike_Nbjets  = -99999;
		EstMin_TTlike_Nbjets = 99999; EstMax_TTlike_Nbjets = -99999;
		for(int j=0;j<3;j++){
			if(EstMin_Wlike_Nbjets  < EstMinMax_Wlike[j][i][0])  EstMin_Wlike_Nbjets  = EstMinMax_Wlike[j][i][0];
			if(EstMin_TTlike_Nbjets < EstMinMax_TTlike[j][i][0]) EstMin_TTlike_Nbjets = EstMinMax_TTlike[j][i][0];
			if(EstMax_Wlike_Nbjets  < EstMinMax_Wlike[j][i][1])  EstMax_Wlike_Nbjets  = EstMinMax_Wlike[j][i][1];
			if(EstMax_TTlike_Nbjets < EstMinMax_TTlike[j][i][1]) EstMax_TTlike_Nbjets = EstMinMax_TTlike[j][i][1];
		}
		sprintf(name,"Nw_%d_bjets_PseudoExp",i);
		hNw_Nbjets[i] = new TH1F(name,"",250,EstMin_Wlike_Nbjets,EstMax_Wlike_Nbjets);
		hNw_Nbjets[i]->SetStats();
		hNw_Nbjets[i]->GetXaxis()->SetTitle("#Delta(N^{exp}-N^{est})");
		sprintf(name,"Ntt_%d_bjets_PseudoExp",i);
		hNtt_Nbjets[i] = new TH1F(name,"",250,EstMin_TTlike_Nbjets,EstMax_TTlike_Nbjets);
		hNtt_Nbjets[i]->SetStats();
		hNtt_Nbjets[i]->GetXaxis()->SetTitle("#Delta(N^{exp}-N^{est})");

		EstMin_Wlike_Incl  += EstMin_Wlike_Nbjets;
		EstMax_Wlike_Incl  += EstMax_Wlike_Nbjets;
		EstMin_TTlike_Incl += EstMin_TTlike_Nbjets;
		EstMax_TTlike_Incl += EstMax_TTlike_Nbjets;
	}

	sprintf(name,"Nw_Incl_PseudoExp");
	hNw_Incl  = new TH1F(name,"",250,EstMin_Wlike_Incl,EstMax_Wlike_Incl);
	hNw_Incl->SetStats();
	hNw_Incl->GetXaxis()->SetTitle("#Delta(N^{exp}-N^{est})");
	sprintf(name,"Ntt_Incl_PseudoExp");
	hNtt_Incl = new TH1F(name,"",250,EstMin_TTlike_Incl,EstMax_TTlike_Incl);
	hNtt_Incl->SetStats();
	hNtt_Incl->GetXaxis()->SetTitle("#Delta(N^{exp}-N^{est})");

	for(int i=0; i<NbOfPseudoExp_; i++){
		PseudoExpResults_WlikeIncl  = 0;
		PseudoExpResults_TTlikeIncl = 0;
		for(int j=0;j<3;j++){
			hEudsc[j]    ->Fill(PseudoExpResults_Effudsc[i][j]);
			hEb[j]       ->Fill(PseudoExpResults_Effb[i][j]);
			hNw_Njets[j] ->Fill(PseudoExpResults[i][j][0]);
			hNtt_Njets[j]->Fill(PseudoExpResults[i][j][1]);
			hNtotal[j]   ->Fill(PseudoExpResults[i][j][2]);
			
			hRandN0bjets[j][0]->Fill(PseudoExpResults_Wlike[i][j][0]);
			hRandN1bjets[j][0]->Fill(PseudoExpResults_Wlike[i][j][1]);
			hRandN2bjets[j][0]->Fill(PseudoExpResults_Wlike[i][j][2]);
			hRandN3bjets[j][0]->Fill(PseudoExpResults_Wlike[i][j][3]);

			hRandN0bjets[j][1]->Fill(PseudoExpResults_TTlike[i][j][0]);
			hRandN1bjets[j][1]->Fill(PseudoExpResults_TTlike[i][j][1]);
			hRandN2bjets[j][1]->Fill(PseudoExpResults_TTlike[i][j][2]);
			hRandN3bjets[j][1]->Fill(PseudoExpResults_TTlike[i][j][3]);
		}
		for(int j=0;j<4;j++){
			hNw_Nbjets[j] ->Fill(PseudoExpResults_Wlike[i][0][j] +PseudoExpResults_Wlike[i][1][j] +PseudoExpResults_Wlike[i][2][j]);
			hNtt_Nbjets[j]->Fill(PseudoExpResults_TTlike[i][0][j]+PseudoExpResults_TTlike[i][1][j]+PseudoExpResults_TTlike[i][2][j]);
			PseudoExpResults_WlikeIncl  += PseudoExpResults_Wlike[i][0][j] +PseudoExpResults_Wlike[i][1][j] +PseudoExpResults_Wlike[i][2][j];
			PseudoExpResults_TTlikeIncl += PseudoExpResults_TTlike[i][0][j]+PseudoExpResults_TTlike[i][1][j]+PseudoExpResults_TTlike[i][2][j];
		}
		hNw_Incl ->Fill(PseudoExpResults_WlikeIncl);
		hNtt_Incl->Fill(PseudoExpResults_TTlikeIncl);
	}
	wjEst_->FillInputs(InitValues[0],InitValues[1],InitValues[2],InitValues[3],multijets_est);
}

void WJetEstPseudoExp::GetPseudoExperiments(double*** randN0bjets, double*** randN1bjets, double*** randN2bjets, double*** randN3bjets){
	randN0bjets = randN0bjets_;
	randN1bjets = randN1bjets_;
	randN2bjets = randN2bjets_;
	randN3bjets = randN3bjets_;
}


double WJetEstPseudoExp::GetNWBias(bool doFit) const{
	double bias = 0;
	for(int i=0;i<3;i++){
		bias+=pow(GetNWBias(i,doFit),2);
	}
	return sqrt(bias);
}
double WJetEstPseudoExp::GetNTTBias(bool doFit) const{
	double bias = 0;
	for(int i=0;i<3;i++){
		bias+=pow(GetNTTBias(i,doFit),2);
	}
	return sqrt(bias);
}

double WJetEstPseudoExp::GetNWBias(int ijets, bool doFit) const{
	if(ijets>=0 && ijets<3 && hNw_Njets[ijets] ){
		if(doFit){
			if(hNw_Njets[ijets]->GetFunction("gaus")==0) hNw_Njets[ijets]->Fit("gaus");
			return hNw_Njets[ijets]->GetFunction("gaus")->GetParameter(1);
		}
		return hNw_Njets[ijets]->GetMean();
	}
	return -9999.;
}
double WJetEstPseudoExp::GetNTTBias(int ijets, bool doFit) const{
	if(ijets>=0 && ijets<3 && hNtt_Njets[ijets] ){
		if(doFit){
			if(hNtt_Njets[ijets]->GetFunction("gaus")==0) hNtt_Njets[ijets]->Fit("gaus");
			return hNtt_Njets[ijets]->GetFunction("gaus")->GetParameter(1);
		}
		return hNtt_Njets[ijets]->GetMean();
	}
	return -9999.;
}


double WJetEstPseudoExp::GetNWError(bool doFit) const{
	double error = 0;
	for(int i=0;i<3;i++){
		error+=pow(GetNWError(i,doFit),2);
	}
	return sqrt(error);
}
double WJetEstPseudoExp::GetNTTError(bool doFit) const{
	double error = 0;
	for(int i=0;i<3;i++){
		error+=pow(GetNTTError(i,doFit),2);
	}
	return sqrt(error);
}

double WJetEstPseudoExp::GetNWError(int ijets, bool doFit) const{
	if(ijets>=0 && ijets<3 && hNw_Njets[ijets] && hNw_Njets[ijets]->GetRMS()>0){
		if(doFit){
			if(hNw_Njets[ijets]->GetFunction("gaus")==0) hNw_Njets[ijets]->Fit("gaus");
			return hNw_Njets[ijets]->GetFunction("gaus")->GetParameter(2);
		}
		return hNw_Njets[ijets]->GetRMS();
	}
	return -9999.;
}
double WJetEstPseudoExp::GetNTTError(int ijets, bool doFit) const{
	if(ijets>=0 && ijets<3 && hNtt_Njets[ijets] && hNtt_Njets[ijets]->GetRMS()>0){
		if(doFit){
			if(hNtt_Njets[ijets]->GetFunction("gaus")==0) hNtt_Njets[ijets]->Fit("gaus");
			return hNtt_Njets[ijets]->GetFunction("gaus")->GetParameter(2);
		}
		return hNtt_Njets[ijets]->GetRMS();
	}
	return -9999.;
}

TH1F* WJetEstPseudoExp::Get_hNw(int njets, int nbjets){
	if(nbjets==0) return hRandN0bjets[njets][0];
	if(nbjets==1) return hRandN1bjets[njets][0];
	if(nbjets==2) return hRandN2bjets[njets][0];
	if(nbjets>=3) return hRandN3bjets[njets][0];
	else return NULL;
}
TH1F* WJetEstPseudoExp::Get_hNw(int njets){
	if(njets<3) return hNw_Njets[njets];
	else        return hNw_Njets[2];
}

TH1F* WJetEstPseudoExp::Get_hNw(){
	return hNw_Incl;
}

TH1F* WJetEstPseudoExp::Get_hNtt(int njets, int nbjets){
	if(nbjets==0) return hRandN0bjets[njets][1];
	if(nbjets==1) return hRandN1bjets[njets][1];
	if(nbjets==2) return hRandN2bjets[njets][1];
	if(nbjets>=3) return hRandN3bjets[njets][1];
	else return NULL;
}
TH1F* WJetEstPseudoExp::Get_hNtt(int njets){
	if(njets<3) return hNtt_Njets[njets];
	else        return hNtt_Njets[2];
}

TH1F* WJetEstPseudoExp::Get_hNtt(){
	return hNtt_Incl;
}

void WJetEstPseudoExp::Write(TFile* file, string label){
	file->cd();
	string dirname = "WJetEstPseudoExp"+label;
	file->mkdir(dirname.c_str());
	file->cd(dirname.c_str());
	for(int i=0;i<3;i++){

		hEudsc[i]    ->Write();
		hEb[i]       ->Write();
		hNw_Njets[i] ->Write();
		hNtt_Njets[i]->Write();
		hNtotal[i]->Write();

		for(int j=0; j<2; j++){
			hRandN0bjets[i][j]->Write();
			hRandN1bjets[i][j]->Write();
			hRandN2bjets[i][j]->Write();
			hRandN3bjets[i][j]->Write();
		}

	}
	for(int i=0;i<4;i++){
		hNw_Nbjets[i] ->Write();
		hNtt_Nbjets[i]->Write();
	}
	hNw_Incl ->Write();
	hNtt_Incl->Write();
}

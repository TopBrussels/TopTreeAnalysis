#include "../interface/WJetEstimation.h"

// WJetEstimation::WJetEstimation(){
// 	NbOfDatasets = 0;
// 	It _ = 40;
// 	//ini jets
// 	NofJetsBins = 3;
// 	Njets_ = new int[NofJetsBins];
// 	Njets_[0] = 4; Njets_[1] = 5; Njets_[2] = 6;
// 	
// 	//init histos
// 	histos = new TH1F**[NofJetsBins];
// 	char name[100];char title[100];
// 	for(int i=0;i<NofJetsBins;i++){
// 		//histos[i] = new TH1F*[4];
// 		//for(unsigned int j=0;j<4;j++) histos[i][j]=0;
// 		histos[i]    = new TH1F*[4];
// 		sprintf(name,"Eudsc_%d_jets",Njets_[i]);
// 		sprintf(title,"Eudsc vs iterations for %d jets",Njets_[i]);
// 		histos[i][0] = new TH1F(name,title,It_,1,It_+1);
// 		sprintf(name,"Eb_%d_jets",Njets_[i]);
// 		sprintf(title,"Eb vs iterations for %d jets",Njets_[i]);
// 		histos[i][1] = new TH1F(name,title,It_,1,It_+1);
// 		sprintf(name,"Ntt_%d_jets",Njets_[i]);
// 		sprintf(title,"Ntt vs iterations for %d jets",Njets_[i]);
// 		histos[i][2] = new TH1F(name,title,It_,1,It_+1);
// 		sprintf(name,"Nw_%d_jets",Njets_[i]);
// 		sprintf(title,"Nw vs iterations for %d jets",Njets_[i]);
// 		histos[i][3] = new TH1F(name,title,It_,1,It_+1);
// 	}
// 
// 	RescaledTTlikeEstimation = 0;
// 	RescaledWlikeEstimation = 0;
// 	tCanva_RescaledTTlikeEstimation = 0;
// 	tCanva_RescaledWlikeEstimation = 0;
//         
// 	//init summary histos
// 	hsNbjets_MC  = new THStack*[NofJetsBins+1];
// 	hsNbjets_Est = new THStack*[NofJetsBins+1];
// 	hsNjets_MC   = new THStack*[4+1];
// 	hsNjets_Est  = new THStack*[4+1];
// 	
// 	MyLeg = new TLegend(0.6,0.6,0.9,0.9);
// 	
// 	hNbjetsEstSummary = new TH1F**[NofJetsBins+1];
// 	hNbjetsMCSummary  = new TH1F**[NofJetsBins+1];
// 	hNjetsEstSummary  = new TH1F**[4+1];
// 	hNjetsMCSummary   = new TH1F**[4+1];
// 	
// 	tCanva_Nbjets_Summary = new TCanvas*[NofJetsBins+1];
// 	tCanva_Njets_Summary  = new TCanvas*[4+1];
// 
// 	for(int i=0;i<NofJetsBins;i++){
// 		hNbjetsEstSummary[i]    = new TH1F*[2];
// 		sprintf(name,"hNbjetsEstSummary_Wlike_%d_jets",Njets_[i]);
// 		sprintf(title,"W-like events estimation summary for %d jets",Njets_[i]);
// 		hNbjetsEstSummary[i][0] = new TH1F(name,title,4,0,4);
// 		sprintf(name,"hNbjetsEstSummary_TTlike_%d_jets",Njets_[i]);
// 		sprintf(title,"TT-like events estimation summary for %d jets",Njets_[i]);
// 		hNbjetsEstSummary[i][1] = new TH1F(name,title,4,0,4);
// 		hNbjetsMCSummary[i]     = new TH1F*[2];
// 		sprintf(name,"hNbjetsMCSummary_Wlike_%d_jets",Njets_[i]);
// 		sprintf(title,"W-like events MC prediction summary for %d jets",Njets_[i]);
// 		hNbjetsMCSummary[i][0]  = new TH1F(name,title,4,0,4);
// 		sprintf(name,"hNbjetsMCSummary_TTlike_%d_jets",Njets_[i]);
// 		sprintf(title,"TT-like events MC prediction summary for %d jets",Njets_[i]);
// 		hNbjetsMCSummary[i][1]  = new TH1F(name,title,4,0,4);
// 		sprintf(name,"hNbjetsSummary_%d_jets",Njets_[i]);
// 		sprintf(title,"MC prediction/EStimation summary for %d jets",Njets_[i]);
// 		tCanva_Nbjets_Summary[i] = new TCanvas(name,title,-1);
// 	}
// 
// 	hNbjetsEstSummary[NofJetsBins]    = new TH1F*[2];
// 	sprintf(name ,"hNbjetsEstSummary_Wlike_Inclusive");
// 	sprintf(title,"W-like events estimation summary (inclusive)");
// 	hNbjetsEstSummary[NofJetsBins][0] = new TH1F(name,title,4,0,4);
// 	sprintf(name ,"hNbjetsEstSummary_TTlike_Inclusive");
// 	sprintf(title,"TT-like events estimation summary (inclusive)");
// 	hNbjetsEstSummary[NofJetsBins][1] = new TH1F(name,title,4,0,4);
// 
// 	hNbjetsMCSummary[NofJetsBins]     = new TH1F*[2];
// 	sprintf(name ,"hNbjetsMCSummary_Wlike_Inclusive");
// 	sprintf(title,"W-like events MC prediction summary (inclusive)");
// 	hNbjetsMCSummary[NofJetsBins][0]  = new TH1F(name,title,4,0,4);
// 	sprintf(name ,"hNbjetsMCSummary_TTlike_Inclusive");
// 	sprintf(title,"TT-like events MC prediction summary (inclusive)");
// 	hNbjetsMCSummary[NofJetsBins][1]  = new TH1F(name,title,4,0,4);
// 	sprintf(name ,"hNbjetsSummary_Inclusive");
// 	sprintf(title,"MC prediction/EStimation summary (inclusive)");
// 	tCanva_Nbjets_Summary[NofJetsBins] = new TCanvas(name,title,-1);
// 
// 	for(int i=0;i<4;i++){
// 		hNjetsEstSummary[i]    = new TH1F*[2];
// 		sprintf(name,"hNjetsEstSummary_Wlike_%d_bjets",i);
// 		sprintf(title,"W-like events estimation summary for %d b-jets",i);
// 		hNjetsEstSummary[i][0] = new TH1F(name,title,3,0,3);
// 		sprintf(name,"hNjetsEstSummary_TTlike_%d_bjets",i);
// 		sprintf(title,"TT-like events estimation summary for %d b-jets",i);
// 		hNjetsEstSummary[i][1] = new TH1F(name,title,3,0,3);
// 		hNjetsMCSummary[i]     = new TH1F*[2];
// 		sprintf(name,"hNjetsMCSummary_Wlike_%d_bjets",i);
// 		sprintf(title,"W-like events MC prediction summary for %d b-jets",i);
// 		hNjetsMCSummary[i][0]  = new TH1F(name,title,3,0,3);
// 		sprintf(name,"hNjetsMCSummary_TTlike_%d_bjets",i);
// 		sprintf(title,"TT-like events MC prediction summary for %d b-jets",i);
// 		hNjetsMCSummary[i][1]  = new TH1F(name,title,3,0,3);
// 		sprintf(name,"hNjetsSummary_%d_bjets",i);
// 		sprintf(title,"MC prediction/Estimation summary for %d b-jets",i);		
//         	tCanva_Njets_Summary[i] = new TCanvas(name,title,-1);
// 	}
// 	hNjetsEstSummary[4]    = new TH1F*[2];
// 	sprintf(name ,"hNjetsEstSummary_Wlike_bInclusive");
// 	sprintf(title,"W-like events estimation summary (b-inclusive)");
// 	hNjetsEstSummary[4][0] = new TH1F(name,title,3,0,3);
// 	sprintf(name ,"hNjetsEstSummary_TTlike_bInclusive");
// 	sprintf(title,"TT-like events estimation summary (b-inclusive)");
// 	hNjetsEstSummary[4][1] = new TH1F(name,title,3,0,3);
// 	hNjetsMCSummary[4]     = new TH1F*[2];
// 	sprintf(name ,"hNjetsMCSummary_Wlike_bInclusive");
// 	sprintf(title,"W-like events MC prediction summary (b-inclusive)");
// 	hNjetsMCSummary[4][0]  = new TH1F(name,title,3,0,3);
// 	sprintf(name ,"hNjetsMCSummary_TTlike_bInclusive");
// 	sprintf(title,"TT-like events MC prediction summary (b-inclusive)");
// 	hNjetsMCSummary[4][1]  = new TH1F(name,title,3,0,3);
// 	sprintf(name ,"hNjetsSummary_bInclusive");
// 	sprintf(title,"MC prediction/Estimation summary (b-inclusive)");
//        	tCanva_Njets_Summary[4] = new TCanvas(name,title,-1);
// 
// 	//init multijets estimated events
// 	MultiJet_Est_N0_ = new double[NofJetsBins];
//         MultiJet_Est_N1_ = new double[NofJetsBins];
//         MultiJet_Est_N2_ = new double[NofJetsBins];
//         MultiJet_Est_N3_ = new double[NofJetsBins];
// 	for(int i=0;i<NofJetsBins;i++){
// 		MultiJet_Est_N0_[i] = 0;
//         	MultiJet_Est_N1_[i] = 0;
//         	MultiJet_Est_N2_[i] = 0;
//         	MultiJet_Est_N3_[i] = 0;
// 	}
// 
// 	//init nof events per btag bin
//         N0bjets_ = new double[NofJetsBins];
//         N1bjets_ = new double[NofJetsBins];
//         N2bjets_ = new double[NofJetsBins];
//         N3bjets_ = new double[NofJetsBins];
// 	for(int i=0;i<NofJetsBins;i++){
// 		N0bjets_[i] = 0;
//         	N1bjets_[i] = 0;
//         	N2bjets_[i] = 0;
//         	N3bjets_[i] = 0;
// 	}
// 	//per datasets & bins
// 	N0_ = new double*[NofJetsBins];
// 	N1_ = new double*[NofJetsBins];
// 	N2_ = new double*[NofJetsBins];
// 	N3_ = new double*[NofJetsBins];
// 
// 	//init estimators
//         eudsc_ = new double[NofJetsBins];
//         eb_ = new double[NofJetsBins];
//         Ntt_ = new double[NofJetsBins];
//         Nw_ = new double[NofJetsBins];
// 	statisticalWLikeError_= new double*[NofJetsBins];
// 	systematicWLikeError_ = new double*[NofJetsBins];
// 	statisticalTTLikeError_= new double*[NofJetsBins];
// 	systematicTTLikeError_ = new double*[NofJetsBins];
// 	for(int i=0;i<NofJetsBins;i++){
// 		eudsc_[i] = 0;
// 		eb_[i] = 0;
// 		Ntt_[i] = 0;
// 		Nw_[i] = 0;
// 		statisticalWLikeError_[i]= new double[5];
// 		systematicWLikeError_[i] = new double[5];
// 		statisticalTTLikeError_[i]= new double[5];
// 		systematicTTLikeError_[i] = new double[5];
// 		for(int j=0;j<5;j++){
// 			statisticalWLikeError_[i][j] = 0;
// 			systematicWLikeError_[i][j]  = 0;
// 			statisticalTTLikeError_[i][j] = 0;
// 			systematicTTLikeError_[i][j]  = 0;
// 		}
// 	}
// 	//init histo for the stability check
// 	EstimationChiSq_Eb_vs_Eudsc = new TH2D*[NofJetsBins];
// 	EstimationChiSq_Ntt_vs_Nw   = new TH2D*[NofJetsBins];
// 	for(int i=0;i<NofJetsBins;i++){
// 		EstimationChiSq_Eb_vs_Eudsc[i] = 0;
// 		EstimationChiSq_Ntt_vs_Nw[i] = 0;
// 	}
// }
// 
// WJetEstimation::WJetEstimation(int Iterations){
// 	NbOfDatasets = 0;
// 	It_ = Iterations;
// 	//ini jets
// 	NofJetsBins = 3;
// 	Njets_ = new int[NofJetsBins];
// 	Njets_[0] = 4; Njets_[1] = 5; Njets_[2] = 6;
// 	
// 	//init histos
// 	histos = new TH1F**[NofJetsBins];
// 	char name[100];char title[100];
// 	for(int i=0;i<NofJetsBins;i++){
// 		histos[i]    = new TH1F*[4];
// 		sprintf(name,"Eudsc_%d_jets",Njets_[i]);
// 		sprintf(title,"Eudsc vs iterations for %d jets",Njets_[i]);
// 		histos[i][0] = new TH1F(name,title,Iterations,1,Iterations+1);
// 		sprintf(name,"Eb_%d_jets",Njets_[i]);
// 		sprintf(title,"Eb vs iterations for %d jets",Njets_[i]);
// 		histos[i][1] = new TH1F(name,title,Iterations,1,Iterations+1);
// 		sprintf(name,"Ntt_%d_jets",Njets_[i]);
// 		sprintf(title,"Ntt vs iterations for %d jets",Njets_[i]);
// 		histos[i][2] = new TH1F(name,title,Iterations,1,Iterations+1);
// 		sprintf(name,"Nw_%d_jets",Njets_[i]);
// 		sprintf(title,"Nw vs iterations for %d jets",Njets_[i]);
// 		histos[i][3] = new TH1F(name,title,Iterations,1,Iterations+1);
// 	}
//         
// 	RescaledTTlikeEstimation = 0;
// 	RescaledWlikeEstimation = 0;
// 	tCanva_RescaledTTlikeEstimation = 0;
// 	tCanva_RescaledWlikeEstimation = 0;
//         
// 	//init summary histos
// 	hsNbjets_MC  = new THStack*[NofJetsBins+1];
// 	hsNbjets_Est = new THStack*[NofJetsBins+1];
// 	hsNjets_MC   = new THStack*[4+1];
// 	hsNjets_Est  = new THStack*[4+1];
// 	
// 	MyLeg = new TLegend(0.6,0.6,0.9,0.9);
// 	
// 	hNbjetsEstSummary = new TH1F**[NofJetsBins+1];
// 	hNbjetsMCSummary  = new TH1F**[NofJetsBins+1];
// 	hNjetsEstSummary  = new TH1F**[4+1];
// 	hNjetsMCSummary   = new TH1F**[4+1];
// 	
// 	tCanva_Nbjets_Summary = new TCanvas*[NofJetsBins+1];
// 	tCanva_Njets_Summary  = new TCanvas*[4+1];
// 
// 	for(int i=0;i<NofJetsBins;i++){
// 		hNbjetsEstSummary[i]    = new TH1F*[2];
// 		sprintf(name,"hNbjetsEstSummary_Wlike_%d_jets",Njets_[i]);
// 		sprintf(title,"W-like events estimation summary for %d jets",Njets_[i]);
// 		hNbjetsEstSummary[i][0] = new TH1F(name,title,4,0,4);
// 		sprintf(name,"hNbjetsEstSummary_TTlike_%d_jets",Njets_[i]);
// 		sprintf(title,"TT-like events estimation summary for %d jets",Njets_[i]);
// 		hNbjetsEstSummary[i][1] = new TH1F(name,title,4,0,4);
// 		hNbjetsMCSummary[i]     = new TH1F*[2];
// 		sprintf(name,"hNbjetsMCSummary_Wlike_%d_jets",Njets_[i]);
// 		sprintf(title,"W-like events MC prediction summary for %d jets",Njets_[i]);
// 		hNbjetsMCSummary[i][0]  = new TH1F(name,title,4,0,4);
// 		sprintf(name,"hNbjetsMCSummary_TTlike_%d_jets",Njets_[i]);
// 		sprintf(title,"TT-like events MC prediction summary for %d jets",Njets_[i]);
// 		hNbjetsMCSummary[i][1]  = new TH1F(name,title,4,0,4);
// 		sprintf(name,"hNbjetsSummary_%d_jets",Njets_[i]);
// 		sprintf(title,"MC prediction/EStimation summary for %d jets",Njets_[i]);
// 		tCanva_Nbjets_Summary[i] = new TCanvas(name,title,-1);
// 	}
// 
// 	hNbjetsEstSummary[NofJetsBins]    = new TH1F*[2];
// 	sprintf(name ,"hNbjetsEstSummary_Wlike_Inclusive");
// 	sprintf(title,"W-like events estimation summary (inclusive)");
// 	hNbjetsEstSummary[NofJetsBins][0] = new TH1F(name,title,4,0,4);
// 	sprintf(name ,"hNbjetsEstSummary_TTlike_Inclusive");
// 	sprintf(title,"TT-like events estimation summary (inclusive)");
// 	hNbjetsEstSummary[NofJetsBins][1] = new TH1F(name,title,4,0,4);
// 
// 	hNbjetsMCSummary[NofJetsBins]     = new TH1F*[2];
// 	sprintf(name ,"hNbjetsMCSummary_Wlike_Inclusive");
// 	sprintf(title,"W-like events MC prediction summary (inclusive)");
// 	hNbjetsMCSummary[NofJetsBins][0]  = new TH1F(name,title,4,0,4);
// 	sprintf(name ,"hNbjetsMCSummary_TTlike_Inclusive");
// 	sprintf(title,"TT-like events MC prediction summary (inclusive)");
// 	hNbjetsMCSummary[NofJetsBins][1]  = new TH1F(name,title,4,0,4);
// 	sprintf(name ,"hNbjetsSummary_Inclusive");
// 	sprintf(title,"MC prediction/EStimation summary (inclusive)");
// 	tCanva_Nbjets_Summary[NofJetsBins] = new TCanvas(name,title,-1);
// 
// 	for(int i=0;i<4;i++){
// 		hNjetsEstSummary[i]    = new TH1F*[2];
// 		sprintf(name,"hNjetsEstSummary_Wlike_%d_bjets",i);
// 		sprintf(title,"W-like events estimation summary for %d b-jets",i);
// 		hNjetsEstSummary[i][0] = new TH1F(name,title,3,0,3);
// 		sprintf(name,"hNjetsEstSummary_TTlike_%d_bjets",i);
// 		sprintf(title,"TT-like events estimation summary for %d b-jets",i);
// 		hNjetsEstSummary[i][1] = new TH1F(name,title,3,0,3);
// 		hNjetsMCSummary[i]     = new TH1F*[2];
// 		sprintf(name,"hNjetsMCSummary_Wlike_%d_bjets",i);
// 		sprintf(title,"W-like events MC prediction summary for %d b-jets",i);
// 		hNjetsMCSummary[i][0]  = new TH1F(name,title,3,0,3);
// 		sprintf(name,"hNjetsMCSummary_TTlike_%d_bjets",i);
// 		sprintf(title,"TT-like events MC prediction summary for %d b-jets",i);
// 		hNjetsMCSummary[i][1]  = new TH1F(name,title,3,0,3);
// 		sprintf(name,"hNjetsSummary_%d_bjets",i);
// 		sprintf(title,"MC prediction/Estimation summary for %d b-jets",i);		
//         	tCanva_Njets_Summary[i] = new TCanvas(name,title,-1);
// 	}
// 	hNjetsEstSummary[4]    = new TH1F*[2];
// 	sprintf(name ,"hNjetsEstSummary_Wlike_bInclusive");
// 	sprintf(title,"W-like events estimation summary (b-inclusive)");
// 	hNjetsEstSummary[4][0] = new TH1F(name,title,3,0,3);
// 	sprintf(name ,"hNjetsEstSummary_TTlike_bInclusive");
// 	sprintf(title,"TT-like events estimation summary (b-inclusive)");
// 	hNjetsEstSummary[4][1] = new TH1F(name,title,3,0,3);
// 	hNjetsMCSummary[4]     = new TH1F*[2];
// 	sprintf(name ,"hNjetsMCSummary_Wlike_bInclusive");
// 	sprintf(title,"W-like events MC prediction summary (b-inclusive)");
// 	hNjetsMCSummary[4][0]  = new TH1F(name,title,3,0,3);
// 	sprintf(name ,"hNjetsMCSummary_TTlike_bInclusive");
// 	sprintf(title,"TT-like events MC prediction summary (b-inclusive)");
// 	hNjetsMCSummary[4][1]  = new TH1F(name,title,3,0,3);
// 	sprintf(name ,"hNjetsSummary_bInclusive");
// 	sprintf(title,"MC prediction/Estimation summary (b-inclusive)");
//        	tCanva_Njets_Summary[4] = new TCanvas(name,title,-1);
//         
// 	//init multijets estimated events
// 	MultiJet_Est_N0_ = new double[NofJetsBins];
//         MultiJet_Est_N1_ = new double[NofJetsBins];
//         MultiJet_Est_N2_ = new double[NofJetsBins];
//         MultiJet_Est_N3_ = new double[NofJetsBins];
// 	for(int i=0;i<NofJetsBins;i++){
// 		MultiJet_Est_N0_[i] = 0;
//         	MultiJet_Est_N1_[i] = 0;
//         	MultiJet_Est_N2_[i] = 0;
//         	MultiJet_Est_N3_[i] = 0;
// 	}
// 
// 	//init nof events per btag bin
//         N0bjets_ = new double[NofJetsBins];
//         N1bjets_ = new double[NofJetsBins];
//         N2bjets_ = new double[NofJetsBins];
//         N3bjets_ = new double[NofJetsBins];
// 	for(int i=0;i<NofJetsBins;i++){
// 		N0bjets_[i] = 0;
//         	N1bjets_[i] = 0;
//         	N2bjets_[i] = 0;
//         	N3bjets_[i] = 0;
// 	}
// 	//per datasets & bins
// 	N0_ = new double*[NofJetsBins];
// 	N1_ = new double*[NofJetsBins];
// 	N2_ = new double*[NofJetsBins];
// 	N3_ = new double*[NofJetsBins];
// 
// 	//init estimators
//         eudsc_ = new double[NofJetsBins];
//         eb_ = new double[NofJetsBins];
//         Ntt_ = new double[NofJetsBins];
//         Nw_ = new double[NofJetsBins];
// 	statisticalWLikeError_= new double*[NofJetsBins];
// 	systematicWLikeError_ = new double*[NofJetsBins];
// 	statisticalTTLikeError_= new double*[NofJetsBins];
// 	systematicTTLikeError_ = new double*[NofJetsBins];
// 	for(int i=0;i<NofJetsBins;i++){
// 		eudsc_[i] = 0;
// 		eb_[i] = 0;
// 		Ntt_[i] = 0;
// 		Nw_[i] = 0;
// 		statisticalWLikeError_[i]= new double[5];
// 		systematicWLikeError_[i] = new double[5];
// 		statisticalTTLikeError_[i]= new double[5];
// 		systematicTTLikeError_[i] = new double[5];
// 		for(int j=0;j<5;j++){
// 			statisticalWLikeError_[i][j] = 0;
// 			systematicWLikeError_[i][j]  = 0;
// 			statisticalTTLikeError_[i][j] = 0;
// 			systematicTTLikeError_[i][j]  = 0;
// 		}
// 	}
// 	//init histo for the stability check
// 	EstimationChiSq_Eb_vs_Eudsc = new TH2D*[NofJetsBins];
// 	EstimationChiSq_Ntt_vs_Nw   = new TH2D*[NofJetsBins];
// 	for(int i=0;i<NofJetsBins;i++){
// 		EstimationChiSq_Eb_vs_Eudsc[i] = 0;
// 		EstimationChiSq_Ntt_vs_Nw[i] = 0;
// 	}
// }

WJetEstimation::WJetEstimation(int NofDatasets, vector<int> iDTTLike, vector<int> iDWLike){
	NbOfDatasets = NofDatasets;
	MCdata_ = true;
	It_ = 40;
	iDatasetsTTLike.clear();
	iDatasetsWLike.clear();
	iDatasetsTTLike = iDTTLike;
	iDatasetsWLike = iDWLike;
	
	//ini jets
	NofJetsBins = 3;
	Njets_ = new int[NofJetsBins];
	Njets_[0] = 4; Njets_[1] = 5; Njets_[2] = 6;
	
	//init histos
	histos = new TH1F**[NofJetsBins];
	char name[100];char title[100];
	for(int i=0;i<NofJetsBins;i++){
		//histos[i] = new TH1F*[4];
		//for(unsigned int j=0;j<4;j++) histos[i][j]=0;
		histos[i] = new TH1F*[4];
		sprintf(name,"Eudsc_%d_jets",Njets_[i]);
		sprintf(title,"Eudsc vs iterations for %d jets",Njets_[i]);
		histos[i][0] = new TH1F(name,title,It_,1,It_+1);
		sprintf(name,"Eb_%d_jets",Njets_[i]);
		sprintf(title,"Eb vs iterations for %d jets",Njets_[i]);
		histos[i][1] = new TH1F(name,title,It_,1,It_+1);
		sprintf(name,"Ntt_%d_jets",Njets_[i]);
		sprintf(title,"Ntt vs iterations for %d jets",Njets_[i]);
		histos[i][2] = new TH1F(name,title,It_,1,It_+1);
		sprintf(name,"Nw_%d_jets",Njets_[i]);
		sprintf(title,"Nw vs iterations for %d jets",Njets_[i]);
		histos[i][3] = new TH1F(name,title,It_,1,It_+1);
	}
        
	RescaledTTlikeEstimation = 0;
	RescaledWlikeEstimation = 0;
	tCanva_RescaledTTlikeEstimation = 0;
	tCanva_RescaledWlikeEstimation = 0;
        
	//init summary histos
	hsNbjets_MC  = new THStack*[NofJetsBins+1];
	hsNbjets_Est = new THStack*[NofJetsBins+1];
	hsNjets_MC   = new THStack*[4+1];
	hsNjets_Est  = new THStack*[4+1];
	
	MyLeg = new TLegend(0.6,0.6,0.9,0.9);
	
	hNbjetsEstSummary = new TH1F**[NofJetsBins+1];
	hNbjetsMCSummary  = new TH1F**[NofJetsBins+1];
	hNjetsEstSummary  = new TH1F**[4+1];
	hNjetsMCSummary   = new TH1F**[4+1];
	
	tCanva_Nbjets_Summary = new TCanvas*[NofJetsBins+1];
	tCanva_Njets_Summary  = new TCanvas*[4+1];

	for(int i=0;i<NofJetsBins;i++){
		hNbjetsEstSummary[i]    = new TH1F*[2];
		sprintf(name,"hNbjetsEstSummary_Wlike_%d_jets",Njets_[i]);
		sprintf(title,"W-like events estimation summary for %d jets",Njets_[i]);
		hNbjetsEstSummary[i][0] = new TH1F(name,title,4,0,4);
		sprintf(name,"hNbjetsEstSummary_TTlike_%d_jets",Njets_[i]);
		sprintf(title,"TT-like events estimation summary for %d jets",Njets_[i]);
		hNbjetsEstSummary[i][1] = new TH1F(name,title,4,0,4);
		hNbjetsMCSummary[i]     = new TH1F*[2];
		sprintf(name,"hNbjetsMCSummary_Wlike_%d_jets",Njets_[i]);
		sprintf(title,"W-like events MC prediction summary for %d jets",Njets_[i]);
		hNbjetsMCSummary[i][0]  = new TH1F(name,title,4,0,4);
		sprintf(name,"hNbjetsMCSummary_TTlike_%d_jets",Njets_[i]);
		sprintf(title,"TT-like events MC prediction summary for %d jets",Njets_[i]);
		hNbjetsMCSummary[i][1]  = new TH1F(name,title,4,0,4);
		sprintf(name,"hNbjetsSummary_%d_jets",Njets_[i]);
		sprintf(title,"MC prediction/EStimation summary for %d jets",Njets_[i]);
		tCanva_Nbjets_Summary[i] = new TCanvas(name,title,-1);
	}

	hNbjetsEstSummary[NofJetsBins]    = new TH1F*[2];
	sprintf(name ,"hNbjetsEstSummary_Wlike_Inclusive");
	sprintf(title,"W-like events estimation summary (inclusive)");
	hNbjetsEstSummary[NofJetsBins][0] = new TH1F(name,title,4,0,4);
	sprintf(name ,"hNbjetsEstSummary_TTlike_Inclusive");
	sprintf(title,"TT-like events estimation summary (inclusive)");
	hNbjetsEstSummary[NofJetsBins][1] = new TH1F(name,title,4,0,4);

	hNbjetsMCSummary[NofJetsBins]     = new TH1F*[2];
	sprintf(name ,"hNbjetsMCSummary_Wlike_Inclusive");
	sprintf(title,"W-like events MC prediction summary (inclusive)");
	hNbjetsMCSummary[NofJetsBins][0]  = new TH1F(name,title,4,0,4);
	sprintf(name ,"hNbjetsMCSummary_TTlike_Inclusive");
	sprintf(title,"TT-like events MC prediction summary (inclusive)");
	hNbjetsMCSummary[NofJetsBins][1]  = new TH1F(name,title,4,0,4);
	sprintf(name ,"hNbjetsSummary_Inclusive");
	sprintf(title,"MC prediction/EStimation summary (inclusive)");
	tCanva_Nbjets_Summary[NofJetsBins] = new TCanvas(name,title,-1);

	for(int i=0;i<4;i++){
		hNjetsEstSummary[i]    = new TH1F*[2];
		sprintf(name,"hNjetsEstSummary_Wlike_%d_bjets",i);
		sprintf(title,"W-like events estimation summary for %d b-jets",i);
		hNjetsEstSummary[i][0] = new TH1F(name,title,3,0,3);
		sprintf(name,"hNjetsEstSummary_TTlike_%d_bjets",i);
		sprintf(title,"TT-like events estimation summary for %d b-jets",i);
		hNjetsEstSummary[i][1] = new TH1F(name,title,3,0,3);
		hNjetsMCSummary[i]     = new TH1F*[2];
		sprintf(name,"hNjetsMCSummary_Wlike_%d_bjets",i);
		sprintf(title,"W-like events MC prediction summary for %d b-jets",i);
		hNjetsMCSummary[i][0]  = new TH1F(name,title,3,0,3);
		sprintf(name,"hNjetsMCSummary_TTlike_%d_bjets",i);
		sprintf(title,"TT-like events MC prediction summary for %d b-jets",i);
		hNjetsMCSummary[i][1]  = new TH1F(name,title,3,0,3);
		sprintf(name,"hNjetsSummary_%d_bjets",i);
		sprintf(title,"MC prediction/Estimation summary for %d b-jets",i);		
        	tCanva_Njets_Summary[i] = new TCanvas(name,title,-1);
	}
	hNjetsEstSummary[4]    = new TH1F*[2];
	sprintf(name ,"hNjetsEstSummary_Wlike_bInclusive");
	sprintf(title,"W-like events estimation summary (b-inclusive)");
	hNjetsEstSummary[4][0] = new TH1F(name,title,3,0,3);
	sprintf(name ,"hNjetsEstSummary_TTlike_bInclusive");
	sprintf(title,"TT-like events estimation summary (b-inclusive)");
	hNjetsEstSummary[4][1] = new TH1F(name,title,3,0,3);
	hNjetsMCSummary[4]     = new TH1F*[2];
	sprintf(name ,"hNjetsMCSummary_Wlike_bInclusive");
	sprintf(title,"W-like events MC prediction summary (b-inclusive)");
	hNjetsMCSummary[4][0]  = new TH1F(name,title,3,0,3);
	sprintf(name ,"hNjetsMCSummary_TTlike_bInclusive");
	sprintf(title,"TT-like events MC prediction summary (b-inclusive)");
	hNjetsMCSummary[4][1]  = new TH1F(name,title,3,0,3);
	sprintf(name ,"hNjetsSummary_bInclusive");
	sprintf(title,"MC prediction/Estimation summary (b-inclusive)");
       	tCanva_Njets_Summary[4] = new TCanvas(name,title,-1);

	//init multijets estimated events
	MultiJet_Est_N0_ = new double[NofJetsBins];
        MultiJet_Est_N1_ = new double[NofJetsBins];
        MultiJet_Est_N2_ = new double[NofJetsBins];
        MultiJet_Est_N3_ = new double[NofJetsBins];
	for(int i=0;i<NofJetsBins;i++){
		MultiJet_Est_N0_[i] = 0;
        	MultiJet_Est_N1_[i] = 0;
        	MultiJet_Est_N2_[i] = 0;
        	MultiJet_Est_N3_[i] = 0;
	}

	//init nof events per btag bin
        N0bjets_ = new double[NofJetsBins];
        N1bjets_ = new double[NofJetsBins];
        N2bjets_ = new double[NofJetsBins];
        N3bjets_ = new double[NofJetsBins];
	for(int i=0;i<NofJetsBins;i++){
		N0bjets_[i] = 0;
        	N1bjets_[i] = 0;
        	N2bjets_[i] = 0;
        	N3bjets_[i] = 0;
	}
	//per datasets & bins
	N0_ = new double*[NofJetsBins];
	N1_ = new double*[NofJetsBins];
	N2_ = new double*[NofJetsBins];
	N3_ = new double*[NofJetsBins];
	// Multi-jets/tt+jets/W+jets/Z+jets/st-s/st-t/st-tW
	for(int i=0;i<NofJetsBins;i++){
		N0_[i] = new double[NbOfDatasets];
		N1_[i] = new double[NbOfDatasets];
		N2_[i] = new double[NbOfDatasets];
		N3_[i] = new double[NbOfDatasets];
		for(int j=0;j<NofDatasets;j++){
			N0_[i][j] = 0;
			N1_[i][j] = 0;
			N2_[i][j] = 0;
			N3_[i][j] = 0;
		}
	}

	//init estimators
        eudsc_ = new double[NofJetsBins];
        eb_ = new double[NofJetsBins];
        Ntt_ = new double[NofJetsBins];
        Nw_ = new double[NofJetsBins];
	statisticalWLikeError_= new double*[NofJetsBins];
	systematicWLikeError_ = new double*[NofJetsBins];
	statisticalTTLikeError_= new double*[NofJetsBins];
	systematicTTLikeError_ = new double*[NofJetsBins];
	for(int i=0;i<NofJetsBins;i++){
		eudsc_[i] = 0;
		eb_[i] = 0;
		Ntt_[i] = 0;
		Nw_[i] = 0;
		statisticalWLikeError_[i]= new double[5];
		systematicWLikeError_[i] = new double[5];
		statisticalTTLikeError_[i]= new double[5];
		systematicTTLikeError_[i] = new double[5];
		for(int j=0;j<5;j++){
			statisticalWLikeError_[i][j] = 0;
			systematicWLikeError_[i][j]  = 0;
			statisticalTTLikeError_[i][j] = 0;
			systematicTTLikeError_[i][j]  = 0;
		}
	}
	//init histo for the stability check
	EstimationChiSq_Eb_vs_Eudsc = new TH2D*[NofJetsBins];
	EstimationChiSq_Ntt_vs_Nw   = new TH2D*[NofJetsBins];
	for(int i=0;i<NofJetsBins;i++){
		EstimationChiSq_Eb_vs_Eudsc[i] = 0;
		EstimationChiSq_Ntt_vs_Nw[i] = 0;
	}
}

WJetEstimation::WJetEstimation(int Iterations, int NofDatasets, vector<int> iDTTLike, vector<int> iDWLike){
	NbOfDatasets = NofDatasets;
	MCdata_ = true;
	It_ = Iterations;
	iDatasetsTTLike.clear();
	iDatasetsWLike.clear();
	iDatasetsTTLike = iDTTLike;
	iDatasetsWLike = iDWLike;
	
	//ini jets
	NofJetsBins = 3;
	Njets_ = new int[NofJetsBins];
	Njets_[0] = 4; Njets_[1] = 5; Njets_[2] = 6;
	
	//init histos
	histos = new TH1F**[NofJetsBins];
	char name[100];char title[100];
	for(int i=0;i<NofJetsBins;i++){
		//histos[i] = new TH1F*[4];
		//for(unsigned int j=0;j<4;j++) histos[i][j]=0;
		histos[i] = new TH1F*[4];
		sprintf(name,"Eudsc_%d_jets",Njets_[i]);
		sprintf(title,"Eudsc vs iterations for %d jets",Njets_[i]);
		histos[i][0] = new TH1F(name,title,Iterations,1,Iterations+1);
		sprintf(name,"Eb_%d_jets",Njets_[i]);
		sprintf(title,"Eb vs iterations for %d jets",Njets_[i]);
		histos[i][1] = new TH1F(name,title,Iterations,1,Iterations+1);
		sprintf(name,"Ntt_%d_jets",Njets_[i]);
		sprintf(title,"Ntt vs iterations for %d jets",Njets_[i]);
		histos[i][2] = new TH1F(name,title,Iterations,1,Iterations+1);
		sprintf(name,"Nw_%d_jets",Njets_[i]);
		sprintf(title,"Nw vs iterations for %d jets",Njets_[i]);
		histos[i][3] = new TH1F(name,title,Iterations,1,Iterations+1);
	}
        
	RescaledTTlikeEstimation = 0;
	RescaledWlikeEstimation = 0;
	tCanva_RescaledTTlikeEstimation = 0;
	tCanva_RescaledWlikeEstimation = 0;
        
	//init summary histos
	hsNbjets_MC  = new THStack*[NofJetsBins+1];
	hsNbjets_Est = new THStack*[NofJetsBins+1];
	hsNjets_MC   = new THStack*[4+1];
	hsNjets_Est  = new THStack*[4+1];
	
	MyLeg = new TLegend(0.6,0.6,0.9,0.9);
	
	hNbjetsEstSummary = new TH1F**[NofJetsBins+1];
	hNbjetsMCSummary  = new TH1F**[NofJetsBins+1];
	hNjetsEstSummary  = new TH1F**[4+1];
	hNjetsMCSummary   = new TH1F**[4+1];
	
	tCanva_Nbjets_Summary = new TCanvas*[NofJetsBins+1];
	tCanva_Njets_Summary  = new TCanvas*[4+1];

	for(int i=0;i<NofJetsBins;i++){
		hNbjetsEstSummary[i]    = new TH1F*[2];
		sprintf(name,"hNbjetsEstSummary_Wlike_%d_jets",Njets_[i]);
		sprintf(title,"W-like events estimation summary for %d jets",Njets_[i]);
		hNbjetsEstSummary[i][0] = new TH1F(name,title,4,0,4);
		sprintf(name,"hNbjetsEstSummary_TTlike_%d_jets",Njets_[i]);
		sprintf(title,"TT-like events estimation summary for %d jets",Njets_[i]);
		hNbjetsEstSummary[i][1] = new TH1F(name,title,4,0,4);
		hNbjetsMCSummary[i]     = new TH1F*[2];
		sprintf(name,"hNbjetsMCSummary_Wlike_%d_jets",Njets_[i]);
		sprintf(title,"W-like events MC prediction summary for %d jets",Njets_[i]);
		hNbjetsMCSummary[i][0]  = new TH1F(name,title,4,0,4);
		sprintf(name,"hNbjetsMCSummary_TTlike_%d_jets",Njets_[i]);
		sprintf(title,"TT-like events MC prediction summary for %d jets",Njets_[i]);
		hNbjetsMCSummary[i][1]  = new TH1F(name,title,4,0,4);
		sprintf(name,"hNbjetsSummary_%d_jets",Njets_[i]);
		sprintf(title,"MC prediction/EStimation summary for %d jets",Njets_[i]);
		tCanva_Nbjets_Summary[i] = new TCanvas(name,title,-1);
	}

	hNbjetsEstSummary[NofJetsBins]    = new TH1F*[2];
	sprintf(name ,"hNbjetsEstSummary_Wlike_Inclusive");
	sprintf(title,"W-like events estimation summary (inclusive)");
	hNbjetsEstSummary[NofJetsBins][0] = new TH1F(name,title,4,0,4);
	sprintf(name ,"hNbjetsEstSummary_TTlike_Inclusive");
	sprintf(title,"TT-like events estimation summary (inclusive)");
	hNbjetsEstSummary[NofJetsBins][1] = new TH1F(name,title,4,0,4);

	hNbjetsMCSummary[NofJetsBins]     = new TH1F*[2];
	sprintf(name ,"hNbjetsMCSummary_Wlike_Inclusive");
	sprintf(title,"W-like events MC prediction summary (inclusive)");
	hNbjetsMCSummary[NofJetsBins][0]  = new TH1F(name,title,4,0,4);
	sprintf(name ,"hNbjetsMCSummary_TTlike_Inclusive");
	sprintf(title,"TT-like events MC prediction summary (inclusive)");
	hNbjetsMCSummary[NofJetsBins][1]  = new TH1F(name,title,4,0,4);
	sprintf(name ,"hNbjetsSummary_Inclusive");
	sprintf(title,"MC prediction/EStimation summary (inclusive)");
	tCanva_Nbjets_Summary[NofJetsBins] = new TCanvas(name,title,-1);

	for(int i=0;i<4;i++){
		hNjetsEstSummary[i]    = new TH1F*[2];
		sprintf(name,"hNjetsEstSummary_Wlike_%d_bjets",i);
		sprintf(title,"W-like events estimation summary for %d b-jets",i);
		hNjetsEstSummary[i][0] = new TH1F(name,title,3,0,3);
		sprintf(name,"hNjetsEstSummary_TTlike_%d_bjets",i);
		sprintf(title,"TT-like events estimation summary for %d b-jets",i);
		hNjetsEstSummary[i][1] = new TH1F(name,title,3,0,3);
		hNjetsMCSummary[i]     = new TH1F*[2];
		sprintf(name,"hNjetsMCSummary_Wlike_%d_bjets",i);
		sprintf(title,"W-like events MC prediction summary for %d b-jets",i);
		hNjetsMCSummary[i][0]  = new TH1F(name,title,3,0,3);
		sprintf(name,"hNjetsMCSummary_TTlike_%d_bjets",i);
		sprintf(title,"TT-like events MC prediction summary for %d b-jets",i);
		hNjetsMCSummary[i][1]  = new TH1F(name,title,3,0,3);
		sprintf(name,"hNjetsSummary_%d_bjets",i);
		sprintf(title,"MC prediction/Estimation summary for %d b-jets",i);		
        	tCanva_Njets_Summary[i] = new TCanvas(name,title,-1);
	}
	hNjetsEstSummary[4]    = new TH1F*[2];
	sprintf(name ,"hNjetsEstSummary_Wlike_bInclusive");
	sprintf(title,"W-like events estimation summary (b-inclusive)");
	hNjetsEstSummary[4][0] = new TH1F(name,title,3,0,3);
	sprintf(name ,"hNjetsEstSummary_TTlike_bInclusive");
	sprintf(title,"TT-like events estimation summary (b-inclusive)");
	hNjetsEstSummary[4][1] = new TH1F(name,title,3,0,3);
	hNjetsMCSummary[4]     = new TH1F*[2];
	sprintf(name ,"hNjetsMCSummary_Wlike_bInclusive");
	sprintf(title,"W-like events MC prediction summary (b-inclusive)");
	hNjetsMCSummary[4][0]  = new TH1F(name,title,3,0,3);
	sprintf(name ,"hNjetsMCSummary_TTlike_bInclusive");
	sprintf(title,"TT-like events MC prediction summary (b-inclusive)");
	hNjetsMCSummary[4][1]  = new TH1F(name,title,3,0,3);
	sprintf(name ,"hNjetsSummary_bInclusive");
	sprintf(title,"MC prediction/Estimation summary (b-inclusive)");
       	tCanva_Njets_Summary[4] = new TCanvas(name,title,-1);

	//init multijets estimated events
	MultiJet_Est_N0_ = new double[NofJetsBins];
        MultiJet_Est_N1_ = new double[NofJetsBins];
        MultiJet_Est_N2_ = new double[NofJetsBins];
        MultiJet_Est_N3_ = new double[NofJetsBins];
	for(int i=0;i<NofJetsBins;i++){
		MultiJet_Est_N0_[i] = 0;
        	MultiJet_Est_N1_[i] = 0;
        	MultiJet_Est_N2_[i] = 0;
        	MultiJet_Est_N3_[i] = 0;
	}

	//init nof events per btag bin
        N0bjets_ = new double[NofJetsBins];
        N1bjets_ = new double[NofJetsBins];
        N2bjets_ = new double[NofJetsBins];
        N3bjets_ = new double[NofJetsBins];
	for(int i=0;i<NofJetsBins;i++){
		N0bjets_[i] = 0;
        	N1bjets_[i] = 0;
        	N2bjets_[i] = 0;
        	N3bjets_[i] = 0;
	}
	//per datasets & bins
	N0_ = new double*[NofJetsBins];
	N1_ = new double*[NofJetsBins];
	N2_ = new double*[NofJetsBins];
	N3_ = new double*[NofJetsBins];
	// Multi-jets/tt+jets/W+jets/Z+jets/st-s/st-t/st-tW
	for(int i=0;i<NofJetsBins;i++){
		N0_[i] = new double[NbOfDatasets];
		N1_[i] = new double[NbOfDatasets];
		N2_[i] = new double[NbOfDatasets];
		N3_[i] = new double[NbOfDatasets];
		for(int j=0;j<NofDatasets;j++){
			N0_[i][j] = 0;
			N1_[i][j] = 0;
			N2_[i][j] = 0;
			N3_[i][j] = 0;
		}
	}

	//init estimators
        eudsc_ = new double[NofJetsBins];
        eb_ = new double[NofJetsBins];
        Ntt_ = new double[NofJetsBins];
        Nw_ = new double[NofJetsBins];
	statisticalWLikeError_= new double*[NofJetsBins];
	systematicWLikeError_ = new double*[NofJetsBins];
	statisticalTTLikeError_= new double*[NofJetsBins];
	systematicTTLikeError_ = new double*[NofJetsBins];
	for(int i=0;i<NofJetsBins;i++){
		eudsc_[i] = 0;
		eb_[i] = 0;
		Ntt_[i] = 0;
		Nw_[i] = 0;
		statisticalWLikeError_[i]= new double[5];
		systematicWLikeError_[i] = new double[5];
		statisticalTTLikeError_[i]= new double[5];
		systematicTTLikeError_[i] = new double[5];
		for(int j=0;j<5;j++){
			statisticalWLikeError_[i][j] = 0;
			systematicWLikeError_[i][j]  = 0;
			statisticalTTLikeError_[i][j] = 0;
			systematicTTLikeError_[i][j]  = 0;
		}
	}
	//init histo for the stability check
	EstimationChiSq_Eb_vs_Eudsc = new TH2D*[NofJetsBins];
	EstimationChiSq_Ntt_vs_Nw   = new TH2D*[NofJetsBins];
	for(int i=0;i<NofJetsBins;i++){
		EstimationChiSq_Eb_vs_Eudsc[i] = 0;
		EstimationChiSq_Ntt_vs_Nw[i] = 0;
	}
}

WJetEstimation::WJetEstimation(int Iterations){
	NbOfDatasets = 1;
	MCdata_ = false;
	It_ = Iterations;
	iDatasetsTTLike.clear();
	iDatasetsWLike.clear();
	//iDatasetsTTLike = iDTTLike;
	//iDatasetsWLike = iDWLike;
	
	//ini jets
	NofJetsBins = 3;
	Njets_ = new int[NofJetsBins];
	Njets_[0] = 4; Njets_[1] = 5; Njets_[2] = 6;
	
	//init histos
	histos = new TH1F**[NofJetsBins];
	char name[100];char title[100];
	for(int i=0;i<NofJetsBins;i++){
		//histos[i] = new TH1F*[4];
		//for(unsigned int j=0;j<4;j++) histos[i][j]=0;
		histos[i] = new TH1F*[4];
		sprintf(name,"Eudsc_%d_jets",Njets_[i]);
		sprintf(title,"Eudsc vs iterations for %d jets",Njets_[i]);
		histos[i][0] = new TH1F(name,title,Iterations,1,Iterations+1);
		sprintf(name,"Eb_%d_jets",Njets_[i]);
		sprintf(title,"Eb vs iterations for %d jets",Njets_[i]);
		histos[i][1] = new TH1F(name,title,Iterations,1,Iterations+1);
		sprintf(name,"Ntt_%d_jets",Njets_[i]);
		sprintf(title,"Ntt vs iterations for %d jets",Njets_[i]);
		histos[i][2] = new TH1F(name,title,Iterations,1,Iterations+1);
		sprintf(name,"Nw_%d_jets",Njets_[i]);
		sprintf(title,"Nw vs iterations for %d jets",Njets_[i]);
		histos[i][3] = new TH1F(name,title,Iterations,1,Iterations+1);
	}
        
	RescaledTTlikeEstimation = 0;
	RescaledWlikeEstimation = 0;
	tCanva_RescaledTTlikeEstimation = 0;
	tCanva_RescaledWlikeEstimation = 0;
        
	//init summary histos
	hsNbjets_MC  = new THStack*[NofJetsBins+1];
	hsNbjets_Est = new THStack*[NofJetsBins+1];
	hsNjets_MC   = new THStack*[4+1];
	hsNjets_Est  = new THStack*[4+1];
	
	MyLeg = new TLegend(0.6,0.6,0.9,0.9);
	
	hNbjetsEstSummary = new TH1F**[NofJetsBins+1];
	hNbjetsMCSummary  = new TH1F**[NofJetsBins+1];
	hNjetsEstSummary  = new TH1F**[4+1];
	hNjetsMCSummary   = new TH1F**[4+1];
	
	tCanva_Nbjets_Summary = new TCanvas*[NofJetsBins+1];
	tCanva_Njets_Summary  = new TCanvas*[4+1];

	for(int i=0;i<NofJetsBins;i++){
		hNbjetsEstSummary[i]    = new TH1F*[2];
		sprintf(name,"hNbjetsEstSummary_Wlike_%d_jets",Njets_[i]);
		sprintf(title,"W-like events estimation summary for %d jets",Njets_[i]);
		hNbjetsEstSummary[i][0] = new TH1F(name,title,4,0,4);
		sprintf(name,"hNbjetsEstSummary_TTlike_%d_jets",Njets_[i]);
		sprintf(title,"TT-like events estimation summary for %d jets",Njets_[i]);
		hNbjetsEstSummary[i][1] = new TH1F(name,title,4,0,4);
		hNbjetsMCSummary[i]     = new TH1F*[2];
		sprintf(name,"hNbjetsMCSummary_Wlike_%d_jets",Njets_[i]);
		sprintf(title,"W-like events MC prediction summary for %d jets",Njets_[i]);
		hNbjetsMCSummary[i][0]  = new TH1F(name,title,4,0,4);
		sprintf(name,"hNbjetsMCSummary_TTlike_%d_jets",Njets_[i]);
		sprintf(title,"TT-like events MC prediction summary for %d jets",Njets_[i]);
		hNbjetsMCSummary[i][1]  = new TH1F(name,title,4,0,4);
		sprintf(name,"hNbjetsSummary_%d_jets",Njets_[i]);
		sprintf(title,"MC prediction/EStimation summary for %d jets",Njets_[i]);
		tCanva_Nbjets_Summary[i] = new TCanvas(name,title,-1);
	}

	hNbjetsEstSummary[NofJetsBins]    = new TH1F*[2];
	sprintf(name ,"hNbjetsEstSummary_Wlike_Inclusive");
	sprintf(title,"W-like events estimation summary (inclusive)");
	hNbjetsEstSummary[NofJetsBins][0] = new TH1F(name,title,4,0,4);
	sprintf(name ,"hNbjetsEstSummary_TTlike_Inclusive");
	sprintf(title,"TT-like events estimation summary (inclusive)");
	hNbjetsEstSummary[NofJetsBins][1] = new TH1F(name,title,4,0,4);

	hNbjetsMCSummary[NofJetsBins]     = new TH1F*[2];
	sprintf(name ,"hNbjetsMCSummary_Wlike_Inclusive");
	sprintf(title,"W-like events MC prediction summary (inclusive)");
	hNbjetsMCSummary[NofJetsBins][0]  = new TH1F(name,title,4,0,4);
	sprintf(name ,"hNbjetsMCSummary_TTlike_Inclusive");
	sprintf(title,"TT-like events MC prediction summary (inclusive)");
	hNbjetsMCSummary[NofJetsBins][1]  = new TH1F(name,title,4,0,4);
	sprintf(name ,"hNbjetsSummary_Inclusive");
	sprintf(title,"MC prediction/EStimation summary (inclusive)");
	tCanva_Nbjets_Summary[NofJetsBins] = new TCanvas(name,title,-1);

	for(int i=0;i<4;i++){
		hNjetsEstSummary[i]    = new TH1F*[2];
		sprintf(name,"hNjetsEstSummary_Wlike_%d_bjets",i);
		sprintf(title,"W-like events estimation summary for %d b-jets",i);
		hNjetsEstSummary[i][0] = new TH1F(name,title,3,0,3);
		sprintf(name,"hNjetsEstSummary_TTlike_%d_bjets",i);
		sprintf(title,"TT-like events estimation summary for %d b-jets",i);
		hNjetsEstSummary[i][1] = new TH1F(name,title,3,0,3);
		hNjetsMCSummary[i]     = new TH1F*[2];
		sprintf(name,"hNjetsMCSummary_Wlike_%d_bjets",i);
		sprintf(title,"W-like events MC prediction summary for %d b-jets",i);
		hNjetsMCSummary[i][0]  = new TH1F(name,title,3,0,3);
		sprintf(name,"hNjetsMCSummary_TTlike_%d_bjets",i);
		sprintf(title,"TT-like events MC prediction summary for %d b-jets",i);
		hNjetsMCSummary[i][1]  = new TH1F(name,title,3,0,3);
		sprintf(name,"hNjetsSummary_%d_bjets",i);
		sprintf(title,"MC prediction/Estimation summary for %d b-jets",i);		
        	tCanva_Njets_Summary[i] = new TCanvas(name,title,-1);
	}
	hNjetsEstSummary[4]    = new TH1F*[2];
	sprintf(name ,"hNjetsEstSummary_Wlike_bInclusive");
	sprintf(title,"W-like events estimation summary (b-inclusive)");
	hNjetsEstSummary[4][0] = new TH1F(name,title,3,0,3);
	sprintf(name ,"hNjetsEstSummary_TTlike_bInclusive");
	sprintf(title,"TT-like events estimation summary (b-inclusive)");
	hNjetsEstSummary[4][1] = new TH1F(name,title,3,0,3);
	hNjetsMCSummary[4]     = new TH1F*[2];
	sprintf(name ,"hNjetsMCSummary_Wlike_bInclusive");
	sprintf(title,"W-like events MC prediction summary (b-inclusive)");
	hNjetsMCSummary[4][0]  = new TH1F(name,title,3,0,3);
	sprintf(name ,"hNjetsMCSummary_TTlike_bInclusive");
	sprintf(title,"TT-like events MC prediction summary (b-inclusive)");
	hNjetsMCSummary[4][1]  = new TH1F(name,title,3,0,3);
	sprintf(name ,"hNjetsSummary_bInclusive");
	sprintf(title,"MC prediction/Estimation summary (b-inclusive)");
       	tCanva_Njets_Summary[4] = new TCanvas(name,title,-1);

	//init multijets estimated events
	MultiJet_Est_N0_ = new double[NofJetsBins];
        MultiJet_Est_N1_ = new double[NofJetsBins];
        MultiJet_Est_N2_ = new double[NofJetsBins];
        MultiJet_Est_N3_ = new double[NofJetsBins];
	for(int i=0;i<NofJetsBins;i++){
		MultiJet_Est_N0_[i] = 0;
        	MultiJet_Est_N1_[i] = 0;
        	MultiJet_Est_N2_[i] = 0;
        	MultiJet_Est_N3_[i] = 0;
	}

	//init nof events per btag bin
        N0bjets_ = new double[NofJetsBins];
        N1bjets_ = new double[NofJetsBins];
        N2bjets_ = new double[NofJetsBins];
        N3bjets_ = new double[NofJetsBins];
	for(int i=0;i<NofJetsBins;i++){
		N0bjets_[i] = 0;
        	N1bjets_[i] = 0;
        	N2bjets_[i] = 0;
        	N3bjets_[i] = 0;
	}
	//per datasets & bins
	N0_ = new double*[NofJetsBins];
	N1_ = new double*[NofJetsBins];
	N2_ = new double*[NofJetsBins];
	N3_ = new double*[NofJetsBins];
	// Multi-jets/tt+jets/W+jets/Z+jets/st-s/st-t/st-tW
	for(int i=0;i<NofJetsBins;i++){
		N0_[i] = new double[NbOfDatasets];
		N1_[i] = new double[NbOfDatasets];
		N2_[i] = new double[NbOfDatasets];
		N3_[i] = new double[NbOfDatasets];
		for(int j=0;j<NbOfDatasets;j++){
			N0_[i][j] = 0;
			N1_[i][j] = 0;
			N2_[i][j] = 0;
			N3_[i][j] = 0;
		}
	}

	//init estimators
        eudsc_ = new double[NofJetsBins];
        eb_ = new double[NofJetsBins];
        Ntt_ = new double[NofJetsBins];
        Nw_ = new double[NofJetsBins];
	statisticalWLikeError_= new double*[NofJetsBins];
	systematicWLikeError_ = new double*[NofJetsBins];
	statisticalTTLikeError_= new double*[NofJetsBins];
	systematicTTLikeError_ = new double*[NofJetsBins];
	for(int i=0;i<NofJetsBins;i++){
		eudsc_[i] = 0;
		eb_[i] = 0;
		Ntt_[i] = 0;
		Nw_[i] = 0;
		statisticalWLikeError_[i]= new double[5];
		systematicWLikeError_[i] = new double[5];
		statisticalTTLikeError_[i]= new double[5];
		systematicTTLikeError_[i] = new double[5];
		for(int j=0;j<5;j++){
			statisticalWLikeError_[i][j] = 0;
			systematicWLikeError_[i][j]  = 0;
			statisticalTTLikeError_[i][j] = 0;
			systematicTTLikeError_[i][j]  = 0;
		}
	}
	//init histo for the stability check
	EstimationChiSq_Eb_vs_Eudsc = new TH2D*[NofJetsBins];
	EstimationChiSq_Ntt_vs_Nw   = new TH2D*[NofJetsBins];
	for(int i=0;i<NofJetsBins;i++){
		EstimationChiSq_Eb_vs_Eudsc[i] = 0;
		EstimationChiSq_Ntt_vs_Nw[i] = 0;
	}
}

WJetEstimation::~WJetEstimation(){
	for(int i=0;i<NofJetsBins;i++){
		for(int j=0;j<4;j++){
			delete histos[i][j];
		}
		delete statisticalWLikeError_[i];
		delete systematicWLikeError_[i];
		delete statisticalTTLikeError_[i];
		delete systematicTTLikeError_[i];
		delete[] histos[i];
		delete hNbjetsEstSummary[i][0]; delete hNbjetsEstSummary[i][1];
		delete hNbjetsMCSummary[i][0];  delete hNbjetsMCSummary[i][1];
		delete hsNbjets_MC[i];          delete hsNbjets_Est[i];
		delete tCanva_Nbjets_Summary[i];
		if(EstimationChiSq_Eb_vs_Eudsc[i]) delete EstimationChiSq_Eb_vs_Eudsc[i];
		if(EstimationChiSq_Ntt_vs_Nw[i])   delete EstimationChiSq_Ntt_vs_Nw[i];

	}
	if(EstimationChiSq_Eb_vs_Eudsc) delete[] EstimationChiSq_Eb_vs_Eudsc;
	if(EstimationChiSq_Ntt_vs_Nw)   delete[] EstimationChiSq_Ntt_vs_Nw;
	delete hNbjetsEstSummary[NofJetsBins][0]; delete hNbjetsEstSummary[NofJetsBins][1];
	delete hNbjetsMCSummary[NofJetsBins][0];  delete hNbjetsMCSummary[NofJetsBins][1];

	for(int i=0;i<4;i++){
		delete hNjetsEstSummary[i][0]; delete hNjetsEstSummary[i][1];
		delete hNjetsMCSummary[i][0];  delete hNjetsMCSummary[i][1];
		delete hsNjets_MC[i];          delete hsNjets_Est[i];
		delete tCanva_Njets_Summary[i];
	}
	delete hNjetsEstSummary[4][0]; delete hNjetsEstSummary[4][1];
	delete hNjetsMCSummary[4][0];  delete hNjetsMCSummary[4][1];

	delete[] tCanva_Nbjets_Summary;
	delete[] tCanva_Njets_Summary;
	delete[] hNbjetsEstSummary;
	delete[] hNjetsEstSummary;
	delete[] hNbjetsMCSummary;
	delete[] hNjetsMCSummary;
	delete[] hsNjets_MC;
	delete[] hsNjets_Est;
	delete[] hsNbjets_MC;
	delete[] hsNbjets_Est;	delete[] histos;
	delete[] N1_;
	delete[] N2_;
	delete[] N3_;
	delete[] Njets_;
	delete[] N0bjets_;
	delete[] N1bjets_;
	delete[] N2bjets_;
	delete[] N3bjets_;
	delete[] MultiJet_Est_N0_;
	delete[] MultiJet_Est_N1_;
	delete[] MultiJet_Est_N2_;
	delete[] MultiJet_Est_N3_;
	delete[] eudsc_;
	delete[] eb_;
	delete[] Ntt_ ;
	delete[] Nw_;
	if(RescaledTTlikeEstimation) delete RescaledTTlikeEstimation;
	if(RescaledWlikeEstimation)  delete RescaledWlikeEstimation;
	if(tCanva_RescaledTTlikeEstimation) delete tCanva_RescaledTTlikeEstimation;
	if(tCanva_RescaledWlikeEstimation)  delete tCanva_RescaledWlikeEstimation;
}


WJetEstimation::WJetEstimation(const WJetEstimation& wj){
	NbOfDatasets = wj.NbOfDatasets;
        histos = 0;
	histos = new TH1F**[wj.NofJetsBins];
	//init summary histos
	hsNbjets_MC  = new THStack*[wj.NofJetsBins];
	hsNbjets_Est = new THStack*[wj.NofJetsBins];
	hsNjets_MC   = new THStack*[4];
	hsNjets_Est  = new THStack*[4];
	hNbjetsEstSummary = new TH1F**[wj.NofJetsBins];
	hNbjetsMCSummary  = new TH1F**[wj.NofJetsBins];
	hNjetsEstSummary  = new TH1F**[4];
	hNjetsMCSummary   = new TH1F**[4];

	tCanva_Nbjets_Summary = new TCanvas*[wj.NofJetsBins];
	tCanva_Njets_Summary  = new TCanvas*[4];
	
	for(int i=0;i<(wj.NofJetsBins+1);i++){
		hNbjetsEstSummary[i] = new TH1F*[2];
		hNbjetsEstSummary[i][0] = (TH1F*)wj.hNbjetsEstSummary[i][0]->Clone();
		hNbjetsEstSummary[i][1] = (TH1F*)wj.hNbjetsEstSummary[i][0]->Clone();

		hNbjetsMCSummary[i]  = new TH1F*[2];
		hNbjetsMCSummary[i][0] = (TH1F*)wj.hNbjetsMCSummary[i][0]->Clone();
		hNbjetsMCSummary[i][1] = (TH1F*)wj.hNbjetsMCSummary[i][0]->Clone();

		hsNbjets_MC[i]  = (THStack*)wj.hsNbjets_MC[i]->Clone();
		hsNbjets_Est[i] = (THStack*)wj.hsNbjets_Est[i]->Clone();
	}
	for(int i=0;i<wj.NofJetsBins;i++){
		histos[i] = new TH1F*[4];
		for(unsigned int j=0;j<4;j++) histos[i][j]=(TH1F*)wj.histos[i][j]->Clone();
	}
	for(int i=0;i<(4+1);i++){
		hNjetsEstSummary[i] = new TH1F*[2];
		hNjetsEstSummary[i][0] = (TH1F*)wj.hNjetsEstSummary[i][0]->Clone();
		hNjetsEstSummary[i][1] = (TH1F*)wj.hNjetsEstSummary[i][0]->Clone();

		hNjetsMCSummary[i]  = new TH1F*[2];
		hNjetsMCSummary[i][0] = (TH1F*)wj.hNjetsMCSummary[i][0]->Clone();
		hNjetsMCSummary[i][1] = (TH1F*)wj.hNjetsMCSummary[i][0]->Clone();

		hsNjets_MC[i]  = (THStack*)wj.hsNjets_MC[i]->Clone();
		hsNjets_Est[i] = (THStack*)wj.hsNjets_Est[i]->Clone();		
	}
	RescaledTTlikeEstimation        = (TGraphErrors*)wj.RescaledTTlikeEstimation->Clone();
	RescaledWlikeEstimation         = (TGraphErrors*)wj.RescaledTTlikeEstimation->Clone();
	tCanva_RescaledTTlikeEstimation = (TCanvas*)wj.tCanva_RescaledTTlikeEstimation->Clone();
	tCanva_RescaledWlikeEstimation  = (TCanvas*)wj.tCanva_RescaledWlikeEstimation->Clone();
        
	MyLeg = (TLegend*)wj.MyLeg->Clone();

        iDatasetsTTLike = wj.iDatasetsTTLike;
        iDatasetsWLike = wj.iDatasetsWLike;
        NofJetsBins = wj.NofJetsBins;
        Njets_ = wj.Njets_;
        N0bjets_ =  new double[wj.NofJetsBins];
        N1bjets_ =  new double[wj.NofJetsBins];
        N2bjets_ =  new double[wj.NofJetsBins];
        N3bjets_ =  new double[wj.NofJetsBins];
        MultiJet_Est_N0_ = new double[wj.NofJetsBins];
       	MultiJet_Est_N1_ = new double[wj.NofJetsBins];
       	MultiJet_Est_N2_ = new double[wj.NofJetsBins];
       	MultiJet_Est_N3_ = new double[wj.NofJetsBins];
	eudsc_ = new double[wj.NofJetsBins];
        eb_ = new double[wj.NofJetsBins];
        Ntt_  = new double[wj.NofJetsBins];
        Nw_ = new double[wj.NofJetsBins];
	for(int i=0;i<wj.NofJetsBins;i++){
		N0bjets_[i] = wj.N0bjets_[i];
        	N1bjets_[i] = wj.N1bjets_[i];
        	N2bjets_[i] = wj.N2bjets_[i];
        	N3bjets_[i] = wj.N3bjets_[i];
        	MultiJet_Est_N0_[i] = wj.MultiJet_Est_N0_[i];
       		MultiJet_Est_N1_[i] = wj.MultiJet_Est_N1_[i];
       		MultiJet_Est_N2_[i] = wj.MultiJet_Est_N2_[i];
       		MultiJet_Est_N3_[i] = wj.MultiJet_Est_N3_[i];
		eudsc_[i] = wj.eudsc_[i];
        	eb_[i] = wj.eb_[i];
        	Ntt_ [i] = wj.Ntt_[i];
        	Nw_[i] = wj.Nw_[i];
	}
        N0_ = new double*[wj.NofJetsBins];
        N1_ = new double*[wj.NofJetsBins];
        N2_ = new double*[wj.NofJetsBins];
        N3_ = new double*[wj.NofJetsBins];
        for(int i=0;i<wj.NofJetsBins;i++){
        	N0_[i] = new double[wj.NbOfDatasets];
        	N1_[i] = new double[wj.NbOfDatasets];
        	N2_[i] = new double[wj.NbOfDatasets];
        	N3_[i] = new double[wj.NbOfDatasets];
		for(int j=0;j<wj.NbOfDatasets;j++){
        		N0_[i][j] = wj.N0_[i][j];
        		N1_[i][j] = wj.N1_[i][j];
        		N2_[i][j] = wj.N2_[i][j];
        		N3_[i][j] = wj.N3_[i][j];
        	}
	}
	//init histo for the stability check
	for(int i=0;i<wj.NofJetsBins;i++){
		if(wj.EstimationChiSq_Eb_vs_Eudsc[i]) EstimationChiSq_Eb_vs_Eudsc[i] = wj.EstimationChiSq_Eb_vs_Eudsc[i];
		if(wj.EstimationChiSq_Ntt_vs_Nw[i])   EstimationChiSq_Ntt_vs_Nw[i]   = wj.EstimationChiSq_Ntt_vs_Nw[i];
	}
}

void WJetEstimation::ReScale(float factor){
        for(int i=0;i<NofJetsBins;i++){
		for(int j=0;j<NbOfDatasets;j++){
			N0_[i][j] *= factor;
			N1_[i][j] *= factor;
			N2_[i][j] *= factor;
			N3_[i][j] *= factor;
		}
		N0bjets_[i] *= factor;
		N1bjets_[i] *= factor;
		N2bjets_[i] *= factor;
		N3bjets_[i] *= factor;
		MultiJet_Est_N0_[i] *= factor;
		MultiJet_Est_N1_[i] *= factor;
		MultiJet_Est_N2_[i] *= factor;
		MultiJet_Est_N3_[i] *= factor;
	}
}

// void WJetEstimation::FillInputs( int NofDatasets, double** n0, double** n1, double** n2, double** n3, double** MultiJets_Estimated_Nj, vector<int> iDTTLike, vector<int> iDWLike){
// 	NbOfDatasets = NofDatasets;
// 	iDatasetsTTLike.clear();
// 	iDatasetsWLike.clear();
// 	iDatasetsTTLike = iDTTLike;
// 	iDatasetsWLike = iDWLike;
// 	
// 	// Multi-jets/tt+jets/W+jets/Z+jets/st-s/st-t/st-tW
// 	for(int i=0;i<NofJetsBins;i++){
// 		N0_[i] = new double[NbOfDatasets];
// 		N1_[i] = new double[NbOfDatasets];
// 		N2_[i] = new double[NbOfDatasets];
// 		N3_[i] = new double[NbOfDatasets];
// 		for(int j=0;j<NofDatasets;j++){
// 			N0_[i][j] = n0[i][j];
// 			N1_[i][j] = n1[i][j];
// 			N2_[i][j] = n2[i][j];
// 			N3_[i][j] = n3[i][j];
// 		}
// 	}
// 	
// 	for(int i=0;i<NofJetsBins;i++){
// 		N0bjets_[i]  = 0;
// 		N1bjets_[i]  = 0;
// 		N2bjets_[i]  = 0;
// 		N3bjets_[i]  = 0;
// 		for(int j=0;j<NbOfDatasets;j++)
// 		{
// 			N0bjets_[i]  += N0_[i][j];
// 			N1bjets_[i]  += N1_[i][j];
// 			N2bjets_[i]  += N2_[i][j];
// 			N3bjets_[i]  += N3_[i][j];
// 		}
// 		N0bjets_[i]   += -MultiJets_Estimated_Nj[i][0];
// 		N1bjets_[i]   += -MultiJets_Estimated_Nj[i][1];
// 		N2bjets_[i]   += -MultiJets_Estimated_Nj[i][2];
// 		N3bjets_[i]   += -MultiJets_Estimated_Nj[i][3];
// 	}
// }

void WJetEstimation::FillInputs(double** n0, double** n1, double** n2, double** n3, double** MultiJets_Estimated_Nj){
	
	// Multi-jets/tt+jets/W+jets/Z+jets/st-s/st-t/st-tW
	for(int i=0;i<NofJetsBins;i++){

		N0bjets_[i]  = 0;
		N1bjets_[i]  = 0;
		N2bjets_[i]  = 0;
		N3bjets_[i]  = 0;

		for(int j=0;j<NbOfDatasets;j++){

			N0_[i][j] = n0[i][j];
			N1_[i][j] = n1[i][j];
			N2_[i][j] = n2[i][j];
			N3_[i][j] = n3[i][j];

			N0bjets_[i]  += N0_[i][j];
			N1bjets_[i]  += N1_[i][j];
			N2bjets_[i]  += N2_[i][j];
			N3bjets_[i]  += N3_[i][j];
		}
		N0bjets_[i]   += -MultiJets_Estimated_Nj[i][0];
		N1bjets_[i]   += -MultiJets_Estimated_Nj[i][1];
		N2bjets_[i]   += -MultiJets_Estimated_Nj[i][2];
		N3bjets_[i]   += -MultiJets_Estimated_Nj[i][3];
	}
}
void WJetEstimation::FillInputs(double** n0, double** n1, double** n2, double** n3){
	
	// Multi-jets/tt+jets/W+jets/Z+jets/st-s/st-t/st-tW
	for(int i=0;i<NofJetsBins;i++){

		N0bjets_[i]  = 0;
		N1bjets_[i]  = 0;
		N2bjets_[i]  = 0;
		N3bjets_[i]  = 0;

		for(int j=0;j<NbOfDatasets;j++){

			N0_[i][j] = n0[i][j];
			N1_[i][j] = n1[i][j];
			N2_[i][j] = n2[i][j];
			N3_[i][j] = n3[i][j];

			N0bjets_[i]  += N0_[i][j];
			N1bjets_[i]  += N1_[i][j];
			N2bjets_[i]  += N2_[i][j];
			N3bjets_[i]  += N3_[i][j];
		}
		N0bjets_[i]   += -MultiJet_Est_N0_[i];
		N1bjets_[i]   += -MultiJet_Est_N1_[i];
		N2bjets_[i]   += -MultiJet_Est_N2_[i];
		N3bjets_[i]   += -MultiJet_Est_N3_[i];
	}
}

void WJetEstimation::PrintInputs(){
	cout<<endl;
	cout<<"***********************************************"<<endl;
	cout<<"*****      WJets Estimation: input values    ***"<<endl;
	cout<<"***********************************************"<<endl;
	cout<<"Nof datasets: "<<NbOfDatasets<<endl;
	for(int i=0;i<NofJetsBins;i++) PrintInputs(i);
	cout<<"***********************************************"<<endl<<endl;;
}

void WJetEstimation::PrintInputs(int Njets){
	cout<<"Initial values : N (0/1/2/3 b-jets) = "<<N0bjets_[Njets]<<"/"<<N1bjets_[Njets]<<"/"<<N2bjets_[Njets]<<"/"<<N3bjets_[Njets]<<endl;
	for(int i=0;i<NbOfDatasets;i++){
		cout<<"Dataset "<<i+1<<" : N (0/1/2/3 jets) = "<<N0_[Njets][i]<<"/"<<N1_[Njets][i]<<"/"<<N2_[Njets][i]<<"/"<<N3_[Njets][i]<<endl;
	}
}


void WJetEstimation::Write(TFile* file, string label){
	file->cd();
	string dirname = "WJetEstimation"+label;
	file->mkdir(dirname.c_str());
	file->cd(dirname.c_str());

	char name[100];char title[100];
	
	for(int i=0;i<NofJetsBins;i++){
		for(int j=0;j<4;j++) histos[i][j]->Write();
	}
	for(int i=0;i<(NofJetsBins+1);i++){
		
		// For Monte Carlo prediction 
		MyLeg->Clear();
		tCanva_Nbjets_Summary[i]->cd();
		
		//hNbjetsMCSummary[i][0]->Write();
		//hNbjetsMCSummary[i][1]->Write();

		if(i!=NofJetsBins) sprintf(name,"hNbjetsMCStackSummary_%d_jets",Njets_[i]);
		else               sprintf(name,"hNbjetsMCStackSummary_Inclusive");
		if(i!=NofJetsBins) sprintf(title,"W-like and TT-like MC prediction summary for %d jets",Njets_[i]);
		else               sprintf(title,"W-like and TT-like MC prediction summary (inclusive)");

		hsNbjets_MC[i] = new THStack(name,title);

		hNbjetsMCSummary[i][0]->SetFillStyle(3004);
		hNbjetsMCSummary[i][0]->SetFillColor(kGreen+1);
		hNbjetsMCSummary[i][1]->SetFillStyle(3005);
		hNbjetsMCSummary[i][1]->SetFillColor(kRed+1);

		hsNbjets_MC[i]->Add(hNbjetsMCSummary[i][0]);
		hsNbjets_MC[i]->Add(hNbjetsMCSummary[i][1]);
		
		MyLeg->AddEntry(hNbjetsMCSummary[i][0],"W-like events","f");
		MyLeg->AddEntry(hNbjetsMCSummary[i][1],"TT-like events","f");
		
		hsNbjets_MC[i]->Draw();
		hsNbjets_MC[i]->GetXaxis()->SetTitle("Nb of b-jets");

		//hNbjetsEstSummary[i][0]->Write();
		//hNbjetsEstSummary[i][1]->Write();

		if(i!=NofJetsBins) sprintf(name,"hNbjetsEstStackSummary_%d_jets",Njets_[i]);
		else               sprintf(name,"hNbjetsEstStackSummary_Inclusive");
		if(i!=NofJetsBins) sprintf(title,"W-like and TT-like estimation summary for %d jets",Njets_[i]);
		else               sprintf(title,"W-like and TT-like estimation summary (inclusive)");

		hsNbjets_Est[i] = new THStack(name,title);

		hNbjetsEstSummary[i][0]->SetMarkerColor(kYellow+1);
		hNbjetsEstSummary[i][1]->SetMarkerColor(kBlue+1);
		hNbjetsEstSummary[i][0]->SetMarkerStyle(22);
		hNbjetsEstSummary[i][1]->SetMarkerStyle(23);
		hNbjetsEstSummary[i][0]->SetMarkerSize(1.5);
		hNbjetsEstSummary[i][1]->SetMarkerSize(1.5);

		hsNbjets_Est[i]->Add(hNbjetsEstSummary[i][0]);
		hsNbjets_Est[i]->Add(hNbjetsEstSummary[i][1]);
		
		MyLeg->AddEntry(hNbjetsEstSummary[i][0],"W-like estimation","p");
		MyLeg->AddEntry(hNbjetsEstSummary[i][1],"TT-like estimation","p");
		
		//hsNbjets_Est[i]->Write();
		//MyLeg->Write();
		
		//tCanva_Nbjets_Summary[i]->cd();
		//hsNbjets_MC[i]->Draw();
		//hsNbjets_MC[i]->GetXaxis()->SetTitle("Nb of b-jets");
		hsNbjets_Est[i]->Draw("PEsame");
		hsNbjets_Est[i]->GetXaxis()->SetTitle("Nb of b-jets");
		MyLeg->Draw("same");

		tCanva_Nbjets_Summary[i]->Update();
		cout<<"Writing summary histograms for "<<i+4<<" jets"<<endl;
		tCanva_Nbjets_Summary[i]->Write();
	}
	for(int i=0;i<(4+1);i++){

		// For Monte Carlo prediction 
		MyLeg->Clear();
		tCanva_Njets_Summary[i]->cd();
		
		//hNjetsMCSummary[i][0] ->Write();
		//hNjetsMCSummary[i][1] ->Write();

		if(i!=4) sprintf(name,"hNjetsMCStackSummary_%d_bjets",i);
		else     sprintf(name,"hNjetsMCStackSummary_bInclusive");
		if(i!=4) sprintf(title,"W-like and TT-like MC prediction summary for %d b-jets",i);
		else     sprintf(title,"W-like and TT-like MC prediction summary (b-inclusive)");
		hsNjets_MC[i] = new THStack(name,title);

		hNjetsMCSummary[i][0]->SetFillStyle(3004);
		hNjetsMCSummary[i][0]->SetFillColor(kGreen+1);
		hNjetsMCSummary[i][1]->SetFillStyle(3005);
		hNjetsMCSummary[i][1]->SetFillColor(kRed+1);

		hsNjets_MC[i]->Add(hNjetsMCSummary[i][0]);
		hsNjets_MC[i]->Add(hNjetsMCSummary[i][1]);
		
		MyLeg->AddEntry(hNjetsMCSummary[i][0],"W-like events","f");
		MyLeg->AddEntry(hNjetsMCSummary[i][1],"TT-like events","f");
		
		hsNjets_MC[i]->Draw();
		hsNjets_MC[i]->GetXaxis()->SetTitle("Nb of jets");

		//hNjetsEstSummary[i][0]->Write();
		//hNjetsEstSummary[i][1]->Write();

		if(i!=4) sprintf(name,"hNjetsEstStackSummary_%d_bjets",i);
		else     sprintf(name,"hNjetsEstStackSummary_bInclusive");
		if(i!=4) sprintf(title,"W-like and TT-like estimation summary for %d b-jets",i);
		else     sprintf(title,"W-like and TT-like estimation summary (b-inclusive)");
		hsNjets_Est[i] = new THStack(name,title);

		//sprintf(name,"hNjetsEstStackLegend_%d_bjets",i);
		//MyLeg->SetName(name);

		hNjetsEstSummary[i][0]->SetMarkerColor(kYellow+1);
		hNjetsEstSummary[i][1]->SetMarkerColor(kBlue+1);
		hNjetsEstSummary[i][0]->SetMarkerStyle(22);
		hNjetsEstSummary[i][1]->SetMarkerStyle(23);
		hNjetsEstSummary[i][0]->SetMarkerSize(1.5);
		hNjetsEstSummary[i][1]->SetMarkerSize(1.5);

		hsNjets_Est[i]->Add(hNjetsEstSummary[i][0]);
		hsNjets_Est[i]->Add(hNjetsEstSummary[i][1]);
		
		MyLeg->AddEntry(hNjetsEstSummary[i][0],"W-like estimation","p");
		MyLeg->AddEntry(hNjetsEstSummary[i][1],"TT-like estimation","p");
		
		//hsNjets_Est[i]->Write();
		//MyLeg->Write();

		//tCanva_Njets_Summary[i]->cd();
		//hsNjets_MC[i]->Draw();
		//hsNjets_MC[i]->GetXaxis()->SetTitle("Nb of jets");
		hsNjets_Est[i]->Draw("PE SAME");
		hsNjets_Est[i]->GetXaxis()->SetTitle("Nb of jets");
		MyLeg->Draw("same");

		tCanva_Njets_Summary[i]->Update();
		tCanva_Njets_Summary[i]->Write();
		cout<<"Writing summary histograms for "<<i<<" b-jets"<<endl;
	}

	if(RescaledTTlikeEstimation)  RescaledTTlikeEstimation->Write();
	if(RescaledWlikeEstimation)   RescaledWlikeEstimation->Write();
	if(tCanva_RescaledTTlikeEstimation) tCanva_RescaledTTlikeEstimation->Write();
	if(tCanva_RescaledWlikeEstimation)  tCanva_RescaledWlikeEstimation->Write();
	for(int i=0;i<NofJetsBins;i++){
		if(EstimationChiSq_Eb_vs_Eudsc[i]) EstimationChiSq_Eb_vs_Eudsc[i]->Write();
		if(EstimationChiSq_Ntt_vs_Nw[i])   EstimationChiSq_Ntt_vs_Nw[i]->Write();
	}
}
							


// void WJetEstimation::Estimation(int Iterations, bool Print, bool DoFillSummary)
// {
// 	//Initialisation of histos
// 	for(int i=0;i<NofJetsBins;i++){
// 		char name[100];char title[100];
// 		sprintf(name,"Eudsc_%d_jets",Njets_[i]);
// 		sprintf(title,"Eudsc vs iterations for %d jets",Njets_[i]);
// 		histos[i][0] = new TH1F(name,title,Iterations,1,Iterations+1);
// 		sprintf(name,"Eb_%d_jets",Njets_[i]);
// 		sprintf(title,"Eb vs iterations for %d jets",Njets_[i]);
// 		histos[i][1] = new TH1F(name,title,Iterations,1,Iterations+1);
// 		sprintf(name,"Ntt_%d_jets",Njets_[i]);
// 		sprintf(title,"Ntt vs iterations for %d jets",Njets_[i]);
// 		histos[i][2] = new TH1F(name,title,Iterations,1,Iterations+1);
// 		sprintf(name,"Nw_%d_jets",Njets_[i]);
// 		sprintf(title,"Nw vs iterations for %d jets",Njets_[i]);
// 		histos[i][3] = new TH1F(name,title,Iterations,1,Iterations+1);
// 		
// 		Initialisation(Print, i);
// 		Iteration(Iterations, Print, i);
// 	}
// 	if(DoFillSummary) FillSummaryHistos();
// }
void WJetEstimation::Estimation(bool Print, bool DoFillSummary)
{
//	int It = 40;
	for(int i=0;i<NofJetsBins;i++){
		Initialisation(Print, i);
		Iteration(It_,Print, i);
	}
	if(DoFillSummary && MCdata_) FillSummaryHistos();
}
void WJetEstimation::FillSummaryHistos(){

	double NWLike[4][NofJetsBins];
	double NWLikeNbjets[4] = {0,0,0,0};
	double NWLikeError = 0;
	double NttLike[4][NofJetsBins];
	double NttLikeNbjets[4] = {0,0,0,0};
	double NttLikeError = 0;
	TString NjetsXLabel[3]  = {"4","5","#geq 6"};
	TString NbjetsXLabel[4] = {"0","1","2","#geq 3"};
	
	for(int Njets=0;Njets<NofJetsBins;Njets++){
		
		for(int i=0;i<4;i++){
			NWLike[i][Njets]=0;
			NttLike[i][Njets]=0;
		}
		
		for(unsigned int i=0;i<iDatasetsWLike.size();i++){
			NWLike[0][Njets] += N0_[Njets][iDatasetsWLike[i]]; NWLikeNbjets[0] += N0_[Njets][iDatasetsWLike[i]]; 
			NWLike[1][Njets] += N1_[Njets][iDatasetsWLike[i]]; NWLikeNbjets[1] += N1_[Njets][iDatasetsWLike[i]];
			NWLike[2][Njets] += N2_[Njets][iDatasetsWLike[i]]; NWLikeNbjets[2] += N2_[Njets][iDatasetsWLike[i]];
			NWLike[3][Njets] += N3_[Njets][iDatasetsWLike[i]]; NWLikeNbjets[3] += N3_[Njets][iDatasetsWLike[i]]; 
		}
		for(unsigned int i=0;i<iDatasetsTTLike.size();i++){
			NttLike[0][Njets] += N0_[Njets][iDatasetsTTLike[i]]; NttLikeNbjets[0] += N0_[Njets][iDatasetsTTLike[i]];
			NttLike[1][Njets] += N1_[Njets][iDatasetsTTLike[i]]; NttLikeNbjets[1] += N1_[Njets][iDatasetsTTLike[i]];
			NttLike[2][Njets] += N2_[Njets][iDatasetsTTLike[i]]; NttLikeNbjets[2] += N2_[Njets][iDatasetsTTLike[i]];
			NttLike[3][Njets] += N3_[Njets][iDatasetsTTLike[i]]; NttLikeNbjets[3] += N3_[Njets][iDatasetsTTLike[i]]; 			
		}
		
		for(int i=0;i<4;i++){
			// For W-like and TT-like estimation
			hNbjetsEstSummary[Njets][0]->SetBinContent(i+1,GetNWEventsEstimatedNbjets( i, Njets));
			if(statisticalWLikeError_[Njets][i]!=0)  hNbjetsEstSummary[Njets][0]->SetBinError(i+1,statisticalWLikeError_[Njets][i]);
			hNbjetsEstSummary[Njets][1]->SetBinContent(i+1,GetNttEventsEstimatedNbjets(i, Njets));
			if(statisticalTTLikeError_[Njets][i]!=0) hNbjetsEstSummary[Njets][1]->SetBinError(i+1,statisticalTTLikeError_[Njets][i]);
			hNbjetsEstSummary[Njets][0]->GetXaxis()->SetBinLabel(i+1,NbjetsXLabel[i]);
			hNbjetsEstSummary[Njets][1]->GetXaxis()->SetBinLabel(i+1,NbjetsXLabel[i]);
			// For W-like and TT-like MC prediction
			hNbjetsMCSummary[Njets][0]->SetBinContent(i+1,NWLike[i][Njets]);
			hNbjetsMCSummary[Njets][1]->SetBinContent(i+1,NttLike[i][Njets]);
			hNbjetsMCSummary[Njets][0]->GetXaxis()->SetBinLabel(i+1,NbjetsXLabel[i]);
			hNbjetsMCSummary[Njets][1]->GetXaxis()->SetBinLabel(i+1,NbjetsXLabel[i]);
		}
		hNbjetsEstSummary[Njets][0]->GetXaxis()->CenterLabels();hNbjetsEstSummary[Njets][0]->GetXaxis()->SetTitle("Nb of b-jets");
		hNbjetsEstSummary[Njets][1]->GetXaxis()->CenterLabels();hNbjetsEstSummary[Njets][1]->GetXaxis()->SetTitle("Nb of b-jets");
		hNbjetsMCSummary[Njets][0]->GetXaxis()->CenterLabels();hNbjetsMCSummary[Njets][0]->GetXaxis()->SetTitle("Nb of b-jets");
		hNbjetsMCSummary[Njets][1]->GetXaxis()->CenterLabels();hNbjetsMCSummary[Njets][1]->GetXaxis()->SetTitle("Nb of b-jets");
	}
	NWLikeError = 0; NttLikeError = 0;
	for(int i=0;i<4;i++){
		// For W-like and TT-like estimation
		hNbjetsEstSummary[NofJetsBins][0]->SetBinContent(i+1,GetNWEventsEstimatedNbjets(i));
		for(int j=0;j<NofJetsBins;j++) NWLikeError += pow(statisticalWLikeError_[j][i],2);
		hNbjetsEstSummary[NofJetsBins][0]->SetBinError(i+1,sqrt(NWLikeError));

		hNbjetsEstSummary[NofJetsBins][1]->SetBinContent(i+1,GetNttEventsEstimatedNbjets(i));
		for(int j=0;j<NofJetsBins;j++) NttLikeError += pow(statisticalTTLikeError_[j][i],2);
		hNbjetsEstSummary[NofJetsBins][1]->SetBinError(i+1,sqrt(NttLikeError));

		hNbjetsEstSummary[NofJetsBins][0]->GetXaxis()->SetBinLabel(i+1,NbjetsXLabel[i]);
		hNbjetsEstSummary[NofJetsBins][1]->GetXaxis()->SetBinLabel(i+1,NbjetsXLabel[i]);
		// For W-like and TT-like MC prediction
		hNbjetsMCSummary[NofJetsBins][0]->SetBinContent(i+1,NWLikeNbjets[i]);
		hNbjetsMCSummary[NofJetsBins][1]->SetBinContent(i+1,NttLikeNbjets[i]);
		hNbjetsMCSummary[NofJetsBins][0]->GetXaxis()->SetBinLabel(i+1,NbjetsXLabel[i]);
		hNbjetsMCSummary[NofJetsBins][1]->GetXaxis()->SetBinLabel(i+1,NbjetsXLabel[i]);
	}
	hNbjetsEstSummary[NofJetsBins][0]->GetXaxis()->CenterLabels();hNbjetsEstSummary[NofJetsBins][0]->GetXaxis()->SetTitle("Nb of b-jets");
	hNbjetsEstSummary[NofJetsBins][1]->GetXaxis()->CenterLabels();hNbjetsEstSummary[NofJetsBins][1]->GetXaxis()->SetTitle("Nb of b-jets");
	hNbjetsMCSummary[NofJetsBins][0]->GetXaxis()->CenterLabels();hNbjetsMCSummary[NofJetsBins][0]->GetXaxis()->SetTitle("Nb of b-jets");
	hNbjetsMCSummary[NofJetsBins][1]->GetXaxis()->CenterLabels();hNbjetsMCSummary[NofJetsBins][1]->GetXaxis()->SetTitle("Nb of b-jets");

	for(int Nbjets=0;Nbjets<4;Nbjets++){
		// For W-like and TT-like estimation
		for(int i=0;i<NofJetsBins;i++){
			// For W-like and TT-like estimation
			hNjetsEstSummary[Nbjets][0]->SetBinContent(i+1,GetNWEventsEstimatedNbjets( Nbjets, i));
			if(statisticalWLikeError_[i][Nbjets]!=0)  hNjetsEstSummary[Nbjets][0]->SetBinError(i+1,statisticalWLikeError_[i][Nbjets]);
			hNjetsEstSummary[Nbjets][1]->SetBinContent(i+1,GetNttEventsEstimatedNbjets(Nbjets, i));
			if(statisticalTTLikeError_[i][Nbjets]!=0) hNjetsEstSummary[Nbjets][1]->SetBinError(i+1,statisticalTTLikeError_[i][Nbjets]);
			hNjetsEstSummary[Nbjets][0]->GetXaxis()->SetBinLabel(i+1,NjetsXLabel[i]);
			hNjetsEstSummary[Nbjets][1]->GetXaxis()->SetBinLabel(i+1,NjetsXLabel[i]);
			// For W-like and TT-like MC prediction
			hNjetsMCSummary[Nbjets][0]->SetBinContent(i+1,NWLike[Nbjets][i]);
			hNjetsMCSummary[Nbjets][1]->SetBinContent(i+1,NttLike[Nbjets][i]);
			hNjetsMCSummary[Nbjets][0]->GetXaxis()->SetBinLabel(i+1,NjetsXLabel[i]);
			hNjetsMCSummary[Nbjets][1]->GetXaxis()->SetBinLabel(i+1,NjetsXLabel[i]);
		}
		hNjetsEstSummary[Nbjets][0]->GetXaxis()->CenterLabels();hNjetsEstSummary[Nbjets][0]->GetXaxis()->SetTitle("Nb of jets");
		hNjetsEstSummary[Nbjets][1]->GetXaxis()->CenterLabels();hNjetsEstSummary[Nbjets][1]->GetXaxis()->SetTitle("Nb of jets");
		hNjetsMCSummary[Nbjets][0] ->GetXaxis()->CenterLabels();hNjetsMCSummary[Nbjets][0]->GetXaxis()->SetTitle("Nb of jets");
		hNjetsMCSummary[Nbjets][1] ->GetXaxis()->CenterLabels();hNjetsMCSummary[Nbjets][1]->GetXaxis()->SetTitle("Nb of jets");
	}
	for(int i=0;i<NofJetsBins;i++){
		// For W-like and TT-like estimation
		hNjetsEstSummary[4][0]->SetBinContent(i+1,GetNWEventsEstimated(i));
		if(statisticalWLikeError_[i][4]!=0)  hNjetsEstSummary[4][0]->SetBinError(i+1,statisticalWLikeError_[i][4]);
		hNjetsEstSummary[4][1]->SetBinContent(i+1,GetNttEventsEstimated(i));
		if(statisticalTTLikeError_[i][4]!=0) hNjetsEstSummary[4][1]->SetBinError(i+1,statisticalTTLikeError_[i][4]);
		hNjetsEstSummary[4][0]->GetXaxis()->SetBinLabel(i+1,NjetsXLabel[i]);
		hNjetsEstSummary[4][1]->GetXaxis()->SetBinLabel(i+1,NjetsXLabel[i]);
		// For W-like and TT-like MC prediction
		hNjetsMCSummary[4][0]->SetBinContent(i+1,GetNWlikeMC(i));
		hNjetsMCSummary[4][1]->SetBinContent(i+1,GetNTTlikeMC(i));
		hNjetsMCSummary[4][0]->GetXaxis()->SetBinLabel(i+1,NjetsXLabel[i]);
		hNjetsMCSummary[4][1]->GetXaxis()->SetBinLabel(i+1,NjetsXLabel[i]);
	}
	hNjetsEstSummary[4][0]->GetXaxis()->CenterLabels(); hNjetsEstSummary[4][0]->GetXaxis()->SetTitle("Nb of jets");
	hNjetsEstSummary[4][1]->GetXaxis()->CenterLabels(); hNjetsEstSummary[4][1]->GetXaxis()->SetTitle("Nb of jets");
	hNjetsMCSummary[4][0] ->GetXaxis()->CenterLabels(); hNjetsMCSummary[4][0] ->GetXaxis()->SetTitle("Nb of jets");
	hNjetsMCSummary[4][1] ->GetXaxis()->CenterLabels(); hNjetsMCSummary[4][1] ->GetXaxis()->SetTitle("Nb of jets");
}
void WJetEstimation::Print4Estimators(){
	for(int i=0;i<NofJetsBins;i++) Print4Estimators(i);
}

void WJetEstimation::Print4Estimators(int Njets){
	if(Njets>=NofJetsBins) return;
	cout<<"For N = "<<Njets_[Njets]<<" jets"<<endl;
	cout<<"E(b) = "<<eb_[Njets]<<"\t E(udsc) = "<<eudsc_[Njets]<<"\t N(ttbar) = "<<Ntt_[Njets]<<"\t N(Wjets) = "<<Nw_[Njets]<<endl;
}

void WJetEstimation::CheckConsistency(){
	for(int i=0;i<NofJetsBins;i++) CheckConsistency(i); 
}

void WJetEstimation::CheckConsistency(int Njets){
	cout<<"***********************************************************"<<endl;
	if(Njets>=NofJetsBins) return;
	cout<<"For N = "<<Njets_[Njets]<<" jets"<<endl;
	cout<<"N0bjet : input:  "<<N0bjets_[Njets]<<" - output: "<<WJetEstimation::N0bjet(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])<<" - diff: "<<(WJetEstimation::N0bjet(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])-N0bjets_[Njets])<<endl;
	cout<<"N1bjet : input:  "<<N1bjets_[Njets]<<" - output: "<<WJetEstimation::N1bjet(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])<<" - diff: "<<(WJetEstimation::N1bjet(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])-N1bjets_[Njets])<<endl;
	cout<<"N2bjet : input:  "<<N2bjets_[Njets]<<" - output: "<<WJetEstimation::N2bjets(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])<<" - diff: "<<(WJetEstimation::N2bjets(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])-N2bjets_[Njets])<<endl;
	cout<<"N3bjet : input:  "<<N3bjets_[Njets]<<" - output: "<<WJetEstimation::N3bjets(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])<<" - diff: "<<(WJetEstimation::N3bjets(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])-N3bjets_[Njets])<<endl;
	cout<<"***********************************************************"<<endl;
}

bool WJetEstimation::CheckEstimation(int Njets){
	bool ValidEstimate = true;
	if(Njets>=NofJetsBins) return ValidEstimate;
	if(GetNttEventsEstimated(Njets) == 0 || GetNWEventsEstimated(Njets) == 0) ValidEstimate = false;
	return ValidEstimate;
}

bool WJetEstimation::CheckEstimation(bool doCheckByNbjet){
	bool ValidEstimate = true;
	if(!doCheckByNbjet){
		for(int i=0;i<NofJetsBins;i++){
			if((GetNttEventsEstimated(i) == 0 && GetNTTlikeMC(i) >1 )
			|| (GetNWEventsEstimated(i)  == 0 && GetNWlikeMC(i)  >1 )) ValidEstimate = false;
		}
	}
	else{
		for(int i=0;i<NofJetsBins;i++){
			for(int j=0;j<3;j++){
				if(GetNttEventsEstimatedNbjets(j,i) == 0 
				|| GetNWEventsEstimatedNbjets(j,i)  == 0) ValidEstimate = false;
			}
		}
	}
	return ValidEstimate;
}

double WJetEstimation::EstimationChiSq(int Njets){
	double ChiSq = 9999;
	if(GetNTTlikeMC(Njets)>0 && GetNWlikeMC(Njets)>0){
		ChiSq = 0.5*((pow(GetNttEventsEstimated(Njets)-GetNTTlikeMC(Njets),2)/GetNTTlikeMC(Njets)) + (pow(GetNWEventsEstimated(Njets)-GetNWlikeMC(Njets),2)/GetNWlikeMC(Njets)) );
	}
	return ChiSq;
}
double WJetEstimation::EstimationChiSq(){
	double ChiSq = 9999;
	for(int i=0;i<NofJetsBins;i++)	ChiSq += EstimationChiSq(i);
	return ChiSq;
}

void WJetEstimation::CheckEstimationStability(double ebRange_Percent, int ebRange_Nbins, double eudscRange_Percent, int eudscRange_Nbins, double NttRange_NSigma, int NttRange_Nbins, double NwRange_NSigma, int NwRange_Nbins, bool Verbose){
	int It = 40;
	bool Print = false;
	double ChiSq[NofJetsBins];

	double eb_Init[NofJetsBins];
	double eb_MinChiSq[NofJetsBins];
	double eudsc_Init[NofJetsBins];
	double eudsc_MinChiSq[NofJetsBins];
	double Ntt_Init[NofJetsBins];
	double Ntt_MinChiSq[NofJetsBins];
	double Nw_Init[NofJetsBins];
	double Nw_MinChiSq[NofJetsBins];

	char name[100];

	for(int i=0;i<NofJetsBins;i++){
		Initialisation(Print, i);
		ChiSq[i] = 9999;
		// Get initial values
		eb_Init[i]    = eb_[i];
		eudsc_Init[i] = eudsc_[i];
		Ntt_Init[i]   = Ntt_[i];
		Nw_Init[i]    = Nw_[i];
	}
	//Modify the initial values according to the users inputs :
	// Range for :
	// eb    = (1-X)eb,(1+X)eb;
	// eudsc = (1-X)eudsc,(1+X)eudsc;
	// Ntt   = (Ntt-nSigma),(Ntt+nSigma); Sigma = sqrt(Ntt)
	// Nw    = (Nw -nSigma),(Nw +nSigma); Sigma = sqrt(Nw)
	cout<<"*****************Stability check*****************"<<endl;
	cout<<"--- Initial/MinChiSq values :"<<endl;
	for(int i=0;i<NofJetsBins;i++){
		sprintf(name,"EstimationChiSq_Eb_vs_Eudsc_%d_jets",Njets_[i]);
		EstimationChiSq_Eb_vs_Eudsc[i] = new TH2D(name,"",ebRange_Nbins,(1-ebRange_Percent)*eb_Init[i],(1+ebRange_Percent)*eb_Init[i],eudscRange_Nbins,(1-eudscRange_Percent)*eudsc_Init[i],(1+eudscRange_Percent)*eudsc_Init[i]);
		sprintf(name,"EstimationChiSq_Ntt_vs_Nw_%d_jets",Njets_[i]);
		EstimationChiSq_Ntt_vs_Nw[i]   = new TH2D(name,"",NttRange_Nbins,Ntt_Init[i]-(NttRange_NSigma*sqrt(Ntt_Init[i])),Ntt_Init[i]+(NttRange_NSigma*sqrt(Ntt_Init[i])),NwRange_Nbins,Nw_Init[i]-(NwRange_NSigma*sqrt(Nw_Init[i])),Nw_Init[i]+(NwRange_NSigma*sqrt(Nw_Init[i])));
		for(int j=0;j<ebRange_Nbins;j++){
			for(int k=0;k<eudscRange_Nbins;k++){
				for(int l=0;l<NttRange_Nbins;l++){
					for(int m=0;m<NwRange_Nbins;m++){
						eb_[i]    = (2*ebRange_Percent*eb_Init[i]/ebRange_Nbins)*j;
						eudsc_[i] = (2*eudscRange_Percent*eudsc_Init[i]/eudscRange_Nbins)*k;
						Ntt_[i]   = (2*NttRange_NSigma*sqrt(Ntt_Init[i])/NttRange_Nbins)*l;
						Nw_[i]    = (2*NwRange_NSigma*sqrt(Nw_Init[i])/NwRange_Nbins)*m;
						Iteration(It,Print,i);
						if(Verbose){
							cout<<"eb/eudsc/Ntt/Nw : "<<eb_[i]<<" / "<<eudsc_[i]<<" / "<<Ntt_[i]<<" / "<<Nw_[i]<<endl;
							cout<<"-- associated ChiSq = "<<EstimationChiSq(i)<<endl;
						}
					}
					if(EstimationChiSq(i)<ChiSq[i]){
						ChiSq[i] = EstimationChiSq(i);
						eb_MinChiSq[i]    = eb_[i];
						eudsc_MinChiSq[i] = eudsc_[i];
						Ntt_MinChiSq[i]   = Ntt_[i];
						Nw_MinChiSq[i]    = Nw_[i];
					}
				}
			}
		}
		for(int j=0;j<ebRange_Nbins;j++){
			for(int k=0;k<eudscRange_Nbins;k++){
				Ntt_[i] = Ntt_MinChiSq[i];
				Nw_[i]  = Nw_MinChiSq[i];
				eb_[i]    = (2*ebRange_Percent*eb_Init[i]/ebRange_Nbins)*j;
				eudsc_[i] = (2*eudscRange_Percent*eudsc_Init[i]/eudscRange_Nbins)*k;
				Iteration(It,Print,i);
				EstimationChiSq_Eb_vs_Eudsc[i]->SetBinContent(j,k,EstimationChiSq(i));
			}
		}
		for(int l=0;l<NttRange_Nbins;l++){
			for(int m=0;m<NwRange_Nbins;m++){
				eb_[i]    = eb_MinChiSq[i];
				eudsc_[i] = eudsc_MinChiSq[i];
				Ntt_[i]   = (2*NttRange_NSigma*sqrt(Ntt_Init[i])/NttRange_Nbins)*l;
				Nw_[i]    = (2*NwRange_NSigma*sqrt(Nw_Init[i])/NwRange_Nbins)*m;
				Iteration(It,Print,i);
				EstimationChiSq_Ntt_vs_Nw[i]->SetBinContent(l,m,EstimationChiSq(i));
			}
		}
		cout<<"----- for N = "<<Njets_[i]<<" jets :"<<endl;
		cout<<"----- MinChiSq = "<<ChiSq[i]<<endl;
		cout<<"------- eb    = "<<eb_Init[i]<<" / "<<eb_MinChiSq[i]<<endl;
		cout<<"------- eudsc = "<<eudsc_Init[i]<<" / "<<eudsc_MinChiSq[i]<<endl;
		cout<<"------- Ntt   = "<<Ntt_Init[i]<<" / "<<Ntt_MinChiSq[i]<<endl;
		cout<<"------- Nw    = "<<Nw_Init[i]<<" / "<<Nw_MinChiSq[i]<<endl;
		cout<<"--------------------"<<endl;
		// Put back intial values and compute the estimates
		eb_[i]    = eb_Init[i];
		eudsc_[i] = eudsc_Init[i];
		Ntt_[i]   = Ntt_Init[i];
		Nw_[i]    = Nw_Init[i];
		Iteration(It,Print,i);		
	}
}

void WJetEstimation::CheckEstimationLinearity(const int &nbOfDatasets, int nbOfRescaleFact, double **RescaleFact, int TTlike_idx, int Wlike_idx){
	
	if(nbOfDatasets != NbOfDatasets){cout<<"*****************Linearity check*****************"<<endl;
					cout<<"The number of rescaled factors is not equal to the number of datasets used in the analysis..."<<endl;
					cout<<"*************************************************************************"<<endl;
					return;
					}

	RescaledTTlikeEstimation = new TGraphErrors(nbOfRescaleFact);
	RescaledWlikeEstimation  = new TGraphErrors(nbOfRescaleFact);
	tCanva_RescaledTTlikeEstimation = new TCanvas("tCanva_RescaledTTlikeEstimation","",-1);
	tCanva_RescaledWlikeEstimation  = new TCanvas("tCanva_RescaledWlikeEstimation", "",-1);

	double **n0_Init = new double*[NofJetsBins];
	double **n0      = new double*[NofJetsBins];
	double **n1_Init = new double*[NofJetsBins];
	double **n1      = new double*[NofJetsBins];
	double **n2_Init = new double*[NofJetsBins];
	double **n2      = new double*[NofJetsBins];
	double **n3_Init = new double*[NofJetsBins];
	double **n3      = new double*[NofJetsBins];

	// Make a copy of the original number of events
	for(int i=0;i<NofJetsBins;i++){

		n0_Init[i] = new double[nbOfDatasets];
		n0[i]      = new double[nbOfDatasets];
		n1_Init[i] = new double[nbOfDatasets];
		n1[i]      = new double[nbOfDatasets];
		n2_Init[i] = new double[nbOfDatasets];
		n2[i]      = new double[nbOfDatasets];
		n3_Init[i] = new double[nbOfDatasets];
		n3[i]      = new double[nbOfDatasets];

		for(int j=0;j<nbOfDatasets;j++) {
			n0_Init[i][j] = N0_[i][j];
			n1_Init[i][j] = N1_[i][j];
			n2_Init[i][j] = N2_[i][j];
			n3_Init[i][j] = N3_[i][j];
		}
	}
	
	cout<<"*****************Linearity check : results *****************"<<endl;
	cout<<" -- Number of rescaled factors :"<<nbOfRescaleFact<<endl;
	cout<<" -- List of rescaled factors for each dataset :"<<endl;
	for(int i=0;i<nbOfRescaleFact;i++){
		cout<<" --- ";
		for(int j=0;j<NofJetsBins;j++){
			for(int k=0;k<nbOfDatasets;k++){
			//Multiply the original numbers of events (per datasets) by a given factor "RescalFact"
			n0[j][k] = n0_Init[j][k]*RescaleFact[i][k];
			n1[j][k] = n1_Init[j][k]*RescaleFact[i][k];
			n2[j][k] = n2_Init[j][k]*RescaleFact[i][k];
			n3[j][k] = n3_Init[j][k]*RescaleFact[i][k];
			if(j==0)cout<<RescaleFact[i][k]<<"/";
			if(j==0 && k==nbOfDatasets-1) cout<<endl;
			}
		}
		//Estimate the number of tt-like and w-like events with the modified input numbers
		this->FillInputs(n0,n1,n2,n3);
		this->Estimation(false,false);

		RescaledTTlikeEstimation->SetPoint(i,100*RescaleFact[i][TTlike_idx],(double)this->GetNttEventsEstimated());
		RescaledTTlikeEstimation->SetPointError(i,0,RescaleFact[i][TTlike_idx]*((double)this->GetStatisticalTTLikeError()));
		RescaledWlikeEstimation->SetPoint(i,100*RescaleFact[i][TTlike_idx],(double)this->GetNWEventsEstimated());
		RescaledWlikeEstimation->SetPointError(i,0, RescaleFact[i][TTlike_idx]*((double)this->GetStatisticalWLikeError()));

		cout<<"-- Results :"<<endl;
		cout<<"--- Wlike estimate  = "<<this->GetNWEventsEstimated() <<"+/- "<<RescaleFact[i][TTlike_idx]*((double)this->GetStatisticalWLikeError())<<endl;
		cout<<"--- TTlike estimate = "<<this->GetNttEventsEstimated()<<"+/- "<<RescaleFact[i][TTlike_idx]*((double)this->GetStatisticalTTLikeError())<<endl;
		cout<<"------------"<<endl;
	}

	RescaledTTlikeEstimation ->SetNameTitle("RescaledTTlikeEstimation","");
	RescaledWlikeEstimation  ->SetNameTitle("RescaledWlikeEstimation","");

	tCanva_RescaledTTlikeEstimation->cd();
	RescaledTTlikeEstimation->Draw("AC*");
	RescaledTTlikeEstimation->GetXaxis()->SetTitle("Rescaling factor (%)");
	RescaledTTlikeEstimation->GetYaxis()->SetTitle("Estimated nb of TT-like events");

	tCanva_RescaledWlikeEstimation->cd();
	RescaledWlikeEstimation->Draw("AC*");
	RescaledWlikeEstimation->GetXaxis()->SetTitle("Rescaling factor (%)");
	RescaledWlikeEstimation->GetYaxis()->SetTitle("Estimated nb of W-like events");

	// Put the number of events back to their original values
	this->FillInputs(n0_Init,n1_Init,n2_Init,n3_Init);
	this->Estimation(false,false);

	for(int i=0;i<NofJetsBins;i++){
		delete n0_Init[i];
		delete n0[i];
		delete n1_Init[i];
		delete n1[i];
		delete n2_Init[i];
		delete n2[i];
		delete n3_Init[i];
		delete n3[i];
	}
	delete[] n0_Init;
	delete[] n0;
	delete[] n1_Init;
	delete[] n1;
	delete[] n2_Init;
	delete[] n2;
	delete[] n3_Init;
	delete[] n3;
}

void WJetEstimation::PrintResults(){
	cout<<endl;
	cout<<"***********************************************"<<endl;
	cout<<"*****      WJets Estimation Results         ***"<<endl;
	cout<<"***********************************************"<<endl;
				
	for(int i=0;i<NofJetsBins;i++) PrintResults(i);

	cout<<"***********************************************************"<<endl;
	cout<<"************        ALL INCLUSIF             **************"<<endl; 
	cout<<"Estimated eudsc    : "<<GetEffudscEstimated()<<endl;
	cout<<"Estimated eb       : "<<GetEffbEstimated()<<endl;
	double NttLike=0, NttLike_0=0, NttLike_1=0, NttLike_2=0, NttLike_3=0;
	for(int Njets=0;Njets<NofJetsBins;Njets++){
		for(unsigned int i=0;i<iDatasetsTTLike.size();i++){
			NttLike+= N0_[Njets][iDatasetsTTLike[i]];  NttLike+= N1_[Njets][iDatasetsTTLike[i]];  NttLike+= N2_[Njets][iDatasetsTTLike[i]]; NttLike+= N3_[Njets][iDatasetsTTLike[i]]; 
			NttLike_0+= N0_[Njets][iDatasetsTTLike[i]];  NttLike_1+= N1_[Njets][iDatasetsTTLike[i]];  NttLike_2+= N2_[Njets][iDatasetsTTLike[i]]; NttLike_3+= N3_[Njets][iDatasetsTTLike[i]]; 			}
		}
	double NWLike=0, NWLike_0=0, NWLike_1=0, NWLike_2=0, NWLike_3=0;
	for(int Njets=0;Njets<NofJetsBins;Njets++){
		for(unsigned int i=0;i<iDatasetsWLike.size();i++){
			NWLike+= N0_[Njets][iDatasetsWLike[i]]; NWLike+= N1_[Njets][iDatasetsWLike[i]]; NWLike+= N2_[Njets][iDatasetsWLike[i]]; NWLike+= N3_[Njets][iDatasetsWLike[i]]; 
			NWLike_0+= N0_[Njets][iDatasetsWLike[i]]; NWLike_1+= N1_[Njets][iDatasetsWLike[i]]; NWLike_2+= N2_[Njets][iDatasetsWLike[i]]; NWLike_3+= N3_[Njets][iDatasetsWLike[i]]; 
		}
	}
	cout<<"Estimated Ntt-like : "<<GetNttEventsEstimated()<<" / Predicted value from MC : "<< NttLike <<endl;
	double Ntt_0bjet=0, Ntt_1bjet=0, Ntt_2bjets=0, Ntt_3bjets=0 ;
	for(int Njets=0;Njets<NofJetsBins;Njets++){
		Ntt_0bjet+= WJetEstimation::Ntt_0bjet( Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		Ntt_1bjet+= WJetEstimation::Ntt_1bjet( Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		Ntt_2bjets+= WJetEstimation::Ntt_2bjets( Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		Ntt_3bjets+= WJetEstimation::Ntt_3bjets( Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
	}
	cout<<"- 0 b-jet          : "<<Ntt_0bjet<<" / "<< NttLike_0 <<endl;
	cout<<"- 1 b-jet          : "<<Ntt_1bjet<<" / "<< NttLike_1 <<endl;
	cout<<"- 2 b-jets         : "<<Ntt_2bjets<<" / "<< NttLike_2 <<endl;
	cout<<"- 3 b-jets         : "<<Ntt_3bjets<<" / "<< NttLike_3 <<endl;
	cout<<"Estimated Nw-like  : "<<GetNWEventsEstimated()<<" / Predicted value from MC : "<< NWLike <<endl;
	double Nw_0bjet=0, Nw_1bjet=0, Nw_2bjets=0, Nw_3bjets=0;
	for(int Njets=0;Njets<NofJetsBins;Njets++){
		Nw_0bjet+= WJetEstimation::Nw_0bjet( Nw_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		Nw_1bjet+= WJetEstimation::Nw_1bjet( Nw_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		Nw_2bjets+= WJetEstimation::Nw_2bjets( Nw_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		Nw_3bjets+= WJetEstimation::Nw_3bjets( Nw_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
	}
	cout<<"- 0 b-jet          : "<<Nw_0bjet<<" / "<< NWLike_0 <<endl;
	cout<<"- 1 b-jet          : "<<Nw_1bjet<<" / "<< NWLike_1 <<endl;
	cout<<"- 2 b-jets         : "<<Nw_2bjets<<" / "<< NWLike_2 <<endl;
	cout<<"- 3 b-jets         : "<<Nw_3bjets<<" / "<< NWLike_3 <<endl;
	cout<<"***********************************************************"<<endl<<endl;;
}

void WJetEstimation::PrintResults(int Njets){

		cout<<"***********************************************************"<<endl;
		cout<<"For N = "<<Njets_[Njets]<<" jets"<<endl;
		cout<<"Estimated eudsc    : "<<eudsc_[Njets]<<endl;
		cout<<"Estimated eb       : "<<eb_[Njets]<<endl;
		double NttLike = 0;
		double NttLike_0 = 0;
		double NttLike_1 = 0;
		double NttLike_2 = 0;
		double NttLike_3 = 0;
		for(unsigned int i=0;i<iDatasetsTTLike.size();i++){
			NttLike+= N0_[Njets][iDatasetsTTLike[i]];  NttLike+= N1_[Njets][iDatasetsTTLike[i]];  NttLike+= N2_[Njets][iDatasetsTTLike[i]]; NttLike+= N3_[Njets][iDatasetsTTLike[i]]; 
			NttLike_0+= N0_[Njets][iDatasetsTTLike[i]];  NttLike_1+= N1_[Njets][iDatasetsTTLike[i]];  NttLike_2+= N2_[Njets][iDatasetsTTLike[i]]; NttLike_3+= N3_[Njets][iDatasetsTTLike[i]]; 
			
		}
		double NWLike = 0;
		double NWLike_0 = 0;
		double NWLike_1 = 0;
		double NWLike_2 = 0;
		double NWLike_3 = 0;
		for(unsigned int i=0;i<iDatasetsWLike.size();i++){
			NWLike+= N0_[Njets][iDatasetsWLike[i]]; NWLike+= N1_[Njets][iDatasetsWLike[i]]; NWLike+= N2_[Njets][iDatasetsWLike[i]]; NWLike+= N3_[Njets][iDatasetsWLike[i]]; 
			NWLike_0+= N0_[Njets][iDatasetsWLike[i]]; NWLike_1+= N1_[Njets][iDatasetsWLike[i]]; NWLike_2+= N2_[Njets][iDatasetsWLike[i]]; NWLike_3+= N3_[Njets][iDatasetsWLike[i]]; 
		}
		cout<<"Estimated Ntt-like : "<<Ntt_[Njets]<<" / Predicted value from MC : "<< NttLike <<endl;
		cout<<"- 0 b-jet          : "<<Ntt_0bjet( Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])<<" / "<< NttLike_0 <<endl;
		cout<<"- 1 b-jet          : "<<Ntt_1bjet( Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])<<" / "<< NttLike_1 <<endl;
		cout<<"- 2 b-jets         : "<<Ntt_2bjets(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])<<" / "<< NttLike_2 <<endl;
		cout<<"- 3 b-jets         : "<<Ntt_3bjets(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])<<" / "<< NttLike_3 <<endl;
		cout<<"Estimated Nw-like  : "<<Nw_[Njets]<<" / Predicted value from MC : "<< NWLike <<endl;
		cout<<"- 0 b-jet          : "<<Nw_0bjet( Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])<<" / "<< NWLike_0 <<endl;
		cout<<"- 1 b-jet          : "<<Nw_1bjet( Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])<<" / "<< NWLike_1 <<endl;
		cout<<"- 2 b-jets         : "<<Nw_2bjets(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])<<" / "<< NWLike_2 <<endl;
		cout<<"- 3 b-jets         : "<<Nw_3bjets(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])<<" / "<< NWLike_3 <<endl;
		cout<<"***********************************************************"<<endl;
}

void WJetEstimation::PrintResults_Latex(){

	double NttLike[4]   = {0,0,0};
	double NttLike_0[3] = {0,0,0};
	double NttLike_1[3] = {0,0,0};
	double NttLike_2[3] = {0,0,0};
	double NttLike_3[3] = {0,0,0};
	for(unsigned int i=0;i<iDatasetsTTLike.size();i++){
		for(int j=0; j<3; j++){
			NttLike[0]  += N0_[j][iDatasetsTTLike[i]]; NttLike[1]  += N1_[j][iDatasetsTTLike[i]]; NttLike[2]  += N2_[j][iDatasetsTTLike[i]]; NttLike[3]  += N3_[j][iDatasetsTTLike[i]]; 
			NttLike_0[j]+= N0_[j][iDatasetsTTLike[i]]; NttLike_1[j]+= N1_[j][iDatasetsTTLike[i]]; NttLike_2[j]+= N2_[j][iDatasetsTTLike[i]]; NttLike_3[j]+= N3_[j][iDatasetsTTLike[i]]; 
		}
	}
		
	double NWLike[4]   = {0,0,0};
	double NWLike_0[3] = {0,0,0};
	double NWLike_1[3] = {0,0,0};
	double NWLike_2[3] = {0,0,0};
	double NWLike_3[3] = {0,0,0};
	for(unsigned int i=0;i<iDatasetsWLike.size();i++){
		for(int j=0; j<3; j++){
			NWLike[0]  += N0_[j][iDatasetsWLike[i]]; NWLike[1]  += N1_[j][iDatasetsWLike[i]]; NWLike[2]  += N2_[j][iDatasetsWLike[i]]; NWLike[3]  += N3_[j][iDatasetsWLike[i]]; 
			NWLike_0[j]+= N0_[j][iDatasetsWLike[i]]; NWLike_1[j]+= N1_[j][iDatasetsWLike[i]]; NWLike_2[j]+= N2_[j][iDatasetsWLike[i]]; NWLike_3[j]+= N3_[j][iDatasetsWLike[i]];
		}        
	}
		
	double Ntt_0bjet=0, Ntt_1bjet=0, Ntt_2bjets=0, Ntt_3bjets=0 ;
	for(int Njets=0;Njets<NofJetsBins;Njets++){
		Ntt_0bjet += WJetEstimation::Ntt_0bjet( Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		Ntt_1bjet += WJetEstimation::Ntt_1bjet( Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		Ntt_2bjets+= WJetEstimation::Ntt_2bjets(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		Ntt_3bjets+= WJetEstimation::Ntt_3bjets(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
	}
	double Nw_0bjet=0, Nw_1bjet=0, Nw_2bjets=0, Nw_3bjets=0;
	for(int Njets=0;Njets<NofJetsBins;Njets++){
		Nw_0bjet += WJetEstimation::Nw_0bjet( Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		Nw_1bjet += WJetEstimation::Nw_1bjet( Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		Nw_2bjets+= WJetEstimation::Nw_2bjets(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		Nw_3bjets+= WJetEstimation::Nw_3bjets(Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
	}

        cout<<"\\begin{table}"<<endl;
        cout<<"	\\centering"<<endl;
        cout<<"		\\begin{tabular}{l|cccc}"<<endl;
        cout<<"			\\hline"<<endl;
        cout<<" \\multicolumn{5}{l}{Estimation / Monte Carlo prediction : }\\\\"<<endl;
        cout<<"			             & $4$ jets & $5$ jets & $\\geq 6$ jets & Inclusive \\\\ "<<endl;
        cout<<"\\hline"<<endl;
        cout<<" Estimated $\\epsilon_b$      &"<<setprecision(3)<<eb_[0]<<"&"<<eb_[1]<<"&"<<eb_[2]<<"&"<<GetEffbEstimated()<<" \\\\"<<endl;
        cout<<"	Estimated $\\epsilon_{eudsc}$&"<<setprecision(4)<<eudsc_[0]<<"&"<<eudsc_[1]<<"&"<<eudsc_[2]<<"&"<<GetEffudscEstimated()<<" \\\\"<<endl;
	cout.precision(1);        
	cout<<"\\hline"<<endl;
        cout<<" \\multicolumn{5}{l}{$t\\bar{t}+jets$ and single top : }\\\\"<<endl;
        cout<<"\\hline"<<endl;
        cout<<" $0$ b-jet                    &"<<WJetEstimation::Ntt_0bjet( Ntt_[0], Nw_[0], eb_[0], eudsc_[0], Njets_[0])<<" / "<< NttLike_0[0] <<"&"<<WJetEstimation::Ntt_0bjet( Ntt_[1], Nw_[1], eb_[1], eudsc_[1], Njets_[1])<<" / "<< NttLike_0[1] <<"&"<<WJetEstimation::Ntt_0bjet( Ntt_[2], Nw_[2], eb_[2], eudsc_[2], Njets_[2])<<" / "<< NttLike_0[2] <<"&"<<Ntt_0bjet<<" / "<< NttLike[0]<<" \\\\"<<endl;
        cout<<" $1$ b-jet                    &"<<WJetEstimation::Ntt_1bjet( Ntt_[0], Nw_[0], eb_[0], eudsc_[0], Njets_[0])<<" / "<< NttLike_1[0] <<"&"<<WJetEstimation::Ntt_1bjet( Ntt_[1], Nw_[1], eb_[1], eudsc_[1], Njets_[1])<<" / "<< NttLike_1[1] <<"&"<<WJetEstimation::Ntt_1bjet( Ntt_[2], Nw_[2], eb_[2], eudsc_[2], Njets_[2])<<" / "<< NttLike_1[2] <<"&"<<Ntt_1bjet<<" / "<< NttLike[1]<<" \\\\"<<endl;
        cout<<" $2$ b-jets                   &"<<WJetEstimation::Ntt_2bjets(Ntt_[0], Nw_[0], eb_[0], eudsc_[0], Njets_[0])<<" / "<< NttLike_2[0] <<"&"<<WJetEstimation::Ntt_2bjets(Ntt_[1], Nw_[1], eb_[1], eudsc_[1], Njets_[1])<<" / "<< NttLike_2[1] <<"&"<<WJetEstimation::Ntt_2bjets(Ntt_[2], Nw_[2], eb_[2], eudsc_[2], Njets_[2])<<" / "<< NttLike_2[2] <<"&"<<Ntt_2bjets<<" / "<<NttLike[2]<<" \\\\"<<endl;
        cout<<" $3$ b-jets                   &"<<WJetEstimation::Ntt_3bjets(Ntt_[0], Nw_[0], eb_[0], eudsc_[0], Njets_[0])<<" / "<< NttLike_3[0] <<"&"<<WJetEstimation::Ntt_3bjets(Ntt_[1], Nw_[1], eb_[1], eudsc_[1], Njets_[1])<<" / "<< NttLike_3[1] <<"&"<<WJetEstimation::Ntt_3bjets(Ntt_[2], Nw_[2], eb_[2], eudsc_[2], Njets_[2])<<" / "<< NttLike_3[2] <<"&"<<Ntt_3bjets<<" / "<<NttLike[3]<<" \\\\"<<endl;
        cout<<"\\hline"<<endl;
        cout<<" Inclusive                    &"<<GetNttEventsEstimated(0)<<" / "<<GetNTTlikeMC(0)<<"&"<<GetNttEventsEstimated(1)<<" / "<<GetNTTlikeMC(1)<<"&"<<GetNttEventsEstimated(2)<<" / "<<GetNTTlikeMC(2)<<"&"<<GetNttEventsEstimated()<<" / "<<GetNTTlikeMC()<<" \\\\"<<endl;
        cout<<"\\hline"<<endl;
        cout<<"\\hline"<<endl;
        cout<<" \\multicolumn{5}{l}{Multijet and $W/Z+jets$ : }\\\\"<<endl;
        cout<<"\\hline"<<endl;
        cout<<" $0$ b-jet                    &"<<WJetEstimation::Nw_0bjet( Ntt_[0], Nw_[0], eb_[0], eudsc_[0], Njets_[0])<<" / "<< NWLike_0[0] <<"&"<<WJetEstimation::Nw_0bjet( Ntt_[1], Nw_[1], eb_[1], eudsc_[1], Njets_[1])<<" / "<< NWLike_0[1] <<"&"<<WJetEstimation::Nw_0bjet( Ntt_[2], Nw_[2], eb_[2], eudsc_[2], Njets_[2])<<" / "<< NWLike_0[2] <<"&"<<Nw_0bjet<<" / "<< NWLike[0]<<" \\\\"<<endl;
        cout<<" $1$ b-jet                    &"<<WJetEstimation::Nw_1bjet( Ntt_[0], Nw_[0], eb_[0], eudsc_[0], Njets_[0])<<" / "<< NWLike_1[0] <<"&"<<WJetEstimation::Nw_1bjet( Ntt_[1], Nw_[1], eb_[1], eudsc_[1], Njets_[1])<<" / "<< NWLike_1[1] <<"&"<<WJetEstimation::Nw_1bjet( Ntt_[2], Nw_[2], eb_[2], eudsc_[2], Njets_[2])<<" / "<< NWLike_1[2] <<"&"<<Nw_1bjet<<" / "<< NWLike[1]<<" \\\\"<<endl;
        cout<<" $2$ b-jets                   &"<<WJetEstimation::Nw_2bjets(Ntt_[0], Nw_[0], eb_[0], eudsc_[0], Njets_[0])<<" / "<< NWLike_2[0] <<"&"<<WJetEstimation::Nw_2bjets(Ntt_[1], Nw_[1], eb_[1], eudsc_[1], Njets_[1])<<" / "<< NWLike_2[1] <<"&"<<WJetEstimation::Nw_2bjets(Ntt_[2], Nw_[2], eb_[2], eudsc_[2], Njets_[2])<<" / "<< NWLike_2[2] <<"&"<<Nw_2bjets<<" / "<<NWLike[2]<<" \\\\"<<endl;
        cout<<" $3$ b-jets                   &"<<WJetEstimation::Nw_3bjets(Ntt_[0], Nw_[0], eb_[0], eudsc_[0], Njets_[0])<<" / "<< NWLike_3[0] <<"&"<<WJetEstimation::Nw_3bjets(Ntt_[1], Nw_[1], eb_[1], eudsc_[1], Njets_[1])<<" / "<< NWLike_3[1] <<"&"<<WJetEstimation::Nw_3bjets(Ntt_[2], Nw_[2], eb_[2], eudsc_[2], Njets_[2])<<" / "<< NWLike_3[2] <<"&"<<Nw_3bjets<<" / "<<NWLike[3]<<" \\\\"<<endl;
        cout<<"\\hline"<<endl;
        cout<<" Inclusive                    &"<<GetNWEventsEstimated(0)<<" / "<<GetNWlikeMC(0)<<"&"<<GetNWEventsEstimated(1)<<" / "<<GetNWlikeMC(1)<<"&"<<GetNWEventsEstimated(2)<<" / "<<GetNWlikeMC(2)<<"&"<<GetNWEventsEstimated()<<" / "<<GetNWlikeMC()<<" \\\\"<<endl;
        cout<<"\\hline"<<endl;
        cout<<"\\hline"<<endl;
        cout<<"		\\end{tabular}"<<endl;
        cout<<"	\\caption{}"<<endl;
        cout<<"	\\label{tab:}"<<endl;
        cout<<"\\end{table}"<<endl;
	cout<<resetiosflags (ios::left);

}

double WJetEstimation::GetEffudscEstimated() const {double a = 0; for(int i=0;i<NofJetsBins;i++) a+= GetEffudscEstimated(i)*NofEvents(i); return a/NofEvents();}

double WJetEstimation::GetEffudscEstimated(int Njets) const { return (Njets<NofJetsBins? eudsc_[Njets]: -9999.);}

double WJetEstimation::GetEffbEstimated() const {double a = 0; for(int i=0;i<NofJetsBins;i++) a+= GetEffbEstimated(i)*NofEvents(i); return a/NofEvents();}

double WJetEstimation::GetEffbEstimated(int Njets) const { return (Njets<NofJetsBins? eb_[Njets]: -9999.);}

double WJetEstimation::GetNWEventsEstimated() const {double a = 0; for(int i=0;i<NofJetsBins;i++) a+= GetNWEventsEstimated(i); return a;}

double WJetEstimation::GetNWEventsEstimated(int Njets) const { return (Njets<NofJetsBins? Nw_[Njets]: -9999.);}

double WJetEstimation::GetNttEventsEstimated() const {double a = 0; for(int i=0;i<NofJetsBins;i++) a+= GetNttEventsEstimated(i); return a;}

double WJetEstimation::GetNttEventsEstimated(int Njets) const { return (Njets<NofJetsBins? Ntt_[Njets]: -9999.);} 

double WJetEstimation::GetNtotalEventsEstimated() const {double a = 0; for(int i=0;i<NofJetsBins;i++) a+= GetNtotalEventsEstimated(i); return a;}

double WJetEstimation::GetNtotalEventsEstimated(int Njets) const { return (Njets<NofJetsBins? GetNttEventsEstimated(Njets)+GetNWEventsEstimated(Njets) : -9999.);} 

double WJetEstimation::GetNWEventsRatioEstimatedNbjets(int nbjets) const{
	if(GetNWEventsEstimatedNbjets(nbjets)==-9999. || NofEventsNbjets(nbjets)==-9999. || NofEventsNbjets(nbjets)==0.) return -9999.;
	return GetNWEventsEstimatedNbjets(nbjets)/NofEventsNbjets(nbjets);
}

double WJetEstimation::GetNWEventsEstimatedNbjets(int nbjets) const{
	double a = 0;
	for(int i=0;i<NofJetsBins;i++) a+=GetNWEventsEstimatedNbjets(nbjets,i);
	return a;
}

double WJetEstimation::GetNWEventsEstimatedNbjets(int nbjets, int njets) const{
	if(nbjets<0)  return -9999.;
	if(nbjets==0) return Nw_0bjet( Ntt_[njets],Nw_[njets],eb_[njets],eudsc_[njets],Njets_[njets]);
	if(nbjets==1) return Nw_1bjet( Ntt_[njets],Nw_[njets],eb_[njets],eudsc_[njets],Njets_[njets]);
	if(nbjets==2) return Nw_2bjets(Ntt_[njets],Nw_[njets],eb_[njets],eudsc_[njets],Njets_[njets]);
                      return Nw_3bjets(Ntt_[njets],Nw_[njets],eb_[njets],eudsc_[njets],Njets_[njets]);
}

double WJetEstimation::GetNttEventsRatioEstimatedNbjets(int nbjets) const{
	if(GetNttEventsEstimatedNbjets(nbjets)==-9999. || NofEventsNbjets(nbjets)==-9999. || NofEventsNbjets(nbjets)==0.) return -9999.;
	return GetNttEventsEstimatedNbjets(nbjets)/NofEventsNbjets(nbjets);
}
			     
double WJetEstimation::GetNttEventsEstimatedNbjets(int nbjets) const{
	double a = 0;
	for(int i=0;i<NofJetsBins;i++) a+=GetNttEventsEstimatedNbjets(nbjets,i);
	return a;
}

double WJetEstimation::GetNttEventsEstimatedNbjets(int nbjets, int njets) const{
	if(nbjets<0)  return -9999.;
	if(nbjets==0) return Ntt_0bjet( Ntt_[njets],Nw_[njets],eb_[njets],eudsc_[njets],Njets_[njets]);
	if(nbjets==1) return Ntt_1bjet( Ntt_[njets],Nw_[njets],eb_[njets],eudsc_[njets],Njets_[njets]);
	if(nbjets==2) return Ntt_2bjets(Ntt_[njets],Nw_[njets],eb_[njets],eudsc_[njets],Njets_[njets]);
                      return Ntt_3bjets(Ntt_[njets],Nw_[njets],eb_[njets],eudsc_[njets],Njets_[njets]);
	
}
double WJetEstimation::GetNTTlikeMC() const{
	double NttLike = 0;
	for(unsigned int i=0;i<iDatasetsTTLike.size();i++){
		for(int j=0; j<3; j++){
			NttLike  += N0_[j][iDatasetsTTLike[i]];
			NttLike  += N1_[j][iDatasetsTTLike[i]];
			NttLike  += N2_[j][iDatasetsTTLike[i]];
			NttLike  += N3_[j][iDatasetsTTLike[i]]; 
		}
	}
	return NttLike;
}
double WJetEstimation::GetNTTlikeMC(int njets) const{
	double NttLike = 0;
	for(unsigned int i=0;i<iDatasetsTTLike.size();i++){
		NttLike  += N0_[njets][iDatasetsTTLike[i]];
		NttLike  += N1_[njets][iDatasetsTTLike[i]];
		NttLike  += N2_[njets][iDatasetsTTLike[i]];
		NttLike  += N3_[njets][iDatasetsTTLike[i]]; 
	}
	return NttLike;
}
double WJetEstimation::GetNWlikeMC() const{
	double NwLike = 0;
	for(unsigned int i=0;i<iDatasetsWLike.size();i++){
		for(int j=0; j<3; j++){
			NwLike  += N0_[j][iDatasetsWLike[i]];
			NwLike  += N1_[j][iDatasetsWLike[i]];
			NwLike  += N2_[j][iDatasetsWLike[i]];
			NwLike  += N3_[j][iDatasetsWLike[i]]; 
		}
	}
	return NwLike;
}
double WJetEstimation::GetNWlikeMC(int njets) const{
	double NwLike = 0;
	for(unsigned int i=0;i<iDatasetsWLike.size();i++){
		NwLike  += N0_[njets][iDatasetsWLike[i]];
		NwLike  += N1_[njets][iDatasetsWLike[i]];
		NwLike  += N2_[njets][iDatasetsWLike[i]];
		NwLike  += N3_[njets][iDatasetsWLike[i]]; 
	}
	return NwLike;
}

double WJetEstimation::GetStatisticalWLikeError()  const {
	double error = 0;
	for(int i=0;i<NofJetsBins;i++) error += pow(statisticalWLikeError_[i][4],2);
	return sqrt(error);
}

double WJetEstimation::GetSystematicWLikeError()   const {
	double error = 0;
	for(int i=0;i<NofJetsBins;i++) error += pow(systematicWLikeError_[i][4],2);
	return sqrt(error);
}

double WJetEstimation::GetStatisticalTTLikeError() const {
	double error = 0;
	for(int i=0;i<NofJetsBins;i++) error += pow(statisticalTTLikeError_[i][4],2);
	return sqrt(error);
}

double WJetEstimation::GetSystematicTTLikeError()  const {
	double error = 0;
	for(int i=0;i<NofJetsBins;i++) error += pow(systematicTTLikeError_[i][4],2);
	return sqrt(error);
}

void WJetEstimation::Initialisation(bool &Print){
	for(int i=0;i<NofJetsBins;i++) Initialisation(Print, i);
}

void WJetEstimation::Initialisation(bool &Print, int Njets)
{	
	if(Print) cout<<"For N = "<<Njets_[Njets]<<" jets"<<endl;
	if(Print) cout<<"Total nb of events (0/1/2/3 b-jets) = "<<N0bjets_[Njets]+N1bjets_[Njets]+N2bjets_[Njets]+N3bjets_[Njets]<<"("<<N0bjets_[Njets]<<"/"<<N1bjets_[Njets]<<"/"<<N2bjets_[Njets]<<"/"<<N3bjets_[Njets]<<")"<<endl;
	double R21 = (N1bjets_[Njets] !=0 ? N2bjets_[Njets]/N1bjets_[Njets] : 0);
	eb_[Njets]  = (2*R21)/(1+2*R21);
	if(Print) cout<<"Estimated b-tag efficiency : "<<eb_[Njets]<<endl;
	
	if(eb_[Njets] > 0) Ntt_[Njets] = N2bjets_[Njets]/(eb_[Njets]*eb_[Njets]);
	else               Ntt_[Njets] = 0;
	if(Print) cout<<"Estimated nb of tt+jets event : "<<Ntt_[Njets]<<endl;

	if(N0bjets_[Njets]+N1bjets_[Njets]+N2bjets_[Njets]+N3bjets_[Njets]>=Ntt_[Njets]) Nw_[Njets] = N0bjets_[Njets]+N1bjets_[Njets]+N2bjets_[Njets]+N3bjets_[Njets]-Ntt_[Njets];
	else Nw_[Njets] = 0;
	if(Print) cout<<"Estimated nb of w+jets event : "<<Nw_[Njets]<<endl;
	
	eudsc_[Njets] = 0;
	double Delta = 0, r1 =0, r2 = 0;
	Delta = pow(eb_[Njets],4)+4*(eb_[Njets]*(1-eb_[Njets])*N3bjets_[Njets]/(2*Ntt_[Njets]));
	if(Print) cout<<"Delta = "<<Delta<<endl;
	r1 = (-eb_[Njets]*eb_[Njets]-sqrt(Delta)) /(2*eb_[Njets]*(1-eb_[Njets]));
	r2 = (-eb_[Njets]*eb_[Njets]+sqrt(Delta))/(2*eb_[Njets]*(1-eb_[Njets]));
	if(Print) cout<<"R1 = "<<r1<<endl;
	if(Print) cout<<"R2 = "<<r2<<endl;
	
	//cout<<"edusc_[Njets] "<<eudsc_[Njets]<<endl;
	if(r1<0) eudsc_[Njets] = r2;
	else if(r2<0) eudsc_[Njets] = r1;
	else (r1<r2 ? eudsc_[Njets] = r1 : eudsc_[Njets] = r2) ;
	
	//double ChiSq = GetGlobalEffudscEstimate(N3bjets_[Njets], N2bjets_[Njets], N1bjets_[Njets],N0bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets],Print);

	if(Print){
		cout<<"Estimated mis-tag efficiency : "<<eudsc_[Njets]<<endl;//"/ associated ChiSq ="<<ChiSq<<endl;
		Print4Estimators(Njets);
	}
}

void WJetEstimation::Iteration(int &It, bool &Print)
{
	for(int i=0;i<NofJetsBins;i++) Iteration(It, Print, i);
}

void WJetEstimation::Iteration(int &It, bool &Print, int Njets) 
{
	if(Njets>=NofJetsBins) return;
	double Chi2_eb = 0, Chi2_eudsc = 0;//, eb_old = 0, eudsc_old = 0, Ntt_old = 0, Nw_old = 0;
	bool verbose = false;
	bool UseGlobalMethod = false;
/*
	// First equation : eudsc from the "3 b-jets" bins
	histos[Njets][0]->Fill(1,fabs(eudsc_[Njets]-eudsc_fromN3bjets(N3bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])));
	eudsc_[Njets] = eudsc_fromN3bjets(N3bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
	//eudsc_[Njets] = GetGlobalEffudscEstimate(N3bjets_[Njets], N2bjets_[Njets], N1bjets_[Njets],N0bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], Njets_[Njets],true);

	// Second equation : eb from the "2 b-jets" bins
	histos[Njets][1]->Fill(1,fabs(eb_[Njets]-eb_fromN2bjets(N2bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])));
	eb_[Njets]    = eb_fromN2bjets(N2bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
	//eb_[Njets] = GetGlobalEffbEstimate(N3bjets_[Njets], N2bjets_[Njets], N1bjets_[Njets],N0bjets_[Njets], Ntt_[Njets], Nw_[Njets], eudsc_[Njets], Njets_[Njets],true);
*/

	// First equation : eb from the "2 b-jets" bins
	histos[Njets][0]->Fill(1,fabs(eb_[Njets]-eb_fromN2bjets(N2bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])));
	if(!UseGlobalMethod) eb_[Njets]    = eb_fromN2bjets(N2bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
	else Chi2_eb = GetGlobalEffbEstimate(N3bjets_[Njets], N2bjets_[Njets], N1bjets_[Njets],N0bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets],verbose);

	// Second equation : eudsc from the "3 b-jets" bins
	histos[Njets][1]->Fill(1,fabs(eudsc_[Njets]-eudsc_fromN3bjets(N3bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])));
	if(!UseGlobalMethod) eudsc_[Njets] = eudsc_fromN3bjets(N3bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
	else Chi2_eudsc = GetGlobalEffudscEstimate(N3bjets_[Njets], N2bjets_[Njets], N1bjets_[Njets],N0bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets],verbose);

	// Third equation : Nb of tt+jets events from the "1 b-jet" bins
	histos[Njets][2]->Fill(1,fabs(Ntt_[Njets]-Ntt_fromN1bjet(N1bjets_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])));
	Ntt_[Njets]   = Ntt_fromN1bjet(N1bjets_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
	
	// Fourth equation : Nb of W+jets events from the "0 b-jets" bins
	histos[Njets][3]->Fill(1,fabs(Nw_[Njets]-Nw_fromN0bjet(N0bjets_[Njets], Ntt_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])));
	Nw_[Njets]    = Nw_fromN0bjet(N0bjets_[Njets], Ntt_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
	
	if(Print){
		cout<<"0th iteration"<<endl;
		Print4Estimators(Njets);
		CheckConsistency(Njets);
	}

	for(int i = 0;i<It;i++)
	{

		if(Print) cout<<i+1<<"th iteration"<<endl;

		// The following equations are solved sequentially
		
		// First equation : eb from the "2 b-jets" bins
		histos[Njets][0]->Fill(i+2,fabs(eb_[Njets]-eb_fromN2bjets(N2bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])));
		//eb_old = eb_[Njets];
		if(!UseGlobalMethod) eb_[Njets] = eb_fromN2bjets(N2bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		else                    Chi2_eb = GetGlobalEffbEstimate(N3bjets_[Njets], N2bjets_[Njets], N1bjets_[Njets],N0bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets],verbose);
		//std::cout<<"eb_[Njets] = "<<(double)eb_[Njets]<<std::endl;
		
		// Second equation : eudsc from the "3 b-jets" bins
		histos[Njets][1]->Fill(i+2,fabs(eudsc_[Njets]-eudsc_fromN3bjets(N3bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])));
		//eudsc_old = eudsc_[Njets];
		if(!UseGlobalMethod) eudsc_[Njets] = eudsc_fromN3bjets(N3bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		else                    Chi2_eudsc = GetGlobalEffudscEstimate(N3bjets_[Njets], N2bjets_[Njets], N1bjets_[Njets],N0bjets_[Njets], Ntt_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets],verbose);
		//std::cout<<"eudsc_[Njets] = "<<(double)eudsc_[Njets]<<std::endl;
		
		// Third equation : Nb of tt+jets events from the "1 b-jet" bins
		histos[Njets][2]->Fill(i+2,fabs(Ntt_[Njets]-Ntt_fromN1bjet(N1bjets_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])));
		//Ntt_old = Ntt_[Njets];
		Ntt_[Njets] = Ntt_fromN1bjet(N1bjets_[Njets], Nw_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
		
		// Fourth equation : Nb of W+jets events from the "0 b-jets" bins
		histos[Njets][3]->Fill(i+2,fabs(Nw_[Njets]-Nw_fromN0bjet(N0bjets_[Njets], Ntt_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets])));
		//Nw_old = Nw_[Njets];
		Nw_[Njets] = Nw_fromN0bjet(N0bjets_[Njets], Ntt_[Njets], eb_[Njets], eudsc_[Njets], Njets_[Njets]);
/*
		if(Print) Print4Estimators(Njets);
		if(Print) CheckConsistency(Njets);
		
		if(fabs(eb_old   - eb_[Njets])   <0.00001) break;
		if(fabs(eudsc_old- eudsc_[Njets])<0.00001) break;
		if(fabs(Ntt_old- Ntt_[Njets])<0.1) break;
		if(fabs(Nw_old- Nw_[Njets])<0.1) break;
*/
	}
	if(Print) Print4Estimators(Njets);
	if(Print) CheckConsistency(Njets);
}

double WJetEstimation::N0bjet(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const
{
	double Ntt_0bjet = (1-eb)*(1-eb)*pow((1-eudsc),n-2)*Ntt;
	double Nw_0bjet  = pow((1-eudsc),n)*Nw;
	if(Ntt_0bjet<0) Ntt_0bjet = 0;
	if(Nw_0bjet<0)  Nw_0bjet  = 0;
	
	if(Ntt_0bjet+Nw_0bjet>=0) return (Ntt_0bjet+Nw_0bjet);
	else return 0;
}

double WJetEstimation::Ntt_0bjet(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const
{
	double Ntt_0bjet = (1-eb)*(1-eb)*pow((1-eudsc),n-2)*Ntt;
	if(Ntt_0bjet<0) Ntt_0bjet = 0;
	return (Ntt_0bjet) ;
}

double WJetEstimation::Nw_0bjet(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const
{
	double Nw_0bjet  = pow((1-eudsc),n)*Nw;	
	if(Nw_0bjet<0)  Nw_0bjet  = 0;
	return (Nw_0bjet) ;
}

double WJetEstimation::N1bjet(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const
{
	double Ntt_1bjet = (2*eb*(1-eb)*pow((1-eudsc),n-2)+(n-2)*pow((1-eb),2)*eudsc*pow((1-eudsc),n-3))*Ntt;
	double Nw_1bjet  = n*eudsc*pow((1-eudsc),n-1)*Nw;
	if(Ntt_1bjet<0) Ntt_1bjet = 0;
	if(Nw_1bjet<0)  Nw_1bjet  = 0;
	
	if(Ntt_1bjet+Nw_1bjet>=0) return (Ntt_1bjet+Nw_1bjet);
	else return 0;
}

double WJetEstimation::Ntt_1bjet(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const
{
	double Ntt_1bjet = (2*eb*(1-eb)*pow((1-eudsc),n-2)+(n-2)*pow((1-eb),2)*eudsc*pow((1-eudsc),n-3))*Ntt;
	if(Ntt_1bjet<0) Ntt_1bjet = 0;
	return (Ntt_1bjet) ;
}

double WJetEstimation::Nw_1bjet(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const
{
	double Nw_1bjet  = n*eudsc*pow((1-eudsc),n-1)*Nw;
	if(Nw_1bjet<0)  Nw_1bjet  = 0;
	return (Nw_1bjet) ;
}

double WJetEstimation::N2bjets(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const
{
	double Ntt_2bjets = (eb*eb*pow((1-eudsc),n-2)+2*eb*(1-eb)*(n-2)*eudsc*pow((1-eudsc),n-3)+pow((1-eb),2)*((n-2)*(n-3)/2)*eudsc*eudsc*pow((1-eudsc),n-4))*Ntt;
	double Nw_2bjets  = ((n*(n-1)/2)*eudsc*eudsc*pow((1-eudsc),n-2))*Nw;
	if(Ntt_2bjets<0) Ntt_2bjets = 0;
	if(Nw_2bjets<0)  Nw_2bjets  = 0;
	
	if(Ntt_2bjets+Nw_2bjets>=0) return (Ntt_2bjets+Nw_2bjets);
	else return 0;
}

double WJetEstimation::Ntt_2bjets(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const
{
	double Ntt_2bjets = (eb*eb*pow((1-eudsc),n-2)+2*eb*(1-eb)*(n-2)*eudsc*pow((1-eudsc),n-3)+pow((1-eb),2)*((n-2)*(n-3)/2)*eudsc*eudsc*pow((1-eudsc),n-4))*Ntt;
	if(Ntt_2bjets<0) Ntt_2bjets = 0;
	return (Ntt_2bjets) ;
}

double WJetEstimation::Nw_2bjets(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const
{
	double Nw_2bjets  = ((n*(n-1)/2)*eudsc*eudsc*pow((1-eudsc),n-2))*Nw;
	if(Nw_2bjets<0)  Nw_2bjets  = 0;
	return (Nw_2bjets) ;
}

double WJetEstimation::N3bjets(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const
{
        double Ntt_3bjets = (eb*eb*(n-2)*eudsc*pow((1-eudsc),n-3)+2*eb*(1-eb)*((n-2)*(n-3)/2)*pow(eudsc,2)*pow((1-eudsc),n-4)+(n>4 ? pow((1-eb),2)*((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow((1-eudsc),n-5) : 0 ))*Ntt;
	double Nw_3bjets  = ((n*(n-1)*(n-2)/6)*pow((eudsc),3)*pow((1-eudsc),n-3))*Nw;
	if(Ntt_3bjets<0) Ntt_3bjets = 0;
	if(Nw_3bjets<0)  Nw_3bjets  = 0;
	
	if(Ntt_3bjets+Nw_3bjets>=0) return (Ntt_3bjets+Nw_3bjets);
	else return 0;
			
}

double WJetEstimation::Ntt_3bjets(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const
{
	double Ntt_3bjets = (eb*eb*(n-2)*eudsc*pow((1-eudsc),n-3)+2*eb*(1-eb)*((n-2)*(n-3)/2)*pow(eudsc,2)*pow((1-eudsc),n-4)+(n>4 ? pow((1-eb),2)*((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow((1-eudsc),n-5) : 0 ))*Ntt;
	if(Ntt_3bjets<0) Ntt_3bjets = 0;
	return (Ntt_3bjets) ;
}

double WJetEstimation::Nw_3bjets(double &Ntt, double &Nw, double &eb, double &eudsc, int &n) const
{
	double Nw_3bjets  = ((n*(n-1)*(n-2)/6)*pow((eudsc),3)*pow((1-eudsc),n-3))*Nw;
	if(Nw_3bjets<0)  Nw_3bjets  = 0;
	return (Nw_3bjets) ;
}

double WJetEstimation::eb_fromN2bjets(double &N2bjets, double &Ntt, double &Nw, double &eb_old, double &eudsc, int &n) const
{
/*
        double alpha = 2*(n-2)*eudsc*pow((1-eudsc),n-3);
        double beta  = ((n-2)*(n-3)/2)*eudsc*eudsc*pow((1-eudsc),n-4);

        double a = pow((1-eudsc),n-2)-alpha+beta;
        double b = alpha-2*beta;
        double c = beta+(n*(n-1)/2)*eudsc*eudsc*pow((1-eudsc),n-2)*Nw/Ntt-N2bjets/Ntt;
*/
       	vector<double> roots;
	double a = (pow(1-eudsc,n-2) - 2*(n-2)*eudsc*pow(1-eudsc,n-3) + ((n-2)*(n-3)/2)*pow(eudsc,2)*pow(1-eudsc,n-4))*Ntt;
       	double b = (2*(n-2)*eudsc*pow(1-eudsc,n-3) - 2*((n-2)*(n-3)/2)*pow(eudsc,2)*pow(1-eudsc,n-4))*Ntt;
       	double c = (((n-2)*(n-3)/2)*pow(eudsc,2)*pow(1-eudsc,n-4))*Ntt + ((n)*(n-1)/2)*pow(eudsc,2)*pow(1-eudsc,n-2)*Nw - N2bjets;

        if(a!= 0 && b!=0){
		ROOT::Math::Polynomial poly(a,b,c);
        	roots = poly.FindRealRoots();
	}
	
        if(roots.size()<1) return eb_old;
        else if(roots.size() == 1)
        {
                if(roots[0]<0) return eb_old;
                else return roots[0];
        }
        else if(roots.size() == 2)
        {
                if(roots[1]<0)         return eb_old;
                if(roots[0]<0)
                {
                        if(roots[1]<0) return eb_old;
                        if(roots[1]<1) return roots[1];
                        else           return eb_old;
                }
                if(roots[0]>1)         return eb_old;
                else if(roots[1]>1)    return roots[0];
                else if(fabs(roots[0]-eb_old)<fabs(roots[1]-eb_old)) return roots[0];
                else return roots[1];
        }
        else return eb_old;
}

double WJetEstimation::eudsc_fromN3bjets(double &N3bjets, double &Ntt, double &Nw, double &eb, double &eudsc_old, int &n) const
{    
	if(N3bjets<= 0 || Ntt<=0 || Nw<=0 ) return eudsc_old;
	double eudsc =0;
	vector<double> roots;
	if(n==4)/////// Equation holds only for n = 4  ///////////////
        {
                double a4 = -4*Nw;
                double a3 =  4*Nw;
                double a2 =      (-4*eb*eb+2*eb)*Ntt;
                double a1 =      ( 2*eb*eb     )*Ntt;
                double a0 = -N3bjets;
                double par[5] = {a0,a1,a2,a3,a4};

                ROOT::Math::Polynomial poly(4);
                poly.SetParameters(par);
		roots = poly.FindRealRoots();
                
		for(unsigned int i = 0 ; i < roots.size(); i++)
                {
                        if(roots[i]<0 || roots[i]>1) continue;
                        if(fabs(eudsc_old-roots[i])<fabs(eudsc_old-eudsc)) eudsc = roots[i];
                }
                
		if(eudsc != 0) return eudsc;
                else return eudsc_old;
        }
        else if(n==5)/////// Equation holds only for n = 5  ///////////////
        {
                double a5 = 10*Nw;
                double a4 =-20*Nw;
                double a3 = 10*Nw+( 10*eb*eb-8*eb+1)*Ntt;
                double a2 =       (-12*eb*eb+6*eb  )*Ntt;
                double a1 =       (  3*eb*eb       )*Ntt;
                double a0 = -N3bjets;
                double par[6] = {a0,a1,a2,a3,a4,a5};

                ROOT::Math::Polynomial poly(5);
                poly.SetParameters(par);

                vector<double> roots = poly.FindRealRoots();
                for(unsigned int i = 0 ; i < roots.size(); i++)
                {
                        if(roots[i]<0 || roots[i]>1) continue;
                        if(fabs(eudsc_old-roots[i])<fabs(eudsc_old-eudsc)) eudsc = roots[i];
                }
                if(eudsc != 0) return eudsc;
                else return eudsc_old;
        }
        else /////// Equation holds only for n = 6 but is used for n>6 as well  ///////////////
        {
                double a6 =-20*Nw;
                double a5 = 60*Nw;
                double a4 =-60*Nw+(-20*eb*eb+20*eb-4)*Ntt;
                double a3 = 20*Nw+( 40*eb*eb-32*eb+4)*Ntt;
                double a2 =       (-24*eb*eb+12*eb  )*Ntt;
                double a1 =       (  4*eb*eb        )*Ntt;
                double a0 = -N3bjets;
                double par[7] = {a0,a1,a2,a3,a4,a5,a6};

                ROOT::Math::Polynomial poly(6);
                poly.SetParameters(par);

                vector<double> roots = poly.FindRealRoots();
                for(unsigned int i = 0 ; i < roots.size(); i++)
                {
                        if(roots[i]<0 || roots[i]>1) continue;
                        if(fabs(eudsc_old-roots[i])<fabs(eudsc_old-eudsc)) eudsc = roots[i];
                }
                if(eudsc != 0) return eudsc;
                else return eudsc_old;
        }
}

double WJetEstimation::Ntt_fromN2bjets(double &N2bjets, double &Nw, double &eb, double &eudsc, int &n) const
{
	double Nw_2bjets  = ((n*(n-1)/2)*eudsc*eudsc*pow((1-eudsc),n-2))*Nw;
	double Ntt = 0;
	if(eb!=0 && eudsc!=0) Ntt = (N2bjets-Nw_2bjets)/(eb*eb*pow((1-eudsc),n-2)+2*eb*(1-eb)*(n-2)*eudsc*pow((1-eudsc),n-3)+pow((1-eb),2)*((n-2)*(n-3)/2)*eudsc*eudsc*pow((1-eudsc),n-4));
	
	return (Ntt>=0 ? Ntt : 0);
}

double WJetEstimation::Ntt_fromN1bjet(double &N1bjet, double &Nw, double &eb, double &eudsc, int &n) const
{
	double Nw_1bjet  = (n*eudsc*pow((1-eudsc),n-1))*Nw;
	double Ntt = (N1bjet-Nw_1bjet)/(2*eb*(1-eb)*pow((1-eudsc),n-2)+(1-eb)*(1-eb)*(n-2)*eudsc*pow((1-eudsc),n-3));
	
	return (Ntt>=0 ? Ntt : 0);
}

double WJetEstimation::Nw_fromNtotal(double &Ntotal, double &Ntt) const
{
	return (Ntotal - Ntt);
}

double WJetEstimation::Nw_fromN0bjet(double &N0bjet, double &Ntt, double &eb, double &eudsc, int &n) const
{
	double Ntt_0bjet = (1-eb)*(1-eb)*pow((1-eudsc),n-2)*Ntt;
	double Nw = (N0bjet-Ntt_0bjet)/(pow((1-eudsc),n));
	
	return (Nw>=0 ? Nw : 0);
}


double WJetEstimation::NofEventsNbjets(int nbjets) const { 
	if(nbjets<0) return -9999.;
	double a = 0;
	if(nbjets==0) {for(int i=0;i<NofJetsBins;i++) a+=N0bjets_[i]; return a;}
	if(nbjets==1) {for(int i=0;i<NofJetsBins;i++) a+=N1bjets_[i]; return a;}
	if(nbjets==2) {for(int i=0;i<NofJetsBins;i++) a+=N2bjets_[i]; return a;}
	for(int i=0;i<NofJetsBins;i++) a+=N3bjets_[i]; return a;
}

double WJetEstimation::NofEvents(int njets, int nbjets) const {
	if(njets<NofJetsBins)  return -9999.;
	if(nbjets==0) return N0bjets_[njets];
	if(nbjets==1) return N1bjets_[njets];
	if(nbjets==2) return N2bjets_[njets];
	if(nbjets==3) return N3bjets_[njets];
	return -9999.;
};
												
vector<double> WJetEstimation::eudsc_cand_fromN3bjets(double &N3bjets, double &Ntt, double &Nw, double &eb, int &n) const
{    
	vector<double> roots;
	if(N3bjets <=0 || Ntt <=0 || Nw <=0 ) return roots;

	if(n==4)/////// Equation holds only for n = 4  ///////////////
        {
                double a4 = -4*Nw;
                double a3 =  4*Nw;
                double a2 =      (-4*eb*eb+2*eb)*Ntt;
                double a1 =      ( 2*eb*eb     )*Ntt;
                double a0 = -N3bjets;
                double par[5] = {a0,a1,a2,a3,a4};

                ROOT::Math::Polynomial poly(4);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
        }
        else if(n==5)/////// Equation holds only for n = 5  ///////////////
        {
                double a5 = 10*Nw;
                double a4 =-20*Nw;
                double a3 = 10*Nw+( 10*eb*eb-8*eb+1)*Ntt;
                double a2 =       (-12*eb*eb+6*eb  )*Ntt;
                double a1 =       (  3*eb*eb       )*Ntt;
                double a0 = -N3bjets;
                double par[6] = {a0,a1,a2,a3,a4,a5};

                ROOT::Math::Polynomial poly(5);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
        }
        else /////// Equation holds only for n = 6 but is used for n>6 as well  ///////////////
        {
                double a6 =-20*Nw;
                double a5 = 60*Nw;
                double a4 =-60*Nw+(-20*eb*eb+20*eb-4)*Ntt;
                double a3 = 20*Nw+( 40*eb*eb-32*eb+4)*Ntt;
                double a2 =       (-24*eb*eb+12*eb  )*Ntt;
                double a1 =       (  4*eb*eb        )*Ntt;
                double a0 = -N3bjets;
                double par[7] = {a0,a1,a2,a3,a4,a5,a6};

                ROOT::Math::Polynomial poly(6);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
        }
        return roots;
}
vector<double> WJetEstimation::eudsc_cand_fromN2bjets(double &N2bjets, double &Ntt, double &Nw, double &eb, int &n) const
{    
     vector<double> roots;
     if(N2bjets<=0 || Ntt<=0 || Nw<=0) return roots;

     if(n==4)/////// Equation holds only for n = 4  ///////////////
     {
                double a4 =  6*Nw;
                double a3 =-12*Nw;
                double a2 =  6*Nw+(   eb*eb-4*eb*(1-eb)+pow(1-eb,2))*Ntt;
                double a1 =       (-2*eb*eb+4*eb*(1-eb)            )*Ntt;
                double a0 =       (   eb*eb                        )*Ntt-N2bjets;
                double par[5] = {a0,a1,a2,a3,a4};

                ROOT::Math::Polynomial poly(4);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
     }
     else if(n==5)/////// Equation holds only for n = 5  ///////////////
     {
                double a5 =-10*Nw;
                double a4 = 30*Nw;
                double a3 =-30*Nw+(  -eb*eb+ 6*eb*(1-eb)-3*pow(1-eb,2))*Ntt;
                double a2 = 10*Nw+( 3*eb*eb-12*eb*(1-eb)+3*pow(1-eb,2))*Ntt;
                double a1 =       (-3*eb*eb+ 6*eb*(1-eb)              )*Ntt;
                double a0 =       (   eb*eb                           )*Ntt-N2bjets;
                double par[6] = {a0,a1,a2,a3,a4,a5};

                ROOT::Math::Polynomial poly(5);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
     }
     else /////// Equation holds only for n = 6 but is used for n>6 as well  ///////////////
     {
                double a6 = 15*Nw;
                double a5 =-60*Nw;
                double a4 = 90*Nw+(   eb*eb- 8*eb*(1-eb)+ 6*pow(1-eb,2))*Ntt;
                double a3 =-60*Nw+(-4*eb*eb+24*eb*(1-eb)-12*pow(1-eb,2))*Ntt;
                double a2 = 15*Nw+( 6*eb*eb-24*eb*(1-eb)+ 6*pow(1-eb,2))*Ntt;
                double a1 =       (-4*eb*eb+ 8*eb*(1-eb)               )*Ntt;
                double a0 =       (   eb*eb                            )*Ntt-N2bjets;
                double par[7] = {a0,a1,a2,a3,a4,a5,a6};

                ROOT::Math::Polynomial poly(6);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
     }
     return roots;
}

vector<double> WJetEstimation::eudsc_cand_fromN1bjet(double &N1bjet, double &Ntt, double &Nw, double &eb, int &n) const
{    
     vector<double> roots;
     if(N1bjet <=0 || Ntt<=0 || Nw<=0) return roots;

     if(n==4)/////// Equation holds only for n = 4  ///////////////
     {
                double a4 = -4*Nw;
                double a3 = 12*Nw;
                double a2 =-12*Nw+(-4*eb*eb+6*eb-2)*Ntt;
                double a1 =  4*Nw+( 6*eb*eb-8*eb+2)*Ntt;
                double a0 = -N1bjet+2*eb*(1-eb)*Ntt;
                double par[5] = {a0,a1,a2,a3,a4};

                ROOT::Math::Polynomial poly(4);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
     }
     else if(n==5)/////// Equation holds only for n = 5  ///////////////
     {
                double a5 =  5*Nw;
                double a4 =-20*Nw;
                double a3 = 30*Nw+(-2*eb*(1-eb)+3*pow(1-eb,2))*Ntt;
                double a2 =-20*Nw+( 6*eb*(1-eb)-6*pow(1-eb,2))*Ntt;
                double a1 =  5*Nw+(-6*eb*(1-eb)+3*pow(1-eb,2))*Ntt;
                double a0 =       ( 2*eb*(1-eb)              )*Ntt-N1bjet;
                double par[6] = {a0,a1,a2,a3,a4,a5};

                ROOT::Math::Polynomial poly(5);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
     }
     else /////// Equation holds only for n = 6 but is used for n>6 as well  ///////////////
     {
                double a6 = 15*Nw;
                double a5 =-60*Nw;
                double a4 = 90*Nw+(   eb*eb- 8*eb*(1-eb)+ 6*pow(1-eb,2))*Ntt;
                double a3 =-60*Nw+(-4*eb*eb+24*eb*(1-eb)-12*pow(1-eb,2))*Ntt;
                double a2 = 15*Nw+( 6*eb*eb-24*eb*(1-eb)+ 6*pow(1-eb,2))*Ntt;
                double a1 =       (-4*eb*eb+ 8*eb*(1-eb)               )*Ntt;
                double a0 =       (   eb*eb                            )*Ntt-N1bjet;
                double par[7] = {a0,a1,a2,a3,a4,a5,a6};

                ROOT::Math::Polynomial poly(6);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
     }
     return roots;
}

vector<double> WJetEstimation::eudsc_cand_fromN0bjet(double &N0bjet, double &Ntt, double &Nw, double &eb, int &n) const
{    
     vector<double> roots;
     if(N0bjet <= 0 || Ntt <=0 || Nw <=0) return roots;

     if(n==4)/////// Equation holds only for n = 4  ///////////////
     {
                double a4 =    Nw;
                double a3 = -4*Nw;
                double a2 =  6*Nw+1*pow(1-eb,2)*Ntt;
                double a1 = -4*Nw-2*pow(1-eb,2)*Ntt;
                double a0 =    Nw+1*pow(1-eb,2)*Ntt-N0bjet;
                double par[5] = {a0,a1,a2,a3,a4};

                ROOT::Math::Polynomial poly(4);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
     }
     if(n==5)/////// Equation holds only for n = 5  ///////////////
     {
                double a5 =    Nw;
                double a4 =  5*Nw;
                double a3 =-10*Nw-1*pow(1-eb,2)*Ntt;
                double a2 = 10*Nw+3*pow(1-eb,2)*Ntt;
                double a1 = -5*Nw-3*pow(1-eb,2)*Ntt;
                double a0 =    Nw+1*pow(1-eb,2)*Ntt-N0bjet;
                double par[6] = {a0,a1,a2,a3,a4,a5};

                ROOT::Math::Polynomial poly(5);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
     }
     else /////// Equation holds only for n = 6 but is used for n>6 as well  ///////////////
     {
                double a6 =    Nw;
                double a5 = -6*Nw;
                double a4 = 15*Nw+1*pow(1-eb,2)*Ntt;
                double a3 =-20*Nw-4*pow(1-eb,2)*Ntt;
                double a2 = 15*Nw+6*pow(1-eb,2)*Ntt;
                double a1 = -6*Nw-4*pow(1-eb,2)*Ntt;
                double a0 =    Nw+1*pow(1-eb,2)*Ntt-N0bjet;
                double par[7] = {a0,a1,a2,a3,a4,a5,a6};

                ROOT::Math::Polynomial poly(6);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
     }
     return roots;
}

vector<double> WJetEstimation::eb_cand_fromN3bjets(double &N3bjets, double &Ntt, double &Nw, double &eudsc, int &n) const
{
	vector<double> roots;
	
       	double a2 = ((n-2)*eudsc*pow(1-eudsc,n-3) - (n-2)*(n-3)*eudsc*eudsc*pow(1-eudsc,n-4) + ((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow(1-eudsc,n-5))*Ntt;
       	double a1 = ((n-2)*(n-3)*eudsc*eudsc*pow(1-eudsc,n-4) - 2*((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow(1-eudsc,n-5))*Ntt;
       	double a0 = ((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow(1-eudsc,n-5)*Ntt + ((n)*(n-1)*(n-2)/6)*pow(eudsc,3)*pow(1-eudsc,n-3)*Nw - N3bjets;
	double par[3] = {a0,a1,a2};
	
       	ROOT::Math::Polynomial poly(2);
	poly.SetParameters(par);
       	
	if(N3bjets>0) roots = poly.FindRealRoots();

	return roots;
}

vector<double> WJetEstimation::eb_cand_fromN2bjets(double &N2bjets, double &Ntt, double &Nw, double &eudsc, int &n) const
{
	vector<double> roots;

       	double a2 = (pow(1-eudsc,n-2) - 2*(n-2)*eudsc*pow(1-eudsc,n-3) + ((n-2)*(n-3)/2)*pow(eudsc,2)*pow(1-eudsc,n-4))*Ntt;
       	double a1 = (2*(n-2)*eudsc*pow(1-eudsc,n-3) - 2*((n-2)*(n-3)/2)*pow(eudsc,2)*pow(1-eudsc,n-4))*Ntt;
       	double a0 = (((n-2)*(n-3)/2)*pow(eudsc,2)*pow(1-eudsc,n-4))*Ntt + ((n)*(n-1)/2)*pow(eudsc,2)*pow(1-eudsc,n-2)*Nw - N2bjets;
	double par[3] = {a0,a1,a2};

       	ROOT::Math::Polynomial poly(2);
	poly.SetParameters(par);

       	if(N2bjets>0) roots = poly.FindRealRoots();

	return roots;
}

vector<double> WJetEstimation::eb_cand_fromN1bjet(double &N1bjet, double &Ntt, double &Nw, double &eudsc, int &n) const
{
	vector<double> roots;

       	double a2 = (-2*pow(1-eudsc,n-2) +   (n-2)*eudsc*pow(1-eudsc,n-3))*Ntt;
       	double a1 = ( 2*pow(1-eudsc,n-2) - 2*(n-2)*eudsc*pow(1-eudsc,n-3))*Ntt;
       	double a0 = (n-2)*eudsc*pow(1-eudsc,n-3)*Ntt + n*eudsc*pow(1-eudsc,n-1)*Nw - N1bjet;
	double par[3] = {a0,a1,a2};

       	ROOT::Math::Polynomial poly(2);
	poly.SetParameters(par);

       	if(N1bjet>0) roots = poly.FindRealRoots();

	return roots;
}

vector<double> WJetEstimation::eb_cand_fromN0bjet(double &N0bjet, double &Ntt, double &Nw, double &eudsc, int &n) const
{
	vector<double> roots;

       	double a2 =    pow(1-eudsc,n-2)*Ntt;
       	double a1 = -2*pow(1-eudsc,n-2)*Ntt;
       	double a0 =    pow(1-eudsc,n-2)*Ntt + pow(1-eudsc,n)*Nw - N0bjet;
	double par[3] = {a0,a1,a2};

       	ROOT::Math::Polynomial poly(2);
	poly.SetParameters(par);

       	if(N0bjet>0) roots = poly.FindRealRoots();

	return roots;
}

double WJetEstimation::GetGlobalEffudscEstimate(double &n3bjets, double &n2bjets, double &n1bjet, double &n0bjet, double &Ntt, double &Nw, double &eb, double &eudsc, int &n, bool &Print) const
{
       map<double,double> eudsc_cand;
       double N[4] = {0,0,0,0};
       if(Print)std::cout<<"--Filling vectors with eudsc candidates"<<std::endl;
       vector<double> eudsc_cand_N3bjets; if(n3bjets >= 0) {eudsc_cand_N3bjets = eudsc_cand_fromN3bjets(n3bjets, Ntt, Nw, eb, n); N[3]= 1;}
       vector<double> eudsc_cand_N2bjets; if(n2bjets >= 0) {eudsc_cand_N2bjets = eudsc_cand_fromN2bjets(n2bjets, Ntt, Nw, eb, n); N[2]= 1;}
       vector<double> eudsc_cand_N1bjet ; if(n1bjet  >= 0) {eudsc_cand_N1bjet  = eudsc_cand_fromN1bjet( n1bjet , Ntt, Nw, eb, n); N[1]= 1;}
       vector<double> eudsc_cand_N0bjet ; if(n0bjet  >= 0) {eudsc_cand_N0bjet  = eudsc_cand_fromN0bjet( n0bjet , Ntt, Nw, eb, n); N[0]= 1;}
       
       double ChiSq[4]  = {999999,999999,999999,999999};
       double ChiSq_tmp =  0;
       int idx = -1;
       if(Print){
       		std::cout<<"--Calculating eudsc candidate from N3bjets"<<std::endl;
       		std::cout<<"--- eudsc candidate size : "<<eudsc_cand_N3bjets.size()<<std::endl;
       }
       for(unsigned int i = 0 ; i < eudsc_cand_N3bjets.size(); ++i)
       {
       		if(eudsc_cand_N3bjets[i]<0 || eudsc_cand_N3bjets[i]>1) continue;
       		if(Print){
			std::cout<<"---- eudsc candidate = "<<eudsc_cand_N3bjets[i]<<std::endl;
			std::cout<<"---- n3bjets/ N3bjets = "<<n3bjets<<"/"<<N3bjets(Ntt, Nw, eb, eudsc_cand_N3bjets[i], n)<<std::endl;
			std::cout<<"---- n2bjets/ N2bjets = "<<n2bjets<<"/"<<N2bjets(Ntt, Nw, eb, eudsc_cand_N3bjets[i], n)<<std::endl;
			std::cout<<"---- n1bjet / N1bjet  = "<<n1bjet<<"/"<< N1bjet( Ntt, Nw, eb, eudsc_cand_N3bjets[i], n)<<std::endl;
			std::cout<<"---- n0bjet / N0bjet  = "<<n0bjet<<"/"<< N0bjet( Ntt, Nw, eb, eudsc_cand_N3bjets[i], n)<<std::endl;
		}

                ChiSq_tmp = (n2bjets>0 ? pow(n2bjets - (double)N2bjets(Ntt, Nw, eb, eudsc_cand_N3bjets[i], n),2)/n2bjets : pow(n2bjets - (double)N2bjets(Ntt, Nw, eb, eudsc_cand_N3bjets[i], n),2))
                	  + (n1bjet >0 ? pow(n1bjet  - (double)N1bjet( Ntt, Nw, eb, eudsc_cand_N3bjets[i], n),2)/n1bjet  : pow(n1bjet  - (double)N1bjet( Ntt, Nw, eb, eudsc_cand_N3bjets[i], n),2))
                          + (n0bjet >0 ? pow(n0bjet  - (double)N0bjet( Ntt, Nw, eb, eudsc_cand_N3bjets[i], n),2)/n0bjet  : pow(n0bjet  - (double)N0bjet( Ntt, Nw, eb, eudsc_cand_N3bjets[i], n),2));

		ChiSq_tmp /= (N[2]+N[1]+N[0]!=0 ? N[2]+N[1]+N[0] : 1);
		
		if(Print)std::cout<<"---- ChiSq_tmp = "<<ChiSq_tmp<<std::endl;
                if(ChiSq_tmp<ChiSq[0])
                {
                	ChiSq[0] = ChiSq_tmp;
                	idx = i;
                }
       }
       if(idx != -1) eudsc_cand[ChiSq[0]] = eudsc_cand_N3bjets[idx];
       
       ChiSq_tmp =  0;
       idx = -1;
       if(Print){
       		std::cout<<"--Calculating eudsc candidate from N2bjets"<<std::endl;
       		std::cout<<"--- eudsc candidate size : "<<eudsc_cand_N2bjets.size()<<std::endl;
       }
       for(unsigned int i = 0 ; i < eudsc_cand_N2bjets.size(); ++i)
       {
       		if(eudsc_cand_N2bjets[i]<0 || eudsc_cand_N2bjets[i]>1) continue;
       		if(Print){
			std::cout<<"---- eudsc candidate = "<<eudsc_cand_N2bjets[i]<<std::endl;
			std::cout<<"---- n3bjets/ N3bjets = "<<n3bjets<<"/"<<N3bjets(Ntt, Nw, eb, eudsc_cand_N2bjets[i], n)<<std::endl;
			std::cout<<"---- n2bjets/ N2bjets = "<<n2bjets<<"/"<<N2bjets(Ntt, Nw, eb, eudsc_cand_N2bjets[i], n)<<std::endl;
			std::cout<<"---- n1bjet / N1bjet  = "<<n1bjet<<"/"<< N1bjet( Ntt, Nw, eb, eudsc_cand_N2bjets[i], n)<<std::endl;
			std::cout<<"---- n0bjet / N0bjet  = "<<n0bjet<<"/"<< N0bjet( Ntt, Nw, eb, eudsc_cand_N2bjets[i], n)<<std::endl;
		}
		ChiSq_tmp = (n3bjets>0 ? pow(n3bjets - (double)N3bjets(Ntt, Nw, eb, eudsc_cand_N2bjets[i], n),2)/n3bjets : pow(n3bjets - (double)N3bjets(Ntt, Nw, eb, eudsc_cand_N2bjets[i], n),2))
			  + (n1bjet >0 ? pow(n1bjet  - (double)N1bjet( Ntt, Nw, eb, eudsc_cand_N2bjets[i], n),2)/n1bjet  : pow(n1bjet  - (double)N1bjet( Ntt, Nw, eb, eudsc_cand_N2bjets[i], n),2))
			  + (n0bjet >0 ? pow(n0bjet  - (double)N0bjet( Ntt, Nw, eb, eudsc_cand_N2bjets[i], n),2)/n0bjet  : pow(n0bjet  - (double)N0bjet( Ntt, Nw, eb, eudsc_cand_N2bjets[i], n),2));

		ChiSq_tmp /= (N[3]+N[1]+N[0]!=0 ? N[3]+N[1]+N[0] : 1);

		if(Print)std::cout<<"---- ChiSq_tmp = "<<ChiSq_tmp<<std::endl;
		if(ChiSq_tmp<ChiSq[1])
		{
			ChiSq[1] = ChiSq_tmp;
			idx = i;
		}
       }
       if(idx != -1) eudsc_cand[ChiSq[1]] = eudsc_cand_N2bjets[idx];
       
       ChiSq_tmp =  0;
       idx = -1;
       if(Print){
       		std::cout<<"--Calculating eudsc candidate from N1bjet"<<std::endl;
       		std::cout<<"--- eudsc candidate size : "<<eudsc_cand_N1bjet.size()<<std::endl;
       }
       for(unsigned int i = 0 ; i < eudsc_cand_N1bjet.size(); i++)
       {
       		if(eudsc_cand_N1bjet[i]<0 || eudsc_cand_N1bjet[i]>1) continue;
       		if(Print){
			std::cout<<"---- eudsc candidate = "<<eudsc_cand_N1bjet[i]<<std::endl;
			std::cout<<"---- n3bjets/ N3bjets = "<<n3bjets<<"/"<<N3bjets(Ntt, Nw, eb, eudsc_cand_N1bjet[i], n)<<std::endl;
			std::cout<<"---- n2bjets/ N2bjets = "<<n2bjets<<"/"<<N2bjets(Ntt, Nw, eb, eudsc_cand_N1bjet[i], n)<<std::endl;
			std::cout<<"---- n1bjet / N1bjet  = "<<n1bjet<<"/"<< N1bjet( Ntt, Nw, eb, eudsc_cand_N1bjet[i], n)<<std::endl;
			std::cout<<"---- n0bjet / N0bjet  = "<<n0bjet<<"/"<< N0bjet( Ntt, Nw, eb, eudsc_cand_N1bjet[i], n)<<std::endl;
		}
		ChiSq_tmp = (n3bjets>0 ? pow(n3bjets - (double)N3bjets(Ntt, Nw, eb, eudsc_cand_N1bjet[i], n),2)/n3bjets : pow(n3bjets - (double)N3bjets(Ntt, Nw, eb, eudsc_cand_N1bjet[i], n),2))
			  + (n2bjets>0 ? pow(n2bjets - (double)N2bjets(Ntt, Nw, eb, eudsc_cand_N1bjet[i], n),2)/n2bjets : pow(n2bjets - (double)N2bjets(Ntt, Nw, eb, eudsc_cand_N1bjet[i], n),2))
			  + (n0bjet >0 ? pow(n0bjet  - (double)N0bjet( Ntt, Nw, eb, eudsc_cand_N1bjet[i], n),2)/n0bjet  : pow(n0bjet  - (double)N0bjet( Ntt, Nw, eb, eudsc_cand_N1bjet[i], n),2));

		ChiSq_tmp /= (N[3]+N[2]+N[0]!=0 ? N[3]+N[2]+N[0] : 1);

		if(Print)std::cout<<"---- ChiSq_tmp = "<<ChiSq_tmp<<std::endl;
		if(ChiSq_tmp<ChiSq[2])
		{
			ChiSq[2] = ChiSq_tmp;
			idx = i;
		}
       }
       if(idx != -1) eudsc_cand[ChiSq[2]] = eudsc_cand_N1bjet[idx];
       
       ChiSq_tmp =  0;
       idx = -1;
       if(Print){
       		std::cout<<"--Calculating eudsc candidate from N0bjet"<<std::endl;
       		std::cout<<"--- eudsc candidate size : "<<eudsc_cand_N0bjet.size()<<std::endl;
       }
       for(unsigned int i = 0 ; i < eudsc_cand_N0bjet.size(); i++)
       {
       		if(eudsc_cand_N0bjet[i]<0 || eudsc_cand_N0bjet[i]>1) continue;
       		if(Print){
			std::cout<<"---- eudsc candidate = "<<eudsc_cand_N0bjet[i]<<std::endl;
			std::cout<<"---- n3bjets/ N3bjets = "<<n3bjets<<"/"<<N3bjets(Ntt, Nw, eb, eudsc_cand_N0bjet[i], n)<<std::endl;
			std::cout<<"---- n2bjets/ N2bjets = "<<n2bjets<<"/"<<N2bjets(Ntt, Nw, eb, eudsc_cand_N0bjet[i], n)<<std::endl;
			std::cout<<"---- n1bjet / N1bjet  = "<<n1bjet<<"/"<< N1bjet( Ntt, Nw, eb, eudsc_cand_N0bjet[i], n)<<std::endl;
			std::cout<<"---- n0bjet / N0bjet  = "<<n0bjet<<"/"<< N0bjet( Ntt, Nw, eb, eudsc_cand_N0bjet[i], n)<<std::endl;
                }
		ChiSq_tmp = (n3bjets>0 ? pow(n3bjets - (double)N3bjets(Ntt, Nw, eb, eudsc_cand_N0bjet[i], n),2)/n3bjets : pow(n3bjets - (double)N3bjets(Ntt, Nw, eb, eudsc_cand_N0bjet[i], n),2))
			  + (n2bjets>0 ? pow(n2bjets - (double)N2bjets(Ntt, Nw, eb, eudsc_cand_N0bjet[i], n),2)/n2bjets : pow(n2bjets - (double)N2bjets(Ntt, Nw, eb, eudsc_cand_N0bjet[i], n),2))
			  + (n1bjet >0 ? pow(n1bjet  - (double)N1bjet( Ntt, Nw, eb, eudsc_cand_N0bjet[i], n),2)/n1bjet  : pow(n1bjet  - (double)N1bjet( Ntt, Nw, eb, eudsc_cand_N0bjet[i], n),2));

		ChiSq_tmp /= (N[3]+N[2]+N[1]!=0 ? N[3]+N[2]+N[1] : 1);

		if(Print)std::cout<<"---- ChiSq_tmp = "<<ChiSq_tmp<<std::endl;
		if(ChiSq_tmp<ChiSq[3])
		{
			ChiSq[3] = ChiSq_tmp;
			idx = i;
		}
       }
       if(idx != -1) eudsc_cand[ChiSq[3]] = eudsc_cand_N0bjet[idx];
       
       ChiSq_tmp = 999999;
       std::map<double,double>::iterator it;
       if(Print)std::cout<<"--Calculating eudsc candidate with lowest X2"<<std::endl;
       for ( it=eudsc_cand.begin() ; it !=eudsc_cand.end(); ++it )
       {
                    if(it->first<ChiSq_tmp) ChiSq_tmp = it->first;
       }
       eudsc = (eudsc_cand.find(ChiSq_tmp) != eudsc_cand.end() ? eudsc_cand[ChiSq_tmp] : eudsc);
       if(Print)std::cout<<"--- eudsc = "<<eudsc<<std::endl;
       if(eudsc!=0) return ChiSq_tmp;
       else         return -1;
}

double WJetEstimation::GetGlobalEffbEstimate(double &n3bjets, double &n2bjets, double &n1bjet, double &n0bjet, double &Ntt, double &Nw, double &eb, double &eudsc, int &n, bool &Print) const
{
       map<double,double> eb_cand;
       double N[4] = {0,0,0,0};
       if(Print)std::cout<<"--Filling vectors with eb candidates"<<std::endl;

       vector<double> eb_cand_N3bjets; if(n3bjets >= 0) {eb_cand_N3bjets = eb_cand_fromN3bjets(n3bjets, Ntt, Nw, eudsc, n); N[3]= 1;}
       vector<double> eb_cand_N2bjets; if(n2bjets >= 0) {eb_cand_N2bjets = eb_cand_fromN2bjets(n2bjets, Ntt, Nw, eudsc, n); N[2]= 1;}
       vector<double> eb_cand_N1bjet ; if(n1bjet  >= 0) {eb_cand_N1bjet  = eb_cand_fromN1bjet( n1bjet , Ntt, Nw, eudsc, n); N[1]= 1;}
       vector<double> eb_cand_N0bjet ; if(n0bjet  >= 0) {eb_cand_N0bjet  = eb_cand_fromN0bjet( n0bjet , Ntt, Nw, eudsc, n); N[0]= 1;}
       
       double ChiSq[4]  = {999999,999999,999999,999999};
       double ChiSq_tmp =  0;
       int idx = -1;
       if(Print){
       		std::cout<<"--Calculating eb candidate from N3bjets"<<std::endl;
       		std::cout<<"--- eb candidate size : "<<eb_cand_N3bjets.size()<<std::endl;
       }
       for(unsigned int i = 0 ; i < eb_cand_N3bjets.size(); ++i)
       {
       		if(eb_cand_N3bjets[i]<0 || eb_cand_N3bjets[i]>1) continue;
       		if(Print){
			std::cout<<"---- eb candidate = "<<eb_cand_N3bjets[i]<<std::endl;
			std::cout<<"---- n3bjets/ N3bjets = "<<n3bjets<<"/"<<N3bjets(Ntt, Nw, eb_cand_N3bjets[i], eudsc, n)<<std::endl;
			std::cout<<"---- n2bjets/ N2bjets = "<<n2bjets<<"/"<<N2bjets(Ntt, Nw, eb_cand_N3bjets[i], eudsc, n)<<std::endl;
			std::cout<<"---- n1bjet / N1bjet  = "<<n1bjet<<"/"<< N1bjet( Ntt, Nw, eb_cand_N3bjets[i], eudsc, n)<<std::endl;
			std::cout<<"---- n0bjet / N0bjet  = "<<n0bjet<<"/"<< N0bjet( Ntt, Nw, eb_cand_N3bjets[i], eudsc, n)<<std::endl;
		}

                ChiSq_tmp = (n2bjets>0 ? pow(n2bjets - (double)N2bjets(Ntt, Nw, eb_cand_N3bjets[i], eudsc, n),2)/n2bjets : pow(n2bjets - (double)N2bjets(Ntt, Nw, eb_cand_N3bjets[i], eudsc, n),2))
                	  + (n1bjet >0 ? pow(n1bjet  - (double)N1bjet( Ntt, Nw, eb_cand_N3bjets[i], eudsc, n),2)/n1bjet  : pow(n1bjet  - (double)N1bjet( Ntt, Nw, eb_cand_N3bjets[i], eudsc, n),2))
                          + (n0bjet >0 ? pow(n0bjet  - (double)N0bjet( Ntt, Nw, eb_cand_N3bjets[i], eudsc, n),2)/n0bjet  : pow(n0bjet  - (double)N0bjet( Ntt, Nw, eb_cand_N3bjets[i], eudsc, n),2));
		
		ChiSq_tmp /= (N[2]+N[1]+N[0]!=0 ? N[2]+N[1]+N[0] : 1);
		
		if(Print)std::cout<<"---- ChiSq_tmp = "<<ChiSq_tmp<<std::endl;
                if(ChiSq_tmp<ChiSq[0])
                {
                	ChiSq[0] = ChiSq_tmp;
                	idx = i;
                }
       }
       if(idx != -1) eb_cand[ChiSq[0]] = eb_cand_N3bjets[idx];
       
       ChiSq_tmp =  0;
       idx = -1;
       if(Print){
       		std::cout<<"--Calculating eb candidate from N2bjets"<<std::endl;
       		std::cout<<"--- eb candidate size : "<<eb_cand_N2bjets.size()<<std::endl;
       }
       for(unsigned int i = 0 ; i < eb_cand_N2bjets.size(); ++i)
       {
       		if(eb_cand_N2bjets[i]<0 || eb_cand_N2bjets[i]>1) continue;
       		if(Print){
			std::cout<<"---- eb candidate = "<<eb_cand_N2bjets[i]<<std::endl;
			std::cout<<"---- n3bjets/ N3bjets = "<<n3bjets<<"/"<<N3bjets(Ntt, Nw, eb_cand_N2bjets[i], eudsc, n)<<std::endl;
			std::cout<<"---- n2bjets/ N2bjets = "<<n2bjets<<"/"<<N2bjets(Ntt, Nw, eb_cand_N2bjets[i], eudsc, n)<<std::endl;
			std::cout<<"---- n1bjet / N1bjet  = "<<n1bjet<<"/"<< N1bjet( Ntt, Nw, eb_cand_N2bjets[i], eudsc, n)<<std::endl;
			std::cout<<"---- n0bjet / N0bjet  = "<<n0bjet<<"/"<< N0bjet( Ntt, Nw, eb_cand_N2bjets[i], eudsc, n)<<std::endl;
		}

		ChiSq_tmp = (n3bjets>0 ? pow(n3bjets - (double)N3bjets(Ntt, Nw, eb_cand_N2bjets[i], eudsc, n),2)/n3bjets : pow(n3bjets - (double)N3bjets(Ntt, Nw, eb_cand_N2bjets[i], eudsc, n),2))
			  + (n1bjet >0 ? pow(n1bjet  - (double)N1bjet( Ntt, Nw, eb_cand_N2bjets[i], eudsc, n),2)/n1bjet  : pow(n1bjet  - (double)N1bjet( Ntt, Nw, eb_cand_N2bjets[i], eudsc, n),2))
			  + (n0bjet >0 ? pow(n0bjet  - (double)N0bjet( Ntt, Nw, eb_cand_N2bjets[i], eudsc, n),2)/n0bjet  : pow(n0bjet  - (double)N0bjet( Ntt, Nw, eb_cand_N2bjets[i], eudsc, n),2));

		ChiSq_tmp /= (N[3]+N[1]+N[0]!=0 ? N[3]+N[1]+N[0] : 1);

		if(Print)std::cout<<"---- ChiSq_tmp = "<<ChiSq_tmp<<std::endl;
		if(ChiSq_tmp<ChiSq[1])
		{
			ChiSq[1] = ChiSq_tmp;
			idx = i;
		}
       }
       if(idx != -1) eb_cand[ChiSq[1]] = eb_cand_N2bjets[idx];
       
       ChiSq_tmp =  0;
       idx = -1;
       if(Print){
       		std::cout<<"--Calculating eb candidate from N1bjet"<<std::endl;
       		std::cout<<"--- eb candidate size : "<<eb_cand_N1bjet.size()<<std::endl;
       }
       for(unsigned int i = 0 ; i < eb_cand_N1bjet.size(); i++)
       {
       		if(eb_cand_N1bjet[i]<0 || eb_cand_N1bjet[i]>1) continue;
       		if(Print){
			std::cout<<"---- eb candidate = "<<eb_cand_N1bjet[i]<<std::endl;
			std::cout<<"---- n3bjets/ N3bjets = "<<n3bjets<<"/"<<N3bjets(Ntt, Nw, eb_cand_N1bjet[i], eudsc, n)<<std::endl;
			std::cout<<"---- n2bjets/ N2bjets = "<<n2bjets<<"/"<<N2bjets(Ntt, Nw, eb_cand_N1bjet[i], eudsc, n)<<std::endl;
			std::cout<<"---- n1bjet / N1bjet  = "<<n1bjet<<"/"<< N1bjet( Ntt, Nw, eb_cand_N1bjet[i], eudsc, n)<<std::endl;
			std::cout<<"---- n0bjet / N0bjet  = "<<n0bjet<<"/"<< N0bjet( Ntt, Nw, eb_cand_N1bjet[i], eudsc, n)<<std::endl;
		}
		ChiSq_tmp = (n3bjets>0 ? pow(n3bjets - (double)N3bjets(Ntt, Nw, eb_cand_N1bjet[i], eudsc, n),2)/n3bjets : pow(n3bjets - (double)N3bjets(Ntt, Nw, eb_cand_N1bjet[i], eudsc, n),2))
			  + (n2bjets>0 ? pow(n2bjets - (double)N2bjets(Ntt, Nw, eb_cand_N1bjet[i], eudsc, n),2)/n2bjets : pow(n2bjets - (double)N2bjets(Ntt, Nw, eb_cand_N1bjet[i], eudsc, n),2))
			  + (n0bjet >0 ? pow(n0bjet  - (double)N0bjet( Ntt, Nw, eb_cand_N1bjet[i], eudsc, n),2)/n0bjet  : pow(n0bjet  - (double)N0bjet( Ntt, Nw, eb_cand_N1bjet[i], eudsc, n),2));

		ChiSq_tmp /= (N[3]+N[2]+N[0]!=0 ? N[3]+N[2]+N[0] : 1);

		if(Print)std::cout<<"---- ChiSq_tmp = "<<ChiSq_tmp<<std::endl;
		if(ChiSq_tmp<ChiSq[2])
		{
			ChiSq[2] = ChiSq_tmp;
			idx = i;
		}
       }
       if(idx != -1) eb_cand[ChiSq[2]] = eb_cand_N1bjet[idx];
       
       ChiSq_tmp =  0;
       idx = -1;
       if(Print){
       		std::cout<<"--Calculating eb candidate from N0bjet"<<std::endl;
       		std::cout<<"--- eb candidate size : "<<eb_cand_N0bjet.size()<<std::endl;
       }
       for(unsigned int i = 0 ; i < eb_cand_N0bjet.size(); i++)
       {
       		if(eb_cand_N0bjet[i]<0 || eb_cand_N0bjet[i]>1) continue;
       		if(Print){
			std::cout<<"---- eb candidate = "<<eb_cand_N0bjet[i]<<std::endl;
			std::cout<<"---- n3bjets/ N3bjets = "<<n3bjets<<"/"<<N3bjets(Ntt, Nw, eb_cand_N0bjet[i], eudsc, n)<<std::endl;
			std::cout<<"---- n2bjets/ N2bjets = "<<n2bjets<<"/"<<N2bjets(Ntt, Nw, eb_cand_N0bjet[i], eudsc, n)<<std::endl;
			std::cout<<"---- n1bjet / N1bjet  = "<<n1bjet<<"/"<< N1bjet( Ntt, Nw, eb_cand_N0bjet[i], eudsc, n)<<std::endl;
			std::cout<<"---- n0bjet / N0bjet  = "<<n0bjet<<"/"<< N0bjet( Ntt, Nw, eb_cand_N0bjet[i], eudsc, n)<<std::endl;
                }
		ChiSq_tmp = (n3bjets>0 ? pow(n3bjets - (double)N3bjets(Ntt, Nw, eb_cand_N0bjet[i], eudsc, n),2)/n3bjets : pow(n3bjets - (double)N3bjets(Ntt, Nw, eb_cand_N0bjet[i], eudsc, n),2))
			  + (n2bjets>0 ? pow(n2bjets - (double)N2bjets(Ntt, Nw, eb_cand_N0bjet[i], eudsc, n),2)/n2bjets : pow(n2bjets - (double)N2bjets(Ntt, Nw, eb_cand_N0bjet[i], eudsc, n),2))
			  + (n1bjet >0 ? pow(n1bjet  - (double)N1bjet( Ntt, Nw, eb_cand_N0bjet[i], eudsc, n),2)/n1bjet  : pow(n1bjet  - (double)N1bjet( Ntt, Nw, eb_cand_N0bjet[i], eudsc, n),2));

		ChiSq_tmp /= (N[3]+N[2]+N[1]!=0 ? N[3]+N[2]+N[1] : 1);

		if(Print)std::cout<<"---- ChiSq_tmp = "<<ChiSq_tmp<<std::endl;
		if(ChiSq_tmp<ChiSq[3])
		{
			ChiSq[3] = ChiSq_tmp;
			idx = i;
		}
       }
       if(idx != -1) eb_cand[ChiSq[3]] = eb_cand_N0bjet[idx];
       
       ChiSq_tmp = 999999;
       std::map<double,double>::iterator it;
       if(Print)std::cout<<"--Calculating eb candidate with lowest X2"<<std::endl;
       for ( it=eb_cand.begin() ; it !=eb_cand.end(); ++it )
       {
                    if(it->first<ChiSq_tmp) ChiSq_tmp = it->first;
       }
       eb = (eb_cand.find(ChiSq_tmp) != eb_cand.end() ? eb_cand[ChiSq_tmp] : eb);
       if(Print)std::cout<<"--- eb = "<<eb<<std::endl;
       if(eb!=0) return ChiSq_tmp;
       else      return -1;
}

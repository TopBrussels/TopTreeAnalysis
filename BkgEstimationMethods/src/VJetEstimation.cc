#include "../interface/VJetEstimation.h"

/**______________________________________________________________________________________________________________________________________________________________________________*/
VJetEstimation::VJetEstimation(int NofBtagWorkingPoint, float* BtagWorkingPoint, int NofJets, int NofJetBins, int NofDatasets, vector<int> iDTTLike, vector<int> iDVLike, vector<int> iDVbLike){
	cout<<"Object from the class VJetEstimation being instantiated"<<endl;
	NbOfDatasets_ = NofDatasets;
	NbOfBtagWorkingPoint_ = NofBtagWorkingPoint;
	BtagWorkingPoint_ = new float[NbOfBtagWorkingPoint_];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++) BtagWorkingPoint_[i] = BtagWorkingPoint[i];
	MCdata_ = true;
	iDatasetsTTLike_.clear();
	iDatasetsVLike_.clear();
	iDatasetsTTLike_ = iDTTLike;
	iDatasetsVLike_  = iDVLike;
	iDatasetsVbLike_ = iDVbLike;
	
	NbOfJetsBins_ = NofJetBins;
	Njets_ = new int[NbOfJetsBins_+1];
	for(int i=0;i<NbOfJetsBins_+1;i++) Njets_[i] = NofJets+i;
	NbOfBJetsBins_ = 4;

	init_ = 0;//
	minValue_ = new double[NbOfJetsBins_];
	for(int i=0;i<NbOfJetsBins_;i++) minValue_[i] = 0;
	
	RescaledTTLikeEstimation = 0;
	RescaledVLikeEstimation = 0;
	RescaledVbLikeEstimation = 0;
	tCanva_RescaledTTLikeEstimation = 0;
	tCanva_RescaledVLikeEstimation = 0;
	tCanva_RescaledVbLikeEstimation = 0;
        
	//init summary histos
	hsNbjets_MC  = new THStack**[NbOfBtagWorkingPoint_];
	hsNbjets_Est = new THStack**[NbOfBtagWorkingPoint_];
	hsNjets_MC   = new THStack**[NbOfBtagWorkingPoint_];
	hsNjets_Est  = new THStack**[NbOfBtagWorkingPoint_];
	MyLeg = new TLegend(0.6,0.6,0.9,0.9);
	
	hNbjets_mc_       = new TH1F***[NbOfBtagWorkingPoint_];
	hNbjets_pdf_mc_   = new TH1F***[NbOfBtagWorkingPoint_];
	hNbjets_pdf_est_  = new TH1F***[NbOfBtagWorkingPoint_];
	
	hNbjetsEstSummary = new TH1F***[NbOfBtagWorkingPoint_];
	hNbjetsMCSummary  = new TH1F***[NbOfBtagWorkingPoint_];
	hNjetsEstSummary  = new TH1F***[NbOfBtagWorkingPoint_];
	hNjetsMCSummary   = new TH1F***[NbOfBtagWorkingPoint_];
	
	tCanva_Nbjets_Summary = new TCanvas**[NbOfBtagWorkingPoint_];
	tCanva_Njets_Summary  = new TCanvas**[NbOfBtagWorkingPoint_];
	
	char name[100];
	char title[100];

	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		hsNbjets_MC[i]  = new THStack*[NbOfJetsBins_+1];
		hsNbjets_Est[i] = new THStack*[NbOfJetsBins_+1];
		hsNjets_MC[i]   = new THStack*[NbOfBJetsBins_+1];
		hsNjets_Est[i]  = new THStack*[NbOfBJetsBins_+1];
	
		hNbjets_mc_[i]       = new TH1F**[NbOfJetsBins_];
		hNbjets_pdf_mc_[i]   = new TH1F**[NbOfJetsBins_];
		hNbjets_pdf_est_[i]  = new TH1F**[NbOfJetsBins_];
		
		hNbjetsEstSummary[i] = new TH1F**[NbOfJetsBins_+1];
		hNbjetsMCSummary[i]  = new TH1F**[NbOfJetsBins_+1];
		hNjetsEstSummary[i]  = new TH1F**[NbOfBJetsBins_+1];
		hNjetsMCSummary[i]   = new TH1F**[NbOfBJetsBins_+1];
		
		tCanva_Nbjets_Summary[i] = new TCanvas*[NbOfJetsBins_+1];
		tCanva_Njets_Summary[i]  = new TCanvas*[NbOfBJetsBins_+1];

		for(int j=0;j<NbOfJetsBins_;j++){
		// Booking histograms for VLike,VbLike and TTlike predictions/estimations
		// for different jet multiplicity

			// estimation
			hsNbjets_MC[i][j]          = 0;
			hsNbjets_Est[i][j]         = 0;

			hNbjets_mc_[i][j]          = new TH1F*[3];
			sprintf(name,"hNbjets_mc_VLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_mc_[i][j][0]       = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_mc_[i][j][0]->Sumw2();
			sprintf(name,"hNbjets_mc_VbLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_mc_[i][j][1]       = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_mc_[i][j][1]->Sumw2();
			sprintf(name,"hNbjets_mc_TTLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_mc_[i][j][2]       = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_mc_[i][j][2]->Sumw2();

			hNbjets_pdf_mc_[i][j]      = new TH1F*[3];
			sprintf(name,"hNbjets_pdf_mc_VLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_mc_[i][j][0]   = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_mc_[i][j][0]->Sumw2();
			sprintf(name,"hNbjets_pdf_mc_VbLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_mc_[i][j][1]   = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_mc_[i][j][1]->Sumw2();
			sprintf(name,"hNbjets_pdf_mc_TTLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_mc_[i][j][2]   = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_mc_[i][j][2]->Sumw2();

			hNbjets_pdf_est_[i][j]     = new TH1F*[3];
			sprintf(name,"hNbjets_pdf_est_VLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_est_[i][j][0]  = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_est_[i][j][0]->Sumw2();
			sprintf(name,"hNbjets_pdf_est_VbLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_est_[i][j][1]  = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_est_[i][j][1]->Sumw2();
			sprintf(name,"hNbjets_pdf_est_TTLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_est_[i][j][2]  = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_est_[i][j][2]->Sumw2();

			hNbjetsEstSummary[i][j]    = new TH1F*[3];
			sprintf(name,"hNbjetsEstSummary_VLike_%d_jets_wp_%d",Njets_[j],i);
			//sprintf(title,"W-like events estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsEstSummary[i][j][0] = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsEstSummary_VbLike_%d_jets_wp_%d",Njets_[j],i);
			//sprintf(title,"Wb-like events estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsEstSummary[i][j][1] = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsEstSummary_TTlike_%d_jets_wp_%d",Njets_[j],i);
			//sprintf(title,"TT-like events estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsEstSummary[i][j][2] = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);

			// MC prediction
			hNbjetsMCSummary[i][j]     = new TH1F*[3];
			sprintf(name,"hNbjetsMCSummary_VLike_%d_jets_wp_%d",Njets_[j],i);
			//sprintf(title,"W-like events MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsMCSummary[i][j][0]  = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsMCSummary_VbLike_%d_jets_wp_%d",Njets_[j],i);
			//sprintf(title,"Wb-like events MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsMCSummary[i][j][1]  = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsMCSummary_TTlike_%d_jets_wp_%d",Njets_[j],i);
			//sprintf(title,"TT-like events MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsMCSummary[i][j][2]  = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);

			// MC prediction vs estimation
			sprintf(name,"hNbjetsSummary_%d_jets_wp_%d",Njets_[j],i);
			sprintf(title,"MC prediction/Estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			tCanva_Nbjets_Summary[i][j] = new TCanvas(name,title,600,600);
		}

		// Estimation summary (inclusive)
		hsNbjets_MC[i][NbOfBJetsBins_]         = 0;
		hsNbjets_Est[i][NbOfBJetsBins_]        = 0;
		hNbjetsEstSummary[i][NbOfJetsBins_]    = new TH1F*[3];
		sprintf(name ,"hNbjetsEstSummary_VLike_Inclusive_wp_%d",i);
		sprintf(title,"W-like events estimation summary (inclusive, WP nr%d)",i);
		hNbjetsEstSummary[i][NbOfJetsBins_][0] = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsEstSummary_VbLike_Inclusive_wp_%d",i);
		sprintf(title,"Wb-like events estimation summary (inclusive, WP nr%d)",i);
		hNbjetsEstSummary[i][NbOfJetsBins_][1] = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsEstSummary_TTlike_Inclusive_wp_%d",i);
		sprintf(title,"TT-like events estimation summary (inclusive, WP nr%d)",i);
		hNbjetsEstSummary[i][NbOfJetsBins_][2] = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);

		// MC prediction summary (inclusive)
		hNbjetsMCSummary[i][NbOfJetsBins_]     = new TH1F*[3];
		sprintf(name ,"hNbjetsMCSummary_VLike_Inclusive_wp_%d",i);
		sprintf(title,"W-like events MC prediction summary (inclusive, WP nr%d)",i);
		hNbjetsMCSummary[i][NbOfJetsBins_][0]  = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsMCSummary_VbLike_Inclusive_wp_%d",i);
		sprintf(title,"Wb-like events MC prediction summary (inclusive, WP nr%d)",i);
		hNbjetsMCSummary[i][NbOfJetsBins_][1]  = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsMCSummary_TTlike_Inclusive_wp_%d",i);
		sprintf(title,"TT-like events MC prediction summary (inclusive, WP nr%d)",i);
		hNbjetsMCSummary[i][NbOfJetsBins_][2]  = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);

		// MC prediction/EStimation summary (inclusive)
		sprintf(name ,"hNbjetsSummary_Inclusive_wp_%d",i);
		sprintf(title,"MC prediction/Estimation summary (inclusive, WP nr%d)",i);
		tCanva_Nbjets_Summary[i][NbOfJetsBins_] = new TCanvas(name,title,600,600);

		for(int j=0;j<NbOfBJetsBins_;j++){
		// Booking histograms for VLike,VbLike and TTlike predictions/estimations
		// for different b-jet multiplicity

			// estimation
			hsNjets_MC[i][j]          = 0;
			hsNjets_Est[i][j]         = 0;
			hNjetsEstSummary[i][j]    = new TH1F*[3];
			sprintf(name,"hNjetsEstSummary_VLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"W-like events estimation summary for %d b-jets (WP nr%d)",j,i);
			hNjetsEstSummary[i][j][0] = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsEstSummary_VbLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"Wb-like events estimation summary for %d b-jets (WP nr%d)",j,i);
			hNjetsEstSummary[i][j][1] = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsEstSummary_TTlike_%d_bjets_wp_%d",j,i);
			sprintf(title,"TT-like events estimation summary for %d b-jets (WP nr%d)",j,i);
			hNjetsEstSummary[i][j][2] = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);

			// MC prediction
			hNjetsMCSummary[i][j]     = new TH1F*[3];
			sprintf(name,"hNjetsMCSummary_VLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"W-like events MC prediction summary for %d b-jets (WP nr%d)",j,i);
			hNjetsMCSummary[i][j][0]  = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsMCSummary_VbLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"Wb-like events MC prediction summary for %d b-jets (WP nr%d)",j,i);
			hNjetsMCSummary[i][j][1]  = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsMCSummary_TTlike_%d_bjets_wp_%d",j,i);
			sprintf(title,"TT-like events MC prediction summary for %d b-jets (WP nr%d)",j,i);
			hNjetsMCSummary[i][j][2]  = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);

			// MC prediction vs estimation
			sprintf(name,"hNjetsSummary_%d_bjets_wp_%d",j,i);
			sprintf(title,"MC prediction/Estimation summary for %d b-jets (WP nr%d)",j,i);		
        		tCanva_Njets_Summary[i][j] = new TCanvas(name,title,600,600);
		}

		// Estimation summary (b-inclusive)
		hsNjets_MC[i][NbOfBJetsBins_]          = 0;
		hsNjets_Est[i][NbOfBJetsBins_]         = 0;
		hNjetsEstSummary[i][NbOfBJetsBins_]    = new TH1F*[3];
		sprintf(name ,"hNjetsEstSummary_VLike_bInclusive_wp_%d",i);
		sprintf(title,"W-like events estimation summary (b-inclusive, WP nr%d)",i);
		hNjetsEstSummary[i][NbOfBJetsBins_][0] = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsEstSummary_VbLike_bInclusive_wp_%d",i);
		sprintf(title,"Wb-like events estimation summary (b-inclusive, WP nr%d)",i);
		hNjetsEstSummary[i][NbOfBJetsBins_][1] = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsEstSummary_TTlike_bInclusive_wp_%d",i);
		sprintf(title,"TT-like events estimation summary (b-inclusive, WP nr%d)",i);
		hNjetsEstSummary[i][NbOfBJetsBins_][2] = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);

		// MC prediction summary (b-inclusive)
		hNjetsMCSummary[i][NbOfBJetsBins_]     = new TH1F*[3];
		sprintf(name ,"hNjetsMCSummary_VLike_bInclusive_wp_%d",i);
		sprintf(title,"W-like events MC prediction summary (b-inclusive, WP nr%d)",i);
		hNjetsMCSummary[i][NbOfBJetsBins_][0]  = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsMCSummary_VbLike_bInclusive_wp_%d",i);
		sprintf(title,"Wb-like events MC prediction summary (b-inclusive, WP nr%d)",i);
		hNjetsMCSummary[i][NbOfBJetsBins_][1]  = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsMCSummary_TTlike_bInclusive_wp_%d",i);
		sprintf(title,"TT-like events MC prediction summary (b-inclusive, WP nr%d)",i);
		hNjetsMCSummary[i][NbOfBJetsBins_][2]  = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
	
		// MC prediction/Estimation summary (b-inclusive)
		sprintf(name ,"hNjetsSummary_bInclusive_wp_%d",i);
		sprintf(title,"MC prediction/Estimation summary (b-inclusive, WP nr%d)",i);
       		tCanva_Njets_Summary[i][NbOfBJetsBins_] = new TCanvas(name,title,600,600);
	}
	cout<<" -- Summary histograms correctly instantiated"<<endl;

	//init efficiency estimators
        eudsc_        = new double*[NbOfBtagWorkingPoint_];
        eudsc_err_    = new double*[NbOfBtagWorkingPoint_];
        eudsc_mc_     = new double*[NbOfBtagWorkingPoint_];
        eudsc_err_mc_ = new double*[NbOfBtagWorkingPoint_];
        euds_         = new double*[NbOfBtagWorkingPoint_];
        euds_err_     = new double*[NbOfBtagWorkingPoint_];
        euds_mc_      = new double*[NbOfBtagWorkingPoint_];
        euds_err_mc_  = new double*[NbOfBtagWorkingPoint_];
        eb_           = new double*[NbOfBtagWorkingPoint_];
        eb_err_       = new double*[NbOfBtagWorkingPoint_];
        eb_mc_        = new double*[NbOfBtagWorkingPoint_];
        eb_err_mc_    = new double*[NbOfBtagWorkingPoint_];
	//init containers
	N_            = new double***[NbOfBtagWorkingPoint_];
	N_err_        = new double***[NbOfBtagWorkingPoint_];
	Nbjets_       = new double** [NbOfBtagWorkingPoint_];
	MultiJet_Est_ = new double** [NbOfBtagWorkingPoint_];
	//init error estimators
	NvVar_        = new double**[NbOfBtagWorkingPoint_];
	NvBias_       = new double**[NbOfBtagWorkingPoint_];
	NvbVar_       = new double**[NbOfBtagWorkingPoint_];
	NvbBias_      = new double**[NbOfBtagWorkingPoint_];
	NttVar_       = new double**[NbOfBtagWorkingPoint_];
	NttBias_      = new double**[NbOfBtagWorkingPoint_];

	for(int i=0;i<NbOfBtagWorkingPoint_;i++){

	        eudsc_[i]        = new double[NbOfJetsBins_];
	        eudsc_err_[i]    = new double[NbOfJetsBins_];
	        eudsc_mc_[i]     = new double[NbOfJetsBins_];
	        eudsc_err_mc_[i] = new double[NbOfJetsBins_];
	        euds_[i]         = new double[NbOfJetsBins_];
	        euds_err_[i]     = new double[NbOfJetsBins_];
	        euds_mc_[i]      = new double[NbOfJetsBins_];
	        euds_err_mc_[i]  = new double[NbOfJetsBins_];
	        eb_[i]           = new double[NbOfJetsBins_];
	        eb_err_[i]       = new double[NbOfJetsBins_];
	        eb_mc_[i]        = new double[NbOfJetsBins_];
	        eb_err_mc_[i]    = new double[NbOfJetsBins_];

		N_[i]            = new double**[NbOfJetsBins_];
		N_err_[i]        = new double**[NbOfJetsBins_];
		Nbjets_[i]       = new double* [NbOfJetsBins_];
		MultiJet_Est_[i] = new double* [NbOfJetsBins_];

		for(int j=0;j<NbOfJetsBins_;j++){

			eudsc_[i][j]       = 0;
			eudsc_err_[i][j]   = 0;
			eudsc_mc_[i][j]    = 0;
			eudsc_err_mc_[i][j]= 0;
			euds_[i][j]        = 0;
			euds_err_[i][j]    = 0;
			euds_mc_[i][j]     = 0;
			euds_err_mc_[i][j] = 0;
			eb_[i][j]          = 0;
			eb_err_[i][j]      = 0;
			eb_mc_[i][j]       = 0;
			eb_err_mc_[i][j]   = 0;

			N_[i][j]            = new double*[NbOfBJetsBins_];
			N_err_[i][j]        = new double*[NbOfBJetsBins_];
			Nbjets_[i][j]       = new double [NbOfBJetsBins_];
			MultiJet_Est_[i][j] = new double [NbOfBJetsBins_];

			for(int k=0;k<NbOfBJetsBins_;k++){
				N_[i][j][k]     = new double[NbOfDatasets_];
				N_err_[i][j][k] = new double[NbOfDatasets_];
				Nbjets_[i][j][k]       = 0;
				MultiJet_Est_[i][j][k] = 0;

				for(int l=0;l<NbOfDatasets_;l++){
					N_[i][j][k][l]     = 0;
					N_err_[i][j][k][l] = 0;
				}
			}
		}

		NvVar_[i]   = new double*[NbOfJetsBins_+1];
		NvBias_[i]  = new double*[NbOfJetsBins_+1];
		NvbVar_[i]  = new double*[NbOfJetsBins_+1];
		NvbBias_[i] = new double*[NbOfJetsBins_+1];
		NttVar_[i]  = new double*[NbOfJetsBins_+1];
		NttBias_[i] = new double*[NbOfJetsBins_+1];
		for(int j=0;j<NbOfJetsBins_+1;j++){
			NvVar_[i][j]        = new double[NbOfBJetsBins_+1];
			NvBias_[i][j]       = new double[NbOfBJetsBins_+1];
			NvbVar_[i][j]       = new double[NbOfBJetsBins_+1];
			NvbBias_[i][j]      = new double[NbOfBJetsBins_+1];
			NttVar_[i][j]       = new double[NbOfBJetsBins_+1];
			NttBias_[i][j]      = new double[NbOfBJetsBins_+1];
			for(int k=0;k<NbOfBJetsBins_+1;k++){
				NvVar_[i][j][k]   = 0;
				NvBias_[i][j][k]  = 0;
				NvbVar_[i][j][k]  = 0;
				NvbBias_[i][j][k] = 0;
				NttVar_[i][j][k]  = 0;
				NttBias_[i][j][k] = 0;
			}
		}
	}

	//init histograms used to calculate the b/mis-tagging efficiencies.
	hNbOfBGenJets_                                 = new TH1F**[NbOfDatasets_];
	hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_     = new TH3F* [NbOfDatasets_];
	hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_     = new TH3F* [NbOfDatasets_];
	hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_   = new TH3F* [NbOfDatasets_];
	hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_  = new TH3F* [NbOfDatasets_];

	hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_    = new TH3F**[NbOfDatasets_];
	hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_    = new TH3F**[NbOfDatasets_];
	hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_  = new TH3F**[NbOfDatasets_];
	hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_ = new TH3F**[NbOfDatasets_];

	bTagEff_vs_Njets_   = new TGraphAsymmErrors**[NbOfDatasets_];
	cTagEff_vs_Njets_   = new TGraphAsymmErrors**[NbOfDatasets_];
	udsTagEff_vs_Njets_ = new TGraphAsymmErrors**[NbOfDatasets_];
	misTagEff_vs_Njets_ = new TGraphAsymmErrors**[NbOfDatasets_];

	const int JetPtNbins = 4;             double JetPtBins[JetPtNbins+1]   = {30,50,80,120,250};
	const int JetEtaNbins= 4;             double JetEtaBins[JetEtaNbins+1] = {0,0.4,0.8,1.3,2.4};
	const int JetNbins   = NbOfJetsBins_; double JetNBins[JetNbins+1]; for(int i=0;i<JetNbins+1;i++) JetNBins[i] = Njets_[0]+i;

	for(int i=0;i<NbOfDatasets_;i++){

		hNbOfBGenJets_[i] = new TH1F*[NbOfJetsBins_];
		for(int j=0; j<NbOfJetsBins_;j++){
			sprintf(name ,"hNbOfBGenJets_%d_%d",i,j);
			hNbOfBGenJets_[i][j] = new TH1F(name,"",4,0,4);
			hNbOfBGenJets_[i][j] ->Sumw2();
			hNbOfBGenJets_[i][j] ->GetXaxis()->CenterLabels();
			hNbOfBGenJets_[i][j] ->GetXaxis()->SetTitle("Nb. of b-quark jets");
		}
		sprintf(name ,"hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->SetTitle("Nb. of jets");
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->CenterLabels();

		sprintf(name ,"hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->SetTitle("Nb. of jets");
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->CenterLabels();

		sprintf(name ,"hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->SetTitle("Nb. of jets");
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->CenterLabels();

		sprintf(name ,"hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->SetTitle("Nb. of jets");
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->CenterLabels();

		bTagEff_vs_Njets_[i]   = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];
		cTagEff_vs_Njets_[i]   = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];
		udsTagEff_vs_Njets_[i] = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];
		misTagEff_vs_Njets_[i] = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];

		hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i]    = new TH3F*[NbOfBtagWorkingPoint_];
		hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i]    = new TH3F*[NbOfBtagWorkingPoint_];
		hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i]  = new TH3F*[NbOfBtagWorkingPoint_];
		hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = new TH3F*[NbOfBtagWorkingPoint_];

		for(int j=0; j<NbOfBtagWorkingPoint_;j++){
			sprintf(name ,"hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->SetTitle("Nb. of jets");
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->CenterLabels();

			sprintf(name ,"hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->SetTitle("Nb. of jets");
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->CenterLabels();

			sprintf(name ,"hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->SetTitle("Nb. of jets");
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->CenterLabels();

			sprintf(name ,"hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->SetTitle("Nb. of jets");
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->CenterLabels();

			bTagEff_vs_Njets_[i][j]   = 0;
			cTagEff_vs_Njets_[i][j]   = 0;
			udsTagEff_vs_Njets_[i][j] = 0;
			misTagEff_vs_Njets_[i][j] = 0;
		}
	}

	//init efficiency estimators
        ebq_   = new double[NbOfJetsBins_];
        e0bq_  = new double[NbOfJetsBins_];
        e1bq_  = new double[NbOfJetsBins_];
        e2bq_  = new double[NbOfJetsBins_];
	//init Ntt/Nvb/Nv estimators
        Ntt_          = new double[NbOfJetsBins_];
        Ntt_err_      = new double[NbOfJetsBins_];
        Ntt_err_up_   = new double[NbOfJetsBins_];
        Ntt_err_down_ = new double[NbOfJetsBins_];
        Nv_           = new double[NbOfJetsBins_];
        Nv_err_       = new double[NbOfJetsBins_];
        Nv_err_up_    = new double[NbOfJetsBins_];
        Nv_err_down_  = new double[NbOfJetsBins_];
        Nvb_          = new double[NbOfJetsBins_];
        Nvb_err_      = new double[NbOfJetsBins_];
	for(int i=0;i<NbOfJetsBins_;i++){
		ebq_[i]      = 0;
		e0bq_[i]     = 0;
		e1bq_[i]     = 0;
		e2bq_[i]     = 0;
		Ntt_[i]          = 0;
		Ntt_err_[i]      = 0;
		Ntt_err_up_[i]   = 0;
		Ntt_err_down_[i] = 0;
		Nv_[i]           = 0;
		Nv_err_[i]       = 0;
		Nv_err_up_[i]    = 0;
		Nv_err_down_[i]  = 0;
		Nvb_[i]          = 0;
		Nvb_err_[i]      = 0;
	}

	condProb_ = new double****[NbOfBtagWorkingPoint_-1];
	for(int i=0;i<NbOfBtagWorkingPoint_-1;i++){
		condProb_[i] = new double***[NbOfJetsBins_];
		for(int j=0;j<NbOfJetsBins_;j++){
			condProb_[i][j] = new double**[NbOfDatasets_];
			for(int k=0;k<NbOfDatasets_;k++){
				condProb_[i][j][k] = new double*[NbOfBJetsBins_];
				for(int l=0;l<NbOfBJetsBins_;l++){
					condProb_[i][j][k][l] = new double[NbOfBJetsBins_];
					for(int m=0;m<NbOfBJetsBins_;m++) condProb_[i][j][k][l][m]=0;
				}
			}
		}
	}

	cout<<" -- Stability check histograms correctly instantiated"<<endl;
	//init histo for the stability check
	hScanMin_Eb    = new TGraph**[NbOfBtagWorkingPoint_];
	hScanMin_Eudsc = new TGraph**[NbOfBtagWorkingPoint_];
	hScanMin_Euds  = new TGraph**[NbOfBtagWorkingPoint_];
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		hScanMin_Eb[i]    = new TGraph*[NbOfJetsBins_];
		hScanMin_Eudsc[i] = new TGraph*[NbOfJetsBins_];
		hScanMin_Euds[i]  = new TGraph*[NbOfJetsBins_];
			for(int j=0;j<NbOfJetsBins_;j++){
			hScanMin_Eb[i][j]    = 0;
			hScanMin_Eudsc[i][j] = 0;
			hScanMin_Euds[i][j]  = 0;
		}
	}
	hScanMin_Ntt = new TGraph*[NbOfJetsBins_];
	hScanMin_Nv  = new TGraph*[NbOfJetsBins_];
	hScanMin_Nvb = new TGraph*[NbOfJetsBins_];
	for(int i=0;i<NbOfJetsBins_;i++){
		hScanMin_Ntt[i] = 0;
		hScanMin_Nv[i]  = 0;
		hScanMin_Nvb[i] = 0;
	}
	cout<<"Object from the class VJetEstimation correctly instantiated"<<endl;
}

/**________________________________________________________________________________________________________________*/
VJetEstimation::VJetEstimation(int NofBtagWorkingPoint, float* BtagWorkingPoint, int NofJets, int NofJetBins, double** EffXbq, int NofDatasets, vector<int> iDTTLike, vector<int> iDVLike, vector<int> iDVbLike){
	cout<<"Object from the class VJetEstimation being instantiated"<<endl;
	NbOfDatasets_ = NofDatasets;
	NbOfBtagWorkingPoint_ = NofBtagWorkingPoint;
	BtagWorkingPoint_ = new float[NbOfBtagWorkingPoint_];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++) BtagWorkingPoint_[i] = BtagWorkingPoint[i];
	MCdata_ = true;
	iDatasetsTTLike_.clear();
	iDatasetsVLike_.clear();
	iDatasetsTTLike_ = iDTTLike;
	iDatasetsVLike_  = iDVLike;
	iDatasetsVbLike_ = iDVbLike;
	
	NbOfJetsBins_ = NofJetBins;
	Njets_ = new int[NbOfJetsBins_+1];
	for(int i=0;i<NbOfJetsBins_+1;i++) Njets_[i] = NofJets+i;
	NbOfBJetsBins_ = 4;

	init_ = 0;//
	minValue_ = new double[NbOfJetsBins_];
	for(int i=0;i<NbOfJetsBins_;i++) minValue_[i] = 0;
	
	RescaledTTLikeEstimation = 0;
	RescaledVLikeEstimation = 0;
	RescaledVbLikeEstimation = 0;
	tCanva_RescaledTTLikeEstimation = 0;
	tCanva_RescaledVLikeEstimation = 0;
	tCanva_RescaledVbLikeEstimation = 0;
        
	//init summary histos
	hsNbjets_MC  = new THStack**[NbOfBtagWorkingPoint_];
	hsNbjets_Est = new THStack**[NbOfBtagWorkingPoint_];
	hsNjets_MC   = new THStack**[NbOfBtagWorkingPoint_];
	hsNjets_Est  = new THStack**[NbOfBtagWorkingPoint_];
	MyLeg = new TLegend(0.6,0.6,0.9,0.9);
	
	hNbjets_mc_       = new TH1F***[NbOfBtagWorkingPoint_];
	hNbjets_pdf_mc_   = new TH1F***[NbOfBtagWorkingPoint_];
	hNbjets_pdf_est_  = new TH1F***[NbOfBtagWorkingPoint_];
	
	hNbjetsEstSummary = new TH1F***[NbOfBtagWorkingPoint_];
	hNbjetsMCSummary  = new TH1F***[NbOfBtagWorkingPoint_];
	hNjetsEstSummary  = new TH1F***[NbOfBtagWorkingPoint_];
	hNjetsMCSummary   = new TH1F***[NbOfBtagWorkingPoint_];
	
	tCanva_Nbjets_Summary = new TCanvas**[NbOfBtagWorkingPoint_];
	tCanva_Njets_Summary  = new TCanvas**[NbOfBtagWorkingPoint_];
	
	char name[100];
	char title[100];

	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		hsNbjets_MC[i]  = new THStack*[NbOfJetsBins_+1];
		hsNbjets_Est[i] = new THStack*[NbOfJetsBins_+1];
		hsNjets_MC[i]   = new THStack*[NbOfBJetsBins_+1];
		hsNjets_Est[i]  = new THStack*[NbOfBJetsBins_+1];
	
		hNbjets_mc_[i]       = new TH1F**[NbOfJetsBins_];
		hNbjets_pdf_mc_[i]   = new TH1F**[NbOfJetsBins_];
		hNbjets_pdf_est_[i]  = new TH1F**[NbOfJetsBins_];
		
		hNbjetsEstSummary[i] = new TH1F**[NbOfJetsBins_+1];
		hNbjetsMCSummary[i]  = new TH1F**[NbOfJetsBins_+1];
		hNjetsEstSummary[i]  = new TH1F**[NbOfBJetsBins_+1];
		hNjetsMCSummary[i]   = new TH1F**[NbOfBJetsBins_+1];
		
		tCanva_Nbjets_Summary[i] = new TCanvas*[NbOfJetsBins_+1];
		tCanva_Njets_Summary[i]  = new TCanvas*[NbOfBJetsBins_+1];

		for(int j=0;j<NbOfJetsBins_;j++){
		// Booking histograms for VLike,VbLike and TTlike predictions/estimations
		// for different jet multiplicity

			// estimation
			hsNbjets_MC[i][j]          = 0;
			hsNbjets_Est[i][j]         = 0;

			hNbjets_mc_[i][j]          = new TH1F*[3];
			sprintf(name,"hNbjets_mc_VLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_mc_[i][j][0]       = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_mc_[i][j][0]->Sumw2();
			sprintf(name,"hNbjets_mc_VbLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_mc_[i][j][1]       = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_mc_[i][j][1]->Sumw2();
			sprintf(name,"hNbjets_mc_TTLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_mc_[i][j][2]       = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_mc_[i][j][2]->Sumw2();

			hNbjets_pdf_mc_[i][j]      = new TH1F*[3];
			sprintf(name,"hNbjets_pdf_mc_VLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_mc_[i][j][0]   = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_mc_[i][j][0]->Sumw2();
			sprintf(name,"hNbjets_pdf_mc_VbLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_mc_[i][j][1]   = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_mc_[i][j][1]->Sumw2();
			sprintf(name,"hNbjets_pdf_mc_TTLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_mc_[i][j][2]   = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_mc_[i][j][2]->Sumw2();

			hNbjets_pdf_est_[i][j]     = new TH1F*[3];
			sprintf(name,"hNbjets_pdf_est_VLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_est_[i][j][0]  = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_est_[i][j][0]->Sumw2();
			sprintf(name,"hNbjets_pdf_est_VbLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_est_[i][j][1]  = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_est_[i][j][1]->Sumw2();
			sprintf(name,"hNbjets_pdf_est_TTLike_%d_wp_%d_jets",i,Njets_[j]);
			hNbjets_pdf_est_[i][j][2]  = new TH1F(name,"",NbOfBJetsBins_,0,NbOfBJetsBins_);hNbjets_pdf_est_[i][j][2]->Sumw2();

			hNbjetsEstSummary[i][j]    = new TH1F*[3];
			sprintf(name,"hNbjetsEstSummary_VLike_%d_jets_wp_%d",Njets_[j],i);
			sprintf(title,"W-like events estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsEstSummary[i][j][0] = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsEstSummary_VbLike_%d_jets_wp_%d",Njets_[j],i);
			sprintf(title,"Wb-like events estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsEstSummary[i][j][1] = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsEstSummary_TTlike_%d_jets_wp_%d",Njets_[j],i);
			sprintf(title,"TT-like events estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsEstSummary[i][j][2] = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);

			// MC prediction
			hNbjetsMCSummary[i][j]     = new TH1F*[3];
			sprintf(name,"hNbjetsMCSummary_VLike_%d_jets_wp_%d",Njets_[j],i);
			sprintf(title,"W-like events MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsMCSummary[i][j][0]  = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsMCSummary_VbLike_%d_jets_wp_%d",Njets_[j],i);
			sprintf(title,"Wb-like events MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsMCSummary[i][j][1]  = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
			sprintf(name,"hNbjetsMCSummary_TTlike_%d_jets_wp_%d",Njets_[j],i);
			sprintf(title,"TT-like events MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
			hNbjetsMCSummary[i][j][2]  = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);

			// MC prediction vs estimation
			sprintf(name,"hNbjetsSummary_%d_jets_wp_%d",Njets_[j],i);
			sprintf(title,"MC prediction/Estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			tCanva_Nbjets_Summary[i][j] = new TCanvas(name,title,600,600);
		}

		hsNbjets_MC[i][NbOfJetsBins_]          = 0;
		hsNbjets_Est[i][NbOfJetsBins_]         = 0;
		// Estimation summary (inclusive)
		hNbjetsEstSummary[i][NbOfJetsBins_]    = new TH1F*[3];
		sprintf(name ,"hNbjetsEstSummary_VLike_Inclusive_wp_%d",i);
		sprintf(title,"W-like events estimation summary (inclusive, WP nr%d)",i);
		hNbjetsEstSummary[i][NbOfJetsBins_][0] = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsEstSummary_VbLike_Inclusive_wp_%d",i);
		sprintf(title,"Wb-like events estimation summary (inclusive, WP nr%d)",i);
		hNbjetsEstSummary[i][NbOfJetsBins_][1] = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsEstSummary_TTlike_Inclusive_wp_%d",i);
		sprintf(title,"TT-like events estimation summary (inclusive, WP nr%d)",i);
		hNbjetsEstSummary[i][NbOfJetsBins_][2] = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);

		// MC prediction summary (inclusive)
		hNbjetsMCSummary[i][NbOfJetsBins_]     = new TH1F*[3];
		sprintf(name ,"hNbjetsMCSummary_VLike_Inclusive_wp_%d",i);
		sprintf(title,"W-like events MC prediction summary (inclusive, WP nr%d)",i);
		hNbjetsMCSummary[i][NbOfJetsBins_][0]  = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsMCSummary_VbLike_Inclusive_wp_%d",i);
		sprintf(title,"Wb-like events MC prediction summary (inclusive, WP nr%d)",i);
		hNbjetsMCSummary[i][NbOfJetsBins_][1]  = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);
		sprintf(name ,"hNbjetsMCSummary_TTlike_Inclusive_wp_%d",i);
		sprintf(title,"TT-like events MC prediction summary (inclusive, WP nr%d)",i);
		hNbjetsMCSummary[i][NbOfJetsBins_][2]  = new TH1F(name,title,NbOfBJetsBins_,0,NbOfBJetsBins_);

		// MC prediction/EStimation summary (inclusive)
		sprintf(name ,"hNbjetsSummary_Inclusive_wp_%d",i);
		sprintf(title,"MC prediction/Estimation summary (inclusive, WP nr%d)",i);
		tCanva_Nbjets_Summary[i][NbOfJetsBins_] = new TCanvas(name,title,600,600);

		for(int j=0;j<NbOfBJetsBins_;j++){
		// Booking histograms for VLike,VbLike and TTlike predictions/estimations
		// for different b-jet multiplicity

			// estimation
			hsNjets_MC[i][j]          = 0;
			hsNjets_Est[i][j]         = 0;
			hNjetsEstSummary[i][j]    = new TH1F*[3];
			sprintf(name,"hNjetsEstSummary_VLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"W-like events estimation summary for %d b-jets (WP nr%d)",j,i);
			hNjetsEstSummary[i][j][0] = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsEstSummary_VbLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"Wb-like events estimation summary for %d b-jets (WP nr%d)",j,i);
			hNjetsEstSummary[i][j][1] = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsEstSummary_TTlike_%d_bjets_wp_%d",j,i);
			sprintf(title,"TT-like events estimation summary for %d b-jets (WP nr%d)",j,i);
			hNjetsEstSummary[i][j][2] = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);

			// MC prediction
			hNjetsMCSummary[i][j]     = new TH1F*[3];
			sprintf(name,"hNjetsMCSummary_VLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"W-like events MC prediction summary for %d b-jets (WP nr%d)",j,i);
			hNjetsMCSummary[i][j][0]  = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsMCSummary_VbLike_%d_bjets_wp_%d",j,i);
			sprintf(title,"Wb-like events MC prediction summary for %d b-jets (WP nr%d)",j,i);
			hNjetsMCSummary[i][j][1]  = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
			sprintf(name,"hNjetsMCSummary_TTlike_%d_bjets_wp_%d",j,i);
			sprintf(title,"TT-like events MC prediction summary for %d b-jets (WP nr%d)",j,i);
			hNjetsMCSummary[i][j][2]  = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);

			// MC prediction vs estimation
			sprintf(name,"hNjetsSummary_%d_bjets_wp_%d",j,i);
			sprintf(title,"MC prediction/Estimation summary for %d b-jets (WP nr%d)",j,i);		
        		tCanva_Njets_Summary[i][j] = new TCanvas(name,title,600,600);
		}

		// Estimation summary (b-inclusive)
		hsNjets_MC[i][NbOfBJetsBins_]          = 0;
		hsNjets_Est[i][NbOfBJetsBins_]         = 0;
		hNjetsEstSummary[i][NbOfBJetsBins_]    = new TH1F*[3];
		sprintf(name ,"hNjetsEstSummary_VLike_bInclusive_wp_%d",i);
		sprintf(title,"W-like events estimation summary (b-inclusive, WP nr%d)",i);
		hNjetsEstSummary[i][NbOfBJetsBins_][0] = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsEstSummary_VbLike_bInclusive_wp_%d",i);
		sprintf(title,"Wb-like events estimation summary (b-inclusive, WP nr%d)",i);
		hNjetsEstSummary[i][NbOfBJetsBins_][1] = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsEstSummary_TTlike_bInclusive_wp_%d",i);
		sprintf(title,"TT-like events estimation summary (b-inclusive, WP nr%d)",i);
		hNjetsEstSummary[i][NbOfBJetsBins_][2] = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);

		// MC prediction summary (b-inclusive)
		hNjetsMCSummary[i][NbOfBJetsBins_]     = new TH1F*[3];
		sprintf(name ,"hNjetsMCSummary_VLike_bInclusive_wp_%d",i);
		sprintf(title,"W-like events MC prediction summary (b-inclusive, WP nr%d)",i);
		hNjetsMCSummary[i][NbOfBJetsBins_][0]  = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsMCSummary_VbLike_bInclusive_wp_%d",i);
		sprintf(title,"Wb-like events MC prediction summary (b-inclusive, WP nr%d)",i);
		hNjetsMCSummary[i][NbOfBJetsBins_][1]  = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
		sprintf(name ,"hNjetsMCSummary_TTlike_bInclusive_wp_%d",i);
		sprintf(title,"TT-like events MC prediction summary (b-inclusive, WP nr%d)",i);
		hNjetsMCSummary[i][NbOfBJetsBins_][2]  = new TH1F(name,title,NbOfJetsBins_,0,NbOfJetsBins_);
	
		// MC prediction/Estimation summary (b-inclusive)
		sprintf(name ,"hNjetsSummary_bInclusive_wp_%d",i);
		sprintf(title,"MC prediction/Estimation summary (b-inclusive, WP nr%d)",i);
       		tCanva_Njets_Summary[i][NbOfBJetsBins_] = new TCanvas(name,title,600,600);
	}
	cout<<" -- Summary histograms correctly instantiated"<<endl;

	//init efficiency estimators
        eudsc_        = new double*[NbOfBtagWorkingPoint_];
        eudsc_err_    = new double*[NbOfBtagWorkingPoint_];
        eudsc_mc_     = new double*[NbOfBtagWorkingPoint_];
        eudsc_err_mc_ = new double*[NbOfBtagWorkingPoint_];
        euds_         = new double*[NbOfBtagWorkingPoint_];
        euds_err_     = new double*[NbOfBtagWorkingPoint_];
        euds_mc_      = new double*[NbOfBtagWorkingPoint_];
        euds_err_mc_  = new double*[NbOfBtagWorkingPoint_];
        eb_           = new double*[NbOfBtagWorkingPoint_];
        eb_err_       = new double*[NbOfBtagWorkingPoint_];
        eb_mc_        = new double*[NbOfBtagWorkingPoint_];
        eb_err_mc_    = new double*[NbOfBtagWorkingPoint_];
	//init containers
	N_            = new double***[NbOfBtagWorkingPoint_];
	N_err_        = new double***[NbOfBtagWorkingPoint_];
	Nbjets_       = new double** [NbOfBtagWorkingPoint_];
	MultiJet_Est_ = new double** [NbOfBtagWorkingPoint_];
	//init error estimators
	NvVar_        = new double**[NbOfBtagWorkingPoint_];
	NvBias_       = new double**[NbOfBtagWorkingPoint_];
	NvbVar_       = new double**[NbOfBtagWorkingPoint_];
	NvbBias_      = new double**[NbOfBtagWorkingPoint_];
	NttVar_       = new double**[NbOfBtagWorkingPoint_];
	NttBias_      = new double**[NbOfBtagWorkingPoint_];

	for(int i=0;i<NbOfBtagWorkingPoint_;i++){

	        eudsc_[i]        = new double[NbOfJetsBins_];
	        eudsc_err_[i]    = new double[NbOfJetsBins_];
	        eudsc_mc_[i]     = new double[NbOfJetsBins_];
	        eudsc_err_mc_[i] = new double[NbOfJetsBins_];
	        euds_[i]         = new double[NbOfJetsBins_];
	        euds_err_[i]     = new double[NbOfJetsBins_];
	        euds_mc_[i]      = new double[NbOfJetsBins_];
	        euds_err_mc_[i]  = new double[NbOfJetsBins_];
	        eb_[i]           = new double[NbOfJetsBins_];
	        eb_err_[i]       = new double[NbOfJetsBins_];
	        eb_mc_[i]        = new double[NbOfJetsBins_];
	        eb_err_mc_[i]    = new double[NbOfJetsBins_];

		N_[i]            = new double**[NbOfJetsBins_];
		N_err_[i]        = new double**[NbOfJetsBins_];
		Nbjets_[i]       = new double* [NbOfJetsBins_];
		MultiJet_Est_[i] = new double* [NbOfJetsBins_];
		for(int j=0;j<NbOfJetsBins_;j++){

			eudsc_[i][j]        = 0;
			eudsc_err_[i][j]    = 0;
			eudsc_mc_[i][j]     = 0;
			eudsc_err_mc_[i][j] = 0;
			euds_[i][j]         = 0;
			euds_err_[i][j]     = 0;
			euds_mc_[i][j]      = 0;
			euds_err_mc_[i][j]  = 0;
			eb_[i][j]           = 0;
			eb_err_[i][j]       = 0;
			eb_mc_[i][j]        = 0;
			eb_err_mc_[i][j]    = 0;

			N_[i][j]            = new double*[NbOfBJetsBins_];
			N_err_[i][j]        = new double*[NbOfBJetsBins_];
			Nbjets_[i][j]       = new double [NbOfBJetsBins_];
			MultiJet_Est_[i][j] = new double [NbOfBJetsBins_];

			for(int k=0;k<NbOfBJetsBins_;k++){
				N_[i][j][k]     = new double[NbOfDatasets_];
				N_err_[i][j][k] = new double[NbOfDatasets_];
				Nbjets_[i][j][k]       = 0;
				MultiJet_Est_[i][j][k] = 0;

				for(int l=0;l<NbOfDatasets_;l++){
					N_[i][j][k][l] = 0;
					N_err_[i][j][k][l] = 0;
				}
			}
		}

		NvVar_[i]   = new double*[NbOfJetsBins_+1];
		NvBias_[i]  = new double*[NbOfJetsBins_+1];
		NvbVar_[i]  = new double*[NbOfJetsBins_+1];
		NvbBias_[i] = new double*[NbOfJetsBins_+1];
		NttVar_[i]  = new double*[NbOfJetsBins_+1];
		NttBias_[i] = new double*[NbOfJetsBins_+1];
		for(int j=0;j<NbOfJetsBins_+1;j++){
			NvVar_[i][j]        = new double[NbOfBJetsBins_+1];
			NvBias_[i][j]       = new double[NbOfBJetsBins_+1];
			NvbVar_[i][j]       = new double[NbOfBJetsBins_+1];
			NvbBias_[i][j]      = new double[NbOfBJetsBins_+1];
			NttVar_[i][j]       = new double[NbOfBJetsBins_+1];
			NttBias_[i][j]      = new double[NbOfBJetsBins_+1];
			for(int k=0;k<NbOfBJetsBins_+1;k++){
				NvVar_[i][j][k]   = 0;
				NvBias_[i][j][k]  = 0;
				NvbVar_[i][j][k]  = 0;
				NvbBias_[i][j][k] = 0;
				NttVar_[i][j][k]  = 0;
				NttBias_[i][j][k] = 0;
			}
		}
	}

	//init histograms used to calculate the b/mis-tagging efficiencies.
	hNbOfBGenJets_                                 = new TH1F**[NbOfDatasets_];
	hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_     = new TH3F* [NbOfDatasets_];
	hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_     = new TH3F* [NbOfDatasets_];
	hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_   = new TH3F* [NbOfDatasets_];
	hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_  = new TH3F* [NbOfDatasets_];

	hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_    = new TH3F**[NbOfDatasets_];
	hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_    = new TH3F**[NbOfDatasets_];
	hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_  = new TH3F**[NbOfDatasets_];
	hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_ = new TH3F**[NbOfDatasets_];

	bTagEff_vs_Njets_   = new TGraphAsymmErrors**[NbOfDatasets_];
	cTagEff_vs_Njets_   = new TGraphAsymmErrors**[NbOfDatasets_];
	udsTagEff_vs_Njets_ = new TGraphAsymmErrors**[NbOfDatasets_];
	misTagEff_vs_Njets_ = new TGraphAsymmErrors**[NbOfDatasets_];

	const int JetPtNbins = 4;             double JetPtBins[JetPtNbins+1]   = {30,50,80,120,250};
	const int JetEtaNbins= 4;             double JetEtaBins[JetEtaNbins+1] = {0,0.4,0.8,1.3,2.4};
	const int JetNbins   = NbOfJetsBins_; double JetNBins[JetNbins+1]; for(int i=0;i<JetNbins+1;i++) JetNBins[i] = Njets_[0]+i;

	for(int i=0;i<NbOfDatasets_;i++){

		hNbOfBGenJets_[i] = new TH1F*[NbOfJetsBins_];
		for(int j=0; j<NbOfJetsBins_;j++){
			sprintf(name ,"hNbOfBGenJets_%d_%d",i,j);
			hNbOfBGenJets_[i][j] = new TH1F(name,"",4,0,4);
			hNbOfBGenJets_[i][j] ->Sumw2();
			hNbOfBGenJets_[i][j] ->GetXaxis()->CenterLabels();
			hNbOfBGenJets_[i][j] ->GetXaxis()->SetTitle("Nb. of b-quark jets");
		}
		sprintf(name ,"hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->SetTitle("Nb. of jets");
		hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->CenterLabels();

		sprintf(name ,"hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->SetTitle("Nb. of jets");
		hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->CenterLabels();

		sprintf(name ,"hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->SetTitle("Nb. of jets");
		hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->CenterLabels();

		sprintf(name ,"hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_%d",i);
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetYaxis()->SetTitle("Jet |#eta|");
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->SetTitle("Nb. of jets");
		hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i]->GetZaxis()->CenterLabels();

		bTagEff_vs_Njets_[i]   = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];
		cTagEff_vs_Njets_[i]   = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];
		udsTagEff_vs_Njets_[i] = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];
		misTagEff_vs_Njets_[i] = new TGraphAsymmErrors*[NbOfBtagWorkingPoint_];

		hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i]    = new TH3F*[NbOfBtagWorkingPoint_];
		hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i]    = new TH3F*[NbOfBtagWorkingPoint_];
		hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i]  = new TH3F*[NbOfBtagWorkingPoint_];
		hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i] = new TH3F*[NbOfBtagWorkingPoint_];

		for(int j=0; j<NbOfBtagWorkingPoint_;j++){
			sprintf(name ,"hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->SetTitle("Nb. of jets");
			hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->CenterLabels();

			sprintf(name ,"hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->SetTitle("Nb. of jets");
			hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->CenterLabels();

			sprintf(name ,"hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->SetTitle("Nb. of jets");
			hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->CenterLabels();

			sprintf(name ,"hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_%d_wp_%d",i,j);
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j] = new TH3F(name,"",JetPtNbins,JetPtBins,JetEtaNbins,JetEtaBins,JetNbins,JetNBins);
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetYaxis()->SetTitle("Jet |#eta|");
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->SetTitle("Nb. of jets");
			hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j]->GetZaxis()->CenterLabels();

			bTagEff_vs_Njets_[i][j]   = 0;
			cTagEff_vs_Njets_[i][j]   = 0;
			udsTagEff_vs_Njets_[i][j] = 0;
			misTagEff_vs_Njets_[i][j] = 0;
		}
	}

	//init efficiency estimators
        ebq_   = new double[NbOfJetsBins_];
        e0bq_  = new double[NbOfJetsBins_];//EffXbq[0];
        e1bq_  = new double[NbOfJetsBins_];//EffXbq[1];
        e2bq_  = new double[NbOfJetsBins_];//EffXbq[2];
	//init Ntt/Nvb/Nv estimators
        Ntt_          = new double[NbOfJetsBins_];
        Ntt_err_      = new double[NbOfJetsBins_];
        Ntt_err_up_   = new double[NbOfJetsBins_];
        Ntt_err_down_ = new double[NbOfJetsBins_];
        Nv_           = new double[NbOfJetsBins_];
        Nv_err_       = new double[NbOfJetsBins_];
        Nv_err_up_    = new double[NbOfJetsBins_];
        Nv_err_down_  = new double[NbOfJetsBins_];
        Nvb_          = new double[NbOfJetsBins_];
        Nvb_err_      = new double[NbOfJetsBins_];
	for(int i=0;i<NbOfJetsBins_;i++){
		ebq_[i]      = 0;
	        e0bq_[i]     = EffXbq[i][0];
	        e1bq_[i]     = EffXbq[i][1];
	        e2bq_[i]     = EffXbq[i][2];
		Ntt_[i]          = 0;
		Ntt_err_[i]      = 0;
		Ntt_err_up_[i]   = 0;
		Ntt_err_down_[i] = 0;
		Nv_[i]           = 0;
		Nv_err_[i]       = 0;
		Nv_err_up_[i]    = 0;
		Nv_err_down_[i]  = 0;
		Nvb_[i]          = 0;
		Nvb_err_[i]      = 0;
	}

	condProb_ = new double****[NbOfBtagWorkingPoint_-1];
	for(int i=0;i<NbOfBtagWorkingPoint_-1;i++){
		condProb_[i] = new double***[NbOfJetsBins_];
		for(int j=0;j<NbOfJetsBins_;j++){
			condProb_[i][j] = new double**[NbOfDatasets_];
			for(int k=0;k<NbOfDatasets_;k++){
				condProb_[i][j][k] = new double*[NbOfBJetsBins_];
				for(int l=0;l<NbOfBJetsBins_;l++){
					condProb_[i][j][k][l] = new double[NbOfBJetsBins_];
					for(int m=0;m<NbOfBJetsBins_;m++) condProb_[i][j][k][l][m]=0;
				}
			}
		}
	}

	cout<<" -- Stability check histograms correctly instantiated"<<endl;
	//init histo for the stability check
	hScanMin_Eb    = new TGraph**[NbOfBtagWorkingPoint_];
	hScanMin_Eudsc = new TGraph**[NbOfBtagWorkingPoint_];
	hScanMin_Euds  = new TGraph**[NbOfBtagWorkingPoint_];
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		hScanMin_Eb[i]    = new TGraph*[NbOfJetsBins_];
		hScanMin_Eudsc[i] = new TGraph*[NbOfJetsBins_];
		hScanMin_Euds[i]  = new TGraph*[NbOfJetsBins_];
			for(int j=0;j<NbOfJetsBins_;j++){
			hScanMin_Eb[i][j]    = 0;
			hScanMin_Eudsc[i][j] = 0;
			hScanMin_Euds[i][j]  = 0;
		}
	}
	hScanMin_Ntt = new TGraph*[NbOfJetsBins_];
	hScanMin_Nv  = new TGraph*[NbOfJetsBins_];
	hScanMin_Nvb = new TGraph*[NbOfJetsBins_];
	for(int i=0;i<NbOfJetsBins_;i++){
		hScanMin_Ntt[i] = 0;
		hScanMin_Nv[i]  = 0;
		hScanMin_Nvb[i] = 0;
	}
	cout<<"Object from the class VJetEstimation correctly instantiated"<<endl;
}
/**________________________________________________________________________________________________________________*/
//VJetEstimation::VJetEstimation(int Iterations){
//}

/**________________________________________________________________________________________________________________*/
VJetEstimation::~VJetEstimation(){
	delete BtagWorkingPoint_;

	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		for(int j=0;j<NbOfJetsBins_;j++){
			for(int k=0;k<NbOfBJetsBins_;k++){
				delete [] N_[i][j][k];
				delete [] N_err_[i][j][k];
			}
			delete [] N_[i][j];
			delete [] N_err_[i][j];
			delete [] Nbjets_[i][j];
		}
		delete [] N_[i];
		delete [] N_err_[i];
		delete [] Nbjets_[i];
	}
	delete [] N_;
	delete [] N_err_;
	delete [] Nbjets_;

	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		delete [] eudsc_[i];
		delete [] euds_[i];
		delete [] eb_[i];

		for(int j=0;j<NbOfJetsBins_;j++){
			for(int k=0;k<3;k++){
				delete hNbjetsEstSummary[i][j][k];
				delete hNbjetsMCSummary[i][j][k];
				delete hNbjets_mc_[i][j][k];
				delete hNbjets_pdf_mc_[i][j][k];
				delete hNbjets_pdf_est_[i][j][k];
			}
			delete [] hNbjetsEstSummary[i][j];
			delete [] hNbjetsMCSummary[i][j];
			delete [] hNbjets_mc_[i][j];
			delete [] hNbjets_pdf_mc_[i][j];
			delete [] hNbjets_pdf_est_[i][j];
			delete hsNbjets_MC[i][j];
			delete hsNbjets_Est[i][j];
			delete tCanva_Nbjets_Summary[i][j];

			delete hScanMin_Eb[i][j];
			delete hScanMin_Eudsc[i][j];
			delete hScanMin_Euds[i][j];
		}
		delete [] hNbjetsEstSummary[i][NbOfJetsBins_];
		delete [] hNbjetsMCSummary[i][NbOfJetsBins_];
		delete hsNbjets_MC[i][NbOfJetsBins_];
		delete hsNbjets_Est[i][NbOfJetsBins_];
		delete tCanva_Nbjets_Summary[i][NbOfJetsBins_];

		for(int j=0;j<NbOfBJetsBins_;j++){
			for(int k=0;k<3;k++){
				delete hNjetsEstSummary[i][j][k];
				delete hNjetsMCSummary[i][j][k];
			}
			delete [] hNjetsEstSummary[i][j];
			delete [] hNjetsMCSummary[i][j];
			delete hsNjets_MC[i][j];
			delete hsNjets_Est[i][j];
			delete tCanva_Njets_Summary[i][j];
		}
		delete [] hNjetsEstSummary[i][NbOfBJetsBins_];
		delete [] hNjetsMCSummary[i][NbOfBJetsBins_];
		delete hsNjets_MC[i][NbOfBJetsBins_];
		delete hsNjets_Est[i][NbOfBJetsBins_];
		delete tCanva_Njets_Summary[i][NbOfBJetsBins_];

		delete [] hNbjetsEstSummary[i];
		delete [] hNbjetsMCSummary[i];
		delete [] hNbjets_mc_[i];
		delete [] hNbjets_pdf_mc_[i];
		delete [] hNbjets_pdf_est_[i];
		delete [] hNjetsEstSummary[i];
		delete [] hNjetsMCSummary[i];
		delete [] hsNbjets_MC[i];
		delete [] hsNbjets_Est[i];
		delete [] hsNjets_MC[i];
		delete [] hsNjets_Est[i];
		delete [] tCanva_Nbjets_Summary[i];
		delete [] tCanva_Njets_Summary[i];

		delete [] hScanMin_Eb[i];
		delete [] hScanMin_Eudsc[i];
		delete [] hScanMin_Euds[i];
	}

	for(int i=0;i<NbOfJetsBins_;i++){
		delete hScanMin_Ntt[i];
		delete hScanMin_Nv[i];
		delete hScanMin_Nvb[i];
	}
	delete [] hNbjetsEstSummary;
	delete [] hNbjetsMCSummary;
	delete [] hNbjets_mc_;
	delete [] hNbjets_pdf_mc_;
	delete [] hNbjets_pdf_est_;
	delete [] hNjetsEstSummary;
	delete [] hNjetsMCSummary;
	delete [] hsNbjets_MC;
	delete [] hsNbjets_Est;
	delete [] hsNjets_MC;
	delete [] hsNjets_Est;
	delete [] tCanva_Nbjets_Summary;
	delete [] tCanva_Njets_Summary;
	delete [] hScanMin_Eb;
	delete [] hScanMin_Eudsc;
	delete [] hScanMin_Euds;

	delete [] eudsc_;
	delete [] euds_;
	delete [] eb_;
	delete [] ebq_;
	delete [] Ntt_ ;
	delete [] Nv_;
	delete [] Nvb_;
	delete [] Njets_;
	delete [] MultiJet_Est_;

	for(int i=0;i<NbOfBtagWorkingPoint_-1;i++){
		for(int j=0;j<NbOfJetsBins_;j++){
			for(int k=0;k<NbOfDatasets_;k++){
				for(int l=0;l<NbOfBJetsBins_;l++){
					for(int n=0;n<NbOfBJetsBins_;n++) cout<<condProb_[i][j][k][l][n]<<endl;
					delete [] condProb_[i][j][k][l];
					cout<<"Here I am : "<<i<<"/"<<j<<"/"<<k<<"/"<<l<<endl;
				}
				delete [] condProb_[i][j][k];
				cout<<"Here I loop"<<endl;
			}
			delete [] condProb_[i][j];
			cout<<"Here I loop again"<<endl;
		}
		delete [] condProb_[i];
		cout<<"Here I loop again and again"<<endl;
	}
	delete [] condProb_;
	cout<<"Here I should be"<<endl;

	for(int i=0;i<NbOfDatasets_;i++){
		for(int j=0;j<NbOfJetsBins_;j++){
			delete hNbOfBGenJets_[i][j];
		}
		for(int j=0;j<NbOfBtagWorkingPoint_;j++){
			delete hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j];
			delete hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j];
			delete hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j];
			delete hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i][j];
			delete bTagEff_vs_Njets_[i][j];
			delete cTagEff_vs_Njets_[i][j];
			delete udsTagEff_vs_Njets_[i][j];
			delete misTagEff_vs_Njets_[i][j];
		}
		delete [] hNbOfBGenJets_[i];
		delete [] hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
		delete [] hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
		delete [] hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
		delete [] hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
		delete [] bTagEff_vs_Njets_[i];
		delete [] cTagEff_vs_Njets_[i];
		delete [] udsTagEff_vs_Njets_[i];
		delete [] misTagEff_vs_Njets_[i];
		delete hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
		delete hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
		delete hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
		delete hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[i];
	}
	delete [] hNbOfBGenJets_;
	delete [] hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;
	delete [] hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;
	delete [] hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;
	delete [] hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_;
	delete [] bTagEff_vs_Njets_;
	delete [] cTagEff_vs_Njets_;
	delete [] udsTagEff_vs_Njets_;
	delete [] misTagEff_vs_Njets_;
	delete [] hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_;
	delete [] hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_;
	delete [] hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_;
	delete [] hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_;
	
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		for(int j=0;j<NbOfJetsBins_+1;j++){
			delete [] NvVar_[i][j];
			delete [] NvBias_[i][j];
			delete [] NvbVar_[i][j];
			delete [] NvbBias_[i][j];
			delete [] NttVar_[i][j];
			delete [] NttBias_[i][j];
		}
		delete [] NvVar_[i];
		delete [] NvBias_[i];
		delete [] NvbVar_[i];
		delete [] NvbBias_[i];
		delete [] NttVar_[i];
		delete [] NttBias_[i];
	}
	delete [] NvVar_;
	delete [] NvBias_;
	delete [] NvbVar_;
	delete [] NvbBias_;
	delete [] NttVar_;
	delete [] NttBias_;

	if(RescaledTTLikeEstimation) delete RescaledTTLikeEstimation;
	if(RescaledVLikeEstimation)  delete RescaledVLikeEstimation;
	if(RescaledVbLikeEstimation) delete RescaledVbLikeEstimation;
	if(tCanva_RescaledTTLikeEstimation) delete tCanva_RescaledTTLikeEstimation;
	if(tCanva_RescaledVLikeEstimation)  delete tCanva_RescaledVLikeEstimation;
	if(tCanva_RescaledVbLikeEstimation) delete tCanva_RescaledVbLikeEstimation;

}

/**________________________________________________________________________________________________________________*/
//copy constructor
VJetEstimation::VJetEstimation(const VJetEstimation& vjet){}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::FillInputs(double**** n, double*** MultiJets_Estimated_Nj){
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		for(int j=0;j<NbOfJetsBins_;j++){
			for(int k=0;k<NbOfBJetsBins_;k++){
				Nbjets_[i][j][k]  = 0;
				for(int l=0;l<NbOfDatasets_;l++){
					N_[i][j][k][l] = n[i][j][k][l];
					Nbjets_[i][j][k]  += N_[i][j][k][l];
				}
				Nbjets_[i][j][k]   += -MultiJets_Estimated_Nj[i][j][k];
			}
		}
	}
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::FillInputs(double**** n){
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		for(int j=0;j<NbOfJetsBins_;j++){
			for(int k=0;k<NbOfBJetsBins_;k++){
				Nbjets_[i][j][k]  = 0;
				for(int l=0;l<NbOfDatasets_;l++){
					N_[i][j][k][l] = n[i][j][k][l];
					Nbjets_[i][j][k]  += N_[i][j][k][l];
				}
			}
		}
	}
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::FillInputs(vector<TRootJet*> &SelectedJets, int idx, int btagAlgo){
	double btagDisc = 0;
	
	unsigned int nbofjetsIdx  = ((SelectedJets.size()-Njets_[0])<(unsigned int)(NbOfJetsBins_-1) ? (SelectedJets.size()-Njets_[0]) : (NbOfJetsBins_-1));
	
	int nbofbjetsIdx[NbOfBtagWorkingPoint_];
	for(int i=0;i<NbOfBtagWorkingPoint_;i++) nbofbjetsIdx[i]=0;

	std::vector<TRootJet*> SelectedBGenJets;

	for(unsigned int i=0;i<SelectedJets.size();i++){
       	        if(     btagAlgo == 0) btagDisc = SelectedJets[i]->btag_trackCountingHighEffBJetTags();
       	        else if(btagAlgo == 1) btagDisc = SelectedJets[i]->btag_trackCountingHighPurBJetTags();
       	        else if(btagAlgo == 2) btagDisc = SelectedJets[i]->btag_jetProbabilityBJetTags();
       	        else if(btagAlgo == 3) btagDisc = SelectedJets[i]->btag_jetBProbabilityBJetTags();
       	        else if(btagAlgo == 4) btagDisc = SelectedJets[i]->btag_simpleSecondaryVertexHighEffBJetTags();
       	        else if(btagAlgo == 5) btagDisc = SelectedJets[i]->btag_simpleSecondaryVertexHighPurBJetTags();
	       	else if(btagAlgo == 6) btagDisc = SelectedJets[i]->btag_combinedSecondaryVertexBJetTags();
		else                   btagDisc = -999999;
			
		for(int j=0;j<NbOfBtagWorkingPoint_;j++){
			if(btagDisc>BtagWorkingPoint_[j]) nbofbjetsIdx[j]++;
		}

		if(fabs(SelectedJets[i]->partonFlavour()) > 6) continue;

		if(fabs(SelectedJets[i]->partonFlavour()) == 5)	     hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx]);
		else if(fabs(SelectedJets[i]->partonFlavour()) == 4) hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx]);
		else                                                 hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx]);

		if(fabs(SelectedJets[i]->partonFlavour()) != 5)	     hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx]);
		else SelectedBGenJets.push_back(SelectedJets[i]);

		for(int j=0;j<NbOfBtagWorkingPoint_;j++){
			if(btagDisc>BtagWorkingPoint_[j]){
				if(fabs(SelectedJets[i]->partonFlavour()) == 5)      hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][j]->Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx]);
				else if(fabs(SelectedJets[i]->partonFlavour()) == 4) hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][j]->Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx]);
				else                                                 hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][j]->Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx]);

				if(fabs(SelectedJets[i]->partonFlavour()) != 5)      hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][j]->Fill(SelectedJets[i]->Pt(),fabs(SelectedJets[i]->Eta()),Njets_[nbofjetsIdx]);
			}
	        }
	}
	for(int i=0;i<NbOfBtagWorkingPoint_;i++) N_[i][nbofjetsIdx][(nbofbjetsIdx[i]<NbOfBJetsBins_ ? nbofbjetsIdx[i] : NbOfBJetsBins_-1)][idx]++;
	for(int i=0;i<NbOfBtagWorkingPoint_-1;i++) condProb_[i][nbofjetsIdx][idx][(nbofbjetsIdx[i+1]<NbOfBJetsBins_ ? nbofbjetsIdx[i+1] : NbOfBJetsBins_-1)][(nbofbjetsIdx[i]<NbOfBJetsBins_ ? nbofbjetsIdx[i] : NbOfBJetsBins_-1)]++;
	hNbOfBGenJets_[idx][nbofjetsIdx]->Fill(SelectedBGenJets.size());

}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::FillInputs(vector<TRootJet> &SelectedJets, unsigned int idx, int btagAlgo){
	double btagDisc = 0;
	
	unsigned int nbofjetsIdx  = ((SelectedJets.size()-Njets_[0])<(unsigned int)(NbOfJetsBins_-1) ? (SelectedJets.size()-Njets_[0]) : (NbOfJetsBins_-1));
	
	int nbofbjetsIdx[NbOfBtagWorkingPoint_];
	for(int i=0;i<NbOfBtagWorkingPoint_;i++) nbofbjetsIdx[i]=0;
	
	std::vector<TRootJet> SelectedBGenJets;

	for(unsigned int i=0;i<SelectedJets.size();i++){
       	        if(     btagAlgo == 0) btagDisc = SelectedJets[i].btag_trackCountingHighEffBJetTags();
       	        else if(btagAlgo == 1) btagDisc = SelectedJets[i].btag_trackCountingHighPurBJetTags();
       	        else if(btagAlgo == 2) btagDisc = SelectedJets[i].btag_jetProbabilityBJetTags();
       	        else if(btagAlgo == 3) btagDisc = SelectedJets[i].btag_jetBProbabilityBJetTags();
       	        else if(btagAlgo == 4) btagDisc = SelectedJets[i].btag_simpleSecondaryVertexHighEffBJetTags();
       	        else if(btagAlgo == 5) btagDisc = SelectedJets[i].btag_simpleSecondaryVertexHighPurBJetTags();
	       	else if(btagAlgo == 6) btagDisc = SelectedJets[i].btag_combinedSecondaryVertexBJetTags();
		else                   btagDisc = -999999;
			
		for(int j=0;j<NbOfBtagWorkingPoint_;j++){
			if(btagDisc>BtagWorkingPoint_[j]) nbofbjetsIdx[j]++;
		}

		if(fabs(SelectedJets[i].partonFlavour()) > 6) continue;

		if(fabs(SelectedJets[i].partonFlavour()) == 5)	     hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Fill(SelectedJets[i].Pt(),fabs(SelectedJets[i].Eta()),Njets_[nbofjetsIdx]);
		else if(fabs(SelectedJets[i].partonFlavour()) == 4)  hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Fill(SelectedJets[i].Pt(),fabs(SelectedJets[i].Eta()),Njets_[nbofjetsIdx]);
		else                                                 hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Fill(SelectedJets[i].Pt(),fabs(SelectedJets[i].Eta()),Njets_[nbofjetsIdx]);

		if(fabs(SelectedJets[i].partonFlavour()) != 5)	     hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Fill(SelectedJets[i].Pt(),fabs(SelectedJets[i].Eta()),Njets_[nbofjetsIdx]);
		else SelectedBGenJets.push_back(SelectedJets[i]);

		for(int j=0;j<NbOfBtagWorkingPoint_;j++){
			if(btagDisc>BtagWorkingPoint_[j]){
				if(fabs(SelectedJets[i].partonFlavour()) == 5)      hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][j]->Fill(SelectedJets[i].Pt(),fabs(SelectedJets[i].Eta()),Njets_[nbofjetsIdx]);
				else if(fabs(SelectedJets[i].partonFlavour()) == 4) hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][j]->Fill(SelectedJets[i].Pt(),fabs(SelectedJets[i].Eta()),Njets_[nbofjetsIdx]);
				else                                                hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][j]->Fill(SelectedJets[i].Pt(),fabs(SelectedJets[i].Eta()),Njets_[nbofjetsIdx]);

				if(fabs(SelectedJets[i].partonFlavour()) != 5)      hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][j]->Fill(SelectedJets[i].Pt(),fabs(SelectedJets[i].Eta()),Njets_[nbofjetsIdx]);
			}
	        }
	}
	for(int i=0;i<NbOfBtagWorkingPoint_;i++) N_[i][nbofjetsIdx][(nbofbjetsIdx[i]<NbOfBJetsBins_ ? nbofbjetsIdx[i] : NbOfBJetsBins_-1)][idx]++;
	for(int i=0;i<NbOfBtagWorkingPoint_-1;i++) condProb_[i][nbofjetsIdx][idx][(nbofbjetsIdx[i+1]<NbOfBJetsBins_ ? nbofbjetsIdx[i+1] : NbOfBJetsBins_-1)][(nbofbjetsIdx[i]<NbOfBJetsBins_ ? nbofbjetsIdx[i] : NbOfBJetsBins_-1)]++;
	hNbOfBGenJets_[idx][nbofjetsIdx]->Fill(SelectedBGenJets.size());

}

/**________________________________________________________________________________________________________________*/
// void VJetEstimation::ReScaleInputs(double factor){
// 	for(int i=0;i<NbOfBtagWorkingPoint_;i++) {
// 	        for(int j=0;j<NbOfJetsBins_;j++){
// 			for(int k=0;k<NbOfBJetsBins_;k++){
// 				Nbjets_[i][j][k] *= factor;
// 				MultiJet_Est_[i][j][k] *= factor;
// 				for(int l=0;l<NbOfDatasets_;l++){
// 					N_err_[i][j][k][l]  = sqrt(N_[i][j][k][l]);
// 					N_err_[i][j][k][l] *= factor;
// 					N_[i][j][k][l]    *= factor;
// 				}
// 			}
// 		}
// 	}
// }

/**________________________________________________________________________________________________________________*/
float VJetEstimation::WilsonScoreIntervalHigh(float Non, float Ntot)
{
	double T = (Ntot>0 ? 1/Ntot : 0);
	double p_hat = (Ntot>0 && Non>=0 && Ntot>=Non ? Non/Ntot : 0);
	double Int_High = ((p_hat+(T/2))/(1+T))+(sqrt(p_hat*(1-p_hat)*T+pow(T/2,2))/(1+T));
	return Int_High;
}

/**________________________________________________________________________________________________________________*/
float VJetEstimation::WilsonScoreIntervalLow(float Non, float Ntot)
{
	double T = (Ntot>0 ? 1/Ntot : 0);
	double p_hat = (Ntot>0 && Non>=0 && Ntot>=Non ? Non/Ntot : 0);
	double Int_Low = ((p_hat+(T/2))/(1+T))-(sqrt(p_hat*(1-p_hat)*T+pow(T/2,2))/(1+T));
	return Int_Low;
}

/**________________________________________________________________________________________________________________*/
float VJetEstimation::WilsonScoreIntervalMean(float Non, float Ntot)
{
	double Err_High = WilsonScoreIntervalHigh(Non, Ntot)-Non;
	double Err_Low  = Non-WilsonScoreIntervalLow(Non, Ntot);
	return (Err_High+Err_Low)/2;
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::ReScaleInputs(int idx, float Ntot, double factor){
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		for(int j=0;j<NbOfJetsBins_;j++){
			for(int k=0;k<NbOfBJetsBins_;k++){
				N_err_[i][j][k][idx]  = Ntot*WilsonScoreIntervalMean(N_[i][j][k][idx],Ntot);
				N_err_[i][j][k][idx] *= factor;
				N_[i][j][k][idx]     *= factor;
				//cout<<"N_["<<i<<"]["<<j<<"]["<<k<<"]["<<idx<<"] = "<<N_[i][j][k][idx]<<endl;
				//cout<<"N_err_["<<i<<"]["<<j<<"]["<<k<<"]["<<idx<<"] = "<<N_err_[i][j][k][idx]<<endl;
				//cout<<"Ntot = "<<Ntot<<endl;
			}
		}
	}	
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::SumOverAllInputs(){
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		for(int j=0;j<NbOfJetsBins_;j++){
			for(int k=0;k<NbOfBJetsBins_;k++){
				Nbjets_[i][j][k] = 0; // Reset value
				for(int l=0;l<NbOfDatasets_;l++) Nbjets_[i][j][k] += N_[i][j][k][l]; //Loop over all the datasets and sum their weighted contributions.
			}
		}
	}	
}
/**________________________________________________________________________________________________________________*/
void VJetEstimation::ComputeEffbqFromMC(int idx){
	for(int i=0;i<NbOfJetsBins_;i++){
		e0bq_[i] = hNbOfBGenJets_[idx][i]->GetBinContent(1)/hNbOfBGenJets_[idx][i]->Integral();
		e1bq_[i] = hNbOfBGenJets_[idx][i]->GetBinContent(2)/hNbOfBGenJets_[idx][i]->Integral();
		e2bq_[i] = hNbOfBGenJets_[idx][i]->GetBinContent(3)/hNbOfBGenJets_[idx][i]->Integral();
	}
}
/**________________________________________________________________________________________________________________*/
void VJetEstimation::ComputeEffFromMC(){
	double x=0;
	char name[100];

	for(int idx=0;idx<NbOfDatasets_;idx++){
		for(int i=0;i<NbOfBtagWorkingPoint_;i++){
			if(hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Project3D("z")->GetEntries()>0){
				//cout<<"Nb of events = "<<hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Project3D("z")->GetEntries()<<endl;
				sprintf(name,"bTagEff_vs_Njets_%d_wp_%d",idx,i);
				bTagEff_vs_Njets_[idx][i] = new TGraphAsymmErrors(NbOfJetsBins_);
				bTagEff_vs_Njets_[idx][i]->SetName(name);
				bTagEff_vs_Njets_[idx][i]->BayesDivide((TH1D*)hNbOfBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][i]->Project3D("z"),(TH1D*)hNbOfBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Project3D("z"),"");
				if(bTagEff_vs_Njets_[idx][i]->GetXaxis()) bTagEff_vs_Njets_[idx][i]->GetXaxis()->SetTitle("Nb. of jets");
			}
			if(hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Project3D("z")->GetEntries()>0){
				//cout<<"Nb of events = "<<hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Project3D("z")->GetEntries()<<endl;
				sprintf(name,"cTagEff_vs_Njets_%d_wp_%d",idx,i);
				cTagEff_vs_Njets_[idx][i] = new TGraphAsymmErrors(NbOfJetsBins_);
				cTagEff_vs_Njets_[idx][i]->SetName(name);
				cTagEff_vs_Njets_[idx][i]->BayesDivide((TH1D*)hNbOfCGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][i]->Project3D("z"),(TH1D*)hNbOfCGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Project3D("z"),"");
				if(cTagEff_vs_Njets_[idx][i]->GetXaxis()) cTagEff_vs_Njets_[idx][i]->GetXaxis()->SetTitle("Nb. of jets");
			}
			if(hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Project3D("z")->GetEntries()>0){
				//cout<<"Nb of events = "<<hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Project3D("z")->GetEntries()<<endl;
				sprintf(name,"udsTagEff_vs_Njets_%d_wp_%d",idx,i);
				udsTagEff_vs_Njets_[idx][i] = new TGraphAsymmErrors(NbOfJetsBins_);
				udsTagEff_vs_Njets_[idx][i]->SetName(name);
				udsTagEff_vs_Njets_[idx][i]->BayesDivide((TH1D*)hNbOfUDSGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][i]->Project3D("z"),(TH1D*)hNbOfUDSGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Project3D("z"),"");
				if(udsTagEff_vs_Njets_[idx][i]->GetXaxis()) udsTagEff_vs_Njets_[idx][i]->GetXaxis()->SetTitle("Nb. of jets");
			}
			if(hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Project3D("z")->GetEntries()>0){
				//cout<<"Nb of events = "<<hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Project3D("z")->GetEntries()<<endl;
				sprintf(name,"misTagEff_vs_Njets_%d_wp_%d",idx,i);
				misTagEff_vs_Njets_[idx][i] = new TGraphAsymmErrors(NbOfJetsBins_);
				misTagEff_vs_Njets_[idx][i]->SetName(name);
				misTagEff_vs_Njets_[idx][i]->BayesDivide((TH1D*)hNbOfNonBGenBJets_vs_JetPt_vs_JetEta_vs_Njets_[idx][i]->Project3D("z"),(TH1D*)hNbOfNonBGenJets_vs_JetPt_vs_JetEta_vs_Njets_[idx]->Project3D("z"),"");
				if(misTagEff_vs_Njets_[idx][i]->GetXaxis()) misTagEff_vs_Njets_[idx][i]->GetXaxis()->SetTitle("Nb. of jets");
			}
		}
	}

	double btag_tmp[iDatasetsTTLike_.size()];
	double btag_err_tmp[iDatasetsTTLike_.size()];
	double mistag_tmp[NbOfDatasets_];
	double mistag_err_tmp[NbOfDatasets_];

	for(int j=0;j<NbOfBtagWorkingPoint_;j++){
		for(int k=0;k<NbOfJetsBins_;k++){
			if( GetPredNtt(j,k) == 0 || GetPredNv(j,k) == 0) continue;
			eb_mc_[j][k]    = 0; eb_err_mc_[j][k] = 0;
			eudsc_mc_[j][k] = 0; eudsc_err_mc_[j][k] = 0;
			euds_mc_[j][k]  = 0; euds_err_mc_[j][k] = 0;
			for(unsigned int i=0;i<iDatasetsTTLike_.size();i++){
				if(bTagEff_vs_Njets_[iDatasetsTTLike_[i]][j]!=0){
					if(bTagEff_vs_Njets_[iDatasetsTTLike_[i]][j]->GetPoint(k,x,btag_tmp[i]) != -1){
						//bTagEff_vs_Njets_[iDatasetsTTLike_[i]][j]->GetPoint(k,x,btag_tmp[i]);
						btag_err_tmp[i]   = bTagEff_vs_Njets_[iDatasetsTTLike_[i]][j]->GetErrorY(k);
						eb_mc_[j][k]     += btag_tmp[i]*GetPredN(iDatasetsTTLike_[i],k)/GetPredNtt(0,k);
						eb_err_mc_[j][k] += pow(btag_err_tmp[i],2)*GetPredN(iDatasetsTTLike_[i],k)/GetPredNtt(0,k);
					}
				}
				if(misTagEff_vs_Njets_[iDatasetsTTLike_[i]][j]!=0){
					if(misTagEff_vs_Njets_[iDatasetsTTLike_[i]][j]->GetPoint(k,x,mistag_tmp[i]) != -1){
						//misTagEff_vs_Njets_[iDatasetsTTLike_[i]][j]->GetPoint(k,x,mistag_tmp[i]);
						mistag_err_tmp[i]    = misTagEff_vs_Njets_[iDatasetsTTLike_[i]][j]->GetErrorY(k);
						eudsc_mc_[j][k]     += mistag_tmp[i]*GetPredN(iDatasetsTTLike_[i],k)/GetPredNtt(0,k);
						eudsc_err_mc_[j][k] += pow(mistag_err_tmp[i],2)*GetPredN(iDatasetsTTLike_[i],k)/GetPredNtt(0,k);
					}
				}
			}
			for(unsigned int i=0;i<iDatasetsVLike_.size();i++){
				if(misTagEff_vs_Njets_[iDatasetsVLike_[i]][j]!=0){
					if(misTagEff_vs_Njets_[iDatasetsVLike_[i]][j]->GetPoint(k,x,mistag_tmp[i]) != -1){
						//misTagEff_vs_Njets_[iDatasetsVLike_[i]][j]->GetPoint(k,x,mistag_tmp[i]);
						mistag_err_tmp[i]   = misTagEff_vs_Njets_[iDatasetsVLike_[i]][j]->GetErrorY(k);
						euds_mc_[j][k]     += mistag_tmp[i]*GetPredN(iDatasetsVLike_[i],k)/GetPredNv(0,k);
						euds_err_mc_[j][k] += pow(mistag_err_tmp[i],2)*GetPredN(iDatasetsVLike_[i],k)/GetPredNv(0,k);
					}
				}
			}
			eb_err_mc_[j][k]    = sqrt(eb_err_mc_[j][k]);
			eudsc_err_mc_[j][k] = sqrt(eudsc_err_mc_[j][k]);
			euds_err_mc_[j][k]  = sqrt(euds_err_mc_[j][k]);
		}
	}
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::ComputeCondProb(bool Verbose){
	double sum = 0;
	for(int idx=0;idx<NbOfDatasets_;idx++){	
		for(int i=0;i<NbOfBtagWorkingPoint_-1;i++){
			for(int j=0;j<NbOfJetsBins_;j++){
				for(int k=0;k<NbOfBJetsBins_;k++){
					if(Verbose) cout<<"Sum of the prob for(i,j,idx,k) = ("<<i<<"/"<<j<<"/"<<idx<<"/"<<k<<") is equal to : ";
					sum = 0;
					for(int l=0;l<NbOfBJetsBins_;l++) sum += condProb_[i][j][idx][k][l];
					for(int l=0;l<NbOfBJetsBins_;l++){
						if(sum>0)condProb_[i][j][idx][k][l] = condProb_[i][j][idx][k][l]/sum;
						else if(l==k) condProb_[i][j][idx][k][l] = 1;
						if(Verbose) cout<<condProb_[i][j][idx][k][l]<<"/";
					}
					if(Verbose) cout<<endl;
				}
			}
		}
	}
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::Write(TFile* file, string label, bool verbose){
	cout<<"--Start writing histograms for the VJetEstimation --"<<endl;
	file->cd();
	string dirname = "VJetEstimation"+label;
	file->mkdir(dirname.c_str());
	file->cd(dirname.c_str());

	char name[100];char title[100];
	
	for(int i=0;i<NbOfJetsBins_;i++){
		if(hScanMin_Ntt[i]!=0) hScanMin_Ntt[i]->Write();
		if(hScanMin_Nv[i]!=0)  hScanMin_Nv[i] ->Write();
		if(hScanMin_Nvb[i]!=0) hScanMin_Nvb[i]->Write();
	}
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		for(int j=0;j<NbOfJetsBins_;j++){
			if(hScanMin_Eb[i][j]!=0)    hScanMin_Eb[i][j]->Write();
			if(hScanMin_Eudsc[i][j]!=0) hScanMin_Eudsc[i][j]->Write();
			if(hScanMin_Euds[i][j]!=0)  hScanMin_Euds[i][j]->Write();
		}
		for(int j=0;j<NbOfJetsBins_;j++){
			for(int k=0;k<3;k++){
				hNbjets_mc_[i][j][k]->Write();
				hNbjets_pdf_mc_[i][j][k]->Write();
				hNbjets_pdf_est_[i][j][k]->Write();
			}
		}
		for(int j=0;j<=NbOfJetsBins_;j++){
			// For Monte Carlo prediction 
			MyLeg->Clear();
			tCanva_Nbjets_Summary[i][j]->cd();

			if(j!=NbOfJetsBins_) sprintf(name,"hNbjetsMCStackSummary_%d_jets_wp_%d",Njets_[j],i);
			else                 sprintf(name,"hNbjetsMCStackSummary_Inclusive_wp_%d",i);
			if(j!=NbOfJetsBins_) sprintf(title,"V-like and tt-like MC prediction summary for %d jets (WP nr%d)",Njets_[j],i);
			else                 sprintf(title,"V-like and tt-like MC prediction summary (inclusive, WP nr%d)",i);

			hsNbjets_MC[i][j] = new THStack(name,title);

			hNbjetsMCSummary[i][j][0]->SetFillStyle(3004);
			hNbjetsMCSummary[i][j][0]->SetFillColor(kGreen+1);
			hNbjetsMCSummary[i][j][1]->SetFillStyle(3005);
			hNbjetsMCSummary[i][j][1]->SetFillColor(kBlue+1);
			hNbjetsMCSummary[i][j][2]->SetFillStyle(3006);
			hNbjetsMCSummary[i][j][2]->SetFillColor(kRed+1);

			hsNbjets_MC[i][j]->Add(hNbjetsMCSummary[i][j][0]);
			//hsNbjets_MC[i][j]->Add(hNbjetsMCSummary[i][j][1]);
			hsNbjets_MC[i][j]->Add(hNbjetsMCSummary[i][j][2]);
		
			MyLeg->AddEntry(hNbjetsMCSummary[i][j][0],"V-like events","f");
			//MyLeg->AddEntry(hNbjetsMCSummary[i][j][1],"Vb-like events","f");
			MyLeg->AddEntry(hNbjetsMCSummary[i][j][2],"#bar{t}t-like events","f");
		
			hsNbjets_MC[i][j]->Draw();
			hsNbjets_MC[i][j]->GetXaxis()->SetTitle("Nb of b-jets");

			if(i!=NbOfJetsBins_) sprintf(name,"hNbjetsEstStackSummary_%d_jets_wp_%d",Njets_[j],i);
			else                 sprintf(name,"hNbjetsEstStackSummary_Inclusive_wp_%d",i);
			if(i!=NbOfJetsBins_) sprintf(title,"V-like and tt-like estimation summary for %d jets (WP nr%d)",Njets_[j],i);
			else                 sprintf(title,"V-like and tt-like estimation summary (inclusive, WP nr%d)",i);

			hsNbjets_Est[i][j] = new THStack(name,title);

			hNbjetsEstSummary[i][j][0]->SetMarkerColor(kYellow+1);
			//hNbjetsEstSummary[i][j][1]->SetMarkerColor(kMagenta+1);
			hNbjetsEstSummary[i][j][2]->SetMarkerColor(kBlue+1);
			hNbjetsEstSummary[i][j][0]->SetMarkerStyle(22);
			//hNbjetsEstSummary[i][j][1]->SetMarkerStyle(23);
			hNbjetsEstSummary[i][j][2]->SetMarkerStyle(24);
			hNbjetsEstSummary[i][j][0]->SetMarkerSize(1.5);
			//hNbjetsEstSummary[i][j][1]->SetMarkerSize(1.5);
			hNbjetsEstSummary[i][j][2]->SetMarkerSize(1.5);
	
			hsNbjets_Est[i][j]->Add(hNbjetsEstSummary[i][j][0]);
			//hsNbjets_Est[i][j]->Add(hNbjetsEstSummary[i][j][1]);
			hsNbjets_Est[i][j]->Add(hNbjetsEstSummary[i][j][2]);
		
			MyLeg->AddEntry(hNbjetsEstSummary[i][j][0],"V-like estimation","p");
			//MyLeg->AddEntry(hNbjetsEstSummary[i][j][1],"Vb-like estimation","p");
			MyLeg->AddEntry(hNbjetsEstSummary[i][j][2],"TT-like estimation","p");
		
			hsNbjets_Est[i][j]->Draw("PEsame");
			hsNbjets_Est[i][j]->GetXaxis()->SetTitle("Nb of b-jets");
			MyLeg->Draw("same");

			tCanva_Nbjets_Summary[i][j]->Update();
			if(verbose)cout<<"Writing summary histograms for "<<j+4<<" jets, wp nr "<<i<<endl;
			tCanva_Nbjets_Summary[i][j]->Write();
		}
		for(int j=0;j<=NbOfBJetsBins_;j++){

		// For Monte Carlo prediction 
		MyLeg->Clear();
		tCanva_Njets_Summary[i][j]->cd();
		
		if(j!=4) sprintf(name,"hNjetsMCStackSummary_%d_bjets_wp_%d",j,i);
		else     sprintf(name,"hNjetsMCStackSummary_bInclusive_wp_%d",i);
		if(j!=4) sprintf(title,"V-like and tt-like MC prediction summary for %d b-jets (WP nr%d)",j,i);
		else     sprintf(title,"V-like and tt-like MC prediction summary (b-inclusive,WP nr%d)",i);
		hsNjets_MC[i][j] = new THStack(name,title);

		hNjetsMCSummary[i][j][0]->SetFillStyle(3004);
		hNjetsMCSummary[i][j][0]->SetFillColor(kGreen+1);
		//hNjetsMCSummary[i][j][1]->SetFillStyle(3005);
		//hNjetsMCSummary[i][j][1]->SetFillColor(kBlue+1);
		hNjetsMCSummary[i][j][2]->SetFillStyle(3006);
		hNjetsMCSummary[i][j][2]->SetFillColor(kRed+1);

		hsNjets_MC[i][j]->Add(hNjetsMCSummary[i][j][0]);
		//hsNjets_MC[i][j]->Add(hNjetsMCSummary[i][j][1]);
		hsNjets_MC[i][j]->Add(hNjetsMCSummary[i][j][2]);
		
		MyLeg->AddEntry(hNjetsMCSummary[i][j][0],"V-like events","f");
		//MyLeg->AddEntry(hNjetsMCSummary[i][j][1],"Vb-like events","f");
		MyLeg->AddEntry(hNjetsMCSummary[i][j][2],"#bar{t}t-like events","f");
		
		hsNjets_MC[i][j]->Draw();
		hsNjets_MC[i][j]->GetXaxis()->SetTitle("Nb of jets");

		//hNjetsEstSummary[i][0]->Write();
		//hNjetsEstSummary[i][1]->Write();

		if(j!=4) sprintf(name,"hNjetsEstStackSummary_%d_bjets_wp_%d",j,i);
		else     sprintf(name,"hNjetsEstStackSummary_bInclusive_wp_%d",i);
		if(j!=4) sprintf(title,"V-like and tt-like estimation summary for %d b-jets (WP nr%d)",j,i);
		else     sprintf(title,"V-like and tt-like estimation summary (b-inclusive,WP nr%d)",i);
		hsNjets_Est[i][j] = new THStack(name,title);

		//sprintf(name,"hNjetsEstStackLegend_%d_bjets",i);
		//MyLeg->SetName(name);

		hNjetsEstSummary[i][j][0]->SetMarkerColor(kYellow+1);
		//hNjetsEstSummary[i][j][1]->SetMarkerColor(kMagenta+1);
		hNjetsEstSummary[i][j][2]->SetMarkerColor(kBlue+1);
		hNjetsEstSummary[i][j][0]->SetMarkerStyle(22);
		//hNjetsEstSummary[i][j][1]->SetMarkerStyle(23);
		hNjetsEstSummary[i][j][2]->SetMarkerStyle(24);
		hNjetsEstSummary[i][j][0]->SetMarkerSize(1.5);
		//hNjetsEstSummary[i][j][1]->SetMarkerSize(1.5);
		hNjetsEstSummary[i][j][2]->SetMarkerSize(1.5);

		hsNjets_Est[i][j]->Add(hNjetsEstSummary[i][j][0]);
		//hsNjets_Est[i][j]->Add(hNjetsEstSummary[i][j][1]);
		hsNjets_Est[i][j]->Add(hNjetsEstSummary[i][j][2]);
		
		MyLeg->AddEntry(hNjetsEstSummary[i][j][0],"V-like estimation","p");
		//MyLeg->AddEntry(hNjetsEstSummary[i][j][1],"Vb-like estimation","p");
		MyLeg->AddEntry(hNjetsEstSummary[i][j][2],"TT-like estimation","p");
		
		hsNjets_Est[i][j]->Draw("PE SAME");
		hsNjets_Est[i][j]->GetXaxis()->SetTitle("Nb of jets");
		MyLeg->Draw("same");

		tCanva_Njets_Summary[i][j]->Update();
		tCanva_Njets_Summary[i][j]->Write();
		if(verbose)cout<<"Writing summary histograms for "<<j<<" b-jets, wp nr "<<i<<endl;
	}
}
	for(int idx=0;idx<NbOfDatasets_;idx++){
		for(int i=0;i<NbOfJetsBins_;i++){
			if(hNbOfBGenJets_[idx][i]==0) continue;
			if(hNbOfBGenJets_[idx][i]->Integral() > 0) hNbOfBGenJets_[idx][i]->Scale(1/hNbOfBGenJets_[idx][i]->Integral());
			hNbOfBGenJets_[idx][i]->Write();
		}
		for(int i=0;i<NbOfBtagWorkingPoint_;i++){
			if(  bTagEff_vs_Njets_[idx][i]!=0)   bTagEff_vs_Njets_[idx][i]->Write();
			if(  cTagEff_vs_Njets_[idx][i]!=0)   cTagEff_vs_Njets_[idx][i]->Write();
			if(udsTagEff_vs_Njets_[idx][i]!=0) udsTagEff_vs_Njets_[idx][i]->Write();
			if(misTagEff_vs_Njets_[idx][i]!=0) misTagEff_vs_Njets_[idx][i]->Write();
		}
	}
	if(RescaledTTLikeEstimation) RescaledTTLikeEstimation->Write();
	if(RescaledVLikeEstimation)  RescaledVLikeEstimation->Write();
	if(RescaledVbLikeEstimation) RescaledVbLikeEstimation->Write();
	if(tCanva_RescaledTTLikeEstimation) tCanva_RescaledTTLikeEstimation->Write();
	if(tCanva_RescaledVLikeEstimation)  tCanva_RescaledVLikeEstimation->Write();
	if(tCanva_RescaledVbLikeEstimation) tCanva_RescaledVbLikeEstimation->Write();
	cout<<"--Done with writing histograms for the VJetEstimation --"<<endl;
}

/**________________________________________________________________________________________________________________*/
/////////////////////////////////
//Main Methods to be executed ...
/////////////////////////////////
/**________________________________________________________________________________________________________________*/
void VJetEstimation::BinnedMaximumLikelihoodEst(int btag_wp_idx, int njets, string &method, string &option, vector<int> &FixedVarIdx, bool doScan, bool Verbose){

	if(init_ == 0) {cout<<"Can't solve the system of equations. Missing initial values..."<<endl; return;}
	ROOT::Math::Minimizer* min = 0;	
	min = ROOT::Math::Factory::CreateMinimizer(method,option);

	if(min==0) return;

	min->SetMaxFunctionCalls(1000000);
	min->SetMaxIterations(100000);
	min->SetPrecision(-1);
	min->SetTolerance(0.00001);
	min->SetErrorDef(0.5);
  
	fFunc = new ROOT::Math::Functor (this,&VJetEstimation::LikelihoodFunct,7);
	// 0th parameter : b-tagging working point index
	// 1th parameter : jet multiplicity
	// 2nd parameter : Ntt-like
	// 3rd parameter : Nv-like
	// 4th parameter : eb
	// 5th parameter : euds
	// 6th parameter : eud
	min->SetFunction(*fFunc);

	const int NbOfPar = 5;
	string wp[3] = {"loose","medium","tight"};
	string VarNames[NbOfPar+1];
	VarNames[0] = "btag_wp_idx";
	VarNames[1] = "njets";
	VarNames[2] = "Ntt-like";
	VarNames[3] = "Nv-like";
	VarNames[4] = "eb (";	 VarNames[4] += wp[btag_wp_idx]; VarNames[4] += " wp)"; 
	VarNames[5] = "eudsc ("; VarNames[5] += wp[btag_wp_idx]; VarNames[5] += " wp)";
	VarNames[6] = "euds (";  VarNames[6] += wp[btag_wp_idx]; VarNames[6] += " wp)";
	
	const double Nfrac   = 4;
	const double Efffrac = 9999;

	if(Verbose){
		cout<<setprecision(4);
		cout<<"Set fixed variable ("<<VarNames[0]<<") to value : "<<btag_wp_idx<<endl;
		cout<<"Set fixed variable ("<<VarNames[1]<<") to value : "<<njets<<endl;
		cout<<"Set variable ("<<VarNames[2]<<") to initial value : "<<init_[njets][0]<<", step : 0.001,["<<init_[njets][0]*(1-Nfrac)<<","<<init_[njets][0]*(1+Nfrac)<<"]"<<endl;
		cout<<"Set variable ("<<VarNames[3]<<") to initial value : "<<init_[njets][1]<<", step : 0.001,["<<init_[njets][1]*(1-Nfrac)<<","<<init_[njets][1]*(1+Nfrac)<<"]"<<endl;
	}
	min->SetFixedVariable(0,VarNames[0],btag_wp_idx);
	min->SetFixedVariable(1,VarNames[1],njets);
	min->SetLimitedVariable(2,VarNames[2],init_[njets][0],0.00001,(init_[njets][0]*(1-Nfrac)<0 ? 0 : init_[njets][0]*(1-Nfrac)),init_[njets][0]*(1+Nfrac));
	min->SetLimitedVariable(3,VarNames[3],init_[njets][1],0.00001,(init_[njets][1]*(1-Nfrac)<0 ? 0 : init_[njets][1]*(1-Nfrac)),init_[njets][1]*(1+Nfrac));

	vector<int>::iterator Idx;

	// b/mis-tagging efficiencies declared as LimitedVariable, comprised  between 0 and 1
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 0 );
	if(Idx == FixedVarIdx.end()){
		if(Verbose) cout<<"Set variable ("<<VarNames[4]<<") to initial value : "<<init_[njets][3*btag_wp_idx+2+0]<<", step : 0.00001"<<endl;
		min->SetLimitedVariable(4,VarNames[4],init_[njets][3*btag_wp_idx+2+0],0.00001,(init_[njets][3*btag_wp_idx+2+0]*(1-Efffrac)<0 ? 0 : init_[njets][3*btag_wp_idx+2+0]*(1-Efffrac)),(init_[njets][3*btag_wp_idx+2+0]*(1+Efffrac)>1 ? 1 : init_[njets][3*btag_wp_idx+2+0]*(1+Efffrac)));
	}
	else{
		if(Verbose) cout<<"Set fixed variable ("<<VarNames[4]<<") to initial value : "<<init_[njets][3*btag_wp_idx+2+0]<<endl;
		min->SetFixedVariable(4,VarNames[4],init_[njets][3*btag_wp_idx+2+0]);
	}
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 1 );
	if(Idx == FixedVarIdx.end()){
		if(Verbose) cout<<"Set variable ("<<VarNames[5]<<") to initial value : "<<init_[njets][3*btag_wp_idx+2+1]<<", step : 0.00001"<<endl;
		min->SetLimitedVariable(5,VarNames[5],init_[njets][3*btag_wp_idx+2+1],0.00001,(init_[njets][3*btag_wp_idx+2+1]*(1-Efffrac)<0 ? 0 : init_[njets][3*btag_wp_idx+2+1]*(1-Efffrac)),(init_[njets][3*btag_wp_idx+2+1]*(1+Efffrac)>1 ? 1 : init_[njets][3*btag_wp_idx+2+1]*(1+Efffrac)));
	}
	else{
		if(Verbose) cout<<"Set fixed variable ("<<VarNames[5]<<") to initial value : "<<init_[njets][3*btag_wp_idx+2+1]<<endl;
		min->SetFixedVariable(5,VarNames[5],init_[njets][3*btag_wp_idx+2+1]);
	}
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 2 );
	if(Idx == FixedVarIdx.end()){
		if(Verbose) cout<<"Set variable ("<<VarNames[6]<<") to initial value : "<<init_[njets][3*btag_wp_idx+2+2]<<", step : 0.00001"<<endl;
		min->SetLimitedVariable(6,VarNames[6],init_[njets][3*btag_wp_idx+2+2],0.00001,(init_[njets][3*btag_wp_idx+2+2]*(1-Efffrac)<0 ? 0 : init_[njets][3*btag_wp_idx+2+2]*(1-Efffrac)),(init_[njets][3*btag_wp_idx+2+2]*(1+Efffrac)>1 ? 1 : init_[njets][3*btag_wp_idx+2+2]*(1+Efffrac)));
	}
	else{
		if(Verbose) cout<<"Set fixed variable ("<<VarNames[6]<<") to initial value : "<<init_[njets][3*btag_wp_idx+2+2]<<endl;
		min->SetFixedVariable(6,VarNames[6],init_[njets][3*btag_wp_idx+2+2]);
	}
	
	min->Minimize();
	min->Hesse();
	//min->Minos();

	const double *param = min->X();
	const double *err  = min->Errors();
	//double errUp = 0;
	//double errLow = 0;
	if(Verbose){
		cout << "Number of (free) parameters = " << min->NFree() << " / " << min->NDim() << endl;
		cout << "Error definition : " << min->ErrorDef() << endl;
		cout << "Are the computed erors valid ? " << (min->IsValidError() ? "Yes" : "No" ) << endl;
		cout << "Status code : " << min->Status() << endl;
		cout << "Minimum value : f = " << min->MinValue() << " is obtained for :" << endl;
		if(param !=0 && err !=0){
			for(int i=0;i<NbOfPar+2;i++){
				//min->GetMinosError(i,errLow,errUp);
				cout << fixed<<setprecision(4);
				cout << "--" << VarNames[i] << " = " << param[i] <<" +/- "<< /*errLow <<"/"<< errUp << "(" << */err[i] /*<< ")" */<< endl;
			}
		}
		cout << "Estimated distance to the minimum (EDM) : " << min->Edm() << endl;
		min->PrintResults();
	}

	minValue_[njets] = min->MinValue();

	unsigned int nsteps = 400;
	double *x = new double[nsteps];
	double *y = new double[nsteps];
	char name[100];

	sprintf(name,"hScanMin_Ntt_%d_jets",njets);
	Ntt_[njets]                = param[2]; if(doScan) {min->Scan(2,nsteps,x,y,param[2]*0.8,param[2]*1.2); hScanMin_Ntt[njets]                = new TGraph(nsteps,x,y);hScanMin_Ntt[njets]->SetNameTitle(name,name);}

	sprintf(name,"hScanMin_Nv_%d_jets",njets);
	Nv_[njets]                 = param[3]; if(doScan) {min->Scan(3,nsteps,x,y,param[3]*0.8,param[3]*1.2); hScanMin_Nv[njets]   	       = new TGraph(nsteps,x,y);hScanMin_Nv[njets] ->SetNameTitle(name,name);}

	sprintf(name,"hScanMin_Eb_%d_wp_%d_jets",btag_wp_idx,njets);
	eb_[btag_wp_idx][njets]    = param[4]; if(doScan) {min->Scan(4,nsteps,x,y,param[4]*0.8,param[4]*1.2); hScanMin_Eb[btag_wp_idx][njets]    = new TGraph(nsteps,x,y);hScanMin_Eb[btag_wp_idx][njets]->SetNameTitle(name,name);}

	sprintf(name,"hScanMin_Eudsc_%d_wp_%d_jets",btag_wp_idx,njets);
	eudsc_[btag_wp_idx][njets] = param[5]; if(doScan) {min->Scan(5,nsteps,x,y,param[5]*0.8,param[5]*1.2); hScanMin_Eudsc[btag_wp_idx][njets] = new TGraph(nsteps,x,y);hScanMin_Eudsc[btag_wp_idx][njets]->SetNameTitle(name,name);}

	sprintf(name,"hScanMin_Euds_%d_wp_%d_jets",btag_wp_idx,njets);
	euds_[btag_wp_idx][njets]  = param[6]; if(doScan) {min->Scan(6,nsteps,x,y,param[6]*0.8,param[6]*1.2); hScanMin_Euds[btag_wp_idx][njets]  = new TGraph(nsteps,x,y);hScanMin_Euds[btag_wp_idx][njets]->SetNameTitle(name,name);}

	delete min;
	if(Verbose) cout<<"-- End of the EquationSolverByChiSqMin method --"<<endl;
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::BinnedMaximumJointWPLikelihoodEst(int njets, string &method, string &option, vector<int> &FixedVarIdx, bool doScan, bool Verbose){

	if(init_ == 0) {cout<<"Can't solve the system of equations. Missing initial values..."<<endl; return;}
	ROOT::Math::Minimizer* min = 0;	
	min = ROOT::Math::Factory::CreateMinimizer(method,option);

	if(min==0) return;

	min->SetMaxFunctionCalls(1000000);
	min->SetMaxIterations(100000);
	min->SetPrecision(-1);
	min->SetTolerance(0.00001);
	min->SetErrorDef(0.5);
  
	fFunc = new ROOT::Math::Functor (this,&VJetEstimation::JointWPLikelihoodFunct,(3+3*NbOfBtagWorkingPoint_));
	// 0th parameter : jet multiplicity
	// 1st parameter : Ntt-like
	// 2nd parameter : Nv-like
	// 3rd parameter : Nvb-like
	// 4/7/10th parameter : eb (working point)
	// 5/8/11th parameter : euds (working point)
	// 6/9/12th parameter : eud (working point)
	min->SetFunction(*fFunc);

	const int NbOfPar = 2+3*NbOfBtagWorkingPoint_;
	string wp[3] = {"loose","medium","tight"};
	string VarNames[NbOfPar+1];
	VarNames[0] = "njets";
	VarNames[1] = "Ntt-like";
	VarNames[2] = "Nv-like";
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		VarNames[3+i*3] = "eb (";    VarNames[3+i*3] += wp[i]; VarNames[3+i*3] += " wp)"; 
		VarNames[4+i*3] = "eudsc ("; VarNames[4+i*3] += wp[i]; VarNames[4+i*3] += " wp)";
		VarNames[5+i*3] = "euds (";  VarNames[5+i*3] += wp[i]; VarNames[5+i*3] += " wp)";
	}

	const double Nfrac   = 4;
	const double Efffrac = 9999;

	if(Verbose){
		cout<<"Set fixed variable ("<<VarNames[0]<<") to value : "<<njets<<endl;
		cout<<"Set limited variable ("<<VarNames[1]<<") to initial value : "<<init_[njets][0]<<", step : 0.001,["<<init_[njets][0]*(1-Nfrac)<<","<<init_[njets][0]*(1+Nfrac)<<"]"<<endl;
		cout<<"Set limited variable ("<<VarNames[2]<<") to initial value : "<<init_[njets][1]<<", step : 0.001,["<<init_[njets][1]*(1-Nfrac)<<","<<init_[njets][1]*(1+Nfrac)<<"]"<<endl;
	}
	min->SetFixedVariable(0,VarNames[0],njets);
	min->SetLimitedVariable(1,VarNames[1],init_[njets][0],0.00001,(init_[njets][0]*(1-Nfrac)<0 ? 0 : init_[njets][0]*(1-Nfrac)),init_[njets][0]*(1+Nfrac));
	min->SetLimitedVariable(2,VarNames[2],init_[njets][1],0.00001,(init_[njets][1]*(1-Nfrac)<0 ? 0 : init_[njets][1]*(1-Nfrac)),init_[njets][1]*(1+Nfrac));

	vector<int>::iterator Idx;

	// b/mis-tagging efficiencies declared as LimitedVariable, comprised  between 0 and 1
	for(int i=0;i<NbOfBtagWorkingPoint_;i++) {
		Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 0 );
		if(Idx == FixedVarIdx.end()){
			if(Verbose) cout<<"Set variable ("<<VarNames[3+3*i]<<") to initial value : "<<init_[njets][2+3*i]<<", step : 0.00001"<<endl;
			min->SetLimitedVariable(3+3*i,VarNames[3+3*i],init_[njets][2+3*i],0.00001,(init_[njets][2+3*i]*(1-Efffrac)<0 ? 0 : init_[njets][2+3*i]*(1-Efffrac)),(init_[njets][2+3*i]*(1+Efffrac)>1 ? 1 : init_[njets][2+3*i]*(1+Efffrac)));
		}
		else{
			if(Verbose) cout<<"Set fixed variable ("<<VarNames[3+3*i]<<") to initial value : "<<init_[njets][2+3*i]<<endl;
			min->SetFixedVariable(3+3*i,VarNames[3+3*i],init_[njets][2+3*i]);
		}
		Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 1 );
		if(Idx == FixedVarIdx.end()){
			if(Verbose) cout<<"Set variable ("<<VarNames[4+3*i]<<") to initial value : "<<init_[njets][3+3*i]<<", step : 0.00001"<<endl;
			min->SetLimitedVariable(4+3*i,VarNames[4+3*i],init_[njets][3+3*i],0.00001,(init_[njets][3+3*i]*(1-Efffrac)<0 ? 0 : init_[njets][3+3*i]*(1-Efffrac)),(init_[njets][3+3*i]*(1+Efffrac)>1 ? 1 : init_[njets][3+3*i]*(1+Efffrac)));
		}
		else{
			if(Verbose) cout<<"Set fixed variable ("<<VarNames[4+3*i]<<") to initial value : "<<init_[njets][3+3*i]<<endl;
			min->SetFixedVariable(4+3*i,VarNames[4+3*i],init_[njets][3+3*i]);
		}
		Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 2 );
		if(Idx == FixedVarIdx.end()){
			if(Verbose) cout<<"Set variable ("<<VarNames[5+3*i]<<") to initial value : "<<init_[njets][4+3*i]<<", step : 0.00001"<<endl;
			min->SetLimitedVariable(5+3*i,VarNames[5+3*i],init_[njets][4+3*i],0.00001,(init_[njets][4+3*i]*(1-Efffrac)<0 ? 0 : init_[njets][4+3*i]*(1-Efffrac)),(init_[njets][4+3*i]*(1+Efffrac)>1 ? 1 : init_[njets][4+3*i]*(1+Efffrac)));
		}
		else{
			if(Verbose) cout<<"Set fixed variable ("<<VarNames[5+3*i]<<") to initial value : "<<init_[njets][4+3*i]<<endl;
			min->SetFixedVariable(5+3*i,VarNames[5+3*i],init_[njets][4+3*i]);
		}
	}
	min->Minimize();
	min->Hesse();
	//min->Minos();

	const double *param = min->X();
	const double *err  = min->Errors();
	//double errUp = 0;
	//double errLow = 0;
	if(Verbose){
		cout << "Number of (free) parameters = " << min->NFree() << " / " << min->NDim() << endl;
		cout << "Error definition : " << min->ErrorDef() << endl;
		cout << "Are the computed erors valid ? " << (min->IsValidError() ? "Yes" : "No" ) << endl;
		cout << "Status code : " << min->Status() << endl;
		cout << "Minimum value : f = " << min->MinValue() << " is obtained for :" << endl;
		if(param !=0 && err !=0){
			for(int i=0;i<NbOfPar+1;i++){
				//min->GetMinosError(i,errLow,errUp);
				cout << fixed<<setprecision(4);
				cout << "--" << VarNames[i] << " = " << param[i] <<" +/- "<< /*errLow <<"/"<< errUp << "(" << */err[i] /*<< ")" */<< endl;
			}
		}
		cout << "Estimated distance to the minimum (EDM) : " << min->Edm() << endl;
		min->PrintResults();
	}

	minValue_[njets] = min->MinValue();

	unsigned int nsteps = 400;
	double *x = new double[nsteps];
	double *y = new double[nsteps];
	char name[100];
	sprintf(name,"hScanMin_Ntt_%d_jets",Njets_[njets]);
	Ntt_[njets] = param[1]; if(doScan) {min->Scan(1,nsteps,x,y,param[1]*0.8,param[1]*1.2); hScanMin_Ntt[njets] = new TGraph(nsteps,x,y);hScanMin_Ntt[njets]->SetNameTitle(name,name);}
	sprintf(name,"hScanMin_Nv_%d_jets",Njets_[njets]);
	Nv_[njets]  = param[2]; if(doScan) {min->Scan(2,nsteps,x,y,param[2]*0.8,param[2]*1.2); hScanMin_Nv[njets]  = new TGraph(nsteps,x,y);hScanMin_Nv[njets] ->SetNameTitle(name,name);}
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		sprintf(name,"hScanMin_Eb_%d_wp_%d_jets",i,Njets_[njets]);
		eb_[i][njets]    = param[3+i*3]; if(doScan) {min->Scan(3+i*3,nsteps,x,y,param[3+i*3]*0.8,param[3+i*3]*1.2); hScanMin_Eb[i][njets]    = new TGraph(nsteps,x,y);hScanMin_Eb[i][njets]->SetNameTitle(name,name);}
		sprintf(name,"hScanMin_Eudsc_%d_wp_%d_jets",i,Njets_[njets]);
		eudsc_[i][njets] = param[4+i*3]; if(doScan) {min->Scan(4+i*3,nsteps,x,y,param[4+i*3]*0.8,param[4+i*3]*1.2); hScanMin_Eudsc[i][njets] = new TGraph(nsteps,x,y);hScanMin_Eudsc[i][njets]->SetNameTitle(name,name);}
		sprintf(name,"hScanMin_Euds_%d_wp_%d_jets",i,Njets_[njets]);
		euds_[i][njets]  = param[5+i*3]; if(doScan) {min->Scan(5+i*3,nsteps,x,y,param[5+i*3]*0.8,param[5+i*3]*1.2); hScanMin_Euds[i][njets]  = new TGraph(nsteps,x,y);hScanMin_Euds[i][njets]->SetNameTitle(name,name);}
	}
	delete min;
	if(Verbose) cout<<"-- End of the EquationSolverByChiSqMin method --"<<endl;
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::BinnedMaximumJointLikelihoodEst(string &method, string &option, bool doScan, bool Verbose){

	if(init_ == 0) {cout<<"Can't solve the system of equations. Missing initial values..."<<endl; return;}
	ROOT::Math::Minimizer* min = 0;	
	min = ROOT::Math::Factory::CreateMinimizer(method,option);

	if(min==0) return;

	min->SetMaxFunctionCalls(1000000);
	min->SetMaxIterations(100000);
	min->SetPrecision(-1);
	min->SetTolerance(0.00001);
	min->SetErrorDef(0.5);
  
	fFunc = new ROOT::Math::Functor (this,&VJetEstimation::JointLikelihoodFunct,(2*NbOfJetsBins_+3*NbOfBtagWorkingPoint_));
	// 1st parameter : Ntt-like
	// 2nd parameter : Nv-like
	// 3rd parameter : Nvb-like
	// 4/7/10th parameter : eb (working point)
	// 5/8/11th parameter : euds (working point)
	// 6/9/12th parameter : eud (working point)
	min->SetFunction(*fFunc);

	const int NbOfPar = 2*NbOfJetsBins_+3*NbOfBtagWorkingPoint_;
	string wp[3] = {"loose","medium","tight"};
	char name[100];
	string VarNames[NbOfPar];
	for(int i=0;i<NbOfJetsBins_;i++){
		sprintf(name,"Ntt-like_%d_jets",Njets_[i]);
		VarNames[0+i*NbOfJetsBins_] = name;
		sprintf(name,"Nv-like_%d_jets",Njets_[i]);
		VarNames[1+i*NbOfJetsBins_] = name;
	}
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		VarNames[0+2*NbOfJetsBins_+i*3] = "eb (";    VarNames[0+2*NbOfJetsBins_+i*3] += wp[i]; VarNames[0+2*NbOfJetsBins_+i*3] += " wp)"; 
		VarNames[1+2*NbOfJetsBins_+i*3] = "eudsc ("; VarNames[1+2*NbOfJetsBins_+i*3] += wp[i]; VarNames[1+2*NbOfJetsBins_+i*3] += " wp)";
		VarNames[2+2*NbOfJetsBins_+i*3] = "euds (";  VarNames[2+2*NbOfJetsBins_+i*3] += wp[i]; VarNames[2+2*NbOfJetsBins_+i*3] += " wp)";
	}

	for(int i=0;i<NbOfJetsBins_;i++){
		if(Verbose) cout<<"Set limited variable ("<<VarNames[0+2*i]<<") to initial value : "<<init_[i][0]<<", step : 0.001,["<<init_[i][0]/5<<","<<init_[i][0]*5<<"]"<<endl;
		min->SetLimitedVariable(0+2*i,VarNames[0+2*i],init_[i][0],0.00001,init_[i][0]/5,init_[i][0]*5);

		if(Verbose) cout<<"Set limited variable ("<<VarNames[1+2*i]<<") to initial value : "<<init_[i][1]<<", step : 0.001,["<<init_[i][1]/5<<","<<init_[i][1]*5<<"]"<<endl;
		min->SetLimitedVariable(1+2*i,VarNames[1+2*i],init_[i][1],0.00001,init_[i][1]/5,init_[i][1]*5);
	}
	// proabilities and b/mis-tagging efficiencies declared as LimitedVariable, comprised  between 0 and 1
	for(int i=0;i<NbOfBtagWorkingPoint_;i++) {
		if(Verbose) cout<<"Set variable ("<<VarNames[0+2*NbOfJetsBins_+i*3]<<") to initial value : "<<init_[0][0+2+i*3]<<", step : 0.001"<<endl;
		min->SetLimitedVariable(0+2*NbOfJetsBins_+i*3,VarNames[0+2*NbOfJetsBins_+i*3],init_[0][0+2*NbOfJetsBins_+i*3],0.00001,0.0001,1);
		if(Verbose) cout<<"Set variable ("<<VarNames[1+2*NbOfJetsBins_+i*3]<<") to initial value : "<<init_[0][1+2+i*3]<<", step : 0.001"<<endl;
		min->SetLimitedVariable(1+2*NbOfJetsBins_+i*3,VarNames[1+2*NbOfJetsBins_+i*3],init_[0][1+2*NbOfJetsBins_+i*3],0.00001,0.0001,1);
		if(Verbose) cout<<"Set variable ("<<VarNames[2+2*NbOfJetsBins_+i*3]<<") to initial value : "<<init_[0][2+2+i*3]<<", step : 0.001"<<endl;
		min->SetLimitedVariable(2+2*NbOfJetsBins_+i*3,VarNames[2+2*NbOfJetsBins_+i*3],init_[0][2+2*NbOfJetsBins_+i*3],0.00001,0.0001,1);
	}

	min->Minimize();
	min->Hesse();
	//min->Minos();

	const double *param = min->X();
	const double *err  = min->Errors();
	//double errUp = 0;
	//double errLow = 0;
	if(Verbose){
		cout << "Number of (free) parameters = " << min->NFree() << " / " << min->NDim() << endl;
		cout << "Error definition : " << min->ErrorDef() << endl;
		cout << "Are the computed erors valid ? " << (min->IsValidError() ? "Yes" : "No" ) << endl;
		cout << "Status code : " << min->Status() << endl;
		cout << "Minimum value : f = " << min->MinValue() << " is obtained for :" << endl;
		for(int i=0;i<NbOfPar;i++){
			//min->GetMinosError(i,errLow,errUp);
			cout << fixed<<setprecision(4);
			cout << "--" << VarNames[i] << " = " << param[i] <<" +/- "<< /*errLow <<"/"<< errUp << "(" << */err[i] /*<< ")" */<< endl;
		}
		cout << "Estimated distance to the minimum (EDM) : " << min->Edm() << endl;
		min->PrintResults();
	}

	for(int i=0;i<NbOfJetsBins_;i++) minValue_[i] = min->MinValue();

	unsigned int nsteps = 400;
	double *x = new double[nsteps];
	double *y = new double[nsteps];
	for(int i=0;i<NbOfJetsBins_;i++){
		sprintf(name,"hScanMin_Ntt_%d_jets",Njets_[i]);
		Ntt_[i] = param[0+2*i]; if(doScan) {min->Scan(0+2*i,nsteps,x,y,param[0+2*i]*0.8,param[0+2*i]*1.2); hScanMin_Ntt[i] = new TGraph(nsteps,x,y);hScanMin_Ntt[i]->SetNameTitle(name,name);}
		sprintf(name,"hScanMin_Nv_%d_jets",Njets_[i]);
		Nv_[i]  = param[1+2*i]; if(doScan) {min->Scan(1+2*i,nsteps,x,y,param[1+2*i]*0.8,param[1+2*i]*1.2); hScanMin_Nv[i]  = new TGraph(nsteps,x,y);hScanMin_Nv[i] ->SetNameTitle(name,name);}
	}
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		for(int j=0;j<NbOfJetsBins_;j++){
			sprintf(name,"hScanMin_Eb_%d_wp_%d_jets",i,Njets_[j]);
			eb_[i][j]    = param[0+2*NbOfJetsBins_+i*3]; if(doScan) {min->Scan(0+2*NbOfJetsBins_+i*3,nsteps,x,y,param[0+2*NbOfJetsBins_+i*3]*0.8,param[0+2*NbOfJetsBins_+i*3]*1.2); hScanMin_Eb[i][j]    = new TGraph(nsteps,x,y);hScanMin_Eb[i][j]->SetNameTitle(name,name);}
			sprintf(name,"hScanMin_Eudsc_%d_wp_%d_jets",i,Njets_[j]);
			eudsc_[i][j] = param[1+2*NbOfJetsBins_+i*3]; if(doScan) {min->Scan(1+2*NbOfJetsBins_+i*3,nsteps,x,y,param[1+2*NbOfJetsBins_+i*3]*0.8,param[1+2*NbOfJetsBins_+i*3]*1.2); hScanMin_Eudsc[i][j] = new TGraph(nsteps,x,y);hScanMin_Eudsc[i][j]->SetNameTitle(name,name);}
			sprintf(name,"hScanMin_Euds_%d_wp_%d_jets",i,Njets_[j]);
			euds_[i][j]  = param[2+2*NbOfJetsBins_+i*3]; if(doScan) {min->Scan(2+2*NbOfJetsBins_+i*3,nsteps,x,y,param[2+2*NbOfJetsBins_+i*3]*0.8,param[2+2*NbOfJetsBins_+i*3]*1.2); hScanMin_Euds[i][j]  = new TGraph(nsteps,x,y);hScanMin_Euds[i][j]->SetNameTitle(name,name);}
		}
	}
	delete min;
	if(Verbose) cout<<"-- End of the EquationSolverByChiSqMin method --"<<endl;
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::UnBinnedMaximumLikelihoodEst(int btag_wp_idx, int njets, vector<int> &FixedVarIdx, bool doMinos, bool doToyMC, bool Verbose){

  // Declare observables

	RooRealVar Ntt("Ntt","N_{t#bar{t}-like}",init_[njets][0],init_[njets][0]/5,init_[njets][0]*5);//GetPredNtt(btag_wp_idx,njets),0,10*GetPredNtt(btag_wp_idx,njets));
	RooRealVar Nv( "Nv", "N_{V-like}",       init_[njets][1],init_[njets][1]/5,init_[njets][1]*5);//GetPredNv( btag_wp_idx,njets),0,10*GetPredNv( btag_wp_idx,njets));

	vector<int>::iterator Idx;

	RooRealVar  eb(    "eb",  "#epsilon_{b-tag}",      init_[njets][3*btag_wp_idx+2],0.0,1.0);
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 0 );
	if(Idx != FixedVarIdx.end()) eb.setConstant(kTRUE) ;
	RooRealVar  eudsc("eudsc","#epsilon_{mis-tag}",    init_[njets][3*btag_wp_idx+3],0.0,1.0);
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 1 );
	if(Idx != FixedVarIdx.end()) eudsc.setConstant(kTRUE) ;
	RooRealVar euds( "euds",  "#epsilon^{,}_{mis-tag}",init_[njets][3*btag_wp_idx+4],0.0,1.0);
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 2 );
	if(Idx != FixedVarIdx.end()) euds.setConstant(kTRUE) ;

  // Declare constants

	RooConstVar n("n","number of selected jets",Njets_[njets]) ;
	RooConstVar e0bq("e0bq","e0bq",e0bq_[njets]);
	RooConstVar e1bq("e1bq","e1bq",e1bq_[njets]);
	RooConstVar e2bq("e2bq","e2bq",e2bq_[njets]);
	//RooConstVar e3bq("e3bq","e3bq",e3bq_[njets]);

  // C o n s t r u c t   a   c a t e g o r y   w i t h   l a b e l s    a n d   i n d e c e s
  // -----------------------------------------------------------------------------------------

	RooCategory nbjets("nbjets","Number of b-jets");
	nbjets.defineType("N0bjet", 0);
	nbjets.defineType("N1bjet", 1);
	nbjets.defineType("N2bjets",2);
	nbjets.defineType("N3bjets",3);

  // C o n s t r u c t   p . d . f 's
  // -------------------------------------------------------------------------------

	RooGenericPdf pbjets_tt("pbjets_tt","pbjets_tt","(nbjets==0)*((1-eb)*(1-eb)*pow((1-eudsc),n-2)*e2bq+(1-eb)*pow((1-eudsc),n-1)*e1bq+pow((1-eudsc),n)*e0bq)+(nbjets==1)*((2*eb*(1-eb)*pow(1-eudsc,n-2)+(1-eb)*(1-eb)*(n-2)*eudsc*pow(1-eudsc,n-3))*e2bq+(eb*pow(1-eudsc,n-1)+(1-eb)*(n-1)*eudsc*pow(1-eudsc,n-2))*e1bq+(n*eudsc*pow(1-eudsc,n-1))*e0bq)+(nbjets==2)*((eb*eb*pow(1-eudsc,n-2)+2*eb*(1-eb)*(n-2)*eudsc*pow(1-eudsc,n-3)+(1-eb)*(1-eb)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4))*e2bq+(eb*(n-1)*eudsc*pow(1-eudsc,n-2)+(1-eb)*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3))*e1bq+((n*(n-1)/2)*eudsc*eudsc*pow(1-eudsc,n-2))*e0bq)+(nbjets==3)*((eb*eb*(n-2)*eudsc*pow(1-eudsc,n-3)+2*eb*(1-eb)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4)+(n>4 ? pow((1-eb),2)*((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow((1-eudsc),n-5) : 0 ))*e2bq+(eb*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3)+(1-eb)*((n-1)*(n-2)*(n-3)/6)*pow(eudsc,3)*pow(1-eudsc,n-4))*e1bq+((n*(n-1)*(n-2)/6)*pow(eudsc,3)*pow(1-eudsc,n-3))*e0bq)",RooArgList(nbjets,n,eb,eudsc,e0bq,e1bq,e2bq));
	RooExtendPdf  pbjets_tt_ext("pbjets_tt_ext","pbjets_tt_xt",pbjets_tt,Ntt);

	RooGenericPdf pbjets_v("pbjets_v","pbjets_v","(nbjets==0)*pow(1-euds,n)+(nbjets==1)*n*euds*pow(1-euds,n-1)+(nbjets==2)*(n*(n-1)/2)*euds*euds*pow(1-euds,n-2)+(nbjets==3)*((n)*(n-1)*(n-2)/6)*pow(euds,3)*pow(1-euds,n-3)",RooArgList(nbjets,n,euds));
	RooExtendPdf  pbjets_v_ext("pbjets_v_ext","pbjets_v_xt",pbjets_v,Nv);

	RooAddPdf model("model","model",RooArgList(pbjets_tt,pbjets_v),RooArgList(Ntt,Nv));
	
  // C r e a t e   d a t a s e t 
  // -------------------------------------------------------------------------------
  // Sample a dataset for tt+jets events

	RooDataSet data("data","data",RooArgSet(nbjets)) ;
        for (int i=0 ; i<GetPredNtotal(btag_wp_idx,njets,0) ; i++) { nbjets.setLabel("N0bjet")  ; data.add(RooArgSet(nbjets));}
        for (int i=0 ; i<GetPredNtotal(btag_wp_idx,njets,1) ; i++) { nbjets.setLabel("N1bjet")  ; data.add(RooArgSet(nbjets));}
        for (int i=0 ; i<GetPredNtotal(btag_wp_idx,njets,2) ; i++) { nbjets.setLabel("N2bjets") ; data.add(RooArgSet(nbjets));}
        for (int i=0 ; i<GetPredNtotal(btag_wp_idx,njets,3) ; i++) { nbjets.setLabel("N3bjets") ; data.add(RooArgSet(nbjets));}

  // F i t   t h e   d a t a   a n d   c o n s t r u c t   t h e   l i k e l i h o o d   f u n c t i o n
  // ----------------------------------------------------------------------------------------------

	//RooFitResult* fit_result = model.fitTo(data,Save(),Extended(1),PrintLevel(-1));//,ConditionalObservables(nbjets),SplitRange(1),Optimize(0));
	RooAbsReal* nll = model.createNLL(data);

	RooMinuit minimizer(*nll);

	minimizer.optimizeConst(kTRUE) ;
	minimizer.setPrintLevel(-1);
	minimizer.setNoWarn();
	
	minimizer.migrad();
	
	minimizer.hesse();
	
	if(doMinos) minimizer.minos();
	
	RooFitResult* fit_result = minimizer.save();
	
	if(Verbose)fit_result->Print("v");

	Ntt_[njets]                = Ntt.getVal();
	Ntt_err_[njets]            = Ntt.getError();
	Ntt_err_up_[njets]         = Ntt.getErrorHi();
	Ntt_err_down_[njets]       = Ntt.getErrorLo();
	Nv_[njets]                 = Nv.getVal();
	Nv_err_[njets]             = Nv.getError();
	Nv_err_up_[njets]          = Nv.getErrorHi();
	Nv_err_down_[njets]        = Nv.getErrorLo();

	minValue_[njets]           = fit_result->minNll();

	eb_[btag_wp_idx][njets]        = eb.getVal();
	eb_err_[btag_wp_idx][njets]    = eb.getError();
	eudsc_[btag_wp_idx][njets]     = eudsc.getVal();
	eudsc_err_[btag_wp_idx][njets] = eudsc.getError();
	euds_[btag_wp_idx][njets]      = euds.getVal();
	euds_err_[btag_wp_idx][njets]  = euds.getError();
	
	delete fit_result;
	delete nll;
	
  // G e n e r a t e   a n d   f i t   e v e n t s
  // ---------------------------------------------

  // Generate and fit 1000 samples of Poisson(nExpected) events
	RooMCStudy* mcstudy = 0;
	RooDataSet* data_mcstudy = 0;
	RooWorkspace *w = 0;

if(doToyMC){
	mcstudy = new RooMCStudy(model,nbjets,Binned(kTRUE),Silence(),Extended(kTRUE),FitOptions(Save(kFALSE),PrintEvalErrors(-1),Extended(kTRUE))) ;
	mcstudy->generateAndFit(2000,kFALSE) ;
	*data_mcstudy = mcstudy->fitParDataSet();

  // C r e a t e   w o r k s p a c e ,   i m p o r t   d a t a   a n d   m o d e l
  // -----------------------------------------------------------------------------

  // Create a new empty workspace
	w = new RooWorkspace("w","workspace") ;

  // Import model and all its components into the workspace
	w->import(model) ;

  // Import data into the workspace
	w->import(data) ;
	if(doToyMC) w->import(*data_mcstudy) ;
  // Print workspace contents
	//w->Print() ;

  // S a v e   w o r k s p a c e   i n   f i l e
  // -------------------------------------------

  // Save the workspace into a ROOT file
	char name[100];
	sprintf(name,"VJetEstimation_RooFit_ToyMC_%d_wp_%d_jets.root",btag_wp_idx,Njets_[njets]);
	w->writeToFile(name) ;
}	
	delete mcstudy;
	delete data_mcstudy;
	delete w;
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::UnBinnedMaximumJointWPLikelihoodEst(int njets, vector<int> &FixedVarIdx, bool doMinos, bool Verbose){

  // Declare observables
	TString wp[3] = {"loose","medium","tight"};
	vector<int>::iterator Idx;

	RooRealVar Ntt("Ntt","N_{t#bar{t}-like}",init_[njets][0],init_[njets][0]/5,init_[njets][0]*5);
	RooRealVar Nv( "Nv", "N_{V-like}",       init_[njets][1],init_[njets][1]/5,init_[njets][1]*5);

	RooRealVar* eb[NbOfBtagWorkingPoint_];
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 0 );
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		eb[i] = new RooRealVar( "eb_"+wp[i],"#epsilon_{b-tag} ("+wp[i]+")", init_[njets][3*i+2],0.0,1.0);
		if(Idx != FixedVarIdx.end()) eb[i]->setConstant(kTRUE);
	}
	RooRealVar* eudsc[NbOfBtagWorkingPoint_];
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 1 );
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		eudsc[i] = new RooRealVar( "eudsc_"+wp[i],"#epsilon_{mis-tag} ("+wp[i]+")", init_[njets][3*i+3],0.0,1.0);
		if(Idx != FixedVarIdx.end()) eudsc[i]->setConstant(kTRUE);
	}
	RooRealVar* euds[NbOfBtagWorkingPoint_];
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 2 );
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		euds[i]  = new RooRealVar( "euds_"+wp[i],"#epsilon^{,}_{mis-tag} ("+wp[i]+")", init_[njets][3*i+4],0.0,1.0);//GetPredEuds( i,njets),0.0,1.0);
		if(Idx != FixedVarIdx.end()) euds[i]->setConstant(kTRUE);
	}

  // Declare constants

	RooConstVar n("n","number of selected jets",Njets_[njets]) ;
	RooConstVar e0bq("e0bq","e0bq",e0bq_[njets]);
	RooConstVar e1bq("e1bq","e1bq",e1bq_[njets]);
	RooConstVar e2bq("e2bq","e2bq",e2bq_[njets]);
	//RooConstVar e3bq("e3bq","e3bq",e3bq_[njets]);

  // C o n s t r u c t   a   c a t e g o r y   w i t h   l a b e l s    a n d   i n d e c e s
  // -----------------------------------------------------------------------------------------

	RooCategory* nbjets[NbOfBtagWorkingPoint_];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		nbjets[i] = new RooCategory("nbjets_"+wp[i],"Number of b-jets");
		nbjets[i] -> defineType("N0bjet", 0);
		nbjets[i] -> defineType("N1bjet", 1);
		nbjets[i] -> defineType("N2bjets",2);
		nbjets[i] -> defineType("N3bjets",3);
	}

  // C o n s t r u c t   p . d . f 's
  // -------------------------------------------------------------------------------

	RooGenericPdf* pbjets_tt[NbOfBtagWorkingPoint_];
	RooExtendPdf*  pbjets_tt_ext[NbOfBtagWorkingPoint_];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		//pbjets_tt[i]     = new RooGenericPdf("pbjets_tt_"+wp[i],"pbjets_tt_"+wp[i],"((*nbjets[i])==0)*((1-(*eb[i]))*(1-(*eb[i]))*pow((1-(*eudsc[i])),n-2)*e2bq+(1-(*eb[i]))*pow((1-(*eudsc[i])),n-1)*e1bq+pow((1-(*eudsc[i])),n)*e0bq)+((*nbjets[i])==1)*((2*(*eb[i])*(1-(*eb[i]))*pow(1-(*eudsc[i]),n-2)+(1-(*eb[i]))*(1-(*eb[i]))*(n-2)*(*eudsc[i])*pow(1-(*eudsc[i]),n-3))*e2bq+((*eb[i])*pow(1-(*eudsc[i]),n-1)+(1-(*eb[i]))*(n-1)*(*eudsc[i])*pow(1-(*eudsc[i]),n-2))*e1bq+(n*(*eudsc[i])*pow(1-(*eudsc[i]),n-1))*e0bq)+((*nbjets[i])==2)*(((*eb[i])*(*eb[i])*pow(1-(*eudsc[i]),n-2)+2*(*eb[i])*(1-(*eb[i]))*(n-2)*(*eudsc[i])*pow(1-(*eudsc[i]),n-3)+(1-(*eb[i]))*(1-(*eb[i]))*((n-2)*(n-3)/2)*(*eudsc[i])*(*eudsc[i])*pow(1-(*eudsc[i]),n-4))*e2bq+((*eb[i])*(n-1)*(*eudsc[i])*pow(1-(*eudsc[i]),n-2)+(1-(*eb[i]))*((n-1)*(n-2)/2)*(*eudsc[i])*(*eudsc[i])*pow(1-(*eudsc[i]),n-3))*e1bq+((n*(n-1)/2)*(*eudsc[i])*(*eudsc[i])*pow(1-(*eudsc[i]),n-2))*e0bq)+((*nbjets[i])==3)*(((*eb[i])*(*eb[i])*(n-2)*(*eudsc[i])*pow(1-(*eudsc[i]),n-3)+2*(*eb[i])*(1-(*eb[i]))*((n-2)*(n-3)/2)*(*eudsc[i])*(*eudsc[i])*pow(1-(*eudsc[i]),n-4)+(n>4 ? pow((1-(*eb[i])),2)*((n-2)*(n-3)*(n-4)/6)*pow((*eudsc[i]),3)*pow((1-(*eudsc[i])),n-5) : 0 ))*e2bq+((*eb[i])*((n-1)*(n-2)/2)*(*eudsc[i])*(*eudsc[i])*pow(1-(*eudsc[i]),n-3)+(1-(*eb[i]))*((n-1)*(n-2)*(n-3)/6)*pow((*eudsc[i]),3)*pow(1-(*eudsc[i]),n-4))*e1bq+((n*(n-1)*(n-2)/6)*pow((*eudsc[i]),3)*pow(1-(*eudsc[i]),n-3))*e0bq)",RooArgList(*nbjets[i],n,*eb[i],*eudsc[i],e0bq,e1bq,e2bq));
		pbjets_tt[i]     = new RooGenericPdf("pbjets_tt_"+wp[i],"pbjets_tt_"+wp[i],"(@0==0)*((1-@1)*(1-@1)*pow((1-@2),@3-2)*@6+(1-@1)*pow((1-@2),@3-1)*@5+pow((1-@2),@3)*@4)+(@0==1)*((2*@1*(1-@1)*pow(1-@2,@3-2)+(1-@1)*(1-@1)*(@3-2)*@2*pow(1-@2,@3-3))*@6+(@1*pow(1-@2,@3-1)+(1-@1)*(@3-1)*@2*pow(1-@2,@3-2))*@5+(@3*@2*pow(1-@2,@3-1))*@4)+(@0==2)*((@1*@1*pow(1-@2,@3-2)+2*@1*(1-@1)*(@3-2)*@2*pow(1-@2,@3-3)+(1-@1)*(1-@1)*((@3-2)*(@3-3)/2)*@2*@2*pow(1-@2,@3-4))*@6+(@1*(@3-1)*@2*pow(1-@2,@3-2)+(1-@1)*((@3-1)*(@3-2)/2)*@2*@2*pow(1-@2,@3-3))*@5+((@3*(@3-1)/2)*@2*@2*pow(1-@2,@3-2))*@4)+(@0==3)*((@1*@1*(@3-2)*@2*pow(1-@2,@3-3)+2*@1*(1-@1)*((@3-2)*(@3-3)/2)*@2*@2*pow(1-@2,@3-4)+(@3>4 ? pow((1-@1),2)*((@3-2)*(@3-3)*(@3-4)/6)*pow(@2,3)*pow((1-@2),@3-5) : 0 ))*@6+(@1*((@3-1)*(@3-2)/2)*@2*@2*pow(1-@2,@3-3)+(1-@1)*((@3-1)*(@3-2)*(@3-3)/6)*pow(@2,3)*pow(1-@2,@3-4))*@5+((@3*(@3-1)*(@3-2)/6)*pow(@2,3)*pow(1-@2,@3-3))*@4)",RooArgList(*nbjets[i],*eb[i],*eudsc[i],n,e0bq,e1bq,e2bq));
		pbjets_tt_ext[i] = new RooExtendPdf( "pbjets_tt_ext_"+wp[i],"pbjets_tt_ext_"+wp[i],*pbjets_tt[i],Ntt);
	}

	RooGenericPdf* pbjets_v[NbOfBtagWorkingPoint_];
	RooExtendPdf*  pbjets_v_ext[NbOfBtagWorkingPoint_];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		//pbjets_v[i]     = new RooGenericPdf("pbjets_v_"+wp[i],"pbjets_v_"+wp[i],"((*nbjets[i])==0)*pow(1-(*euds[i]),n)+((*nbjets[i])==1)*n*(*euds[i])*pow(1-(*euds[i]),n-1)+((*nbjets[i])==2)*(n*(n-1)/2)*(*euds[i])*(*euds[i])*pow(1-(*euds[i]),n-2)+((*nbjets[i])==3)*((n)*(n-1)*(n-2)/6)*pow((*euds[i]),3)*pow(1-(*euds[i]),n-3)",RooArgList(*nbjets[i],*euds[i],n));
		pbjets_v[i]     = new RooGenericPdf("pbjets_v_"+wp[i],"pbjets_v_"+wp[i],"(@0==0)*pow(1-@1,@2)+(@0==1)*@2*@1*pow(1-@1,@2-1)+(@0==2)*(@2*(@2-1)/2)*@1*@1*pow(1-@1,@2-2)+(@0==3)*((@2)*(@2-1)*(@2-2)/6)*pow(@1,3)*pow(1-@1,@2-3)",RooArgList(*nbjets[i],*euds[i],n));
		pbjets_v_ext[i] = new RooExtendPdf("pbjets_v_ext_"+wp[i],"pbjets_v_ext_"+wp[i],*pbjets_v[i],Nv);
	}

	RooAddPdf* model[NbOfBtagWorkingPoint_];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		model[i] = new RooAddPdf("model"+wp[i],"model"+wp[i],RooArgList(*pbjets_tt[i],*pbjets_v[i]),RooArgList(Ntt,Nv));
	}
	
  // C r e a t e   d a t a s e t 
  // -------------------------------------------------------------------------------

	RooDataSet* data[NbOfBtagWorkingPoint_];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		data[i] = new RooDataSet("data"+wp[i],"data"+wp[i],RooArgSet(*nbjets[i])) ;
		for (int j=0 ; j<GetPredNtotal(i,njets,0) ; j++) { nbjets[i]->setLabel("N0bjet")  ; data[i]->add(RooArgSet(*nbjets[i]));}
		for (int j=0 ; j<GetPredNtotal(i,njets,1) ; j++) { nbjets[i]->setLabel("N1bjet")  ; data[i]->add(RooArgSet(*nbjets[i]));}
		for (int j=0 ; j<GetPredNtotal(i,njets,2) ; j++) { nbjets[i]->setLabel("N2bjets") ; data[i]->add(RooArgSet(*nbjets[i]));}
		for (int j=0 ; j<GetPredNtotal(i,njets,3) ; j++) { nbjets[i]->setLabel("N3bjets") ; data[i]->add(RooArgSet(*nbjets[i]));}
	}

  // F i t   t h e   d a t a   a n d   c o n s t r u c t   t h e   l i k e l i h o o d   f u n c t i o n
  // ----------------------------------------------------------------------------------------------

	RooAbsReal* nll[NbOfBtagWorkingPoint_];
	RooArgSet nlls("myNLLs");
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		nll[i] = model[i]->createNLL(*data[i]);
		nlls.add(*nll[i]);
	}
	
	RooAddition* combNLL = new RooAddition("combNLL","combined likelihood function",nlls);
	
	RooMinuit minimizer(*combNLL);
	
	minimizer.optimizeConst(kTRUE) ;
	minimizer.setPrintLevel(-1);
	minimizer.setNoWarn();
	
	minimizer.migrad();
	
	minimizer.hesse();
	
	if(doMinos) minimizer.minos();
	
	RooFitResult* fit_result = minimizer.save();

	if(Verbose)fit_result->Print("v");

	Ntt_[njets]                = Ntt.getVal();
	Ntt_err_[njets]            = Ntt.getError();
	Ntt_err_up_[njets]         = Ntt.getErrorHi();
	Ntt_err_down_[njets]       = Ntt.getErrorLo();
	Nv_[njets]                 = Nv.getVal();
	Nv_err_[njets]             = Nv.getError();
	Nv_err_up_[njets]          = Nv.getErrorHi();
	Nv_err_down_[njets]        = Nv.getErrorLo();
	minValue_[njets]           = fit_result->minNll();
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		eb_[i][njets]        = eb[i]->getVal();
		eb_err_[i][njets]    = eb[i]->getError();
		eudsc_[i][njets]     = eudsc[i]->getVal();
		eudsc_err_[i][njets] = eudsc[i]->getError();
		euds_[i][njets]      = euds[i]->getVal();
		euds_err_[i][njets]  = euds[i]->getError();
	}
	//RooMCStudy* mcstudy = new RooMCStudy(model,nbjets,Binned(kTRUE),Silence(kTRUE),Extended(kTRUE),ProtoData(*data,kTRUE),FitOptions(Save(kFALSE),PrintEvalErrors(-1),Extended(kTRUE))) ;
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		delete eb[i];
		delete eudsc[i];
		delete euds[i];
		delete nbjets[i];
		delete pbjets_tt[i];
		delete pbjets_tt_ext[i];
		delete pbjets_v[i];
		delete pbjets_v_ext[i];
		delete model[i];
		delete data[i];
	}
	delete fit_result;
	delete combNLL;
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::UnBinnedMaximumNjetsLikelihoodEst(int btag_wp_idx, vector<int> &FixedVarIdx, bool doMinos, bool Verbose){

  // Declare observables
	char name[100];
	char title[100];
	RooRealVar *Ntt[NbOfJetsBins_];
	for(int i=0; i<NbOfJetsBins_; i++){
		sprintf(name,"Ntt_%d_jets",i+Njets_[0]);
		sprintf(title,"N_{t#bar{t}-like}^{%d jets}",i+Njets_[0]);
		Ntt[i] = new RooRealVar(name,title,init_[i][0],init_[i][0]/5,init_[i][0]*5);
	}
	RooRealVar *Nv[NbOfJetsBins_];
	for(int i=0; i<NbOfJetsBins_; i++){
		sprintf(name,"Nv_%d_jets",i+Njets_[0]);
		sprintf(title,"N_{V-like}^{%d jets}",i+Njets_[0]);
		Nv[i]  = new RooRealVar(name,title,init_[i][1],init_[i][1]/5,init_[i][1]*5);
	}
	vector<int>::iterator Idx;

	RooRealVar  eb(    "eb",  "#epsilon_{b-tag}",      init_[0][3*btag_wp_idx+2],0.0,1.0);
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 0 );
	if(Idx != FixedVarIdx.end()) eb.setConstant(kTRUE) ;
	RooRealVar  eudsc("eudsc","#epsilon_{mis-tag}",    init_[0][3*btag_wp_idx+3],0.0,1.0);
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 1 );
	if(Idx != FixedVarIdx.end()) eudsc.setConstant(kTRUE) ;
	RooRealVar euds( "euds",  "#epsilon^{,}_{mis-tag}",init_[0][3*btag_wp_idx+4],0.0,1.0);
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 2 );
	if(Idx != FixedVarIdx.end()) euds.setConstant(kTRUE) ;

  // Declare constants
	RooConstVar *n[NbOfJetsBins_];
	for(int i=0; i<NbOfJetsBins_; i++){
		sprintf(name,"n_%d_jets",i+Njets_[0]);
		n[i] = new RooConstVar(name,"",i+Njets_[0]);
	}
	RooConstVar *e0bq[NbOfJetsBins_];
	for(int i=0; i<NbOfJetsBins_; i++){
		sprintf(name,"e0bq_%d_jets",i+Njets_[0]);
		e0bq[i] = new RooConstVar(name,"",e0bq_[i]);
	}
	RooConstVar *e1bq[NbOfJetsBins_];
	for(int i=0; i<NbOfJetsBins_; i++){
		sprintf(name,"e1bq_%d_jets",i+Njets_[0]);
		e1bq[i] = new RooConstVar(name,"e1bq",e1bq_[i]);
	}
	RooConstVar *e2bq[NbOfJetsBins_];
	for(int i=0; i<NbOfJetsBins_; i++){
		sprintf(name,"e2bq_%d_jets",i+Njets_[0]);
		e2bq[i] = new RooConstVar(name,"e2bq",e2bq_[i]);
	}
	//RooConstVar e3bq("e3bq","e3bq",e3bq_[njets]);

  // C o n s t r u c t   a   c a t e g o r y   w i t h   l a b e l s    a n d   i n d e c e s
  // -----------------------------------------------------------------------------------------

	RooCategory* nbjets[NbOfJetsBins_];
	for(int i=0; i<NbOfJetsBins_; i++){
		sprintf(name,"nbjets_%d_jets",i+Njets_[0]);
		nbjets[i] = new RooCategory(name,"Number of b-jets");
		nbjets[i] -> defineType("N0bjet", 0);
		nbjets[i] -> defineType("N1bjet", 1);
		nbjets[i] -> defineType("N2bjets",2);
		nbjets[i] -> defineType("N3bjets",3);
	}

  // C o n s t r u c t   p . d . f 's
  // -------------------------------------------------------------------------------

	RooGenericPdf* pbjets_tt[NbOfJetsBins_];
	RooExtendPdf*  pbjets_tt_ext[NbOfJetsBins_];
	for(int i=0; i<NbOfJetsBins_; i++){
		sprintf(name,"pbjets_tt_%d_jets",i+Njets_[0]);
		pbjets_tt[i]     = new RooGenericPdf(name,"","(@0==0)*((1-@1)*(1-@1)*pow((1-@2),@3-2)*@6+(1-@1)*pow((1-@2),@3-1)*@5+pow((1-@2),@3)*@4)+(@0==1)*((2*@1*(1-@1)*pow(1-@2,@3-2)+(1-@1)*(1-@1)*(@3-2)*@2*pow(1-@2,@3-3))*@6+(@1*pow(1-@2,@3-1)+(1-@1)*(@3-1)*@2*pow(1-@2,@3-2))*@5+(@3*@2*pow(1-@2,@3-1))*@4)+(@0==2)*((@1*@1*pow(1-@2,@3-2)+2*@1*(1-@1)*(@3-2)*@2*pow(1-@2,@3-3)+(1-@1)*(1-@1)*((@3-2)*(@3-3)/2)*@2*@2*pow(1-@2,@3-4))*@6+(@1*(@3-1)*@2*pow(1-@2,@3-2)+(1-@1)*((@3-1)*(@3-2)/2)*@2*@2*pow(1-@2,@3-3))*@5+((@3*(@3-1)/2)*@2*@2*pow(1-@2,@3-2))*@4)+(@0==3)*((@1*@1*(@3-2)*@2*pow(1-@2,@3-3)+2*@1*(1-@1)*((@3-2)*(@3-3)/2)*@2*@2*pow(1-@2,@3-4)+(@3>4 ? pow((1-@1),2)*((@3-2)*(@3-3)*(@3-4)/6)*pow(@2,3)*pow((1-@2),@3-5) : 0 ))*@6+(@1*((@3-1)*(@3-2)/2)*@2*@2*pow(1-@2,@3-3)+(1-@1)*((@3-1)*(@3-2)*(@3-3)/6)*pow(@2,3)*pow(1-@2,@3-4))*@5+((@3*(@3-1)*(@3-2)/6)*pow(@2,3)*pow(1-@2,@3-3))*@4)",RooArgList(*nbjets[i],eb,eudsc,*n[i],*e0bq[i],*e1bq[i],*e2bq[i]));
		sprintf(name,"pbjets_tt_ex_%d_jets",i+Njets_[0]);
		pbjets_tt_ext[i] = new RooExtendPdf( name,"",*pbjets_tt[i],*Ntt[i]);
	}

	RooGenericPdf* pbjets_v[NbOfJetsBins_];
	RooExtendPdf*  pbjets_v_ext[NbOfJetsBins_];
	for(int i=0; i<NbOfJetsBins_; i++){
		sprintf(name,"pbjets_v_%d_jets",i+Njets_[0]);
		pbjets_v[i]     = new RooGenericPdf(name,"","(@0==0)*pow(1-@1,@2)+(@0==1)*@2*@1*pow(1-@1,@2-1)+(@0==2)*(@2*(@2-1)/2)*@1*@1*pow(1-@1,@2-2)+(@0==3)*((@2)*(@2-1)*(@2-2)/6)*pow(@1,3)*pow(1-@1,@2-3)",RooArgList(*nbjets[i],euds,*n[i]));
		sprintf(name,"pbjets_v_ext_%d_jets",i+Njets_[0]);
		pbjets_v_ext[i] = new RooExtendPdf(name,"",*pbjets_v[i],*Nv[i]);
	}

	RooAddPdf* model[NbOfJetsBins_];
	for(int i=0; i<NbOfJetsBins_; i++){
		sprintf(name,"model_%d_jets",i+Njets_[0]);
		model[i] = new RooAddPdf(name,"",RooArgList(*pbjets_tt[i],*pbjets_v[i]),RooArgList(*Ntt[i],*Nv[i]));
	}
	
  // C r e a t e   d a t a s e t 
  // -------------------------------------------------------------------------------
  // Sample a dataset for tt+jets events

	RooDataSet* data[NbOfJetsBins_];
	for(int i=0; i<NbOfJetsBins_; i++){
		sprintf(name,"data_%d_jets",i+Njets_[0]);
		data[i] = new RooDataSet(name,"",RooArgSet(*nbjets[i])) ;
		for (int j=0 ; j<GetPredNtotal(btag_wp_idx,i,0) ; j++) { nbjets[i]->setLabel("N0bjet")  ; data[i]->add(RooArgSet(*nbjets[i]));}
		for (int j=0 ; j<GetPredNtotal(btag_wp_idx,i,1) ; j++) { nbjets[i]->setLabel("N1bjet")  ; data[i]->add(RooArgSet(*nbjets[i]));}
		for (int j=0 ; j<GetPredNtotal(btag_wp_idx,i,2) ; j++) { nbjets[i]->setLabel("N2bjets") ; data[i]->add(RooArgSet(*nbjets[i]));}
		for (int j=0 ; j<GetPredNtotal(btag_wp_idx,i,3) ; j++) { nbjets[i]->setLabel("N3bjets") ; data[i]->add(RooArgSet(*nbjets[i]));}
	}

  // F i t   t h e   d a t a   a n d   c o n s t r u c t   t h e   l i k e l i h o o d   f u n c t i o n
  // ----------------------------------------------------------------------------------------------

	//RooFitResult* fit_result = model.fitTo(data,Save(),Extended(1),PrintLevel(-1));//,ConditionalObservables(nbjets),SplitRange(1),Optimize(0));
	RooAbsReal* nll[NbOfJetsBins_];
	RooArgSet nlls("myNLLs");
	for(int i=0; i<NbOfJetsBins_; i++){
		nll[i] = model[i]->createNLL(*data[i]);
		nlls.add(*nll[i]);
	}
	
	RooAddition* combNLL = new RooAddition("combNLL","combined likelihood function",nlls);
	
	RooMinuit minimizer(*combNLL);
	
	minimizer.optimizeConst(kTRUE) ;
	minimizer.setPrintLevel(-1);
	minimizer.setNoWarn();
	
	minimizer.migrad();
	
	minimizer.hesse();
	
	if(doMinos) minimizer.minos();
	
	RooFitResult* fit_result = minimizer.save();

	if(Verbose)fit_result->Print("v");

	for(int i=0; i<NbOfJetsBins_; i++){
		Ntt_[i]                = Ntt[i]->getVal();
		Ntt_err_[i]            = Ntt[i]->getError();
		Ntt_err_up_[i]         = Ntt[i]->getErrorHi();
		Ntt_err_down_[i]       = Ntt[i]->getErrorLo();
		Nv_[i]                 = Nv[i]->getVal();
		Nv_err_[i]             = Nv[i]->getError();
		Nv_err_up_[i]          = Nv[i]->getErrorHi();
		Nv_err_down_[i]        = Nv[i]->getErrorLo();

		minValue_[i]           = fit_result->minNll();

		eb_[btag_wp_idx][i]        = eb.getVal();
		eb_err_[btag_wp_idx][i]    = eb.getError();
		eudsc_[btag_wp_idx][i]     = eudsc.getVal();
		eudsc_err_[btag_wp_idx][i] = eudsc.getError();
		euds_[btag_wp_idx][i]      = euds.getVal();
		euds_err_[btag_wp_idx][i]  = euds.getError();
	}
	
	delete fit_result;
	delete combNLL;
	
}

/**________________________________________________________________________________________________________________*/
/*
void VJetEstimation::UnBinnedMaximumJointLikelihoodEst(bool Verbose){

  // Declare observables
	TString wp[3] = {"loose","medium","tight"};
	RooRealVar* eb[NbOfBtagWorkingPoint_];
	RooRealVar* eudsc[NbOfBtagWorkingPoint_];
	RooRealVar* euds[NbOfBtagWorkingPoint_];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		eb[i]    = new RooRealVar(   "eb_"+wp[i],"#epsilon_{b-tag} ("+wp[i]+")",      GetPredEb(   i),0.0,1.0);
		eudsc[i] = new RooRealVar("eudsc_"+wp[i],"#epsilon_{mis-tag} ("+wp[i]+")",    GetPredEudsc(i),0.0,1.0);
		euds[i]  = new RooRealVar( "euds_"+wp[i],"#epsilon^{,}_{mis-tag} ("+wp[i]+")",GetPredEuds( i),0.0,1.0);
	}

	RooRealVar* Ntt[NbOfJetsBins_];
	RooRealVar* Nv[NbOfJetsBins_];
	for(int i=0; i<NbOfJetsBins_; i++){
		Ntt[i] = new RooRealVar("Ntt","N_{t#bar{t}-like}",GetPredNtt(0,i),GetPredNtt(0,i)/5,GetPredNtt(0,i)*5);
		Nv[i]  = new RooRealVar( "Nv","N_{V-like}",       GetPredNv( 0,i),GetPredNv(0,i)/5, GetPredNv( 0,i)*5);
	}

  // Declare constants

	RooConstVar n("n","number of selected jets",Njets_[njets]) ;
	RooConstVar e0bq("e0bq","e0bq",e0bq_[njets]);
	RooConstVar e1bq("e1bq","e1bq",e1bq_[njets]);
	RooConstVar e2bq("e2bq","e2bq",e2bq_[njets]);
	//RooConstVar e3bq("e3bq","e3bq",e3bq_[njets]);

  // C o n s t r u c t   a   c a t e g o r y   w i t h   l a b e l s    a n d   i n d e c e s
  // -----------------------------------------------------------------------------------------

	RooCategory* nbjets[NbOfBtagWorkingPoint_];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		nbjets[i] = new RooCategory("nbjets_"+wp[i],"Number of b-jets");
		nbjets[i] -> defineType("N0bjet", 0);
		nbjets[i] -> defineType("N1bjet", 1);
		nbjets[i] -> defineType("N2bjets",2);
		nbjets[i] -> defineType("N3bjets",3);
	}

  // C o n s t r u c t   p . d . f 's
  // -------------------------------------------------------------------------------

	RooGenericPdf* pbjets_tt[NbOfBtagWorkingPoint_];
	RooExtendPdf*  pbjets_tt_ext[NbOfBtagWorkingPoint_];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		//pbjets_tt[i]     = new RooGenericPdf("pbjets_tt_"+wp[i],"pbjets_tt_"+wp[i],"((*nbjets[i])==0)*((1-(*eb[i]))*(1-(*eb[i]))*pow((1-(*eudsc[i])),n-2)*e2bq+(1-(*eb[i]))*pow((1-(*eudsc[i])),n-1)*e1bq+pow((1-(*eudsc[i])),n)*e0bq)+((*nbjets[i])==1)*((2*(*eb[i])*(1-(*eb[i]))*pow(1-(*eudsc[i]),n-2)+(1-(*eb[i]))*(1-(*eb[i]))*(n-2)*(*eudsc[i])*pow(1-(*eudsc[i]),n-3))*e2bq+((*eb[i])*pow(1-(*eudsc[i]),n-1)+(1-(*eb[i]))*(n-1)*(*eudsc[i])*pow(1-(*eudsc[i]),n-2))*e1bq+(n*(*eudsc[i])*pow(1-(*eudsc[i]),n-1))*e0bq)+((*nbjets[i])==2)*(((*eb[i])*(*eb[i])*pow(1-(*eudsc[i]),n-2)+2*(*eb[i])*(1-(*eb[i]))*(n-2)*(*eudsc[i])*pow(1-(*eudsc[i]),n-3)+(1-(*eb[i]))*(1-(*eb[i]))*((n-2)*(n-3)/2)*(*eudsc[i])*(*eudsc[i])*pow(1-(*eudsc[i]),n-4))*e2bq+((*eb[i])*(n-1)*(*eudsc[i])*pow(1-(*eudsc[i]),n-2)+(1-(*eb[i]))*((n-1)*(n-2)/2)*(*eudsc[i])*(*eudsc[i])*pow(1-(*eudsc[i]),n-3))*e1bq+((n*(n-1)/2)*(*eudsc[i])*(*eudsc[i])*pow(1-(*eudsc[i]),n-2))*e0bq)+((*nbjets[i])==3)*(((*eb[i])*(*eb[i])*(n-2)*(*eudsc[i])*pow(1-(*eudsc[i]),n-3)+2*(*eb[i])*(1-(*eb[i]))*((n-2)*(n-3)/2)*(*eudsc[i])*(*eudsc[i])*pow(1-(*eudsc[i]),n-4)+(n>4 ? pow((1-(*eb[i])),2)*((n-2)*(n-3)*(n-4)/6)*pow((*eudsc[i]),3)*pow((1-(*eudsc[i])),n-5) : 0 ))*e2bq+((*eb[i])*((n-1)*(n-2)/2)*(*eudsc[i])*(*eudsc[i])*pow(1-(*eudsc[i]),n-3)+(1-(*eb[i]))*((n-1)*(n-2)*(n-3)/6)*pow((*eudsc[i]),3)*pow(1-(*eudsc[i]),n-4))*e1bq+((n*(n-1)*(n-2)/6)*pow((*eudsc[i]),3)*pow(1-(*eudsc[i]),n-3))*e0bq)",RooArgList(*nbjets[i],n,*eb[i],*eudsc[i],e0bq,e1bq,e2bq));
		pbjets_tt[i]     = new RooGenericPdf("pbjets_tt_"+wp[i],"pbjets_tt_"+wp[i],"(@0==0)*((1-@1)*(1-@1)*pow((1-@2),@3-2)*@6+(1-@1)*pow((1-@2),@3-1)*@5+pow((1-@2),@3)*@4)+(@0==1)*((2*@1*(1-@1)*pow(1-@2,@3-2)+(1-@1)*(1-@1)*(@3-2)*@2*pow(1-@2,@3-3))*@6+(@1*pow(1-@2,@3-1)+(1-@1)*(@3-1)*@2*pow(1-@2,@3-2))*@5+(@3*@2*pow(1-@2,@3-1))*@4)+(@0==2)*((@1*@1*pow(1-@2,@3-2)+2*@1*(1-@1)*(@3-2)*@2*pow(1-@2,@3-3)+(1-@1)*(1-@1)*((@3-2)*(@3-3)/2)*@2*@2*pow(1-@2,@3-4))*@6+(@1*(@3-1)*@2*pow(1-@2,@3-2)+(1-@1)*((@3-1)*(@3-2)/2)*@2*@2*pow(1-@2,@3-3))*@5+((@3*(@3-1)/2)*@2*@2*pow(1-@2,@3-2))*@4)+(@0==3)*((@1*@1*(@3-2)*@2*pow(1-@2,@3-3)+2*@1*(1-@1)*((@3-2)*(@3-3)/2)*@2*@2*pow(1-@2,@3-4)+(@3>4 ? pow((1-@1),2)*((@3-2)*(@3-3)*(@3-4)/6)*pow(@2,3)*pow((1-@2),@3-5) : 0 ))*@6+(@1*((@3-1)*(@3-2)/2)*@2*@2*pow(1-@2,@3-3)+(1-@1)*((@3-1)*(@3-2)*(@3-3)/6)*pow(@2,3)*pow(1-@2,@3-4))*@5+((@3*(@3-1)*(@3-2)/6)*pow(@2,3)*pow(1-@2,@3-3))*@4)",RooArgList(*nbjets[i],*eb[i],*eudsc[i],n,e0bq,e1bq,e2bq));
		pbjets_tt_ext[i] = new RooExtendPdf( "pbjets_tt_ext_"+wp[i],"pbjets_tt_ext_"+wp[i],*pbjets_tt[i],Ntt);
	}

	RooGenericPdf* pbjets_v[NbOfBtagWorkingPoint_];
	RooExtendPdf*  pbjets_v_ext[NbOfBtagWorkingPoint_];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		//pbjets_v[i]     = new RooGenericPdf("pbjets_v_"+wp[i],"pbjets_v_"+wp[i],"((*nbjets[i])==0)*pow(1-(*euds[i]),n)+((*nbjets[i])==1)*n*(*euds[i])*pow(1-(*euds[i]),n-1)+((*nbjets[i])==2)*(n*(n-1)/2)*(*euds[i])*(*euds[i])*pow(1-(*euds[i]),n-2)+((*nbjets[i])==3)*((n)*(n-1)*(n-2)/6)*pow((*euds[i]),3)*pow(1-(*euds[i]),n-3)",RooArgList(*nbjets[i],*euds[i],n));
		pbjets_v[i]     = new RooGenericPdf("pbjets_v_"+wp[i],"pbjets_v_"+wp[i],"(@0==0)*pow(1-@1,@2)+(@0==1)*@2*@1*pow(1-@1,@2-1)+(@0==2)*(@2*(@2-1)/2)*@1*@1*pow(1-@1,@2-2)+(@0==3)*((@2)*(@2-1)*(@2-2)/6)*pow(@1,3)*pow(1-@1,@2-3)",RooArgList(*nbjets[i],*euds[i],n));
		pbjets_v_ext[i] = new RooExtendPdf("pbjets_v_ext_"+wp[i],"pbjets_v_ext_"+wp[i],*pbjets_v[i],Nv);
	}

	RooAddPdf* model[NbOfBtagWorkingPoint_];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		model[i] = new RooAddPdf("model"+wp[i],"model"+wp[i],RooArgList(*pbjets_tt[i],*pbjets_v[i]),RooArgList(Ntt,Nv));
	}
	
  // C r e a t e   d a t a s e t 
  // -------------------------------------------------------------------------------

	RooDataSet* data[NbOfBtagWorkingPoint_];
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		data[i] = new RooDataSet("data"+wp[i],"data"+wp[i],RooArgSet(*nbjets[i])) ;
		for (int j=0 ; j<GetPredNtotal(i,njets,0) ; j++) { nbjets[i]->setLabel("N0bjet")  ; data[i]->add(RooArgSet(*nbjets[i]));}
		for (int j=0 ; j<GetPredNtotal(i,njets,1) ; j++) { nbjets[i]->setLabel("N1bjet")  ; data[i]->add(RooArgSet(*nbjets[i]));}
		for (int j=0 ; j<GetPredNtotal(i,njets,2) ; j++) { nbjets[i]->setLabel("N2bjets") ; data[i]->add(RooArgSet(*nbjets[i]));}
		for (int j=0 ; j<GetPredNtotal(i,njets,3) ; j++) { nbjets[i]->setLabel("N3bjets") ; data[i]->add(RooArgSet(*nbjets[i]));}
	}

  // F i t   t h e   d a t a   a n d   c o n s t r u c t   t h e   l i k e l i h o o d   f u n c t i o n
  // ----------------------------------------------------------------------------------------------

	RooAbsReal* nll[NbOfBtagWorkingPoint_];
	RooArgSet nlls("myNLLs");
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		nll[i] = model[i]->createNLL(*data[i]);
		nlls.add(*nll[i]);
	}
	
	RooAddition* combNLL = new RooAddition("combNLL","combined likelihood function",nlls);
	
	RooMinuit minimizer(*combNLL);
	
	minimizer.setPrintLevel(-1);
	
	minimizer.migrad();
	
	minimizer.hesse();
	
	RooFitResult* fit_result = minimizer.save();

	if(Verbose)fit_result->Print("v");

	Ntt_[njets]                = Ntt.getVal();
	Nv_[njets]                 = Nv.getVal();
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		eb_[i][njets]    = eb[i]->getVal();
		eudsc_[i][njets] = eudsc[i]->getVal();
		euds_[i][njets]  = euds[i]->getVal();
	}
	//RooMCStudy* mcstudy = new RooMCStudy(model,nbjets,Binned(kTRUE),Silence(kTRUE),Extended(kTRUE),ProtoData(*data,kTRUE),FitOptions(Save(kFALSE),PrintEvalErrors(-1),Extended(kTRUE))) ;
	for(int i=0; i<NbOfBtagWorkingPoint_; i++){
		delete eb[i];
		delete eudsc[i];
		delete euds[i];
		delete nbjets[i];
		delete pbjets_tt[i];
		delete pbjets_tt_ext[i];
		delete pbjets_v[i];
		delete pbjets_v_ext[i];
		delete model[i];
		delete data[i];
	}
	delete combNLL;
}
*/

/**________________________________________________________________________________________________________________*/
void VJetEstimation::PrintInputs(){
	for(int i=0;i<NbOfJetsBins_;i++) PrintInputs(i);
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::PrintResults(){
	cout<<endl;
	cout<<"***********************************************"<<endl;
	cout<<"*****      VJets Estimation Results         ***"<<endl;
	cout<<"***********************************************"<<endl;

	for(int i=0;i<NbOfJetsBins_;i++) PrintResults(i);

	cout<<"***********************************************************"<<endl;
	cout<<"************        ALL INCLUSIF             **************"<<endl; 
	cout<<setprecision(4);
	cout<<"Btag working point : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetBtagWorkingPoint(i)<<" || ";} cout<<endl;
	cout<<"Estimated/Predicted eb	 : ";	for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstEb(i)   <<" +/- "<<GetEstEbErr(i)	<<" / "<<GetPredEb(i)   <<" || ";} cout<<endl;
	cout<<"Estimated/Predicted eudsc : ";	for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstEudsc(i)<<" +/- "<<GetEstEudscErr(i)	<<" / "<<GetPredEudsc(i)<<" || ";} cout<<endl;
	cout<<"Estimated/Predicted euds  : ";	for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstEuds(i) <<" +/- "<<GetEstEudsErr(i)	<<" / "<<GetPredEuds(i) <<" || ";} cout<<endl;
	cout<<setprecision(1);
	cout<<"Estimated Ntt-like : "<<GetEstNtt(0)<<" / Predicted value from MC : "<< GetPredNtt(0)<<endl;
	cout<<"- 0 b-jet          : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNtt(i,-1,0)<<" +/- "<<GetEstNttErr(i,-1,0)<<" / "<< GetPredNtt(i,-1,0)<<" +/- "<<GetPredNttErr(i,-1,0) <<" || ";} cout<<endl;
	cout<<"- 1 b-jet          : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNtt(i,-1,1)<<" +/- "<<GetEstNttErr(i,-1,1)<<" / "<< GetPredNtt(i,-1,1)<<" +/- "<<GetPredNttErr(i,-1,1) <<" || ";} cout<<endl;
	cout<<"- 2 b-jets         : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNtt(i,-1,2)<<" +/- "<<GetEstNttErr(i,-1,2)<<" / "<< GetPredNtt(i,-1,2)<<" +/- "<<GetPredNttErr(i,-1,2) <<" || ";} cout<<endl;
	cout<<"- 3 b-jets         : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNtt(i,-1,3)<<" +/- "<<GetEstNttErr(i,-1,3)<<" / "<< GetPredNtt(i,-1,3)<<" +/- "<<GetPredNttErr(i,-1,3) <<" || ";} cout<<endl;
// 	cout<<"Estimated Nvb-like : "<<GetEstNvb(0)<<" / Predicted value from MC : "<< GetPredNvb(0)<<endl;
// 	cout<<"- 0 b-jet          : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNvb(i,-1,0)<<" / "<< GetPredNvb(i,-1,0) <<" || ";} cout<<endl;
// 	cout<<"- 1 b-jet          : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNvb(i,-1,1)<<" / "<< GetPredNvb(i,-1,1) <<" || ";} cout<<endl;
// 	cout<<"- 2 b-jets         : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNvb(i,-1,2)<<" / "<< GetPredNvb(i,-1,2) <<" || ";} cout<<endl;
// 	cout<<"- 3 b-jets         : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNvb(i,-1,3)<<" / "<< GetPredNvb(i,-1,3) <<" || ";} cout<<endl;
	cout<<"Estimated Nv-like  : "<<GetEstNv(0)<<" / Predicted value from MC : "<< GetPredNv(0)<<endl;
	cout<<"- 0 b-jet          : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNv(i,-1,0)<<" +/- "<<GetEstNvErr(i,-1,0)<<" / "<<GetPredNv(i,-1,0)<<" +/- "<<GetPredNvErr(i,-1,0) <<" || ";} cout<<endl;
	cout<<"- 1 b-jet          : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNv(i,-1,1)<<" +/- "<<GetEstNvErr(i,-1,1)<<" / "<<GetPredNv(i,-1,1)<<" +/- "<<GetPredNvErr(i,-1,1) <<" || ";} cout<<endl;
	cout<<"- 2 b-jets         : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNv(i,-1,2)<<" +/- "<<GetEstNvErr(i,-1,2)<<" / "<<GetPredNv(i,-1,2)<<" +/- "<<GetPredNvErr(i,-1,2) <<" || ";} cout<<endl;
	cout<<"- 3 b-jets         : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNv(i,-1,3)<<" +/- "<<GetEstNvErr(i,-1,3)<<" / "<<GetPredNv(i,-1,3)<<" +/- "<<GetPredNvErr(i,-1,3) <<" || ";} cout<<endl;
	cout<<"***********************************************************"<<endl;
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::PrintResults_LatexFormat(ofstream &ofile){
for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {
	ofile<<"Btag working point : "<<GetBtagWorkingPoint(i)<<endl;
        ofile<<"\\begin{table}"<<endl;
        ofile<<"	\\centering"<<endl;
        ofile<<"		\\begin{tabular}{l|";for(int j=0;j<NbOfJetsBins_+1;j++){ofile<<"c";};ofile<<"}"<<endl;
        ofile<<"			\\hline"<<endl;
        ofile<<" \\multicolumn{"<<NbOfJetsBins_+1<<"}{l}{Estimation / Monte Carlo prediction : }\\\\"<<endl;
        ofile<<"			      &";for(int j=0;j<NbOfJetsBins_;j++){ofile<<"$"<<Njets_[j]<<"$ jets & ";}ofile<<"Inclusive \\\\"<<endl;
        ofile<<"\\hline"<<endl;
        ofile<<" $\\epsilon_b$                &$"<<setprecision(3);for(int j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstEb(   i,j)<<"\\pm"<<GetEstEbErr(   i,j)<<" $/$ "<<GetPredEb(   i,j)<<"\\pm"<<GetPredEbErr(   i,j)<<"$&$";}ofile<<GetEstEb(    i)<<"\\pm"<<GetEstEbErr(   i)<<"$/$"<<GetPredEb(   i)<<"\\pm"<<GetPredEbErr(   i)<<"$ \\\\"<<endl;
        ofile<<" $\\epsilon_{eudsc}$          &$"<<setprecision(4);for(int j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstEudsc(i,j)<<"\\pm"<<GetEstEudscErr(i,j)<<" $/$ "<<GetPredEudsc(i,j)<<"\\pm"<<GetPredEudscErr(i,j)<<"$&$";}ofile<<GetEstEudsc( i)<<"\\pm"<<GetEstEudscErr(i)<<"$/$"<<GetPredEudsc(i)<<"\\pm"<<GetPredEudscErr(i)<<"$ \\\\"<<endl;
        ofile<<" $\\epsilon_{euds} $          &$"<<setprecision(4);for(int j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstEuds( i,j)<<"\\pm"<<GetEstEudsErr( i,j)<<" $/$ "<<GetPredEuds( i,j)<<"\\pm"<<GetPredEudsErr( i,j)<<"$&$";}ofile<<GetEstEuds(  i)<<"\\pm"<<GetEstEudsErr( i)<<"$/$"<<GetPredEuds( i)<<"\\pm"<<GetPredEudsErr( i)<<"$ \\\\"<<endl;
	cout.precision(1);        
	ofile<<"\\hline"<<endl;
        ofile<<" \\multicolumn{"<<NbOfJetsBins_+1<<"}{l}{$t\\bar{t}+jets$ and single top : }\\\\"<<endl;
        ofile<<"\\hline"<<endl;
        ofile<<" $0$ b-jet                    &$";for(int j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNtt(i,j,0)<<"\\pm"<<GetEstNttErr(i,j,0)<<" $/$ "<<GetPredNtt(i,j,0)<<"\\pm"<<GetPredNttErr(i,j,0) <<"$&$";}ofile<<GetEstNtt(i,-1,0)<<"\\pm"<<GetEstNttErr(i,-1,0)<<" $/$ "<<GetPredNtt(i,-1,0)<<"\\pm"<<GetPredNttErr(i,-1,0)<<"$ \\\\"<<endl;
        ofile<<" $1$ b-jet                    &$";for(int j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNtt(i,j,1)<<"\\pm"<<GetEstNttErr(i,j,1)<<" $/$ "<<GetPredNtt(i,j,1)<<"\\pm"<<GetPredNttErr(i,j,1) <<"$&$";}ofile<<GetEstNtt(i,-1,1)<<"\\pm"<<GetEstNttErr(i,-1,1)<<" $/$ "<<GetPredNtt(i,-1,1)<<"\\pm"<<GetPredNttErr(i,-1,1)<<"$ \\\\"<<endl;
        ofile<<" $2$ b-jets                   &$";for(int j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNtt(i,j,2)<<"\\pm"<<GetEstNttErr(i,j,2)<<" $/$ "<<GetPredNtt(i,j,2)<<"\\pm"<<GetPredNttErr(i,j,2) <<"$&$";}ofile<<GetEstNtt(i,-1,2)<<"\\pm"<<GetEstNttErr(i,-1,2)<<" $/$ "<<GetPredNtt(i,-1,2)<<"\\pm"<<GetPredNttErr(i,-1,2)<<"$ \\\\"<<endl;
        ofile<<" $3$ b-jets                   &$";for(int j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNtt(i,j,3)<<"\\pm"<<GetEstNttErr(i,j,3)<<" $/$ "<<GetPredNtt(i,j,3)<<"\\pm"<<GetPredNttErr(i,j,3) <<"$&$";}ofile<<GetEstNtt(i,-1,3)<<"\\pm"<<GetEstNttErr(i,-1,3)<<" $/$ "<<GetPredNtt(i,-1,3)<<"\\pm"<<GetPredNttErr(i,-1,3)<<"$ \\\\"<<endl;
        ofile<<"\\hline"<<endl;
        ofile<<" Inclusive                    &$";for(int j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNtt(i,j)  <<"\\pm"<<GetEstNttErr(i,j)  <<" $/$ "<<GetPredNtt(i,j)  <<"\\pm"<<GetPredNttErr(i,j)   <<"$&$";}ofile<<GetEstNtt(i)     <<"\\pm"<<GetEstNttErr(i)     <<" $/$ "<<GetPredNtt(i)     <<"\\pm"<<GetPredNttErr(i)<<"$ \\\\"<<endl;
        ofile<<"\\hline"<<endl;
        ofile<<"\\hline"<<endl;
        ofile<<" \\multicolumn{"<<NbOfJetsBins_+1<<"}{l}{Multijet and $W/Z+jets$ : }\\\\"<<endl;
        ofile<<"\\hline"<<endl;
        ofile<<" $0$ b-jet                    &$";for(int j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNv(i,j,0)<<"\\pm"<<GetEstNvErr(i,j,0)<<" $/$ "<<GetPredNv(i,j,0)<<"\\pm"<<GetPredNvErr(i,j,0) <<"$&$";}ofile<<GetEstNv(i,-1,0)<<"\\pm"<<GetEstNvErr(i,-1,0)<<" $/$ "<<GetPredNv(i,-1,0)<<"\\pm"<<GetPredNvErr(i,-1,0)<<"$ \\\\"<<endl;
        ofile<<" $1$ b-jet                    &$";for(int j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNv(i,j,1)<<"\\pm"<<GetEstNvErr(i,j,1)<<" $/$ "<<GetPredNv(i,j,1)<<"\\pm"<<GetPredNvErr(i,j,1) <<"$&$";}ofile<<GetEstNv(i,-1,1)<<"\\pm"<<GetEstNvErr(i,-1,1)<<" $/$ "<<GetPredNv(i,-1,1)<<"\\pm"<<GetPredNvErr(i,-1,1)<<"$ \\\\"<<endl;
        ofile<<" $2$ b-jets                   &$";for(int j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNv(i,j,2)<<"\\pm"<<GetEstNvErr(i,j,2)<<" $/$ "<<GetPredNv(i,j,2)<<"\\pm"<<GetPredNvErr(i,j,2) <<"$&$";}ofile<<GetEstNv(i,-1,2)<<"\\pm"<<GetEstNvErr(i,-1,2)<<" $/$ "<<GetPredNv(i,-1,2)<<"\\pm"<<GetPredNvErr(i,-1,2)<<"$ \\\\"<<endl;
        ofile<<" $3$ b-jets                   &$";for(int j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNv(i,j,3)<<"\\pm"<<GetEstNvErr(i,j,3)<<" $/$ "<<GetPredNv(i,j,3)<<"\\pm"<<GetPredNvErr(i,j,3) <<"$&$";}ofile<<GetEstNv(i,-1,3)<<"\\pm"<<GetEstNvErr(i,-1,3)<<" $/$ "<<GetPredNv(i,-1,3)<<"\\pm"<<GetPredNvErr(i,-1,3)<<"$ \\\\"<<endl;
        ofile<<"\\hline"<<endl;
        ofile<<" Inclusive                    &$";for(int j=0;j<NbOfJetsBins_;j++) {ofile<<GetEstNv(i,j)  <<"\\pm"<<GetEstNvErr(i,j)  <<" $/$ "<<GetPredNv(i,j)  <<"\\pm"<<GetPredNvErr(i,j)   <<"$&$";}ofile<<GetEstNv(i)     <<"\\pm"<<GetEstNvErr(i)     <<" $/$ "<<GetPredNv(i)     <<"\\pm"<<GetPredNvErr(i)<<"$ \\\\"<<endl;
        ofile<<"\\hline"<<endl;
        ofile<<"\\hline"<<endl;
        ofile<<"		\\end{tabular}"<<endl;
        ofile<<"	\\caption{}"<<endl;
        ofile<<"	\\label{tab:}"<<endl;
        ofile<<"\\end{table}"<<endl;
	}
}

/**________________________________________________________________________________________________________________*/
bool VJetEstimation::CheckEstimation(double threshold, int njets){ return (minValue_[njets]>threshold ? false : true);}

/**________________________________________________________________________________________________________________*/
bool VJetEstimation::CheckEstimation(double threshold){
	for(int i=0;i<NbOfJetsBins_;i++)
		if(minValue_[i]>threshold) return false;
	return true;
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::CheckEstimationLinearity(int NbOfRescaleFact, double **RescaleFact, int Idx){

	RescaledTTLikeEstimation = new TGraphErrors(NbOfRescaleFact);
	RescaledVLikeEstimation  = new TGraphErrors(NbOfRescaleFact);
	tCanva_RescaledTTLikeEstimation = new TCanvas("tCanva_RescaledTTLikeEstimation","",-1);
	tCanva_RescaledVLikeEstimation  = new TCanvas("tCanva_RescaledVLikeEstimation", "",-1);

	vector<int> Indeces;
	double ****n_Init  = new double***[NbOfBtagWorkingPoint_];
	double ****n       = new double***[NbOfBtagWorkingPoint_];
	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		n_Init[i] = new double**[NbOfJetsBins_];
		n[i]      = new double**[NbOfJetsBins_];
		for(int j=0;j<NbOfJetsBins_;j++){
			n_Init[i][j] = new double*[NbOfDatasets_];
			n[i][j]      = new double*[NbOfDatasets_];
			for(int k=0;k<NbOfBJetsBins_;k++){
				n_Init[i][j][k] = new double[NbOfDatasets_];
				n[i][j][k]      = new double[NbOfDatasets_];
				for(int l=0;l<NbOfDatasets_;l++) n_Init[i][j][k][l] = N_[i][j][k][l];
			}			
		}
	}
	cout<<"*****************Linearity check : results *****************"<<endl;
	for(int nb=0;nb<NbOfRescaleFact;nb++){
		cout<<" --- ";
		for(int i=0;i<NbOfJetsBins_;i++){
			for(int j=0;j<NbOfBtagWorkingPoint_;j++){
				for(int k=0;k<NbOfBJetsBins_;k++){
					for(int l=0;l<NbOfDatasets_;l++){
					/// Multiply the original numbers of events (per datasets) 
					/// by a given factor "RescalFact"
					n[i][j][k][l] = n_Init[i][j][k][l]*RescaleFact[nb][l];
					}
				}
			}
		}
		//Estimate the number of tt-like and w-like events with the modified input numbers
		string method = "Minuit2";
		string option = "Combined";

		this->FillInputs(n);
		this->BinnedMaximumLikelihoodEst(method, option, Indeces, false, false);

		RescaledTTLikeEstimation->SetPoint(nb,100*RescaleFact[nb][Idx],(double)this->GetEstNtt(0));
		RescaledTTLikeEstimation->SetPointError(nb,0,RescaleFact[nb][Idx]*((double)this->GetNttVar(0)));
		RescaledVLikeEstimation->SetPoint(nb,100*RescaleFact[nb][Idx],(double)this->GetEstNv(0));
		RescaledVLikeEstimation->SetPointError(nb,0, RescaleFact[nb][Idx]*((double)this->GetNvVar(0)));
		RescaledVbLikeEstimation->SetPoint(nb,100*RescaleFact[nb][Idx],(double)this->GetEstNvb(0));
		RescaledVbLikeEstimation->SetPointError(nb,0, RescaleFact[nb][Idx]*((double)this->GetNvbVar(0)));

		cout<<"-- Results :"<<endl;
		cout<<"--- V-like estimate  = "<<this->GetEstNv(0) <<"+/- "<<RescaleFact[nb][Idx]*((double)this->GetNvVar(0))<<endl;
		cout<<"--- Vb-like estimate = "<<this->GetEstNvb(0)<<"+/- "<<RescaleFact[nb][Idx]*((double)this->GetNvbVar(0))<<endl;
		cout<<"--- tt-like estimate = "<<this->GetEstNtt(0)<<"+/- "<<RescaleFact[nb][Idx]*((double)this->GetNttVar(0))<<endl;
		cout<<"------------"<<endl;
	}

	RescaledTTLikeEstimation ->SetNameTitle("RescaledTTLikeEstimation","");
	RescaledVLikeEstimation  ->SetNameTitle("RescaledVLikeEstimation","");

	tCanva_RescaledTTLikeEstimation->cd();
	RescaledTTLikeEstimation->Draw("AC*");
	RescaledTTLikeEstimation->GetXaxis()->SetTitle("Rescaling factor (%)");
	RescaledTTLikeEstimation->GetYaxis()->SetTitle("Estimated nb of TT-like events");

	tCanva_RescaledVLikeEstimation->cd();
	RescaledVLikeEstimation->Draw("AC*");
	RescaledVLikeEstimation->GetXaxis()->SetTitle("Rescaling factor (%)");
	RescaledVLikeEstimation->GetYaxis()->SetTitle("Estimated nb of W-like events");

	// Put the number of events back to their original values
	this->FillInputs(n_Init);

	for(int i=0;i<NbOfBtagWorkingPoint_;i++){
		for(int j=0;j<NbOfJetsBins_;j++){
			for(int k=0;k<NbOfBJetsBins_;k++){
				delete [] n_Init[i][j][k];
				delete [] n[i][j][k];
			}
			delete [] n_Init[i][j];
			delete [] n[i][j];
		}
		delete[] n_Init[i];
		delete[] n[i];
	}
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredEb(int wp) const{
	double Effb = 0;
	for(int i=0;i<NbOfJetsBins_;i++) Effb += GetPredEb(wp,i)*GetPredNtotal(wp,i);
	Effb /= GetPredNtotal(wp);
	return Effb;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredEb(int wp, int njets) const{
	return (njets<NbOfJetsBins_? eb_mc_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredEbErr(int wp) const{
	double EffbErr = 0;
	for(int i=0;i<NbOfJetsBins_;i++) EffbErr += pow(GetPredEbErr(wp,i)*(GetPredNtotal(wp,i)/GetPredNtotal(wp)),2);
	return sqrt(EffbErr);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredEbErr(int wp, int njets) const{
	return (njets<NbOfJetsBins_? eb_err_mc_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetPredEb(int wp, int njets, double value) const{
	eb_mc_[wp][njets] = value;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredEudsc(int wp) const{
	double Effudsc = 0;
	for(int i=0;i<NbOfJetsBins_;i++) Effudsc += GetPredEudsc(wp,i)*GetPredNtotal(wp,i);
	Effudsc /= GetPredNtotal(wp);
	return Effudsc;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredEudsc(int wp, int njets) const{
	return (njets<NbOfJetsBins_? eudsc_mc_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredEudscErr(int wp) const{
	double EffudscErr = 0;
	for(int i=0;i<NbOfJetsBins_;i++) EffudscErr += pow(GetPredEudscErr(wp,i)*(GetPredNtotal(wp,i)/GetPredNtotal(wp)),2);
	return sqrt(EffudscErr);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredEudscErr(int wp, int njets) const{
	return (njets<NbOfJetsBins_? eudsc_err_mc_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetPredEudsc(int wp, int njets, double value) const{
	eudsc_mc_[wp][njets] = value;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredEuds(int wp) const{
	double Effuds = 0;
	for(int i=0;i<NbOfJetsBins_;i++) Effuds += GetPredEuds(wp,i)*GetPredNtotal(wp,i);
	Effuds /= GetPredNtotal(wp);
	return Effuds;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredEuds(int wp, int njets) const{
	return (njets<NbOfJetsBins_? euds_mc_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredEudsErr(int wp) const{
	double EffudsErr = 0;
	for(int i=0;i<NbOfJetsBins_;i++) EffudsErr += pow(GetPredEudsErr(wp,i)*(GetPredNtotal(wp,i)/GetPredNtotal(wp)),2);
	return sqrt(EffudsErr);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredEudsErr(int wp, int njets) const{
	return (njets<NbOfJetsBins_? euds_err_mc_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
void VJetEstimation::SetPredEuds(int wp, int njets, double value) const{
	euds_mc_[wp][njets] = value;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNv(int wp) const{
	double Nvlike = 0;
	for(unsigned int i=0;i<iDatasetsVLike_.size();i++){
		for(int j=0; j<NbOfJetsBins_; j++){
			for(int k=0; k<NbOfBJetsBins_; k++) Nvlike  += N_[wp][j][k][iDatasetsVLike_[i]];
		}
	}
	return Nvlike;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNv(int wp, int njets) const{
	double Nvlike = 0;
	for(unsigned int i=0;i<iDatasetsVLike_.size();i++){
		for(int j=0; j<NbOfBJetsBins_; j++) Nvlike  += N_[wp][njets][j][iDatasetsVLike_[i]];
	}
	return Nvlike;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNv(int wp, int njets, int nbjets) const{
	double Nvlike = 0;
	if(njets == -1){
		for(unsigned int i=0;i<iDatasetsVLike_.size();i++){
			for(int j=0;j<NbOfJetsBins_;j++) Nvlike += N_[wp][j][nbjets][iDatasetsVLike_[i]];
		}
		return Nvlike;
	}
	else{
		for(unsigned int i=0;i<iDatasetsVLike_.size();i++){
			Nvlike  += N_[wp][njets][nbjets][iDatasetsVLike_[i]];
		}
		return Nvlike;
	}
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNvErr(int wp) const{
	double NvErr = 0;
	for(int i=0;i<NbOfJetsBins_;i++) NvErr += pow(GetPredNvErr(wp,i),2);
	return sqrt(NvErr);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNvErr(int wp, int njets) const{
	double Nvlike = 0;
	for(unsigned int i=0;i<iDatasetsVLike_.size();i++){
		for(int j=0; j<NbOfBJetsBins_; j++) Nvlike  += pow(N_err_[wp][njets][j][iDatasetsVLike_[i]],2);
	}
	return sqrt(Nvlike);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNvErr(int wp, int njets, int nbjets) const{
	double Nvlike = 0;
	if(njets == -1){
		for(unsigned int i=0;i<iDatasetsVLike_.size();i++){
			for(int j=0;j<NbOfJetsBins_;j++) Nvlike += pow(N_err_[wp][j][nbjets][iDatasetsVLike_[i]],2);
		}
		return sqrt(Nvlike);
	}
	else{
		for(unsigned int i=0;i<iDatasetsVLike_.size();i++){
			Nvlike  += pow(N_err_[wp][njets][nbjets][iDatasetsVLike_[i]],2);
		}
		return sqrt(Nvlike);
	}
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNvb(int wp) const{
	double Nvblike = 0;
	for(unsigned int i=0;i<iDatasetsVbLike_.size();i++){
		for(int j=0; j<NbOfJetsBins_; j++){
			for(int k=0; k<NbOfBJetsBins_; k++) Nvblike  += N_[wp][j][k][iDatasetsVbLike_[i]];
		}
	}
	return Nvblike;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNvb(int wp, int njets) const{
	double Nvblike = 0;
	for(unsigned int i=0;i<iDatasetsVbLike_.size();i++){
		for(int j=0; j<NbOfBJetsBins_; j++) Nvblike  += N_[wp][njets][j][iDatasetsVbLike_[i]];
	}
	return Nvblike;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNvb(int wp, int njets, int nbjets) const{
	double Nvblike = 0;
	if(njets == -1){
		for(unsigned int i=0;i<iDatasetsVbLike_.size();i++){
			for(int j=0;j<NbOfJetsBins_;j++) Nvblike += N_[wp][j][nbjets][iDatasetsVbLike_[i]];
		}
		return Nvblike;
	}
	else{
		for(unsigned int i=0;i<iDatasetsVbLike_.size();i++){
			Nvblike  += N_[wp][njets][nbjets][iDatasetsVbLike_[i]];
		}
		return Nvblike;
	}
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNvbErr(int wp) const{
	double NvbErr = 0;
	for(int i=0;i<NbOfJetsBins_;i++) NvbErr += pow(GetPredNvbErr(wp,i),2);
	return sqrt(NvbErr);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNvbErr(int wp, int njets) const{
	double Nvblike = 0;
	for(unsigned int i=0;i<iDatasetsVbLike_.size();i++){
		for(int j=0; j<NbOfBJetsBins_; j++) Nvblike  += pow(N_err_[wp][njets][j][iDatasetsVbLike_[i]],2);
	}
	return sqrt(Nvblike);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNvbErr(int wp, int njets, int nbjets) const{
	double Nvblike = 0;
	if(njets == -1){
		for(unsigned int i=0;i<iDatasetsVbLike_.size();i++){
			for(int j=0;j<NbOfJetsBins_;j++) Nvblike += pow(N_err_[wp][j][nbjets][iDatasetsVbLike_[i]],2);
		}
		return sqrt(Nvblike);
	}
	else{
		for(unsigned int i=0;i<iDatasetsVbLike_.size();i++){
			Nvblike  += pow(N_err_[wp][njets][nbjets][iDatasetsVbLike_[i]],2);
		}
		return sqrt(Nvblike);
	}
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNtt(int wp) const{
	double Nttlike = 0;
	for(unsigned int i=0;i<iDatasetsTTLike_.size();i++){
		for(int j=0; j<NbOfJetsBins_; j++){
			for(int k=0; k<NbOfBJetsBins_; k++) Nttlike  += N_[wp][j][k][iDatasetsTTLike_[i]];
		}
	}
	return Nttlike;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNtt(int wp, int njets) const{
	double Nttlike = 0;
	for(unsigned int i=0;i<iDatasetsTTLike_.size();i++){
		for(int j=0; j<NbOfBJetsBins_; j++) Nttlike  += N_[wp][njets][j][iDatasetsTTLike_[i]];
	}
	return Nttlike;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNtt(int wp, int njets, int nbjets) const{
	double Nttlike = 0;
	if(njets == -1){
		for(unsigned int i=0;i<iDatasetsTTLike_.size();i++){
			for(int j=0;j<NbOfJetsBins_;j++) Nttlike += N_[wp][j][nbjets][iDatasetsTTLike_[i]];
		}
		return Nttlike;
	}
	else{
		for(unsigned int i=0;i<iDatasetsTTLike_.size();i++){
			Nttlike  += N_[wp][njets][nbjets][iDatasetsTTLike_[i]];
		}
		return Nttlike;
	}
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNttErr(int wp) const{
	double NttErr = 0;
	for(int i=0;i<NbOfJetsBins_;i++) NttErr += pow(GetPredNttErr(wp,i),2);
	return sqrt(NttErr);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNttErr(int wp, int njets) const{
	double Nttlike = 0;
	for(unsigned int i=0;i<iDatasetsTTLike_.size();i++){
		for(int j=0; j<NbOfBJetsBins_; j++) Nttlike  += pow(N_err_[wp][njets][j][iDatasetsTTLike_[i]],2);
	}
	return sqrt(Nttlike);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNttErr(int wp, int njets, int nbjets) const{
	double Nttlike = 0;
	if(njets == -1){
		for(unsigned int i=0;i<iDatasetsTTLike_.size();i++){
			for(int j=0;j<NbOfJetsBins_;j++) Nttlike += pow(N_err_[wp][j][nbjets][iDatasetsTTLike_[i]],2);
		}
		return sqrt(Nttlike);
	}
	else{
		for(unsigned int i=0;i<iDatasetsTTLike_.size();i++){
			Nttlike  += pow(N_err_[wp][njets][nbjets][iDatasetsTTLike_[i]],2);
		}
		return sqrt(Nttlike);
	}
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNtotal(int wp) const {
	double NtotalMC = 0;
	for(int i=0;i<NbOfJetsBins_;i++){
			for(int j=0;j<NbOfBJetsBins_;j++) NtotalMC += Nbjets_[wp][i][j];
	}
	return NtotalMC;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNtotal(int wp, int njets) const {
	double NtotalMC = 0;
	for(int i=0;i<NbOfBJetsBins_;i++) NtotalMC += Nbjets_[wp][njets][i];
	return NtotalMC;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredNtotal(int wp, int njets, int nbjets) const{
	return Nbjets_[wp][njets][nbjets];
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredN(int idx) const {
	double N_MC = 0;
	for(int i=0;i<NbOfJetsBins_;i++){
			for(int j=0;j<NbOfBJetsBins_;j++) N_MC += N_[0][i][j][idx];
	}
	return N_MC;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetPredN(int idx, int njets) const {
	double N_MC = 0;
	for(int i=0;i<NbOfBJetsBins_;i++) N_MC += N_[0][njets][i][idx];
	return N_MC;
}

//////////////////////////////////////
// Access to the parameter estimated values
//////////////////////////////////////

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstEb(int wp) const{
	double Effb = 0;
	for(int i=0;i<NbOfJetsBins_;i++) Effb += GetEstEb(wp,i)*GetPredNtotal(wp,i);
	Effb /= GetPredNtotal(wp);
	return Effb;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstEb(int wp, int njets) const{
	return (njets<NbOfJetsBins_? eb_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstEbErr(int wp) const{
	double EffbErr = 0;
	for(int i=0;i<NbOfJetsBins_;i++) EffbErr += pow(GetEstEbErr(wp,i)*(GetPredNtotal(wp,i)/GetPredNtotal(wp)),2);
	return sqrt(EffbErr);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstEbErr(int wp, int njets) const{
	return (njets<NbOfJetsBins_? eb_err_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstEudsc(int wp) const{
	double Effudsc = 0;
	for(int i=0;i<NbOfJetsBins_;i++) Effudsc += GetEstEudsc(wp,i)*GetPredNtotal(wp,i);
	Effudsc /= GetPredNtotal(wp);
	return Effudsc;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstEudsc(int wp, int njets) const{
	return (njets<NbOfJetsBins_? eudsc_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstEudscErr(int wp) const{
	double EffudscErr = 0;
	for(int i=0;i<NbOfJetsBins_;i++) EffudscErr += pow(GetEstEudscErr(wp,i)*(GetPredNtotal(wp,i)/GetPredNtotal(wp)),2);
	return sqrt(EffudscErr);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstEudscErr(int wp, int njets) const{
	return (njets<NbOfJetsBins_? eudsc_err_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstEuds(int wp) const{
	double Effuds = 0;
	for(int i=0;i<NbOfJetsBins_;i++) Effuds += GetEstEuds(wp,i)*GetPredNtotal(wp,i);
	Effuds /= GetPredNtotal(wp);
	return Effuds;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstEuds(int wp, int njets) const{
	return (njets<NbOfJetsBins_? euds_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstEudsErr(int wp) const{
	double EffudsErr = 0;
	for(int i=0;i<NbOfJetsBins_;i++) EffudsErr += pow(GetEstEudsErr(wp,i)*(GetPredNtotal(wp,i)/GetPredNtotal(wp)),2);
	return sqrt(EffudsErr);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstEudsErr(int wp, int njets) const{
	return (njets<NbOfJetsBins_? euds_err_[wp][njets] : -9999.);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNv(int wp) const{
	double Nv = 0;
	for(int i=0;i<NbOfJetsBins_;i++) Nv += GetEstNv(wp,i);
	return Nv;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNv(int wp, int njets) const{
	return (njets<NbOfJetsBins_ ? Nv_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNv(int wp, int njets, int nbjets) const{
	double Nvlike = -9999;
	if(njets>=NbOfJetsBins_) return Nvlike;
	else if(njets == -1){
		Nvlike = 0;
		for(int i=0;i<NbOfJetsBins_;i++) Nvlike += Nv_bjets( Nv_[i], euds_[wp][i], Njets_[i], nbjets);
		return Nvlike;
	}
	else{
		Nvlike = Nv_bjets( Nv_[njets], euds_[wp][njets], Njets_[njets], nbjets);
		return Nvlike;
	}
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNvErr(int wp) const{
	double NvErr = 0;
	for(int i=0;i<NbOfJetsBins_;i++) NvErr += pow(GetEstNvErr(wp,i),2);
	return sqrt(NvErr);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNvErr(int wp, int njets) const{
	return (njets<NbOfJetsBins_ ? Nv_err_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNvErrUp(int wp, int njets) const{
	return (njets<NbOfJetsBins_ ? Nv_err_up_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNvErrDown(int wp, int njets) const{
	return (njets<NbOfJetsBins_ ? Nv_err_down_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNvErr(int wp, int njets, int nbjets) const{
	double NvlikeErr = -9999;
	if(njets>=NbOfJetsBins_) return NvlikeErr;
	else if(njets == -1){
		NvlikeErr = 0;
		for(int i=0;i<NbOfJetsBins_;i++) NvlikeErr += pow(Nv_err_bjets( Nv_[i], Nv_err_[i], euds_[wp][i], euds_err_[wp][i], Njets_[i], nbjets),2);
		return sqrt(NvlikeErr);
	}
	else{
		NvlikeErr = Nv_err_bjets( Nv_[njets], Nv_err_[njets], euds_[wp][njets], euds_err_[wp][njets], Njets_[njets], nbjets);
		return NvlikeErr;
	}
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNvb(int wp) const{
	double Nvb = 0;
	for(int i=0;i<NbOfJetsBins_;i++) Nvb += GetEstNvb(wp,i);
	return Nvb;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNvb(int wp, int njets) const{
	return (njets<NbOfJetsBins_ ? Nvb_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNvb(int wp, int njets, int nbjets) const{
	double Nvblike = -9999;
	if(njets>=NbOfJetsBins_) return Nvblike;
	else if(njets == -1){
		Nvblike = 0;
		for(int i=0;i<NbOfJetsBins_;i++) Nvblike += Nvb_bjets( Nvb_[i], eb_[wp][njets], euds_[wp][i], Njets_[i], nbjets);
		return Nvblike;
	}
	else{
		Nvblike = Nvb_bjets( Nvb_[njets], eb_[wp][njets], euds_[wp][njets], Njets_[njets], nbjets);
		return Nvblike;
	}
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNvbErr(int wp) const{
	double NvbErr = 0;
	for(int i=0;i<NbOfJetsBins_;i++) NvbErr += pow(GetEstNvbErr(wp,i),2);
	return sqrt(NvbErr);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNvbErr(int wp, int njets) const{
	return (njets<NbOfJetsBins_ ? Nvb_err_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNvbErr(int wp, int njets, int nbjets) const{
	return -9999;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNtt(int wp) const{
	double Ntt = 0;
	for(int i=0;i<NbOfJetsBins_;i++) Ntt += GetEstNtt(wp,i);
	return Ntt;
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNtt(int wp, int njets) const{
	return (njets<NbOfJetsBins_ ? Ntt_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNtt(int wp, int njets, int nbjets) const{
	double Nttlike = -9999;
	if(njets>=NbOfJetsBins_) return Nttlike;
	else if(njets == -1){
		Nttlike = 0;
		for(int i=0;i<NbOfJetsBins_;i++) Nttlike += Ntt_bjets( Ntt_[i], eb_[wp][i], eudsc_[wp][i], Njets_[i], nbjets);
		return Nttlike;
	}
	else{
		Nttlike = Ntt_bjets( Ntt_[njets], eb_[wp][njets], eudsc_[wp][njets], Njets_[njets], nbjets);
		return Nttlike;
	}
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNttErr(int wp) const{
	double NttErr = 0;
	for(int i=0;i<NbOfJetsBins_;i++) NttErr += pow(GetEstNttErr(wp,i),2);
	return sqrt(NttErr);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNttErr(int wp, int njets) const{
	return (njets<NbOfJetsBins_ ? Ntt_err_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNttErrUp(int wp, int njets) const{
	return (njets<NbOfJetsBins_ ? Ntt_err_up_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNttErrDown(int wp, int njets) const{
	return (njets<NbOfJetsBins_ ? Ntt_err_down_[njets] : -9999);
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNttErr(int wp, int njets, int nbjets) const{
	double NttlikeErr = -9999;
	if(njets>=NbOfJetsBins_) return NttlikeErr;
	else if(njets == -1){
		NttlikeErr = 0;
		for(int i=0;i<NbOfJetsBins_;i++) NttlikeErr += pow(Ntt_err_bjets( Ntt_[i], Ntt_err_[i], eb_[wp][i], eb_err_[wp][i], eudsc_[wp][i], eudsc_err_[wp][i], Njets_[i], nbjets),2);
		return sqrt(NttlikeErr);
	}
	else{
		NttlikeErr = Ntt_err_bjets( Ntt_[njets], Ntt_err_[njets], eb_[wp][njets], eb_err_[wp][njets], eudsc_[wp][njets], eudsc_err_[wp][njets], Njets_[njets], nbjets);
		return NttlikeErr;
	}
}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNtotal(int wp)                        const{ return GetEstNtt(wp)              + GetEstNv(wp)              + GetEstNvb(wp);}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNtotal(int wp, int njets)             const{ return GetEstNtt(wp,njets)        + GetEstNv(wp,njets)        + GetEstNvb(wp,njets);}

/**________________________________________________________________________________________________________________*/
double VJetEstimation::GetEstNtotal(int wp, int njets, int nbjets) const{ return GetEstNtt(wp,njets,nbjets) + GetEstNv(wp,njets,nbjets) + GetEstNvb(wp,njets,nbjets);}


///////////////////////////////////////////////////////////////////
// Getting and setting statistical and systematic errors
///////////////////////////////////////////////////////////////////

double VJetEstimation::GetNvMSE(int wp)  const{
	double MSE = 0;
	MSE = pow(GetNvVar(wp),2)+pow(GetNvBias(wp),2);
	return sqrt(MSE);
}
double VJetEstimation::GetNvVar(int wp)  const{
	double error = 0;
	for(int i=0;i<NbOfJetsBins_;i++) error += pow(NvVar_[wp][i][4],2);
	return sqrt(error);
}
double VJetEstimation::GetNvBias(int wp)   const{
	double error = 0;
	for(int i=0;i<NbOfJetsBins_;i++) error += pow(NvBias_[wp][i][4],2);
	return sqrt(error);
}
double VJetEstimation::GetNvbVar(int wp)  const{
	double error = 0;
	for(int i=0;i<NbOfJetsBins_;i++) error += pow(NvbVar_[wp][i][4],2);
	return sqrt(error);
}
double VJetEstimation::GetNvbBias(int wp)   const{
	double error = 0;
	for(int i=0;i<NbOfJetsBins_;i++) error += pow(NvbBias_[wp][i][4],2);
	return sqrt(error);
}
double VJetEstimation::GetNttMSE(int wp)  const{
	double MSE = 0;
	MSE = pow(GetNttVar(wp),2)+pow(GetNttBias(wp),2);
	return sqrt(MSE);
}
double VJetEstimation::GetNttVar(int wp) const{
	double error = 0;
	for(int i=0;i<NbOfJetsBins_;i++) error += pow(NttVar_[wp][i][4],2);
	return sqrt(error);
}
double VJetEstimation::GetNttBias(int wp)  const{
	double error = 0;
	for(int i=0;i<NbOfJetsBins_;i++) error += pow(NttBias_[wp][i][4],2);
	return sqrt(error);
}

void VJetEstimation::FillSummaryHistos(){

	TString NbjetsXLabel[4] = {"0","1","2","#geq 3"};
	double N = 1;
	double Nerr = 0;

	for(int i=0;i<NbOfBtagWorkingPoint_;i++) {
		for(int j=0;j<NbOfJetsBins_;j++){
			for(int k=0;k<NbOfBJetsBins_;k++){
			// For V-like and TT-like estimation
				hNbjetsEstSummary[i][j][0]->SetBinContent(k+1,GetEstNv( i,j,k));
				hNbjetsEstSummary[i][j][0]->SetBinError(k+1,GetEstNvErr( i,j,k));
				//hNbjetsEstSummary[i][j][1]->SetBinContent(k+1,GetEstNvb(i,j,k));
				hNbjetsEstSummary[i][j][2]->SetBinContent(k+1,GetEstNtt(i,j,k));
				hNbjetsEstSummary[i][j][2]->SetBinError(k+1,GetEstNttErr(i,j,k));
				
				hNbjetsEstSummary[i][j][0]->GetXaxis()->SetBinLabel(k+1,NbjetsXLabel[k]);
				//hNbjetsEstSummary[i][j][1]->GetXaxis()->SetBinLabel(k+1,NbjetsXLabel[k]);
				hNbjetsEstSummary[i][j][2]->GetXaxis()->SetBinLabel(k+1,NbjetsXLabel[k]);

			// For V-like and TT-like MC prediction
				hNbjetsMCSummary[i][j][0]->SetBinContent(k+1,GetPredNv( i,j,k));
				hNbjetsMCSummary[i][j][0]->SetBinError(k+1,GetPredNvErr( i,j,k));
				//hNbjetsMCSummary[i][j][1]->SetBinContent(k+1,GetPredNvb(i,j,k));
				hNbjetsMCSummary[i][j][2]->SetBinContent(k+1,GetPredNtt(i,j,k));
				hNbjetsMCSummary[i][j][2]->SetBinError(k+1,GetPredNttErr(i,j,k));
				
				hNbjetsMCSummary[i][j][0]->GetXaxis()->SetBinLabel(k+1,NbjetsXLabel[k]);
				//hNbjetsMCSummary[i][j][1]->GetXaxis()->SetBinLabel(k+1,NbjetsXLabel[k]);
				hNbjetsMCSummary[i][j][2]->GetXaxis()->SetBinLabel(k+1,NbjetsXLabel[k]);
				
				hNbjets_mc_[i][j][0]->SetBinContent(k+1,GetPredNv( i,j,k));
				hNbjets_mc_[i][j][0]->SetBinError(k+1,GetPredNvErr( i,j,k));
				//hNbjets_mc_[i][j][1]->SetBinContent(k+1,GetPredNvb(i,j,k));
				//hNbjets_mc_[i][j][1]->SetBinError(k+1,GetPredNvbErr( i,j,k));
				hNbjets_mc_[i][j][2]->SetBinContent(k+1,GetPredNtt(i,j,k));
				hNbjets_mc_[i][j][2]->SetBinError(k+1,GetPredNttErr( i,j,k));
				
				hNbjets_pdf_mc_[i][j][0]->SetBinContent(k+1,Nv_bjets( N, euds_mc_[i][j], Njets_[j], k));
				hNbjets_pdf_mc_[i][j][0]->SetBinError(k+1,Nv_err_bjets( N, Nerr, euds_mc_[i][j], euds_err_mc_[i][j], Njets_[j], k));
				//hNbjets_pdf_mc_[i][j][1]->SetBinContent(k+1,Nvb_bjets( N, eb_mc_[i][j], euds_mc_[i][j],  Njets_[j], k));
				//hNbjets_pdf_mc_[i][j][1]->SetBinError(k+1,Nvb_err_bjets(Nvb_[j], eb_mc_[i][j], euds_mc_[i][j],  Njets_[j], k));
				hNbjets_pdf_mc_[i][j][2]->SetBinContent(k+1,Ntt_bjets(N, eb_mc_[i][j], eudsc_mc_[i][j], Njets_[j], k));
				hNbjets_pdf_mc_[i][j][2]->SetBinError(k+1,Ntt_err_bjets( N, Nerr, eb_mc_[i][j], eb_err_mc_[i][j], eudsc_mc_[i][j], eudsc_err_mc_[i][j], Njets_[j], k));
				
				hNbjets_pdf_est_[i][j][0]->SetBinContent(k+1,GetEstNv( i,j,k));
				hNbjets_pdf_est_[i][j][0]->SetBinError(k+1,GetEstNvErr( i,j,k));
				//hNbjets_pdf_est_[i][j][1]->SetBinContent(k+1,GetEstNvb(i,j,k));
				//hNbjets_pdf_est_[i][j][1]->SetBinError(k+1,GetEstNvbErr(i,j,k));
				hNbjets_pdf_est_[i][j][2]->SetBinContent(k+1,GetEstNtt(i,j,k));
				hNbjets_pdf_est_[i][j][2]->SetBinError(k+1,GetEstNttErr(i,j,k));
			}

			hNbjets_mc_[i][j][0]->Scale((hNbjets_mc_[i][j][0]->Integral() != 0 ? 1/hNbjets_mc_[i][j][0]->Integral() : 1));
			hNbjets_mc_[i][j][1]->Scale((hNbjets_mc_[i][j][1]->Integral() != 0 ? 1/hNbjets_mc_[i][j][1]->Integral() : 1));
			hNbjets_mc_[i][j][2]->Scale((hNbjets_mc_[i][j][2]->Integral() != 0 ? 1/hNbjets_mc_[i][j][2]->Integral() : 1));
			hNbjets_mc_[i][j][0]->GetXaxis()->CenterLabels();hNbjets_mc_[i][j][0]->GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_mc_[i][j][1]->GetXaxis()->CenterLabels();hNbjets_mc_[i][j][1]->GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_mc_[i][j][2]->GetXaxis()->CenterLabels();hNbjets_mc_[i][j][2]->GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_pdf_mc_[i][j][0]->Scale((hNbjets_pdf_mc_[i][j][0]->Integral() != 0 ? 1/hNbjets_pdf_mc_[i][j][0]->Integral() : 1));
			hNbjets_pdf_mc_[i][j][1]->Scale((hNbjets_pdf_mc_[i][j][1]->Integral() != 0 ? 1/hNbjets_pdf_mc_[i][j][1]->Integral() : 1));
			hNbjets_pdf_mc_[i][j][2]->Scale((hNbjets_pdf_mc_[i][j][2]->Integral() != 0 ? 1/hNbjets_pdf_mc_[i][j][2]->Integral() : 1));
			hNbjets_pdf_mc_[i][j][0]->GetXaxis()->CenterLabels();hNbjets_pdf_mc_[i][j][0]->GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_pdf_mc_[i][j][1]->GetXaxis()->CenterLabels();hNbjets_pdf_mc_[i][j][1]->GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_pdf_mc_[i][j][2]->GetXaxis()->CenterLabels();hNbjets_pdf_mc_[i][j][2]->GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_pdf_est_[i][j][0]->Scale((hNbjets_pdf_est_[i][j][0]->Integral() != 0 ? 1/hNbjets_pdf_est_[i][j][0]->Integral() : 1));
			hNbjets_pdf_est_[i][j][1]->Scale((hNbjets_pdf_est_[i][j][1]->Integral() != 0 ? 1/hNbjets_pdf_est_[i][j][1]->Integral() : 1));
			hNbjets_pdf_est_[i][j][2]->Scale((hNbjets_pdf_est_[i][j][2]->Integral() != 0 ? 1/hNbjets_pdf_est_[i][j][2]->Integral() : 1));
			hNbjets_pdf_est_[i][j][0]->GetXaxis()->CenterLabels();hNbjets_pdf_est_[i][j][0]->GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_pdf_est_[i][j][1]->GetXaxis()->CenterLabels();hNbjets_pdf_est_[i][j][1]->GetXaxis()->SetTitle("Nb. of b-tagged jets");
			hNbjets_pdf_est_[i][j][2]->GetXaxis()->CenterLabels();hNbjets_pdf_est_[i][j][2]->GetXaxis()->SetTitle("Nb. of b-tagged jets");

		}
		for(int j=0;j<NbOfBJetsBins_;j++){
		// For V-like and TT-like estimation
			hNbjetsEstSummary[i][NbOfJetsBins_][0]->SetBinContent(j+1,GetEstNv( i,-1,j));
			//hNbjetsEstSummary[i][NbOfJetsBins_][1]->SetBinContent(j+1,GetEstNvb(i,-1,j));
			hNbjetsEstSummary[i][NbOfJetsBins_][2]->SetBinContent(j+1,GetEstNtt(i,-1,j));
			hNbjetsEstSummary[i][NbOfJetsBins_][0]->GetXaxis()->SetBinLabel(j+1,NbjetsXLabel[j]);
			hNbjetsEstSummary[i][NbOfJetsBins_][1]->GetXaxis()->SetBinLabel(j+1,NbjetsXLabel[j]);
			hNbjetsEstSummary[i][NbOfJetsBins_][2]->GetXaxis()->SetBinLabel(j+1,NbjetsXLabel[j]);
			// For V-like and TT-like MC prediction
			hNbjetsMCSummary[i][NbOfJetsBins_][0]->SetBinContent(j+1,GetPredNv( i,-1,j));
			//hNbjetsMCSummary[i][NbOfJetsBins_][1]->SetBinContent(j+1,GetPredNvb(i,-1,j));
			hNbjetsMCSummary[i][NbOfJetsBins_][2]->SetBinContent(j+1,GetPredNtt(i,-1,j));
			hNbjetsMCSummary[i][NbOfJetsBins_][0]->GetXaxis()->SetBinLabel(j+1,NbjetsXLabel[i]);
			//hNbjetsMCSummary[i][NbOfJetsBins_][1]->GetXaxis()->SetBinLabel(j+1,NbjetsXLabel[i]);
			hNbjetsMCSummary[i][NbOfJetsBins_][2]->GetXaxis()->SetBinLabel(j+1,NbjetsXLabel[i]);
		}
	}
}
void VJetEstimation::BckgdSubstraction(vector<MCObsExpectation*> &hists, vector<string> &name, float lumi){
	if((unsigned int)NbOfBtagWorkingPoint_*NbOfJetsBins_ != hists.size()){
		cout<<"Mis-match between NbOfBtagWorkingPoint_*NbOfJetsBins_ and hists.size()"<<endl;
		return;
	}
	for(int i=0;i<NbOfBtagWorkingPoint_;i++) {
		for(int j=0;j<NbOfJetsBins_;j++){
			for(unsigned int k=0;k<name.size();k++){
				if(!hists[i*NbOfJetsBins_+j]->GetHistogram(name[k])) continue;
				for(int l=0;l<NbOfBJetsBins_;l++){
					Nbjets_[i][j][l]  += - (hists[i*NbOfJetsBins_+j]->GetHistogram(name[k])->GetBinContent(l+1))*lumi;
				}
			}
		}
	}
}

void VJetEstimation::PrintInputs(int njets){
	cout<<"For N = "<<Njets_[njets]<<" jets _________________________________________________"<<endl;
	for(int i=0;i<NbOfBtagWorkingPoint_;i++) {
		cout<<"B-tagging working point : "; cout<<this->GetBtagWorkingPoint(i)<<endl;
		cout<<"Input values : Nbjets (0/1/2/3 b-jets) = "; cout<<GetPredNtotal(i,njets)<<" ("<<Nbjets_[i][njets][0]<<"/"<<Nbjets_[i][njets][1]<<"/"<<Nbjets_[i][njets][2]<<"/"<<Nbjets_[i][njets][3]<<")"<<endl;
		cout<<"Predicted eudsc : "<<eudsc_mc_[i][njets]<<" +/- "<<eudsc_err_mc_[i][njets]<<endl;
		cout<<"Predicted euds  : "<<euds_mc_[i][njets]<<" +/- "<<euds_err_mc_[i][njets]<<endl;
		cout<<"Predicted eb    : "<<eb_mc_[i][njets]<<" +/- "<<eb_err_mc_[i][njets]<<endl;
		cout<<"For tt-like datasets : "<<endl;
		for(unsigned int k=0;k<iDatasetsTTLike_.size();k++){
			if(N_[i][njets][0][iDatasetsTTLike_[k]]==0 && N_[i][njets][0][iDatasetsTTLike_[k]]==N_[i][njets][1][iDatasetsTTLike_[k]] && N_[i][njets][0][iDatasetsTTLike_[k]]==N_[i][njets][2][iDatasetsTTLike_[k]] && N_[i][njets][0][iDatasetsTTLike_[k]]==N_[i][njets][3][iDatasetsTTLike_[k]]) continue;
			cout<<"Dataset "<<k+1<<" : N (0/1/2/3 jets) = "<<GetPredN(iDatasetsTTLike_[k],njets)<<" ("<<N_[i][njets][0][iDatasetsTTLike_[k]]<<"/"<<N_[i][njets][1][iDatasetsTTLike_[k]]<<"/"<<N_[i][njets][2][iDatasetsTTLike_[k]]<<"/"<<N_[i][njets][3][iDatasetsTTLike_[k]]<<endl;
		}
		cout<<"Total number : N (0/1/2/3 jets) = "<<GetPredNtt(i,njets)<<" ("<<GetPredNtt(i,njets,0)<<"/"<<GetPredNtt(i,njets,1)<<"/"<<GetPredNtt(i,njets,2)<<"/"<<GetPredNtt(i,njets,3)<<")"<<endl;
		cout<<"For W-like datasets : "<<endl;
		for(unsigned int k=0;k<iDatasetsVLike_.size();k++){
			if(N_[i][njets][0][iDatasetsVLike_[k]]==0 && N_[i][njets][0][iDatasetsVLike_[k]]==N_[i][njets][1][iDatasetsVLike_[k]] && N_[i][njets][0][iDatasetsVLike_[k]]==N_[i][njets][2][iDatasetsVLike_[k]] && N_[i][njets][0][iDatasetsVLike_[k]]==N_[i][njets][3][iDatasetsVLike_[k]]) continue;
			cout<<"Dataset "<<k+1<<" : N (0/1/2/3 jets) = "<<GetPredN(iDatasetsVLike_[k],njets)<<" ("<<N_[i][njets][0][iDatasetsVLike_[k]]<<"/"<<N_[i][njets][1][iDatasetsVLike_[k]]<<"/"<<N_[i][njets][2][iDatasetsVLike_[k]]<<"/"<<N_[i][njets][3][iDatasetsVLike_[k]]<<endl;
		}
		cout<<"Total number : N (0/1/2/3 jets) = "<<GetPredNv(i,njets)<<" ("<<GetPredNv(i,njets,0)<<"/"<<GetPredNv(i,njets,1)<<"/"<<GetPredNv(i,njets,2)<<"/"<<GetPredNv(i,njets,3)<<")"<<endl;
	}
}
void VJetEstimation::PrintResults(int njets){
	cout<<"***********************************************************"<<endl;
	cout<<"For N = "<<Njets_[njets]<<" jets"<<endl;
	cout<<setprecision(4);
	cout<<"e0bq/e1bq/e2bq = : "<<e0bq_[njets]<<"/"<<e1bq_[njets]<<"/"<<e2bq_[njets]<<"/"<<endl;
	cout<<"Btag working point : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetBtagWorkingPoint(i)<<" || ";} cout<<endl;
	cout<<"Estimated/Predicted eb	 : ";	for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstEb(i,njets)	<<" +/- "<<GetEstEbErr(i,njets)		<<" // "<<GetPredEb(i,njets)   <<" +/- "<<GetPredEbErr(i,njets)<<" || ";}    cout<<endl;
	cout<<"Estimated/Predicted eudsc : ";	for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstEudsc(i,njets)<<" +/- "<<GetEstEudscErr(i,njets)	<<" // "<<GetPredEudsc(i,njets)<<" +/- "<<GetPredEudscErr(i,njets)<<" || ";} cout<<endl;
	cout<<"Estimated/Predicted euds  : ";	for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstEuds(i,njets) <<" +/- "<<GetEstEudsErr(i,njets)	<<" // "<<GetPredEuds(i,njets) <<" +/- "<<GetPredEudsErr(i,njets)<<" || ";}  cout<<endl;
	cout<<setprecision(1);
	cout<<"Estimated Ntt-like : ";cout<<GetEstNtt(0,njets)<<" +/- "<<GetEstNttErr(0,njets)<<" / Predicted value from MC : "<<GetPredNtt(0,njets)<<"	+/- "<<GetPredNttErr(0,njets)<<endl;
	cout<<"- 0 b-jet          : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNtt(i, njets, 0)<<" +/- "<<GetEstNttErr(i,njets,0)<<" // "<< GetPredNtt(i, njets, 0)<<" +/- "<<GetPredNttErr(i,njets,0)<<" || ";} cout<<endl;
	cout<<"- 1 b-jet          : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNtt(i, njets, 1)<<" +/- "<<GetEstNttErr(i,njets,1)<<" // "<< GetPredNtt(i, njets, 1)<<" +/- "<<GetPredNttErr(i,njets,1)<<" || ";} cout<<endl;
	cout<<"- 2 b-jets         : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNtt(i, njets, 2)<<" +/- "<<GetEstNttErr(i,njets,2)<<" // "<< GetPredNtt(i, njets, 2)<<" +/- "<<GetPredNttErr(i,njets,2)<<" || ";} cout<<endl;
	cout<<"- 3 b-jets         : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNtt(i, njets, 3)<<" +/- "<<GetEstNttErr(i,njets,3)<<" // "<< GetPredNtt(i, njets, 3)<<" +/- "<<GetPredNttErr(i,njets,3)<<" || ";} cout<<endl;
// 	cout<<"Estimated Nvb-like : ";cout<<Nvb_[njets]<<" / Predicted value from MC : "<<GetPredNvb(0,njets)<<endl;
// 	cout<<"- 0 b-jet          : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNvb(i, njets, 0)<<" // "<< GetPredNvb(i, njets, 0) <<" || ";} cout<<endl;
// 	cout<<"- 1 b-jet          : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNvb(i, njets, 1)<<" // "<< GetPredNvb(i, njets, 1) <<" || ";} cout<<endl;
// 	cout<<"- 2 b-jets         : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNvb(i, njets, 2)<<" // "<< GetPredNvb(i, njets, 2) <<" || ";} cout<<endl;
// 	cout<<"- 3 b-jets         : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNvb(i, njets, 3)<<" // "<< GetPredNvb(i, njets, 3) <<" || ";} cout<<endl;
	cout<<"Estimated Nv-like  : ";cout<<GetEstNv(0,njets)<<" +/- "<<GetEstNvErr(0,njets)<<" / Predicted value from MC : "<<GetPredNv(0,njets)<<"	+/- "<<GetPredNvErr(0,njets)<<endl;
	cout<<"- 0 b-jet          : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNv(i, njets, 0)<<" +/- "<<GetEstNvErr(i,njets,0)<<" // "<< GetPredNv(i, njets, 0)<<" +/- "<<GetPredNvErr(i,njets,0)<<" || ";} cout<<endl;
	cout<<"- 1 b-jet          : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNv(i, njets, 1)<<" +/- "<<GetEstNvErr(i,njets,1)<<" // "<< GetPredNv(i, njets, 1)<<" +/- "<<GetPredNvErr(i,njets,1)<<" || ";} cout<<endl;
	cout<<"- 2 b-jets         : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNv(i, njets, 2)<<" +/- "<<GetEstNvErr(i,njets,2)<<" // "<< GetPredNv(i, njets, 2)<<" +/- "<<GetPredNvErr(i,njets,2)<<" || ";} cout<<endl;
	cout<<"- 3 b-jets         : ";for(int i=0;i<GetNbOfBtagWorkingPoint();i++) {cout<<GetEstNv(i, njets, 3)<<" +/- "<<GetEstNvErr(i,njets,3)<<" // "<< GetPredNv(i, njets, 3)<<" +/- "<<GetPredNvErr(i,njets,3)<<" || ";} cout<<endl;
	cout<<"***********************************************************"<<endl;
}

double VJetEstimation::LikelihoodFunct(const double *xx) const{
	int njets = (int)xx[1];
	double Nvb = 0;
	double L = 0;
	double cst = 0;
	int btag_wp_idx = (int)xx[0];
	for(int j=0;j<NbOfBJetsBins_;j++){
		cst = 0;
		for(int k=1;k<=Nbjets_[btag_wp_idx][njets][j];k++) cst += log(k);
		L += Nbjets_[btag_wp_idx][njets][j]*log(Nbjets((double&)xx[2], (double&)xx[3], Nvb, (double&)xx[4], (double&)xx[5],(double&)xx[6],Njets_[njets],j))-Nbjets((double&)xx[2], (double&)xx[3], Nvb, (double&)xx[4], (double&)xx[5],(double&)xx[6],Njets_[njets],j)-cst;
	}
	return (-L);
}

double VJetEstimation::JointNjetsLikelihoodFunct(const double *xx) const{
	double Nvb = 0;
	double L = 0;
	double cst = 0;
	int btag_wp_idx = (int)xx[0];
	for(int i=0;i<NbOfJetsBins_;i++){
		for(int j=0;j<NbOfBJetsBins_;j++){
			cst = 0;
			for(int k=1;k<=Nbjets_[btag_wp_idx][i][j];k++) cst += log(k);
			L += Nbjets_[btag_wp_idx][i][j]*log(Nbjets((double&)xx[1+2*i], (double&)xx[2+2*i], Nvb, (double&)xx[5], (double&)xx[6],(double&)xx[7],Njets_[i],j))-Nbjets((double&)xx[1+2*i], (double&)xx[2+2*i], Nvb, (double&)xx[5], (double&)xx[6],(double&)xx[7],Njets_[i],j)-cst;
		}
	}
	return (-L);
}

double VJetEstimation::JointWPLikelihoodFunct(const double *xx) const{
	int njets = (int)xx[0];
	double Nvb = 0;
	double L = 0;
	double cst = 0;
	for(int i=0;i<NbOfBtagWorkingPoint_;i++) {
		for(int j=0;j<NbOfBJetsBins_;j++){
			cst = 0;
			for(int k=1;k<=Nbjets_[i][njets][j];k++) cst += log(k);
			L += Nbjets_[i][njets][j]*log(Nbjets((double&)xx[1], (double&)xx[2], Nvb, (double&)xx[3+i*3], (double&)xx[4+i*3],(double&)xx[5+i*3],Njets_[njets],j))-Nbjets((double&)xx[1], (double&)xx[2], Nvb, (double&)xx[3+i*3], (double&)xx[4+i*3],(double&)xx[5+i*3],Njets_[njets],j)-cst;
		}
	}
	return (-L);
}

double VJetEstimation::JointLikelihoodFunct(const double *xx) const{
	double Nvb = 0;
	double L = 0;
	double cst = 0;
	for(int wp=0;wp<NbOfBtagWorkingPoint_;wp++){
		for(int i=0;i<NbOfJetsBins_;i++){
			for(int j=0;j<NbOfBJetsBins_;j++){
				cst = 0;
				for(int k=1;k<=Nbjets_[wp][i][j];k++) cst += log(k);
				L += Nbjets_[wp][i][j]*log(Nbjets((double&)xx[0+i*NbOfJetsBins_], (double&)xx[1+i*NbOfJetsBins_], Nvb, (double&)xx[2*NbOfJetsBins_+wp*3], (double&)xx[2*NbOfJetsBins_+1+wp*3],(double&)xx[2*NbOfJetsBins_+2+wp*3],Njets_[i],j))-Nbjets((double&)xx[0+i*NbOfJetsBins_], (double&)xx[1+i*NbOfJetsBins_], Nvb, (double&)xx[2*NbOfJetsBins_+wp*3], (double&)xx[2*NbOfJetsBins_+1+wp*3],(double&)xx[2*NbOfJetsBins_+2+wp*3],Njets_[i],j)-cst;
			}
		}
	}
	return (-L);
}

double VJetEstimation::MinimizerFunct(const double *xx) const{
	int njets = (int)xx[0];
	double ChiSq = 0;
	double VarNbjets[NbOfBJetsBins_];
	for(int i=0;i<NbOfBtagWorkingPoint_;i++) {
		for(int j=0;j<NbOfBJetsBins_;j++){
			VarNbjets[j] = (Nbjets_[i][njets][j] == 0 ? 1 : Nbjets_[i][njets][j]);
			ChiSq += pow(Nbjets( (double&)xx[1], (double&)xx[2], (double&)xx[3], (double&)xx[4+i*3], (double&)xx[5+i*3], (double&)xx[6+i*3], Njets_[njets],j) - Nbjets_[i][njets][j],2)/VarNbjets[j];
		}
	}
	//ChiSq /= (NbOfBtagWorkingPoint_>1 ? (4*NbOfBtagWorkingPoint_+1)-(2+3*NbOfBtagWorkingPoint_) : 1);
	return ChiSq;
}

///////////////////////////////////////////////////////////////////	
// Methods to calculates the estimators ...  (from equations)
///////////////////////////////////////////////////////////////////
double VJetEstimation::Nbjets(double &Ntt, double &Nv, double &Nvb, double &eb, double &eudsc,  double &euds, int &njets, int &nbjets) const{
	double Ntt_bjets = this->Ntt_bjets(Ntt,eb,eudsc,njets, nbjets);
	double Nv_bjets  = this->Nv_bjets( Nv, euds, njets, nbjets);
	double Nvb_bjets = this->Nvb_bjets(Nvb,eb,euds, njets, nbjets);
	return (Ntt_bjets+Nv_bjets+Nvb_bjets);
}

double VJetEstimation::Ntt_bjets(double &Ntt, double &eb, double &eudsc, int &njets, int &nbjets) const{
	if(nbjets == 0)      return Ntt_0bjet( Ntt, eb, eudsc, njets);
	else if(nbjets == 1) return Ntt_1bjet( Ntt, eb, eudsc, njets);
	else if(nbjets == 2) return Ntt_2bjets(Ntt, eb, eudsc, njets);
	else if(nbjets == 3) return Ntt_3bjets(Ntt, eb, eudsc, njets);
	else                 return -9999;
}

double VJetEstimation::Ntt_err_bjets(double &Ntt, double &Ntt_err, double &eb, double &eb_err, double &eudsc, double &eudsc_err, int &njets, int &nbjets) const{
	if(nbjets == 0)      return Ntt_err_0bjet( Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else if(nbjets == 1) return Ntt_err_1bjet( Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else if(nbjets == 2) return Ntt_err_2bjets(Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else if(nbjets == 3) return Ntt_err_3bjets(Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else                 return -9999;
}

double VJetEstimation::Nvb_bjets(double &Nvb, double &eb, double &euds, int &njets, int &nbjets) const{
	if(nbjets == 0)      return Nvb_0bjet( Nvb, eb, euds, njets);
	else if(nbjets == 1) return Nvb_1bjet( Nvb, eb, euds, njets);
	else if(nbjets == 2) return Nvb_2bjets(Nvb, eb, euds, njets);
	else if(nbjets == 3) return Nvb_3bjets(Nvb, eb, euds, njets);
	else                 return -9999;
}

double VJetEstimation::Nv_bjets(double &Nv, double &euds, int &njets, int &nbjets) const{
	if(nbjets == 0)      return Nv_0bjet( Nv, euds, njets);
	else if(nbjets == 1) return Nv_1bjet( Nv, euds, njets);
	else if(nbjets == 2) return Nv_2bjets(Nv, euds, njets);
	else if(nbjets == 3) return Nv_3bjets(Nv, euds, njets);
	else                 return -9999;
}

double VJetEstimation::Nv_err_bjets(double &Nv, double &Nv_err, double &euds, double &euds_err, int &njets, int &nbjets) const{
	if(nbjets == 0)      return Nv_err_0bjet( Nv, Nv_err, euds, euds_err, njets);
	else if(nbjets == 1) return Nv_err_1bjet( Nv, Nv_err, euds, euds_err, njets);
	else if(nbjets == 2) return Nv_err_2bjets(Nv, Nv_err, euds, euds_err, njets);
	else if(nbjets == 3) return Nv_err_3bjets(Nv, Nv_err, euds, euds_err, njets);
	else                 return -9999;
}

double VJetEstimation::N0bjet(double &Ntt, double &Nv, double &Nvb, double &eb, double &eudsc,  double &euds, int &n) const{
	double Ntt_0bjet = this->Ntt_0bjet(Ntt,eb,eudsc,n);
	double Nv_0bjet  = this->Nv_0bjet( Nv, euds, n);
	double Nvb_0bjet = this->Nvb_0bjet(Nvb,eb,euds, n);
	return (Ntt_0bjet+Nv_0bjet+Nvb_0bjet);
}

double VJetEstimation::Ntt_0bjet(double &Ntt, double &eb, double &eudsc, int &n) const{
	double  Ntt_0bjet = ((1-eb)*(1-eb)*pow((1-eudsc),n-2)*GetTTEff2bq(n-Njets_[0])
			  +         (1-eb)*pow((1-eudsc),n-1)*GetTTEff1bq(n-Njets_[0])
			  +                pow((1-eudsc),n)  *GetTTEff0bq(n-Njets_[0]))*Ntt;
	return (Ntt_0bjet < 0 ? 0 : Ntt_0bjet);
}
double VJetEstimation::Ntt_err_0bjet(double &Ntt, double &Ntt_err, double &eb, double &eb_err, double &eudsc, double &eudsc_err, int &n) const{
	double  Ntt_err_0bjet = pow((2*(eb-1)*pow((1-eudsc),n-2)*GetTTEff2bq(n-Njets_[0])
			            +    (-1)*pow((1-eudsc),n-1)*GetTTEff1bq(n-Njets_[0]))*Ntt*eb_err,2)
			      + pow((pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3)*(1-eb)*(1-eb)*GetTTEff2bq(n-Njets_[0])
			            +pow(-1.,n-1)*(n-1)*pow(eudsc-1,n-2)*       (1-eb)*GetTTEff1bq(n-Njets_[0])
			            +pow(-1.,n  )*(n  )*pow(eudsc-1,n-1)*       (1   )*GetTTEff0bq(n-Njets_[0]))*Ntt*eudsc_err,2)
			      + pow(Ntt_0bjet(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_0bjet < 0 ? 0 : sqrt(Ntt_err_0bjet));
}
double VJetEstimation::Nvb_0bjet(double &Nvb, double &eb, double &euds, int &n) const{
	double  Nvb_0bjet = (1-eb)*pow((1-euds),n-1)*Nvb;
	return (Nvb_0bjet < 0 ? 0 : Nvb_0bjet);
}
double VJetEstimation::Nv_0bjet (double &Nv, double &euds, int &n) const{
	double  Nv_0bjet = pow((1-euds),n)*Nv;
	return (Nv_0bjet < 0 ? 0 : Nv_0bjet);
}
double VJetEstimation::Nv_err_0bjet (double &Nv, double &Nv_err, double &euds, double &euds_err, int &n) const{
	double  Nv_err_0bjet = pow(pow(-1.,n)*n*pow(euds-1,n-1)*Nv*euds_err,2)
			     + pow(pow((1-euds),n)*Nv_err,2);
	return (Nv_err_0bjet < 0 ? 0 : sqrt(Nv_err_0bjet));
}
	
double VJetEstimation::N1bjet   (double &Ntt, double &Nv, double &Nvb, double &eb, double &eudsc, double &euds, int &n) const{
	double Ntt_1bjet = this->Ntt_1bjet(Ntt,eb,eudsc,n);
	double Nv_1bjet  = this->Nv_1bjet( Nv, euds, n);
	double Nvb_1bjet = this->Nvb_1bjet(Nvb,eb,euds, n);
	return (Ntt_1bjet+Nv_1bjet+Nvb_1bjet);
}
double VJetEstimation::Ntt_1bjet(double &Ntt, double &eb, double &eudsc, int &n) const{
	double  Ntt_1bjet =((2*eb*(1-eb)*pow(1-eudsc,n-2)+(1-eb)*(1-eb)*(n-2)*eudsc*pow(1-eudsc,n-3))*GetTTEff2bq(n-Njets_[0])
			  +          (eb*pow(1-eudsc,n-1)+       (1-eb)*(n-1)*eudsc*pow(1-eudsc,n-2))*GetTTEff1bq(n-Njets_[0])
			  +                                                (n*eudsc*pow(1-eudsc,n-1))*GetTTEff0bq(n-Njets_[0]))*Ntt;
	return (Ntt_1bjet < 0 ? 0 : Ntt_1bjet);
}
double VJetEstimation::Ntt_err_1bjet(double &Ntt, double &Ntt_err, double &eb,  double &eb_err, double &eudsc, double &eudsc_err, int &n) const{
	double  Ntt_err_1bjet = pow(((2*(1-2*eb)*pow(1-eudsc,n-2)+2*(eb-1)*(n-2)*eudsc*pow(1-eudsc,n-3))*GetTTEff2bq(n-Njets_[0])
			            +           (pow(1-eudsc,n-1)+    (-1)*(n-1)*eudsc*pow(1-eudsc,n-2))*GetTTEff1bq(n-Njets_[0]))*Ntt*eb_err,2)
			      + pow(((2*eb*(1-eb)*(n-2)*pow(-1.,n-2)*pow(eudsc-1,n-3)+2*(eb-1)*(n-2)*(pow(1-eudsc,n-3)+eudsc*pow(-1.,n-3)*(n-3)*pow(1-eudsc,n-4)))*GetTTEff2bq(n-Njets_[0])
			            +(         eb*(n-1)*pow(-1.,n-1)*pow(eudsc-1,n-2)+  (1-eb)*(n-1)*(pow(1-eudsc,n-2)+eudsc*pow(-1.,n-2)*(n-2)*pow(1-eudsc,n-3)))*GetTTEff1bq(n-Njets_[0])
			            +(                                                        (n  )*(pow(1-eudsc,n-1)+eudsc*pow(-1.,n-1)*(n-1)*pow(1-eudsc,n-2)))*GetTTEff0bq(n-Njets_[0]))*Ntt*eudsc_err,2)
			      + pow(Ntt_1bjet(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_1bjet < 0 ? 0 : sqrt(Ntt_err_1bjet));
}
double VJetEstimation::Nvb_1bjet(double &Nvb, double &eb, double &euds, int &n) const{
	double  Nvb_1bjet = (eb*pow(1-euds,n-1)+(1-eb)*(n-1)*euds*pow(1-euds,n-2))*Nvb;
	return (Nvb_1bjet < 0 ? 0 : Nvb_1bjet);
}
double VJetEstimation::Nv_1bjet(double &Nv, double &euds, int &n) const{
	double  Nv_1bjet = n*euds*pow(1-euds,n-1)*Nv;
	return (Nv_1bjet < 0 ? 0 : Nv_1bjet);
}
double VJetEstimation::Nv_err_1bjet(double &Nv, double &Nv_err, double &euds, double &euds_err, int &n) const{
	double  Nv_err_1bjet = pow(n*(pow(1-euds,n-1)+euds*pow(-1.,n-1)*(n-1)*pow(euds-1,n-2))*Nv*euds_err,2)
			     + pow(Nv_1bjet(Nv,euds,n)*Nv_err/Nv,2);
	return (Nv_err_1bjet < 0 ? 0 : sqrt(Nv_err_1bjet));
}
	
double VJetEstimation::N2bjets(double &Ntt, double &Nv, double &Nvb, double &eb, double &eudsc,  double &euds, int &n) const{
	double  Ntt_2bjets = this->Ntt_2bjets(Ntt,eb,eudsc,n);
	double  Nv_2bjets  = this->Nv_2bjets( Nv ,euds, n);
	double  Nvb_2bjets = this->Nvb_2bjets(Nvb,eb,euds, n);
	return (Ntt_2bjets+Nv_2bjets+Nvb_2bjets);
}
double VJetEstimation::Ntt_2bjets(double &Ntt, double &eb, double &eudsc, int &n) const{
	double  Ntt_2bjets =((eb*eb*pow(1-eudsc,n-2)+2*eb*(1-eb)*(n-2)*eudsc*pow(1-eudsc,n-3)+(1-eb)*(1-eb)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4))*GetTTEff2bq(n-Njets_[0])
			   +                                 (eb*(n-1)*eudsc*pow(1-eudsc,n-2)+       (1-eb)*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3))*GetTTEff1bq(n-Njets_[0])
			   +                                                                                   ((n*(n-1)/2)*eudsc*eudsc*pow(1-eudsc,n-2))*GetTTEff0bq(n-Njets_[0]))*Ntt;
	return (Ntt_2bjets<0 ? 0 : Ntt_2bjets);
}
double VJetEstimation::Ntt_err_2bjets(double &Ntt, double &Ntt_err, double &eb,  double &eb_err, double &eudsc, double &eudsc_err, int &n) const{
	double  Ntt_err_2bjets  = pow(((2*eb*pow(1-eudsc,n-2)+2*(1-2*eb)*(n-2)*eudsc*pow(1-eudsc,n-3)+2*(eb-1)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4))*GetTTEff2bq(n-Njets_[0])
			              +                                 ((n-1)*eudsc*pow(1-eudsc,n-2)+    (-1)*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3))*GetTTEff1bq(n-Njets_[0]))*Ntt*eb_err,2)
				+ pow(((pow(eb,2)*pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3)+2*eb*(1-eb)*(n-2)*(pow(1-eudsc,n-3)+eudsc*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))+pow(1-eb,2)*((n-2)*(n-3)/2)*(2*eudsc*pow(1-eudsc,n-4)+pow(eudsc,2)*pow(-1.,n-4)*(n-4)*pow(eudsc-1,n-5)))*GetTTEff2bq(n-Njets_[0])
				      +                                                      (eb*(n-1)*(pow(1-eudsc,n-2)+eudsc*pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3))+     (1-eb)*((n-1)*(n-2)/2)*(2*eudsc*pow(1-eudsc,n-3)+pow(eudsc,2)*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))))*GetTTEff1bq(n-Njets_[0])
				      +                                                                                                                 			(((n)*(n-1)/2)*(2*eudsc*pow(1-eudsc,n-2)+pow(eudsc,2)*pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3)))*GetTTEff0bq(n-Njets_[0])*Ntt*eudsc_err,2)
				+ pow(Ntt_2bjets(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_2bjets < 0 ? 0 : sqrt(Ntt_err_2bjets));
}

double VJetEstimation::Nvb_2bjets(double &Nvb, double &eb, double &euds, int &n) const{
	double  Nvb_2bjets = (eb*(n-1)*euds*pow(1-euds,n-2)+(1-eb)*((n-1)*(n-2)/2)*euds*euds*pow(1-euds,n-3))*Nvb;
	return (Nvb_2bjets<0 ? 0 : Nvb_2bjets);
}
double VJetEstimation::Nv_2bjets (double &Nv, double &euds, int &n) const{
	double  Nv_2bjets  = ((n*(n-1)/2)*euds*euds*pow((1-euds),n-2))*Nv;
	return (Nv_2bjets<0 ? 0 : Nv_2bjets);
}
double VJetEstimation::Nv_err_2bjets(double &Nv, double &Nv_err, double &euds, double &euds_err, int &n) const{
	double  Nv_err_2bjets = pow((n*(n-1)/2)*(2*euds*pow(1-euds,n-2)+pow(euds,2)*pow(-1.,n-2)*(n-2)*pow(euds-1,n-3))*Nv*euds_err,2)
			      + pow(Nv_2bjets(Nv,euds,n)*Nv_err/Nv,2);
	return (Nv_err_2bjets < 0 ? 0 : sqrt(Nv_err_2bjets));
}
	
double VJetEstimation::N3bjets   (double &Ntt, double &Nv, double &Nvb, double &eb, double &eudsc,  double &euds, int &n) const{
        double Ntt_3bjets = this->Ntt_3bjets(Ntt, eb, eudsc,n);
	double Nv_3bjets  = this->Nv_3bjets( Nv , euds, n);
	double Nvb_3bjets = this->Nvb_3bjets(Nvb, eb, euds, n);
	return (Ntt_3bjets+Nv_3bjets+Nvb_3bjets);
}
double VJetEstimation::Ntt_3bjets(double &Ntt, double &eb, double &eudsc, int &n) const{
	double  Ntt_3bjets =((eb*eb*(n-2)*eudsc*pow(1-eudsc,n-3)+2*eb*(1-eb)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4)+(n>4 ? pow((1-eb),2)*((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow((1-eudsc),n-5) : 0 ))*GetTTEff2bq(n-Njets_[0])
			  +                                              (eb*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3)+ (1-eb)*((n-1)*(n-2)*(n-3)/6)*pow(eudsc,3)*pow(1-eudsc,n-4))*GetTTEff1bq(n-Njets_[0])
			  +                                                                                                          ((n*(n-1)*(n-2)/6)*pow(eudsc,3)*pow(1-eudsc,n-3))*GetTTEff0bq(n-Njets_[0]))*Ntt;
	return (Ntt_3bjets<0 ? 0 : Ntt_3bjets);
}
double VJetEstimation::Ntt_err_3bjets(double &Ntt, double &Ntt_err, double &eb,  double &eb_err, double &eudsc, double &eudsc_err, int &n) const{
	double  Ntt_err_3bjets  = pow(((2*eb*(n-2)*eudsc*pow(1-eudsc,n-3)+2*(1-2*eb)*(n-2)*(n-3)/2*pow(eudsc,2)*pow(1-eudsc,n-4)+2*(eb-1)*((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow(1-eudsc,n-5))*GetTTEff2bq(n-Njets_[0])
			              +                                             ((n-1)*(n-2)/2*pow(eudsc,2)*pow(1-eudsc,n-3)+    (-1)*((n-1)*(n-2)*(n-3)/6)*pow(eudsc,3)*pow(1-eudsc,n-4))*GetTTEff1bq(n-Njets_[0]))*Ntt*eb_err,2)
				+ pow(((pow(eb,2)*(n-2)*(pow(1-eudsc,n-3)+eudsc*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))+2*eb*(1-eb)*((n-2)*(n-3)/2)*(2*eudsc*pow(1-eudsc,n-4)+pow(eudsc,2)*pow(-1.,n-4)*(n-4)*pow(eudsc-1,n-5))+pow(1-eb,2)*((n-2)*(n-3)*(n-4)/6)*(3*pow(eudsc,2)*pow(1-eudsc,n-5)+pow(eudsc,3)*pow(-1.,n-5)*(n-5)*pow(eudsc-1,n-6)))*GetTTEff2bq(n-Njets_[0])
				      +                                                     			            (eb*((n-1)*(n-2)/2)*(2*eudsc*pow(1-eudsc,n-3)+pow(eudsc,2)*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))+     (1-eb)*((n-1)*(n-2)*(n-3)/6)*(3*pow(eudsc,2)*pow(1-eudsc,n-4)+pow(eudsc,3)*pow(-1.,n-4)*(n-4)*pow(eudsc-1,n-5))))*GetTTEff1bq(n-Njets_[0])
				      +                                                                                                                 				      			      		        (((n)*(n-1)*(n-2)/6)*(3*pow(eudsc,2)*pow(1-eudsc,n-3)+pow(eudsc,3)*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4)))*GetTTEff0bq(n-Njets_[0])*Ntt*eudsc_err,2)
				+ pow(Ntt_3bjets(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_3bjets < 0 ? 0 : sqrt(Ntt_err_3bjets));
}
double VJetEstimation::Nvb_3bjets(double &Nvb, double &eb, double &euds, int &n) const{
	double  Nvb_3bjets = (eb*((n-1)*(n-2)/2)*euds*euds*pow(1-euds,n-3) + (1-eb)*((n-1)*(n-2)*(n-3)/6)*pow(euds,3)*pow(1-euds,n-Njets_[0]))*Nvb;
	return (Nvb_3bjets<0 ? 0 : Nvb_3bjets);
}
double VJetEstimation::Nv_3bjets (double &Nv, double &euds, int &n) const{
	double  Nv_3bjets  = ((n*(n-1)*(n-2)/6)*pow((euds),3)*pow((1-euds),n-3))*Nv;
	return (Nv_3bjets<0 ? 0 : Nv_3bjets);
}
double VJetEstimation::Nv_err_3bjets(double &Nv, double &Nv_err, double &euds, double &euds_err, int &n) const{
	double  Nv_err_3bjets = pow((n*(n-1)*(n-2)/6)*(3*pow(euds,2)*pow(1-euds,n-3)+pow(euds,3)*pow(-1.,n-3)*(n-3)*pow(euds-1,n-4))*Nv*euds_err,2)
			      + pow(Nv_3bjets(Nv,euds,n)*Nv_err/Nv,2);
	return (Nv_err_3bjets < 0 ? 0 : sqrt(Nv_err_3bjets));
}

double VJetEstimation::eudsc_fromN3bjets(double &N3bjets, double &Ntt, double &Nv, double &Nvb, double &eb, double &eudsc_old, double &euds, int &n) const{
	double eudsc =0;
	vector<double> roots;
	if(n==4)/////// Equation holds only for n = 4  ///////////////
        {
                double a4 = -4*(pow(1-ebq_[n-Njets_[0]],2)*Ntt);
                double a3 =  4*(pow(1-ebq_[n-Njets_[0]],2)*Ntt)+(-3*eb+(1-eb))*ebq_[n-Njets_[0]]*(1-ebq_[n-Njets_[0]])*Ntt;
                double a2 =                             ( 3*eb        *ebq_[n-Njets_[0]]*(1-ebq_[n-Njets_[0]])+(-2*eb*eb+2*eb*(1-eb))*pow(ebq_[n-Njets_[0]],2))*Ntt;
                double a1 =                                                                   (( 2*eb*eb            )*pow(ebq_[n-Njets_[0]],2))*Ntt;
                double a0 = -N3bjets+Nv_3bjets(Nv,euds,n)+Nvb_3bjets(Nvb,eb,euds,n);
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
                double a5 = 10*(pow(1-ebq_[n-Njets_[0]],2)*Ntt);
                double a4 =-20*(pow(1-ebq_[n-Njets_[0]],2)*Ntt) + (  6*eb-4*(1-eb))*ebq_[n-Njets_[0]]*(1-ebq_[n-Njets_[0]])*Ntt;
                double a3 = 10*(pow(1-ebq_[n-Njets_[0]],2)*Ntt) + (-12*eb+4*(1-eb))*ebq_[n-Njets_[0]]*(1-ebq_[n-Njets_[0]])*Ntt + ( 3*eb*eb-6*eb*(1-eb)+pow((1-eb),2))*pow(ebq_[n-Njets_[0]],2)*Ntt;
                double a2 =                               (  6*eb         )*ebq_[n-Njets_[0]]*(1-ebq_[n-Njets_[0]])*Ntt + (-6*eb*eb+6*eb*(1-eb)              )*pow(ebq_[n-Njets_[0]],2)*Ntt;
                double a1 =                                                                               ( 3*eb*eb                          )*pow(ebq_[n-Njets_[0]],2)*Ntt;
                double a0 = -N3bjets+Nv_3bjets(Nv,euds,n)+Nvb_3bjets(Nvb,eb,euds,n);
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
                double a6 =-20*(pow(1-ebq_[n-Njets_[0]],2)*Ntt);
                double a5 = 60*(pow(1-ebq_[n-Njets_[0]],2)*Ntt)+(-10*eb+10*(1-eb))*ebq_[n-Njets_[0]]*(1-ebq_[n-Njets_[0]])*Ntt;
                double a4 =-60*(pow(1-ebq_[n-Njets_[0]],2)*Ntt)+( 30*eb-20*(1-eb))*ebq_[n-Njets_[0]]*(1-ebq_[n-Njets_[0]])*Ntt+( -4*eb*eb+12*eb*(1-eb)-4*pow(1-eb,2))*pow(ebq_[n-Njets_[0]],2)*Ntt;
                double a3 = 20*(pow(1-ebq_[n-Njets_[0]],2)*Ntt)+(-30*eb+10*(1-eb))*ebq_[n-Njets_[0]]*(1-ebq_[n-Njets_[0]])*Ntt+( 12*eb*eb-24*eb*(1-eb)+4*pow(1-eb,2))*pow(ebq_[n-Njets_[0]],2)*Ntt;
                double a2 =                             ( 10*eb          )*ebq_[n-Njets_[0]]*(1-ebq_[n-Njets_[0]])*Ntt+(-12*eb*eb+12*eb*(1-eb)              )*pow(ebq_[n-Njets_[0]],2)*Ntt;
                double a1 =                                                                            (  4*eb*eb                           )*pow(ebq_[n-Njets_[0]],2)*Ntt;
                double a0 = -N3bjets+Nv_3bjets(Nv,euds,n)+Nvb_3bjets(Nvb,eb,euds,n);
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
double    VJetEstimation::eb_fromN2bjets(double &N2bjets, double &Ntt, double &Nv, double &Nvb, double &eb_old, double &eudsc, double &euds, int &n) const{
       	vector<double> roots;
	double a = (pow(1-eudsc,n-2) - 2*(n-2)*eudsc*pow(1-eudsc,n-3) + ((n-2)*(n-3)/2)*pow(eudsc,2)*pow(1-eudsc,n-4))*pow(ebq_[n-Njets_[0]],2)*Ntt;
       	double b = ((2*(n-2)*eudsc*pow(1-eudsc,n-3) - 2*((n-2)*(n-3)/2)*pow(eudsc,2)*pow(1-eudsc,n-4))*pow(ebq_[n-Njets_[0]],2) + ((n-1)*eudsc*pow(1-eudsc,n-2)-((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3))*ebq_[n-Njets_[0]]*(1-ebq_[n-Njets_[0]]))*Ntt+((n-1)*euds*pow(1-euds,n-2)-((n-1)*(n-2)/2)*euds*euds*pow(1-euds,n-3))*Nvb;
       	double c = (((n-2)*(n-3)/2)*pow(eudsc,2)*pow(1-eudsc,n-4)*pow(ebq_[n-Njets_[0]],2) + ((n-1)*(n-2)/2)*pow(eudsc,2)*pow(1-eudsc,n-3)*ebq_[n-Njets_[0]]*(1-ebq_[n-Njets_[0]]))*Ntt + ((n)*(n-1)/2)*pow(eudsc,2)*pow(1-eudsc,n-2)*(pow(1-ebq_[n-Njets_[0]],2)*Ntt) - N2bjets+Nv_2bjets(Nv,euds,n)+((n-1)*(n-2)/2)*euds*euds*pow(1-euds,n-3)*Nvb;

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
double   VJetEstimation::Ntt_fromN1bjet( double &N1bjet,  double &Nv,  double &Nvb, double &eb, double &eudsc, double &euds, int &n) const{
	double Ntt_Coef  =((2*eb*(1-eb)*pow(1-eudsc,n-2)+pow((1-eb),2)*(n-2)*eudsc*pow(1-eudsc,n-3))*pow(ebq_[n-Njets_[0]],2)
			 +          (eb*pow(1-eudsc,n-1)+       (1-eb)*(n-1)*eudsc*pow(1-eudsc,n-2))*2*ebq_[n-Njets_[0]]*(1-ebq_[n-Njets_[0]])
			 +                                                (n*eudsc*pow(1-eudsc,n-1))*pow(1-ebq_[n-Njets_[0]],2));
	double Nv_1b  = Nv_1bjet( Nv,  euds, n);
	double Nvb_1b = Nvb_1bjet(Nvb, eb, euds, n);
	double Ntt = 0;
	if(Ntt_Coef!=0) Ntt = (N1bjet-Nv_1b-Nvb_1b)/Ntt_Coef;
	return (Ntt>=0 ? Ntt : 0);
}
double   VJetEstimation::Nvb_fromN1bjet( double &N1bjet,  double &Ntt, double &Nv,  double &eb, double &eudsc, double &euds, int &n) const{
	double Nvb_Coef = (eb*pow(1-euds,n-1)+(1-eb)*(n-1)*euds*pow(1-euds,n-2));
	double Ntt_1b   = Ntt_1bjet(Ntt, eb, eudsc, n);
	double Nv_1b    = Nv_1bjet( Nv,  euds,  n);
	double Nvb      = 0;
	if(Nvb_Coef!=0)  Nvb = (N1bjet-Ntt_1b-Nv_1b)/Nvb_Coef;
	return (Nvb>=0 ? Nvb : 0);
}
double   VJetEstimation::Nv_fromN0bjet(double &N0bjet,  double &Ntt, double &Nvb, double &eb, double &eudsc, double &euds, int &n) const{
	double Nv_Coef = (pow(1-eudsc,n));
	double Ntt_0b  = Ntt_0bjet(Ntt, eb, eudsc, n);
	double Nvb_0b  = Nvb_0bjet(Nvb, eb, euds,  n);
	double Nv      = 0;
	if(Nv_Coef!=0)  Nv = (N0bjet-Ntt_0b-Nvb_0b)/Nv_Coef;
	return (Nv>=0 ? Nv : 0);
}
double VJetEstimation::euds_fromN0bjet(double &N0bjet, double &Ntt, double &Nv, double &Nvb, double &eb, double &eudsc, double &euds_old, int &n) const{
	double euds =0;
	vector<double> roots;
	if(n==4)/////// Equation holds only for n = 4  ///////////////
        {
                double a4 =    Nv;
                double a3 = -4*Nv -   (1-eb)*Nvb;
                double a2 =  6*Nv + 3*(1-eb)*Nvb;
                double a1 = -4*Nv - 3*(1-eb)*Nvb;
                double a0 =    Nv +   (1-eb)*Nvb - (N0bjet-Ntt_0bjet(Ntt, eb, eudsc, n));
                double par[5] = {a0,a1,a2,a3,a4};

                ROOT::Math::Polynomial poly(4);
                poly.SetParameters(par);

		roots = poly.FindRealRoots();
		for(unsigned int i = 0 ; i < roots.size(); i++)
                {
                        if(roots[i]<0 || roots[i]>1) continue;
                        if(fabs(euds_old-roots[i])<fabs(euds_old-euds)) euds = roots[i];
                }
                
		if(euds != 0) return euds;
                else return euds_old;
        }
        else if(n==5)/////// Equation holds only for n = 5  ///////////////
        {
                double a5 =    Nv;
                double a4 =  5*Nv +   (1-eb)*Nvb;
                double a3 =-10*Nv - 4*(1-eb)*Nvb;
                double a2 = 10*Nv + 6*(1-eb)*Nvb;
                double a1 = -5*Nv - 4*(1-eb)*Nvb;
                double a0 =    Nv +   (1-eb)*Nvb - (N0bjet-Ntt_0bjet(Ntt, eb, eudsc, n));
                double par[6] = {a0,a1,a2,a3,a4,a5};

                ROOT::Math::Polynomial poly(5);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
                for(unsigned int i = 0 ; i < roots.size(); i++)
                {
                        if(roots[i]<0 || roots[i]>1) continue;
                        if(fabs(euds_old-roots[i])<fabs(euds_old-euds)) euds = roots[i];
                }
                if(euds != 0) return euds;
                else return euds_old;
        }
        else /////// Equation holds only for n = 6 but is used for n>6 as well  ///////////////
        {
                double a6 =    Nv;
                double a5 = -6*Nv -   (1-eb)*Nvb;
                double a4 = 15*Nv + 5*(1-eb)*Nvb;
                double a3 =-20*Nv -10*(1-eb)*Nvb;
                double a2 = 15*Nv +10*(1-eb)*Nvb;
                double a1 = -6*Nv - 5*(1-eb)*Nvb;
                double a0 =    Nv +   (1-eb)*Nvb - (N0bjet-Ntt_0bjet(Ntt, eb, eudsc, n));
                double par[7] = {a0,a1,a2,a3,a4,a5,a6};

                ROOT::Math::Polynomial poly(6);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
                for(unsigned int i = 0 ; i < roots.size(); i++)
                {
                        if(roots[i]<0 || roots[i]>1) continue;
                        if(fabs(euds_old-roots[i])<fabs(euds_old-euds)) euds = roots[i];
                }
                if(euds != 0) return euds;
                else return euds_old;
        }

}

double VJetEstimation::euds_fromN1bjet(double &N1bjet, double &Ntt, double &Nv, double &Nvb, double &eb, double &eudsc, double &euds_old, int &n) const{
	double euds =0;
	vector<double> roots;
	if(n==4)/////// Equation holds only for n = 4  ///////////////
        {
                double a4 = -4*Nv;
                double a3 = 12*Nv + (-1*eb+3*(1-eb))*Nvb;
                double a2 =-12*Nv + ( 3*eb-6*(1-eb))*Nvb;
                double a1 =  4*Nv + (-3*eb+3*(1-eb))*Nvb;
                double a0 =        (( 1*eb         ))*Nvb - (N1bjet-Ntt_1bjet(Ntt, eb, eudsc, n));
                double par[5] = {a0,a1,a2,a3,a4};

                ROOT::Math::Polynomial poly(4);
                poly.SetParameters(par);
		roots = poly.FindRealRoots();
                
		for(unsigned int i = 0 ; i < roots.size(); i++)
                {
                        if(roots[i]<0 || roots[i]>1) continue;
                        if(fabs(euds_old-roots[i])<fabs(euds_old-euds)) euds = roots[i];
                }
                
		if(euds != 0) return euds;
                else return euds_old;
        }
        else if(n==5)/////// Equation holds only for n = 5  ///////////////
        {
                double a5 =  5*Nv;
                double a4 =-20*Nv + (   eb- 4*(1-eb))*Nvb;
                double a3 = 30*Nv + (-4*eb+12*(1-eb))*Nvb;
                double a2 =-20*Nv + ( 6*eb-12*(1-eb))*Nvb;
                double a1 =  5*Nv + (-4*eb+ 4*(1-eb))*Nvb;
                double a0 =        ((   eb         ))*Nvb - (N1bjet-Ntt_1bjet(Ntt, eb, eudsc, n));
                double par[6] = {a0,a1,a2,a3,a4,a5};

                ROOT::Math::Polynomial poly(5);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
                for(unsigned int i = 0 ; i < roots.size(); i++)
                {
                        if(roots[i]<0 || roots[i]>1) continue;
                        if(fabs(euds_old-roots[i])<fabs(euds_old-euds)) euds = roots[i];
                }
                if(euds != 0) return euds;
                else return euds_old;
        }
        else /////// Equation holds only for n = 6 but is used for n>6 as well  ///////////////
        {
                double a6 = -6*Nv;
                double a5 = 30*Nv + ( -1*eb +  5*(1-eb))*Nvb;
                double a4 =-60*Nv + (  5*eb - 20*(1-eb))*Nvb;
                double a3 = 60*Nv + (-10*eb + 30*(1-eb))*Nvb;
                double a2 =-30*Nv + ( 10*eb - 20*(1-eb))*Nvb;
                double a1 = -6*Nv + ( -5*eb +  5*(1-eb))*Nvb;
                double a0 =         (    eb            )*Nvb - (N1bjet-Ntt_1bjet(Ntt, eb, eudsc, n));
                double par[7] = {a0,a1,a2,a3,a4,a5,a6};

                ROOT::Math::Polynomial poly(6);
                poly.SetParameters(par);

                roots = poly.FindRealRoots();
                for(unsigned int i = 0 ; i < roots.size(); i++)
                {
                        if(roots[i]<0 || roots[i]>1) continue;
                        if(fabs(euds_old-roots[i])<fabs(euds_old-euds)) euds = roots[i];
                }
                if(euds != 0) return euds;
                else return euds_old;
        }

}


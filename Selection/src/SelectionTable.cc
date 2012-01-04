#include "../interface/SelectionTable.h"

using namespace std;

  
SelectionTable::SelectionTable(){
	Lumi_ = 100;
	merge_ = false;
}

SelectionTable::SelectionTable(vector<string> listOfCuts, vector<Dataset*> listOfDatasets){
	Lumi_ = 100;
	merge_ = false;
	listOfCuts_ = listOfCuts;
	listOfDatasets_ = listOfDatasets;
	nofEventsRaw_                = new float*[listOfCuts.size()];
	nofEventsRawError_           = new float*[listOfCuts.size()];
	nofEventsRawErrorHigh_       = new float*[listOfCuts.size()];
	nofEventsRawErrorLow_        = new float*[listOfCuts.size()];
	nofEvents_                   = new float*[listOfCuts.size()];
	nofEventsError_              = new float*[listOfCuts.size()];
	nofEventsExpErrorHigh_       = new float*[listOfCuts.size()];
	nofEventsExpErrorLow_        = new float*[listOfCuts.size()];
	nofEventsMcErrorHigh_        = new float*[listOfCuts.size()];
	nofEventsMcErrorLow_         = new float*[listOfCuts.size()];
	cutEfficiency_               = new float*[listOfCuts.size()];
	cutEfficiencyErrorHigh_      = new float*[listOfCuts.size()];
	cutEfficiencyErrorLow_       = new float*[listOfCuts.size()];
	cutEfficiencyError_          = new float*[listOfCuts.size()];
	totalCutEfficiency_          = new float*[listOfCuts.size()];
	totalCutEfficiencyErrorHigh_ = new float*[listOfCuts.size()];
	totalCutEfficiencyErrorLow_  = new float*[listOfCuts.size()];
	totalCutEfficiencyError_     = new float*[listOfCuts.size()];
  for(unsigned int i=0;i<listOfCuts.size();i++){
		nofEventsRaw_[i]                = new float[listOfDatasets.size()];
		nofEventsRawError_[i]           = new float[listOfDatasets.size()];
		nofEventsRawErrorHigh_[i]       = new float[listOfDatasets.size()];
		nofEventsRawErrorLow_[i]        = new float[listOfDatasets.size()];
		nofEvents_[i]                   = new float[listOfDatasets.size()];
		nofEventsError_[i]              = new float[listOfDatasets.size()];
		nofEventsExpErrorHigh_[i]       = new float[listOfDatasets.size()];
		nofEventsExpErrorLow_[i]        = new float[listOfDatasets.size()];
		nofEventsMcErrorHigh_[i]        = new float[listOfDatasets.size()];
		nofEventsMcErrorLow_[i]         = new float[listOfDatasets.size()];
		cutEfficiency_[i]               = new float[listOfDatasets.size()];
		cutEfficiencyError_[i]          = new float[listOfDatasets.size()];
		cutEfficiencyErrorHigh_[i]      = new float[listOfDatasets.size()];
		cutEfficiencyErrorLow_[i]       = new float[listOfDatasets.size()];
		totalCutEfficiency_[i]          = new float[listOfDatasets.size()];
		totalCutEfficiencyErrorHigh_[i] = new float[listOfDatasets.size()];
		totalCutEfficiencyErrorLow_[i]  = new float[listOfDatasets.size()];
		totalCutEfficiencyError_[i]     = new float[listOfDatasets.size()];
		for(unsigned int j=0;j<listOfDatasets.size();j++){
			nofEventsRaw_[i][j]                = 0;
			nofEventsRawErrorHigh_[i][j]       = 0;
			nofEventsRawErrorLow_[i][j]        = 0;
			nofEvents_[i][j]                   = 0;
			nofEventsError_[i][j]              = 0;
			nofEventsExpErrorHigh_[i][j]       = 0;
			nofEventsExpErrorLow_[i][j]        = 0;
			nofEventsMcErrorHigh_[i][j]        = 0;
			nofEventsMcErrorLow_[i][j]         = 0;
			cutEfficiency_[i][j]               = 0;
			cutEfficiencyErrorHigh_[i][j]      = 0;
			cutEfficiencyErrorLow_[i][j]       = 0;
			cutEfficiencyError_[i][j]          = 0;
			totalCutEfficiency_[i][j]          = 0;
			totalCutEfficiencyErrorHigh_[i][j] = 0;
			totalCutEfficiencyErrorLow_[i][j]  = 0;
			totalCutEfficiencyError_[i][j]     = 0;
		}
	}
}

SelectionTable::SelectionTable(const SelectionTable& table){
  Lumi_                        = table.Lumi_;
  listOfCuts_                  = table.listOfCuts_;
  listOfDatasets_              = table.listOfDatasets_;
  nofEventsRaw_                = table.nofEventsRaw_;
  nofEventsRawErrorHigh_       = table.nofEventsRawError_;
  nofEventsRawErrorLow_        = table.nofEventsRawError_;
  nofEvents_                   = table.nofEvents_;
  nofEventsError_              = table.nofEventsError_;
  nofEventsExpErrorHigh_       = table.nofEventsExpErrorHigh_;
  nofEventsExpErrorLow_        = table.nofEventsExpErrorLow_;
  nofEventsMcErrorHigh_        = table.nofEventsMcErrorHigh_;
  nofEventsMcErrorLow_         = table.nofEventsMcErrorLow_;
  cutEfficiency_               = table.cutEfficiency_;
  cutEfficiencyErrorHigh_      = table.cutEfficiencyErrorHigh_;
  cutEfficiencyErrorLow_       = table.cutEfficiencyErrorLow_;
  cutEfficiencyError_          = table.cutEfficiencyError_;
  totalCutEfficiency_          = table.totalCutEfficiency_;
  totalCutEfficiencyErrorHigh_ = table.totalCutEfficiencyError_;
  totalCutEfficiencyErrorLow_  = table.totalCutEfficiencyError_;
  totalCutEfficiencyError_     = table.totalCutEfficiencyError_;
}

SelectionTable::~SelectionTable(){
}

	    

void SelectionTable::Scale(float Lumi){
  float factor = 1;
  for(unsigned int i=0;i<listOfCuts_.size();i++){
		for(unsigned int j=0;j<listOfDatasets_.size();j++){
			factor = listOfDatasets_[j]->NormFactor();
			//nofEventsRaw_[i][j] = nofEventsRaw_[i][j]*Lumi;
			nofEvents_[i][j]      = nofEventsRaw_[i][j]*factor*Lumi;
			//nofEventsError_[i][j] = ErrorCalculator(nofEvents_[i][j],nofEvents_[i][j]/nofEvents_[0][j],1);
			//cutEfficiency_[i][j] = cutEfficiency_[i][j]*Lumi;
			//cutEfficiencyError_[i][j] = cutEfficiencyError_[i][j]*Lumi;
			//totalCutEfficiency_[i][j] = totalCutEfficiency_[i][j]*Lumi;
			//totalCutEfficiencyError_[i][j] = totalCutEfficiencyError_[i][j]*Lumi;
		}
	}
}

void SelectionTable::TableCalculator(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
{
  if(mergeTT || mergeQCD || mergeW || mergeZ || mergeST) merge_ = true;

  //scale to lumi
  Scale(Lumi_);
  bool addedTT= false, addedQCD = false, addedW = false, addedZ = false, addedSingleTop = false;

  //loop on cuts
  for(unsigned int i=0;i<listOfCuts_.size();i++)
  {
	//loop on datasets
  	for(unsigned int j=0;j<listOfDatasets_.size();j++)
  	{
		float factor = listOfDatasets_[j]->NormFactor();
		//Errors on the raw numbers of events (Wald formula)
		nofEventsRawError_[i][j]     = ErrorCalculator(nofEventsRaw_[i][j],nofEventsRaw_[i][j]/nofEventsRaw_[0][j],1);
		// Asymmetric errors on the raw numbers of events (Wilson Score formula).
		nofEventsRawErrorHigh_[i][j] = (int)(nofEventsRaw_[0][j]*WilsonScoreIntervalHigh(nofEventsRaw_[i][j],nofEventsRaw_[0][j])-nofEventsRaw_[i][j]+0.5);
		nofEventsRawErrorLow_[i][j]  = (int)(nofEventsRaw_[i][j]-nofEventsRaw_[0][j]*WilsonScoreIntervalLow( nofEventsRaw_[i][j],nofEventsRaw_[0][j])+0.5);
		// Errors on the rescaled numbers of events. (Wald formula)
  		nofEventsError_[i][j]        = ErrorCalculator(nofEvents_[i][j],nofEvents_[i][j]/nofEvents_[0][j],1);
		// Asymmetric errors on the rescaled numbers of events (Wilson Score formula).
		// - expected errors for a given iLumi :
  		nofEventsExpErrorHigh_[i][j] = nofEvents_[0][j]*WilsonScoreIntervalHigh(nofEvents_[i][j],nofEvents_[0][j])-nofEvents_[i][j];
  		nofEventsExpErrorLow_[i][j]  = nofEvents_[i][j]-nofEvents_[0][j]*WilsonScoreIntervalLow( nofEvents_[i][j],nofEvents_[0][j]);
		// - rescaled errors from MC :
  		nofEventsMcErrorHigh_[i][j]  = nofEventsRawErrorHigh_[i][j]*factor*Lumi_;
  		nofEventsMcErrorLow_[i][j]   = nofEventsRawErrorLow_[i][j]*factor*Lumi_;
   		//push_back
  		if(i>0){
  			cutEfficiency_[i][j] =  (float)(nofEventsRaw_[i][j]/nofEventsRaw_[i-1][j]);
  			cutEfficiencyErrorHigh_[i][j] =  WilsonScoreIntervalHigh(nofEventsRaw_[i][j],nofEventsRaw_[i-1][j])-(nofEventsRaw_[i][j]/nofEventsRaw_[i-1][j]);
  			cutEfficiencyErrorLow_[i][j]  =  (nofEventsRaw_[i][j]/nofEventsRaw_[i-1][j])-WilsonScoreIntervalLow(nofEventsRaw_[i][j],nofEventsRaw_[i-1][j]);
  			cutEfficiencyError_[i][j] =  BinomialeError((float)(nofEventsRaw_[i][j]/nofEventsRaw_[i-1][j]),nofEventsRaw_[i-1][j]);
  		}
  		else {
  			cutEfficiency_[i][j] = 1.;
  			cutEfficiencyError_[i][j] = 0.;
			
  			string datasetName = listOfDatasets_[j]->Name();

  			if(mergeTT && datasetName.find("TT_") == 0) {
  			  if(!addedTT) {
  			    listOfDatasetsMerged_.push_back(new Dataset("TT_Merged","t#bar{t}+jets, merged",true,listOfDatasets_[j]->Color(),listOfDatasets_[j]->LineStyle(),listOfDatasets_[j]->LineWidth(),1.,1.));
  			    addedTT = true;
  			  }
  			}
  			else if(mergeQCD && datasetName.find("QCD_") == 0) {
  			  if(!addedQCD) {
  			    listOfDatasetsMerged_.push_back(new Dataset("QCD_Merged","Multi-jets, merged",true,listOfDatasets_[j]->Color(),listOfDatasets_[j]->LineStyle(),listOfDatasets_[j]->LineWidth(),1.,1.));
  			    addedQCD = true;
  			  }
  			}
  			else if(mergeW && datasetName.find("W_") == 0) {
  			  if(!addedW) {
  			    listOfDatasetsMerged_.push_back(new Dataset("WJets_Merged","W+jets, merged",true,listOfDatasets_[j]->Color(),listOfDatasets_[j]->LineStyle(),listOfDatasets_[j]->LineWidth(),1.,1.));
  			    addedW = true;
  			  }
  			}
  			else if(mergeZ && datasetName.find("Z_") == 0) {
  			  if(!addedZ) {
  			    listOfDatasetsMerged_.push_back(new Dataset("ZJets_Merged","Z+jets, merged",true,listOfDatasets_[j]->Color(),listOfDatasets_[j]->LineStyle(),listOfDatasets_[j]->LineWidth(),1.,1.));
  			    addedZ = true;
  			  }
  			}
  			else if(mergeST && datasetName.find("ST_") == 0) {
  			  if(!addedSingleTop) {
  			    listOfDatasetsMerged_.push_back(new Dataset("ST_Merged","st+jets, merged",true,listOfDatasets_[j]->Color(),listOfDatasets_[j]->LineStyle(),listOfDatasets_[j]->LineWidth(),1.,1.));
  			    addedSingleTop = true;
  			  }
  			}
			else listOfDatasetsMerged_.push_back(listOfDatasets_[j]);
  		}
  		totalCutEfficiency_[i][j] = (float)(nofEventsRaw_[i][j]/nofEventsRaw_[0][j]);
  		totalCutEfficiencyError_[i][j] =  BinomialeError((float)(nofEventsRaw_[i][j]/nofEventsRaw_[0][j]),nofEventsRaw_[0][j]);
  		totalCutEfficiencyErrorHigh_[i][j] =  WilsonScoreIntervalHigh(nofEventsRaw_[i][j],nofEventsRaw_[0][j])-(nofEventsRaw_[i][j]/nofEventsRaw_[0][j]);
  		totalCutEfficiencyErrorLow_[i][j]  =  (nofEventsRaw_[i][j]/nofEventsRaw_[0][j])-WilsonScoreIntervalLow(nofEventsRaw_[i][j],nofEventsRaw_[0][j]);
	}
	if(merge_)
	{
  	  // initialize everything with zero's and then fill it!
  	  nofEventsRawMerged_          = new float*[listOfCuts_.size()];
  	  nofEventsRawErrorMerged_     = new float*[listOfCuts_.size()];
  	  nofEventsRawErrorHighMerged_ = new float*[listOfCuts_.size()];
  	  nofEventsRawErrorLowMerged_  = new float*[listOfCuts_.size()];

  	  nofEventsMerged_             = new float*[listOfCuts_.size()];
  	  nofEventsErrorMerged_        = new float*[listOfCuts_.size()];
  	  nofEventsExpErrorHighMerged_ = new float*[listOfCuts_.size()];
  	  nofEventsExpErrorLowMerged_  = new float*[listOfCuts_.size()];
  	  nofEventsMcErrorHighMerged_  = new float*[listOfCuts_.size()];
  	  nofEventsMcErrorLowMerged_   = new float*[listOfCuts_.size()];

  	  cutEfficiencyMerged_          = new float*[listOfCuts_.size()];
  	  cutEfficiencyErrorMerged_     = new float*[listOfCuts_.size()];
  	  cutEfficiencyErrorHighMerged_ = new float*[listOfCuts_.size()];
  	  cutEfficiencyErrorLowMerged_  = new float*[listOfCuts_.size()];

  	  totalCutEfficiencyMerged_          = new float*[listOfCuts_.size()];
  	  totalCutEfficiencyErrorMerged_     = new float*[listOfCuts_.size()];
  	  totalCutEfficiencyErrorHighMerged_ = new float*[listOfCuts_.size()];
  	  totalCutEfficiencyErrorLowMerged_  = new float*[listOfCuts_.size()];

  	  for(unsigned int i=0;i<listOfCuts_.size();i++){
  	  	nofEventsRawMerged_[i]          = new float[listOfDatasetsMerged_.size()];
  	  	nofEventsRawErrorMerged_[i]     = new float[listOfDatasetsMerged_.size()];
  	  	nofEventsRawErrorHighMerged_[i] = new float[listOfDatasetsMerged_.size()];
  	  	nofEventsRawErrorLowMerged_[i]  = new float[listOfDatasetsMerged_.size()];

  	  	nofEventsMerged_[i]             = new float[listOfDatasetsMerged_.size()];
  	  	nofEventsErrorMerged_[i]        = new float[listOfDatasetsMerged_.size()];
  	  	nofEventsExpErrorHighMerged_[i] = new float[listOfDatasetsMerged_.size()];
  	  	nofEventsExpErrorLowMerged_[i]  = new float[listOfDatasetsMerged_.size()];
  	  	nofEventsMcErrorHighMerged_[i]  = new float[listOfDatasetsMerged_.size()];
  	  	nofEventsMcErrorLowMerged_[i]   = new float[listOfDatasetsMerged_.size()];

  	  	cutEfficiencyMerged_[i]          = new float[listOfDatasetsMerged_.size()];
  	  	cutEfficiencyErrorMerged_[i]     = new float[listOfDatasetsMerged_.size()];
  	  	cutEfficiencyErrorHighMerged_[i] = new float[listOfDatasetsMerged_.size()];
  	  	cutEfficiencyErrorLowMerged_[i]  = new float[listOfDatasetsMerged_.size()];

  	  	totalCutEfficiencyMerged_[i]          = new float[listOfDatasetsMerged_.size()];
  	  	totalCutEfficiencyErrorMerged_[i]     = new float[listOfDatasetsMerged_.size()];
  	  	totalCutEfficiencyErrorHighMerged_[i] = new float[listOfDatasetsMerged_.size()];
  	  	totalCutEfficiencyErrorLowMerged_[i]  = new float[listOfDatasetsMerged_.size()];

  	  	for(unsigned int j=0;j<listOfDatasetsMerged_.size();j++){
  	  		nofEventsRawMerged_[i][j]          = 0;
  	  		nofEventsRawErrorMerged_[i][j]     = 0;
  	  		nofEventsRawErrorHighMerged_[i][j] = 0;
  	  		nofEventsRawErrorLowMerged_[i][j]  = 0;

	      		nofEventsMerged_[i][j]             = 0;
      			nofEventsErrorMerged_[i][j]        = 0;
      			nofEventsExpErrorHighMerged_[i][j] = 0;
      			nofEventsExpErrorLowMerged_[i][j]  = 0;
      			nofEventsMcErrorHighMerged_[i][j]  = 0;
      			nofEventsMcErrorLowMerged_[i][j]   = 0;

      			cutEfficiencyMerged_[i][j]          = 0;
      			cutEfficiencyErrorMerged_[i][j]     = 0;
      			cutEfficiencyErrorHighMerged_[i][j] = 0;
      			cutEfficiencyErrorLowMerged_[i][j]  = 0;

      			totalCutEfficiencyMerged_[i][j]          = 0;
      			totalCutEfficiencyErrorMerged_[i][j]     = 0;
      			totalCutEfficiencyErrorHighMerged_[i][j] = 0;
      			totalCutEfficiencyErrorLowMerged_[i][j]  = 0;

      			for(unsigned int k=0;k<listOfDatasets_.size();k++){
      				string datasetName = listOfDatasets_[k]->Name(), datasetMergedName = listOfDatasetsMerged_[j]->Name();
				if( (mergeTT  && datasetName.find("TT_")  == 0 && datasetMergedName.find("TT_Merged")    == 0) ||
				    (mergeQCD && datasetName.find("QCD_") == 0 && datasetMergedName.find("QCD_Merged")   == 0) ||
				    (mergeW   && datasetName.find("W_")   == 0 && datasetMergedName.find("WJets_Merged") == 0) ||
				    (mergeZ   && datasetName.find("Z_")   == 0 && datasetMergedName.find("ZJets_Merged") == 0) ||
				    (mergeST  && datasetName.find("ST_")  == 0 && datasetMergedName.find("ST_Merged")    == 0) ) {
					nofEventsRawMerged_[i][j]      += nofEventsRaw_[i][k];
					//nofEventsRawErrorMerged_[i][j] = sqrt( pow(nofEventsRawErrorMerged_[i][j],2.) + pow(nofEventsRawError_[i][k],2.);
					nofEventsMerged_[i][j]     += nofEvents_[i][k];
					nofEventsErrorMerged_[i][j] = sqrt( pow(nofEventsErrorMerged_[i][j],2) + pow(nofEventsError_[i][k],2) );
					//nofEventsExpErrorHighMerged_[i][j]
					//nofEventsExpErrorLowMerged_[i][j]
					nofEventsMcErrorHighMerged_[i][j]=sqrt(pow(nofEventsMcErrorHighMerged_[i][j],2)+pow(nofEventsMcErrorHigh_[i][k],2));
					nofEventsMcErrorLowMerged_[i][j] =sqrt(pow(nofEventsMcErrorLowMerged_[i][j],2) +pow(nofEventsMcErrorLow_[i][k],2));
					if(i == 0) {
						cutEfficiencyMerged_[i][j]      = 1.;
  			        		cutEfficiencyErrorMerged_[i][j] = 0.;
      		    			}
      		    			else {
      		      				cutEfficiencyMerged_[i][j]      =  (float)(nofEventsMerged_[i][j]/nofEventsMerged_[i-1][j]);
  			        		cutEfficiencyErrorMerged_[i][j] =  BinomialeError((float)(nofEventsMerged_[i][j]/nofEventsMerged_[i-1][j]),nofEventsMerged_[i-1][j]);
  						cutEfficiencyErrorHighMerged_[i][j] =  WilsonScoreIntervalHigh(nofEventsMerged_[i][j],nofEventsMerged_[i-1][j])-(nofEventsMerged_[i][j]/nofEventsMerged_[i-1][j]);
  						cutEfficiencyErrorLowMerged_[i][j]  =  (nofEventsMerged_[i][j]/nofEventsMerged_[i-1][j])-WilsonScoreIntervalLow(nofEventsMerged_[i][j],nofEventsMerged_[i-1][j]);
      		    			}
      		    			totalCutEfficiencyMerged_[i][j] = (float)(nofEventsMerged_[i][j]/nofEventsMerged_[0][j]);
          				totalCutEfficiencyErrorMerged_[i][j] =  BinomialeError((float)(nofEventsMerged_[i][j]/nofEventsMerged_[0][j]),nofEventsRawMerged_[0][j]);
  					totalCutEfficiencyErrorHigh_[i][j] =  WilsonScoreIntervalHigh(nofEventsMerged_[i][j],nofEventsMerged_[0][j])-(nofEventsMerged_[i][j]/nofEventsMerged_[0][j]);
  					totalCutEfficiencyErrorLow_[i][j]  =  (nofEventsMerged_[i][j]/nofEventsMerged_[0][j])-WilsonScoreIntervalLow(nofEventsMerged_[i][j],nofEventsMerged_[0][j]);      		  		}
      		  		else if(datasetMergedName.compare(datasetName) == 0)
      		  		{
      		    			nofEventsRawMerged_[i][j] = nofEventsRaw_[i][k];
      		    			nofEventsMerged_[i][j] = nofEvents_[i][k];
      		    			nofEventsErrorMerged_[i][j] = nofEventsError_[i][k];
      		    			cutEfficiencyMerged_[i][j] = cutEfficiency_[i][k];
      		    			cutEfficiencyErrorMerged_[i][j] = cutEfficiencyError_[i][k];
      		    			totalCutEfficiencyMerged_[i][j] = totalCutEfficiency_[i][k];
      		    			totalCutEfficiencyErrorMerged_[i][j] = totalCutEfficiencyError_[i][k];
      		  		}
  				nofEventsExpErrorHighMerged_[i][j] = nofEventsMerged_[0][j]*WilsonScoreIntervalHigh(nofEventsMerged_[i][j],nofEventsMerged_[0][j])-nofEventsMerged_[i][j];
  				nofEventsExpErrorLowMerged_[i][j]  = nofEventsMerged_[i][j]-nofEventsMerged_[0][j]*WilsonScoreIntervalLow( nofEventsMerged_[i][j],nofEventsMerged_[0][j]);
      			}
      		}
  	  }
  	}
  }
}


void SelectionTable::Fill(unsigned int DatasetNumber, unsigned int CutNumber, float value){
	if(DatasetNumber>=listOfDatasets_.size()) cerr<<"DataNumber argument wrong !"<<endl;
	if(CutNumber>=listOfCuts_.size()) cerr<<"CutNumber argument wrong !"<<endl;
	nofEventsRaw_[CutNumber][DatasetNumber]+= value;
}

void SelectionTable::Fill( unsigned int DatasetNumber, vector<float> PassTheCuts){
	if(DatasetNumber>=listOfDatasets_.size()) cerr<<"DataNumber argument wrong !"<<endl;
	if(PassTheCuts.size()>listOfCuts_.size()) cerr<<"CutNumber argument wrong !"<<endl;
	for(unsigned int i=0;i<PassTheCuts.size();i++) nofEventsRaw_[i][DatasetNumber]+= PassTheCuts[i];
}

void SelectionTable::Write(string filename, bool WithError){
	ofstream fout(filename.c_str(), ios::out);
	SelectionTable::Write(fout, WithError);
}

void SelectionTable::WriteTable(ofstream& fout, float** listTable_, bool writeMerged){
	fout<<"&";
	unsigned int nDatasets = listOfDatasets_.size();
	if(writeMerged) nDatasets = listOfDatasetsMerged_.size();
	for(unsigned int i=0;i<nDatasets;i++) {
	  if(writeMerged) fout<<"$"<<listOfDatasetsMerged_[i]->Title()<<"$\t";
	  else fout<<"$"<<listOfDatasets_[i]->Title()<<"$\t";
		if(i<nDatasets-1) fout<<"&";
		else fout<<"\\\\"<<endl;
	}
	fout<<endl;
	fout<<"\\hline"<<endl;
	for(unsigned int i=0;i<listOfCuts_.size();i++){
		fout<<listOfCuts_[i]<<"&\t";
		for(unsigned int j=0;j<nDatasets;j++){
			fout<<listTable_[i][j]<<"\t";	
			if(j<nDatasets-1) fout<<"&";
			else fout<<"\\\\"<<endl;
		}
		fout<<endl;
	}
}

void SelectionTable::WriteTable(ofstream& fout, float** listTable_,float** listTableError_, bool writeMerged){
	fout<<"&";
	unsigned int nDatasets = listOfDatasets_.size();
	if(writeMerged) nDatasets = listOfDatasetsMerged_.size();
	for(unsigned int i=0;i<nDatasets;i++) {
		if(writeMerged) fout<<"$"<<listOfDatasetsMerged_[i]->Title()<<"$\t";
		else fout<<"$"<<listOfDatasets_[i]->Title()<<"$\t";
		if(i<nDatasets-1) fout<<"&";
		else fout<<"\\\\"<<endl;
	}
	fout<<endl;
	fout<<"\\hline"<<endl;
	for(unsigned int i=0;i<listOfCuts_.size();i++){
		fout<<listOfCuts_[i]<<"&\t";
		for(unsigned int j=0;j<nDatasets;j++){
			fout<<listTable_[i][j]<<" $\\pm$ "<<listTableError_[i][j]<<"\t";	
			if(j<nDatasets-1) fout<<"&";
			else fout<<"\\\\"<<endl;
		}
		fout<<endl;
	}
}

void SelectionTable::WriteTable(ofstream& fout, float** listTable_,float** listTableErrorHigh_ ,float** listTableErrorLow_, bool writeMerged){
	fout<<"&";
	unsigned int nDatasets = listOfDatasets_.size();
	if(writeMerged) nDatasets = listOfDatasetsMerged_.size();
	for(unsigned int i=0;i<nDatasets;i++) {
		if(writeMerged) fout<<"$"<<listOfDatasetsMerged_[i]->Title()<<"$\t";
	  else fout<<"$"<<listOfDatasets_[i]->Title()<<"$\t";
		if(i<nDatasets-1) fout<<"&";
		else fout<<"\\\\"<<endl;
	}
	fout<<endl;
	fout<<"\\hline"<<endl;
	for(unsigned int i=0;i<listOfCuts_.size();i++){
		fout<<listOfCuts_[i]<<"&\t";
		for(unsigned int j=0;j<nDatasets;j++){
			fout<<"$"<<listTable_[i][j]<<"^{+"<<listTableErrorHigh_[i][j]<<"}_{-"<<listTableErrorLow_[i][j]<<"}$"<<"\t";	
			if(j<nDatasets-1) fout<<"&";
			else fout<<"\\\\"<<endl;
		}
		fout<<endl;
	}
}

void SelectionTable::Write(ofstream& fout, bool WithError){
	if(WithError){
  	////////////////////////////
		fout<<"\\begin{table}"<<endl;
  	fout<<"\\centering"<<endl;
  	fout<<"\\begin{tabular}{|c|";
		for(unsigned int i=0;i<listOfDatasets_.size();i++) fout<<"c";
	fout<<"|}"<<endl;
  	fout<<"\\hline"<<endl;
		WriteTable(fout, nofEventsRaw_, nofEventsRawErrorHigh_, nofEventsRawErrorLow_);
  	fout<<"\\hline"<<endl;
  	fout<<"\\end{tabular}"<<endl;
  	fout<<"\\caption{Number of events}"<<endl;
  	fout<<"\\label{tab:SelectionTable}"<<endl;
  	fout<<"\\end{table}"<<endl;
  	////////////////////////////
		fout<<"\\begin{table}"<<endl;
  	fout<<"\\centering"<<endl;
  	fout<<"\\begin{tabular}{|c|";
		for(unsigned int i=0;i<listOfDatasets_.size();i++) fout<<"c";
	fout<<"|}"<<endl;
  	fout<<"\\hline"<<endl;
		WriteTable(fout, nofEvents_, nofEventsMcErrorHigh_, nofEventsMcErrorLow_);
  	fout<<"\\hline"<<endl;
  	fout<<"\\end{tabular}"<<endl;
  	fout<<"\\caption{Expected numbers of events for $"<<Lumi_<<"~pb^{-1}$ of integrated luminosity as calculated from Monte Carlo simulations }"<<endl;
  	fout<<"\\label{tab:SelectionTable}"<<endl;
  	fout<<"\\end{table}"<<endl;
  	////////////////////////////
		fout<<"\\begin{table}"<<endl;
  	fout<<"\\centering"<<endl;
  	fout<<"\\begin{tabular}{|c|";
		for(unsigned int i=0;i<listOfDatasets_.size();i++) fout<<"c";
	fout<<"|}"<<endl;
  	fout<<"\\hline"<<endl;
		WriteTable(fout, nofEvents_, nofEventsExpErrorHigh_, nofEventsExpErrorLow_);
  	fout<<"\\hline"<<endl;
  	fout<<"\\end{tabular}"<<endl;
  	fout<<"\\caption{Expected numbers of events and errors for $"<<Lumi_<<"~pb^{-1}$ of integrated luminosity }"<<endl;
  	fout<<"\\label{tab:SelectionTable}"<<endl;
  	fout<<"\\end{table}"<<endl;
  	////////////////////////////
		fout<<"\\begin{table}"<<endl;
  	fout<<"\\centering"<<endl;
  	fout<<"\\begin{tabular}{|c|";
		for(unsigned int i=0;i<listOfDatasets_.size();i++) fout<<"c";
		fout<<"|}"<<endl;
  	fout<<"\\hline"<<endl;
		WriteTable(fout, cutEfficiency_, cutEfficiencyErrorHigh_, cutEfficiencyErrorLow_);
  	fout<<"\\hline"<<endl;
  	fout<<"\\end{tabular}"<<endl;
  	fout<<"\\caption{Cut efficiencies}"<<endl;
  	fout<<"\\label{tab:EfficienciesTable}"<<endl;
  	fout<<"\\end{table}"<<endl;
  	////////////////////////////
		fout<<"\\begin{table}"<<endl;
  	fout<<"\\centering"<<endl;
  	fout<<"\\begin{tabular}{|c|";
		for(unsigned int i=0;i<listOfDatasets_.size();i++) fout<<"c";
		fout<<"|}"<<endl;
  	fout<<"\\hline"<<endl;
		WriteTable(fout, cutEfficiency_, cutEfficiencyError_);
  	fout<<"\\hline"<<endl;
  	fout<<"\\end{tabular}"<<endl;
  	fout<<"\\caption{Cut efficiencies}"<<endl;
  	fout<<"\\label{tab:EfficienciesTable}"<<endl;
  	fout<<"\\end{table}"<<endl;
  	////////////////////////////
		fout<<"\\begin{table}"<<endl;
  	fout<<"\\centering"<<endl;
  	fout<<"\\begin{tabular}{|c|";
		for(unsigned int i=0;i<listOfDatasets_.size();i++) fout<<"c";
		fout<<"|}"<<endl;
  	fout<<"\\hline"<<endl;
		WriteTable(fout, totalCutEfficiency_, totalCutEfficiencyErrorHigh_, totalCutEfficiencyErrorLow_);
  	fout<<"\\hline"<<endl;
  	fout<<"\\end{tabular}"<<endl;
  	fout<<"\\caption{Total cut efficiencies}"<<endl;
  	fout<<"\\label{tab:TotalEfficienciesTable}"<<endl;
  	fout<<"\\end{table}"<<endl;
  	////////////////////////////
		fout<<"\\begin{table}"<<endl;
  	fout<<"\\centering"<<endl;
  	fout<<"\\begin{tabular}{|c|";
		for(unsigned int i=0;i<listOfDatasets_.size();i++) fout<<"c";
		fout<<"|}"<<endl;
  	fout<<"\\hline"<<endl;
		WriteTable(fout, totalCutEfficiency_, totalCutEfficiencyError_);
  	fout<<"\\hline"<<endl;
  	fout<<"\\end{tabular}"<<endl;
  	fout<<"\\caption{Total cut efficiencies}"<<endl;
  	fout<<"\\label{tab:TotalEfficienciesTable}"<<endl;
  	fout<<"\\end{table}"<<endl;
  	if(merge_) {
    	////////////////////////////
      fout<<"\\begin{table}"<<endl;
      fout<<"\\centering"<<endl;
    	fout<<"\\begin{tabular}{|c|";
      for(unsigned int i=0;i<listOfDatasetsMerged_.size();i++) fout<<"c";
      fout<<"|}"<<endl;
  	  fout<<"\\hline"<<endl;
  	  WriteTable(fout, nofEventsMerged_, nofEventsErrorMerged_, true);
      fout<<"\\hline"<<endl;
    	fout<<"\\end{tabular}"<<endl;
    	fout<<"\\caption{Number of events, Merged}"<<endl;
    	fout<<"\\label{tab:SelectionTableMerged}"<<endl;
    	fout<<"\\end{table}"<<endl;
    	////////////////////////////
      fout<<"\\begin{table}"<<endl;
      fout<<"\\centering"<<endl;
    	fout<<"\\begin{tabular}{|c|";
      for(unsigned int i=0;i<listOfDatasetsMerged_.size();i++) fout<<"c";
      fout<<"|}"<<endl;
  	  fout<<"\\hline"<<endl;
  	  WriteTable(fout, nofEventsMerged_, nofEventsMcErrorHighMerged_, nofEventsMcErrorLowMerged_, true);
      fout<<"\\hline"<<endl;
    	fout<<"\\end{tabular}"<<endl;
    	fout<<"\\caption{Expected numbers of events for $"<<Lumi_<<"~pb^{-1}$ of integrated luminosity as calculated from Monte Carlo simulations, Merged}"<<endl;
    	fout<<"\\label{tab:SelectionTableMerged}"<<endl;
    	fout<<"\\end{table}"<<endl;
    	////////////////////////////
      fout<<"\\begin{table}"<<endl;
      fout<<"\\centering"<<endl;
    	fout<<"\\begin{tabular}{|c|";
      for(unsigned int i=0;i<listOfDatasetsMerged_.size();i++) fout<<"c";
      fout<<"|}"<<endl;
  	  fout<<"\\hline"<<endl;
  	  WriteTable(fout, nofEventsMerged_, nofEventsExpErrorHighMerged_, nofEventsExpErrorLowMerged_, true);
      fout<<"\\hline"<<endl;
    	fout<<"\\end{tabular}"<<endl;
    	fout<<"\\caption{Expected numbers of events and errors for $"<<Lumi_<<"~pb^{-1}$ of integrated luminosity, Merged}"<<endl;
    	fout<<"\\label{tab:SelectionTableMerged}"<<endl;
    	fout<<"\\end{table}"<<endl;
    	////////////////////////////
      fout<<"\\begin{table}"<<endl;
    	fout<<"\\centering"<<endl;
    	fout<<"\\begin{tabular}{|c|";
	  	for(unsigned int i=0;i<listOfDatasetsMerged_.size();i++) fout<<"c";
	  	fout<<"|}"<<endl;
    	fout<<"\\hline"<<endl;
	  	WriteTable(fout, cutEfficiencyMerged_, cutEfficiencyErrorMerged_, true);
    	fout<<"\\hline"<<endl;
    	fout<<"\\end{tabular}"<<endl;
    	fout<<"\\caption{Cut efficiencies, Merged}"<<endl;
    	fout<<"\\label{tab:EfficienciesTableMerged}"<<endl;
    	fout<<"\\end{table}"<<endl;
    	////////////////////////////
      fout<<"\\begin{table}"<<endl;
    	fout<<"\\centering"<<endl;
    	fout<<"\\begin{tabular}{|c|";
	  	for(unsigned int i=0;i<listOfDatasetsMerged_.size();i++) fout<<"c";
	  	fout<<"|}"<<endl;
    	fout<<"\\hline"<<endl;
	  	WriteTable(fout, cutEfficiencyMerged_, cutEfficiencyErrorHighMerged_, cutEfficiencyErrorLowMerged_, true);
    	fout<<"\\hline"<<endl;
    	fout<<"\\end{tabular}"<<endl;
    	fout<<"\\caption{Cut efficiencies, Merged}"<<endl;
    	fout<<"\\label{tab:EfficienciesTableMerged}"<<endl;
    	fout<<"\\end{table}"<<endl;
    	////////////////////////////
    	fout<<"\\begin{table}"<<endl;
    	fout<<"\\centering"<<endl;
    	fout<<"\\begin{tabular}{|c|";
	  	for(unsigned int i=0;i<listOfDatasetsMerged_.size();i++) fout<<"c";
	  	fout<<"|}"<<endl;
    	fout<<"\\hline"<<endl;
	  	WriteTable(fout, totalCutEfficiencyMerged_, totalCutEfficiencyErrorMerged_, true);
    	fout<<"\\hline"<<endl;
    	fout<<"\\end{tabular}"<<endl;
    	fout<<"\\caption{Total cut efficiencies, Merged}"<<endl;
    	fout<<"\\label{tab:TotalEfficienciesTableMerged}"<<endl;
    	fout<<"\\end{table}"<<endl;
    }
	}
	else{
		fout<<"nof Events"<<endl;
		WriteTable(fout, nofEvents_);
		fout<<"Cut Effiency"<<endl;
		WriteTable(fout, cutEfficiency_);
		fout<<"Total Cut Effiency"<<endl;
		WriteTable(fout, totalCutEfficiency_);
	}
}

float SelectionTable::WilsonScoreIntervalHigh(float Non, float Ntot)
{
	float T = (Ntot>0 ? 1/Ntot : 0);
	float p_hat = (Ntot>0 && Non>=0 && Ntot>=Non ? Non/Ntot : 0);
	float Int_High = ((p_hat+(T/2))/(1+T))+(sqrt(p_hat*(1-p_hat)*T+pow(T/2,2))/(1+T));
	return Int_High;
}

float SelectionTable::WilsonScoreIntervalLow(float Non, float Ntot)
{
	float T = (Ntot>0 ? 1/Ntot : 0);
	float p_hat = (Ntot>0 && Non>=0 && Ntot>=Non ? Non/Ntot : 0);
	float Int_Low = ((p_hat+(T/2))/(1+T))-(sqrt(p_hat*(1-p_hat)*T+pow(T/2,2))/(1+T));
	return Int_Low;
}

    

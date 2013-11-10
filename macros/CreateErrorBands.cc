#include <iostream>
#include <stdio.h>
#include <cmath>

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TList.h"

#include "TSystem.h"
#include "TKey.h"
#include "TString.h"

using namespace std;

int main(unsigned int argc, char *argv[])
{
  cout<<"This macro will calculate the systematic error bands."<<endl;
	cout<<"Note that it makes use of the fixed structure of MultiSamplePlots + some hardcoded strings to recognize systematics and physics-process names!"<<endl;

  string basedir = "/user/gvonsem/VectorLikeQuarkSearch/Git/CMSSW_5_3_6_patch1/src/TopBrussels/TopTreeAnalysis/macros/OutputFiles_VLQAnalyzer_withoutQCDestimation_MSPlotTEST2/";
	string InputfilenameNominal = basedir+"VLQTreeAnalyzer_Nominal.root";
	string InputfilenameJESPlus = basedir+"VLQTreeAnalyzer_JESPlus.root";
	string InputfilenameJESMinus = basedir+"VLQTreeAnalyzer_JESMinus.root";
	string InputfilenameJERPlus = basedir+"VLQTreeAnalyzer_JERPlus.root";
	string InputfilenameJERMinus = basedir+"VLQTreeAnalyzer_JERMinus.root";
	string InputfilenamebTagPlus = basedir+"VLQTreeAnalyzer_bTagPlus.root";
	string InputfilenamebTagMinus = basedir+"VLQTreeAnalyzer_bTagMinus.root";
	string InputfilenamemisTagPlus = basedir+"VLQTreeAnalyzer_misTagPlus.root";
	string InputfilenamemisTagMinus = basedir+"VLQTreeAnalyzer_misTagMinus.root";	
	string InputfilenamePUPlus = basedir+"VLQTreeAnalyzer_PUPlus.root";
	string InputfilenamePUMinus = basedir+"VLQTreeAnalyzer_PUMinus.root";
	
	string Outputpath = "VLQSearch_ErrorBand/";
	string Outputfilename = "ErrorBandFile.root";
	mkdir(Outputpath.c_str(),0777);
	//bool mergeSignal = false;
	
	//systematics you want to consider in the errorband. The effects will be quadratically summed bin by bin (i.e. assumption = uncorrelated systematics, which is in most cases rather unrealistic).
  vector<string> syst;
  //syst.push_back("JES");	
  //syst.push_back("JER");
	//syst.push_back("bTag");
	//syst.push_back("misTag");
  //syst.push_back("PU"); // pile-up
	syst.push_back("top_XS");
  //syst.push_back("w_XS");
  //syst.push_back("z_XS");
  syst.push_back("st_XS");
	//syst.push_back("ww_XS");
  //syst.push_back("wz_XS");
	//syst.push_back("zz_XS");
	//syst.push_back("ttw_XS");
	//syst.push_back("ttz_XS");
	//syst.push_back("ssww_XS");
	//syst.push_back("vvv_XS");
	//syst.push_back("lep_eff");
	syst.push_back("lumi");
	
	
	//numbers are to be revisited
	float rel_unc_top = 0.08;
	float rel_unc_w = 0.30;
	float rel_unc_z = 0.05;
	float rel_unc_ww = 0.15;
	float rel_unc_wz = 0.17;
	float rel_unc_zz = 0.42;
  float rel_unc_ttw = 0.23;
	float rel_unc_ttz = 0.39;
	float rel_unc_ssww = 0.95;
	float rel_unc_st = 0.30;
	float rel_unc_vvv = 0.30;
	//now some overall scale systematics
	float rel_unc_mu_eff 	= 0.03; //WARNING: how to deal with this; and with multiple leptons...
	float rel_unc_el_eff 	= 0.05;
	float rel_unc_lumi = 0.026; //2012 rereco

  TFile* outFile = new TFile((Outputpath+Outputfilename).c_str(),"RECREATE");

	TFile* inFileNominal = new TFile(InputfilenameNominal.c_str(),"READ");
	TIter nextkey(inFileNominal->GetListOfKeys());
	TKey *key;
  //loop over de directories (~the multisampleplots) in the file
  while (key = (TKey*)nextkey()) 
	{			
	  TH1F* sumNominal = 0;
	  TH1F* errorMinus = 0; //lower values of error band
    TH1F* errorPlus = 0; //upper values of error band
				
	  TDirectoryFile* subdir = (TDirectoryFile*) key->ReadObj();
		string subdirname = subdir->GetName();
		//cout<<subdirname<<endl;
		
		if(subdirname.find("MultiSamplePlot")==0)
		{
		  TIter nextplotkey(subdir->GetListOfKeys());
	    TKey *plotkey;
			//loop over histograms (in multisampleplot directory) and obtain total nominal histogram
			bool firstdataset = true;
      while (plotkey = (TKey*)nextplotkey()) 
	    {
				    TH1F* histo = 0;
						if(((string) (plotkey->ReadObj())->ClassName()) == "TH1F")
						{
						   histo = (TH1F*) plotkey->ReadObj();
							 if(!( ((string) histo->GetName()).find("Data")!=string::npos || ((string) histo->GetName()).find("data")!=string::npos || ((string) histo->GetName()).find("DATA")!=string::npos || ((string) histo->GetName()).find("_NP_")!=string::npos ) )
							 {
							    //WARNING: slightly dangerous to search for a string _NP_ in the plot 'name', since it is in principle possible to write this string in the title of a plot without realizing this is not allowed..
							    //cout<<"      histo->GetName() = "<<histo->GetName()<<endl;
									if(firstdataset)
						      {
						        sumNominal = histo;
							      firstdataset = false;
						      }
				          else sumNominal->Add( histo );
							 }
						}
			}
			
			//maybe not needed, but to be sure; add the overflow bin to the last bin, because MSPlot does this as well
		  sumNominal->SetBinContent(sumNominal->GetNbinsX(),sumNominal->GetBinContent(sumNominal->GetNbinsX()) + sumNominal->GetBinContent(sumNominal->GetNbinsX()+1));
		  sumNominal->SetBinContent(sumNominal->GetNbinsX()+1, 0);
					
			//loop over systematics
			for(unsigned int iSyst=0; iSyst<syst.size(); iSyst++)
      {			
			   //cout<<" -> systematic "<<syst[iSyst]<<endl;				 
				 //'up' variations	 
		     TFile* inFilePlus = 0;
				 if(syst[iSyst] == "JES") inFilePlus = new TFile(InputfilenameJESPlus.c_str(),"READ");
		     else if(syst[iSyst] == "JER") inFilePlus = new TFile(InputfilenameJERPlus.c_str(),"READ");
         else if(syst[iSyst] == "bTag") inFilePlus = new TFile(InputfilenamebTagPlus.c_str(),"READ");
		     else if(syst[iSyst] == "misTag") inFilePlus = new TFile(InputfilenamemisTagPlus.c_str(),"READ");
		     else if(syst[iSyst] == "PU") inFilePlus = new TFile(InputfilenamePUPlus.c_str(),"READ");
		     else inFilePlus = new TFile(InputfilenameNominal.c_str(),"READ");
		     TDirectory* subdirPlus = (TDirectory*) inFilePlus->Get(subdirname.c_str());
         TH1F* sumPlus = 0;
				 
				 TIter nextplotkeyPlus(subdirPlus->GetListOfKeys());
	       TKey *plotkeyPlus;
				 bool firstdatasetPlus = true;
				 				 				 
				 if(syst[iSyst] == "JES" || syst[iSyst] == "JER" || syst[iSyst] == "bTag" || syst[iSyst] == "misTag" || syst[iSyst] == "PU") //WARNING: this is some necessary hardcoding, unfortunately...
         {
						//loop over histograms (in multisampleplot directory)
            while (plotkeyPlus = (TKey*)nextplotkeyPlus()) 
	          {
				       TH1F* histoPlus = 0;
						   if(((string) (plotkeyPlus->ReadObj())->ClassName()) == "TH1F")
						   {
						      histoPlus = (TH1F*) plotkeyPlus->ReadObj();
							    if(!( ((string) histoPlus->GetName()).find("Data")!=string::npos || ((string) histoPlus->GetName()).find("data")!=string::npos || ((string) histoPlus->GetName()).find("DATA")!=string::npos || ((string) histoPlus->GetName()).find("_NP_")!=string::npos ) )
							    {
							       //cout<<"      histoPlus->GetName() = "<<histoPlus->GetName()<<", integral = "<<histoPlus->GetIntegral()<<endl;										 
								     if(firstdatasetPlus)
						         {
						           sumPlus = histoPlus;
							         firstdataset = false;
						         }
				             else sumPlus->Add( histoPlus );
							    }
						   }
				    }	 //end loop over histograms	
				 }
				 else
				 {
						//loop over histograms (in multisampleplot directories)
            while (plotkeyPlus = (TKey*)nextplotkeyPlus()) 
	          {
				       TH1F* histoPlus = 0;
						   if(((string) (plotkeyPlus->ReadObj())->ClassName()) == "TH1F")
						   {
						      histoPlus = (TH1F*) plotkeyPlus->ReadObj();
							    if(!( ((string) histoPlus->GetName()).find("Data")!=string::npos || ((string) histoPlus->GetName()).find("data")!=string::npos || ((string) histoPlus->GetName()).find("DATA")!=string::npos || ((string) histoPlus->GetName()).find("_NP_")!=string::npos ) )
							    {
							       //cout<<"      histoPlus->GetName() = "<<histoPlus->GetName()<<endl;
								     if(syst[iSyst]=="top_XS" && ((string) histoPlus->GetName()).find("TTbarJets")!=string::npos) //WARNING: this is some necessary hardcoding, unfortunately...
						         {
										    histoPlus->Scale(1+rel_unc_top);
						         }
										 else if(syst[iSyst] == "w_XS" && ((string) histoPlus->GetName()).find("WJets")!=string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_w);
										 }
										 else if(syst[iSyst] == "z_XS" && ((string) histoPlus->GetName()).find("ZJets")!=string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_z);
										 }
										 else if(syst[iSyst] == "st_XS" && ((string) histoPlus->GetName()).find("ST_")!=string::npos)
										 {
										    histoPlus->Scale(1+rel_unc_st);
										 }
										 else if(syst[iSyst] == "ww_XS" && ((string) histoPlus->GetName()).find("WW")!=string::npos && ((string) histoPlus->GetName()).find("WWW")==string::npos && ((string) histoPlus->GetName()).find("WWZ")==string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_ww);
										 }
										 else if(syst[iSyst] == "wz_XS" && ((string) histoPlus->GetName()).find("WZ")!=string::npos && ((string) histoPlus->GetName()).find("WWZ")==string::npos && ((string) histoPlus->GetName()).find("WZZ")==string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_wz);
										 }
										 else if(syst[iSyst] == "zz_XS" && ((string) histoPlus->GetName()).find("ZZ")!=string::npos && ((string) histoPlus->GetName()).find("WZZ")==string::npos && ((string) histoPlus->GetName()).find("ZZZ")==string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_zz);
										 }
										 else if(syst[iSyst] == "ttw_XS" && ((string) histoPlus->GetName()).find("ttW")!=string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_ttw);
										 }
										 else if(syst[iSyst] == "ttz_XS" && ((string) histoPlus->GetName()).find("ttZ")!=string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_ttz);
										 }
										 else if(syst[iSyst] == "ssww_XS" && ( ((string) histoPlus->GetName()).find("WpWp")!=string::npos || ((string) histoPlus->GetName()).find("WmWm")!=string::npos ) )
			               {
										    histoPlus->Scale(1+rel_unc_ssww);
										 }
										 else if(syst[iSyst] == "vvv_XS" && ( ((string) histoPlus->GetName()).find("WWW")!=string::npos || ((string) histoPlus->GetName()).find("WWZ")!=string::npos || ((string) histoPlus->GetName()).find("WZZ")!=string::npos || ((string) histoPlus->GetName()).find("ZZZ")!=string::npos ) )
			               {
										    histoPlus->Scale(1+rel_unc_vvv);
										 }
										 else if(syst[iSyst] == "lumi")
			               {
										    histoPlus->Scale(1+rel_unc_lumi);
										 }
										 
										 if(firstdatasetPlus)
						         {
						           sumPlus = histoPlus;
							         firstdatasetPlus = false;
						         }
				             else sumPlus->Add( histoPlus );										 
							    }
						   }
				    } //end loop over histograms				 
				 }
				 sumPlus->SetDirectory(0);
         inFilePlus->Close();
         delete inFilePlus;
				 
		     //'minus' variations	 
		     TFile* inFileMinus = 0;
				 if(syst[iSyst] == "JES") inFileMinus = new TFile(InputfilenameJESMinus.c_str(),"READ");
		     else if(syst[iSyst] == "JER") inFileMinus = new TFile(InputfilenameJERMinus.c_str(),"READ");
         else if(syst[iSyst] == "bTag") inFileMinus = new TFile(InputfilenamebTagMinus.c_str(),"READ");
		     else if(syst[iSyst] == "misTag") inFileMinus = new TFile(InputfilenamemisTagMinus.c_str(),"READ");
		     else if(syst[iSyst] == "PU") inFileMinus = new TFile(InputfilenamePUMinus.c_str(),"READ");
		     else inFileMinus = new TFile(InputfilenameNominal.c_str(),"READ");
		     TDirectory* subdirMinus = (TDirectory*) inFileMinus->Get(subdirname.c_str());
         TH1F* sumMinus = 0;
				 
				 TIter nextplotkeyMinus(subdirMinus->GetListOfKeys());
	       TKey *plotkeyMinus;
				 bool firstdatasetMinus = true;
				 
				 if(syst[iSyst] == "JES" || syst[iSyst] == "JER" || syst[iSyst] == "bTag" || syst[iSyst] == "misTag" || syst[iSyst] == "PU")
         {
						//loop over histograms (in multisampleplot directories)
            while (plotkeyMinus = (TKey*)nextplotkeyMinus()) 
	          {
				       TH1F* histoMinus = 0;
						   if(((string) (plotkeyMinus->ReadObj())->ClassName()) == "TH1F")
						   {
						      histoMinus = (TH1F*) plotkeyMinus->ReadObj();
							    if(!( ((string) histoMinus->GetName()).find("Data")!=string::npos || ((string) histoMinus->GetName()).find("data")!=string::npos || ((string) histoMinus->GetName()).find("DATA")!=string::npos || ((string) histoMinus->GetName()).find("_NP_")!=string::npos ) )
							    {
							       //cout<<"      histoMinus->GetName() = "<<histoMinus->GetName()<<endl;
								     if(firstdatasetMinus)
						         {
						           sumMinus = histoMinus;
							         firstdataset = false;
						         }
				             else sumMinus->Add( histoMinus );
							    }
						   }
				    }	
				}
				else
				{
						//loop over histograms (in multisampleplot directories)
            while (plotkeyMinus = (TKey*)nextplotkeyMinus()) 
	          {
				       TH1F* histoMinus = 0;
						   if(((string) (plotkeyMinus->ReadObj())->ClassName()) == "TH1F")
						   {
						      histoMinus = (TH1F*) plotkeyMinus->ReadObj();
							    if(!( ((string) histoMinus->GetName()).find("Data")!=string::npos || ((string) histoMinus->GetName()).find("data")!=string::npos || ((string) histoMinus->GetName()).find("DATA")!=string::npos || ((string) histoMinus->GetName()).find("_NP_")!=string::npos ) )
							    {
							       //cout<<"      histoMinus->GetName() = "<<histoMinus->GetName()<<endl;
								     if(syst[iSyst]=="top_XS" && ((string) histoMinus->GetName()).find("TTbarJets")!=string::npos) //WARNING: this is some necessary hardcoding, unfortunately...
						         {
										    histoMinus->Scale(1-rel_unc_top);
						         }
										 else if(syst[iSyst] == "w_XS" && ((string) histoMinus->GetName()).find("WJets")!=string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_w);
										 }
										 else if(syst[iSyst] == "z_XS" && ((string) histoMinus->GetName()).find("ZJets")!=string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_z);
										 }
										 else if(syst[iSyst] == "st_XS" && ((string) histoMinus->GetName()).find("ST_")!=string::npos)
										 {
										    histoMinus->Scale(1-rel_unc_st);
										 }
										 else if(syst[iSyst] == "ww_XS" && ((string) histoMinus->GetName()).find("WW")!=string::npos && ((string) histoMinus->GetName()).find("WWW")==string::npos && ((string) histoMinus->GetName()).find("WWZ")==string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_ww);
										 }
										 else if(syst[iSyst] == "wz_XS" && ((string) histoMinus->GetName()).find("WZ")!=string::npos && ((string) histoMinus->GetName()).find("WWZ")==string::npos && ((string) histoMinus->GetName()).find("WZZ")==string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_wz);
										 }
										 else if(syst[iSyst] == "zz_XS" && ((string) histoMinus->GetName()).find("ZZ")!=string::npos && ((string) histoMinus->GetName()).find("WZZ")==string::npos && ((string) histoMinus->GetName()).find("ZZZ")==string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_zz);
										 }
										 else if(syst[iSyst] == "ttw_XS" && ((string) histoMinus->GetName()).find("ttW")!=string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_ttw);
										 }
										 else if(syst[iSyst] == "ttz_XS" && ((string) histoMinus->GetName()).find("ttZ")!=string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_ttz);
										 }
										 else if(syst[iSyst] == "ssww_XS" && ( ((string) histoMinus->GetName()).find("WpWp")!=string::npos || ((string) histoMinus->GetName()).find("WmWm")!=string::npos ) )
			               {
										    histoMinus->Scale(1-rel_unc_ssww);
										 }
										 else if(syst[iSyst] == "vvv_XS" && ( ((string) histoMinus->GetName()).find("WWW")!=string::npos || ((string) histoMinus->GetName()).find("WWZ")!=string::npos || ((string) histoMinus->GetName()).find("WZZ")!=string::npos || ((string) histoMinus->GetName()).find("ZZZ")!=string::npos ) )
			               {
										    histoMinus->Scale(1-rel_unc_vvv);
										 }
										 else if(syst[iSyst] == "lumi")
			               {
										    histoMinus->Scale(1-rel_unc_lumi);
										 }
										 
										 
										 if(firstdatasetMinus)
						         {
						           sumMinus = histoMinus;
							         firstdatasetMinus = false;
						         }
				             else sumMinus->Add( histoMinus );
										 
							    }
						   }
				    }				 //end loop over histograms		
				}	 
				sumMinus->SetDirectory(0);
        inFileMinus->Close();
        delete inFileMinus;
		
		
		    //cout << subdirname << "    Nominal: " << sumNominal->Integral(1,sumNominal->GetNbinsX()) << "  Plus: " << sumPlus->Integral(1,sumPlus->GetNbinsX()) << "  Minus: " << sumMinus->Integral(1,sumMinus->GetNbinsX()) << endl;
		 		
		    //maybe not needed, but to be sure; add the overflow bin to the last bin, because MSPlot does this as well
		    sumMinus->SetBinContent(sumMinus->GetNbinsX(),sumMinus->GetBinContent(sumMinus->GetNbinsX()) + sumMinus->GetBinContent(sumMinus->GetNbinsX()+1));
		    sumMinus->SetBinContent(sumMinus->GetNbinsX()+1, 0);
		    sumPlus->SetBinContent(sumPlus->GetNbinsX(),sumPlus->GetBinContent(sumPlus->GetNbinsX()) + sumPlus->GetBinContent(sumPlus->GetNbinsX()+1));
		    sumPlus->SetBinContent(sumPlus->GetNbinsX()+1, 0);
				//cout<<"         sumPlus->GetName() = "<<sumPlus->GetName()<<", integral = "<<sumPlus->Integral(1,sumPlus->GetNbinsX())<<endl; 
				//cout<<"         sumNominal->GetName() = "<<sumNominal->GetName()<<", integral = "<<sumNominal->Integral(1,sumPlus->GetNbinsX())<<endl;	
				//cout<<"         sumMinus->GetName() = "<<sumMinus->GetName()<<", integral = "<<sumMinus->Integral(1,sumPlus->GetNbinsX())<<endl;	
		 
		 
		    if(iSyst == 0) //only the first time in the loop over systematics...
        {
           errorMinus = (TH1F*) sumMinus->Clone();
           errorPlus = (TH1F*) sumPlus->Clone();
        }
        for(unsigned int iBin=0; iBin<=errorMinus->GetNbinsX()+1; iBin++)
        { 
          if(sumMinus->GetBinContent(iBin) < sumNominal->GetBinContent(iBin))
					{
            if(iSyst != 0) errorMinus->SetBinContent(iBin, sumNominal->GetBinContent(iBin) - sqrt( pow(errorMinus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) + pow(sumMinus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) ) );
            //cout<<"errorMinus->GetBinContent("<<iBin<<") = "<<errorMinus->GetBinContent(iBin)<<endl;
					}
					else
					{
            if(iSyst != 0) errorPlus->SetBinContent(iBin, sumNominal->GetBinContent(iBin) + sqrt( pow(errorPlus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) + pow(sumMinus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) ) );
            //cout<<"errorPlus->GetBinContent("<<iBin<<") = "<<errorPlus->GetBinContent(iBin)<<endl;
				  }
			
          if(sumPlus->GetBinContent(iBin) < sumNominal->GetBinContent(iBin))
					{
            if(iSyst != 0) errorMinus->SetBinContent(iBin, sumNominal->GetBinContent(iBin) - sqrt( pow(errorMinus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) + pow(sumPlus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) ) );
            //cout<<"errorMinus->GetBinContent("<<iBin<<") = "<<errorMinus->GetBinContent(iBin)<<endl;
					}
					else
					{
            if(iSyst != 0) errorPlus->SetBinContent(iBin, sumNominal->GetBinContent(iBin) + sqrt( pow(errorPlus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) + pow(sumPlus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) ) );
            //cout<<"errorPlus->GetBinContent("<<iBin<<") = "<<errorPlus->GetBinContent(iBin)<<endl;
					}
				}
		    //cout << subdirname << "    errorPlus: " << errorPlus->Integral(1,errorPlus->GetNbinsX()) << "  errorMinus: " << errorMinus->Integral(1,errorMinus->GetNbinsX()) << endl;
		 
		  } //end loop over systematics
			//cout << "Integrals for plot " << subdirname << " ->   errorPlus: " << errorPlus->Integral(1,errorPlus->GetNbinsX()) << ", nominal: " << sumNominal->Integral(1,sumNominal->GetNbinsX()) << ",  errorMinus: " << errorMinus->Integral(1,errorMinus->GetNbinsX()) << endl;

      outFile->cd();
		  if(outFile->Get(subdirname.c_str())==0)
		    outFile->mkdir(subdirname.c_str());
		  outFile->cd(subdirname.c_str());
		  errorMinus->SetNameTitle("Minus","Minus");
			errorMinus->SetLineColor(1);
      errorMinus->Write();
      errorPlus->SetNameTitle("Plus","Plus");
			errorPlus->SetLineColor(1);
      errorPlus->Write();
	    sumNominal->SetNameTitle("Nominal","Nominal");
			sumNominal->SetLineColor(1);
	    sumNominal->Write();

	  }					
  }
	outFile->Close();
  delete outFile;
	inFileNominal->Close();
	delete inFileNominal;
	
	cout<<"Error band histograms written in output file: "<<Outputpath+Outputfilename<<endl;

}

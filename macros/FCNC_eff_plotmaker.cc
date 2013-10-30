// isis.marina.van.parijs@cern.ch 
// 2013
// This is a program that calculates efficiencies and makes plots 
//  



#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "TCanvas.h"
//#include "setTDRStyle.C"





using namespace std;	//needed for cout and stuff

void FCNC_eff_plotmaker(string channel = "3L")
{

	char rootFileName[20];
	char channelchar[10];
	char plotTitle_eff_total[30];
	char plotName_eff_total[20];
	
	
	if(channel.find("3gamma")!=string::npos)	sprintf(channelchar, "3gamma"); 
	if(channel.find("3L")!=string::npos)		sprintf(channelchar, "3L");	
	if(channel.find("1L3B")!=string::npos)		sprintf(channelchar, "1L3B");	
	if(channel.find("SSdilepton")!=string::npos)	sprintf(channelchar, "SSdilepton"); 
	
	sprintf(plotTitle_eff_total,"The total efficiency for %s channel",channelchar); 
	
	sprintf(plotName_eff_total,"eff_total_%s",channelchar);
	sprintf(rootFileName,"Output_FCNC_selection_%s.root",channelchar);
	
	TFile *fin = TFile::Open(rootFileName);
	cout << "[PROCES]	Declared input rootfile  "<< rootFileName << endl;
	
	TH1F* eff_total = new TH1F(plotName_eff_total, plotTitle_eff_total, 6, -0.5,5.5); 
	eff_total->Sumw2();
	
	TH1F* cutflowtotal; 
	
	cutflowtotal = (TH1F*) fin->Get("cutflow_total");
	
	
	
	
	for(int i=2;i<cutflowtotal->GetSize()-2; i++)
	{
		
		
	 	double denum = cutflowtotal->GetBinContent(2);
		cout << "denum: " << denum << endl;
		double num = cutflowtotal->GetBinContent(i);
		cout << "num: " << num << endl; 
		double eff = num/denum; 
		cout << "Efficiency for bin " << i << " is: " << eff << endl; 
		eff_total->SetBinContent(i,eff);
	
	
	}	
	
	TCanvas *c1 = new TCanvas();
       
	eff_total->Draw("h");
	eff_total->SetMaximum(1.1 );
	eff_total->GetXaxis()->SetBinLabel(2,"initial");
	eff_total->GetYaxis()->SetTitle("eff");
	
		if(channel.find("3L")!=string::npos)		
	{	
		eff_total->GetXaxis()->SetBinLabel(3, "3L");
		c1->SaveAs("FCNC_plot_eff_3L.png");
	}
	if(channel.find("1L3B")!=string::npos)
	{		
		eff_total->GetXaxis()->SetBinLabel(3, "1L");
		eff_total->GetXaxis()->SetBinLabel(4, "3J");
		eff_total->GetXaxis()->SetBinLabel(5, "3B");
		c1->SaveAs("FCNC_plot_eff_1L3B.png");	
	}
	if(channel.find("SSdilepton")!=string::npos)
	{	
		eff_total->GetXaxis()->SetBinLabel(3, "2L");
		eff_total->GetXaxis()->SetBinLabel(4, "SS");
		c1->SaveAs("FCNC_plot_eff_SSdilepton.png");
	}
	
	/////// 3L
	////////////////////
	char plotTitle_eff_3L_ttw[30];
	sprintf(plotTitle_eff_3L_ttw,"The efficiency for %s channel: ttw",channelchar); 
	TH1F* eff_3L_ttw = new TH1F("eff_3L_ttw", plotTitle_eff_total, 6, -0.5,5.5); 
	eff_3L_ttw->Sumw2();
	
	TH1F* cutflow_3L_ttw;
	
	char plotTitle_eff_3L_ttbar[30];
	sprintf(plotTitle_eff_3L_ttbar,"The efficiency for %s channel: ttbar",channelchar); 
	TH1F* eff_3L_ttbar = new TH1F("eff_3L_ttbar", plotTitle_eff_total, 6, -0.5,5.5); 
	eff_3L_ttbar->Sumw2();
	
	TH1F* cutflow_3L_ttbar;
	
	char plotTitle_eff_3L_wz[30];
	sprintf(plotTitle_eff_3L_wz,"The efficiency for %s channel: wz",channelchar); 
	TH1F* eff_3L_wz = new TH1F("eff_3L_wz", plotTitle_eff_total, 6, -0.5,5.5); 
	eff_3L_wz->Sumw2();
	
	TH1F* cutflow_3L_wz;
	
	char plotTitle_eff_3L_zz[30];
	sprintf(plotTitle_eff_3L_zz,"The efficiency for %s channel: zz",channelchar); 
	TH1F* eff_3L_zz = new TH1F("eff_3L_zz", plotTitle_eff_total, 6, -0.5,5.5); 
	eff_3L_zz->Sumw2();
	
	TH1F* cutflow_3L_zz;
	
	char plotTitle_eff_3L_ttz[30];
	sprintf(plotTitle_eff_3L_ttz,"The efficiency for %s channel: ttz",channelchar); 
	TH1F* eff_3L_ttz = new TH1F("eff_3L_ttz", plotTitle_eff_total, 6, -0.5,5.5); 
	eff_3L_ttz->Sumw2();
	
	TH1F* cutflow_3L_ttz;
	
	TCanvas *c3L_ttw = new TCanvas();
	TCanvas *c3L_ttbar = new TCanvas();
	TCanvas *c3L_wz = new TCanvas();
	TCanvas *c3L_zz = new TCanvas();
	TCanvas *c3L_ttz = new TCanvas();
	
	if(channel.find("3L")!=string::npos)
	{
	
	     
	     	cutflow_3L_ttw= (TH1F*) fin->Get("cutflow_ttw");
	
		for(int i=2;i<cutflow_3L_ttw->GetSize()-2; i++)
		{
			double denum = cutflow_3L_ttw->GetBinContent(2);
			//cout << "denum: " << denum << endl;
			double num = cutflow_3L_ttw->GetBinContent(i);
			//cout << "num: " << num << endl; 
			double eff = num/denum; 
			//cout << "Efficiency for bin " << i << " is: " << eff << endl; 
			eff_3L_ttw->SetBinContent(i,eff);

		}
		eff_3L_ttw->GetXaxis()->SetBinLabel(2,"initial");
		eff_3L_ttw->GetYaxis()->SetTitle("eff");
		eff_3L_ttw->GetXaxis()->SetBinLabel(3, "3L");
		c3L_ttw->cd(); 
		eff_3L_ttw->Draw("h");
		c3L_ttw->SaveAs("FCNC_plot_eff_3L_ttw.png");
	
		cutflow_3L_ttbar= (TH1F*) fin->Get("cutflow_ttbar");
	
		for(int i=2;i<cutflow_3L_ttbar->GetSize()-2; i++)
		{
			double denum = cutflow_3L_ttbar->GetBinContent(2);
			//cout << "denum: " << denum << endl;
			double num = cutflow_3L_ttbar->GetBinContent(i);
			//cout << "num: " << num << endl; 
			double eff = num/denum; 
			//cout << "Efficiency for bin " << i << " is: " << eff << endl; 
			eff_3L_ttbar->SetBinContent(i,eff);

		}
		eff_3L_ttbar->GetXaxis()->SetBinLabel(2,"initial");
		eff_3L_ttbar->GetYaxis()->SetTitle("eff");
		eff_3L_ttbar->GetXaxis()->SetBinLabel(3, "3L");
		c3L_ttbar->cd(); 
		eff_3L_ttbar->Draw("h");
		c3L_ttbar->SaveAs("FCNC_plot_eff_3L_ttbar.png");
		
		cutflow_3L_wz= (TH1F*) fin->Get("cutflow_wz");
	
		for(int i=2;i<cutflow_3L_wz->GetSize()-2; i++)
		{
			double denum = cutflow_3L_wz->GetBinContent(2);
			//cout << "denum: " << denum << endl;
			double num = cutflow_3L_wz->GetBinContent(i);
			//cout << "num: " << num << endl; 
			double eff = num/denum; 
			//cout << "Efficiency for bin " << i << " is: " << eff << endl; 
			eff_3L_wz->SetBinContent(i,eff);

		}
		eff_3L_wz->GetXaxis()->SetBinLabel(2,"initial");
		eff_3L_wz->GetYaxis()->SetTitle("eff");
		eff_3L_wz->GetXaxis()->SetBinLabel(3, "3L");
		c3L_wz->cd(); 
		eff_3L_wz->Draw("h");
		c3L_wz->SaveAs("FCNC_plot_eff_3L_wz.png");
		
		cutflow_3L_zz= (TH1F*) fin->Get("cutflow_zz");
	
		for(int i=2;i<cutflow_3L_zz->GetSize()-2; i++)
		{
			double denum = cutflow_3L_zz->GetBinContent(2);
			//cout << "denum: " << denum << endl;
			double num = cutflow_3L_zz->GetBinContent(i);
			//cout << "num: " << num << endl; 
			double eff = num/denum; 
			//cout << "Efficiency for bin " << i << " is: " << eff << endl; 
			eff_3L_zz->SetBinContent(i,eff);

		}
		eff_3L_zz->GetXaxis()->SetBinLabel(2,"initial");
		eff_3L_zz->GetYaxis()->SetTitle("eff");
		eff_3L_zz->GetXaxis()->SetBinLabel(3, "3L");
		c3L_zz->cd(); 
		eff_3L_zz->Draw("h");
		c3L_zz->SaveAs("FCNC_plot_eff_3L_zz.png");
		
		cutflow_3L_ttz= (TH1F*) fin->Get("cutflow_ttz");
	
		for(int i=2;i<cutflow_3L_ttz->GetSize()-2; i++)
		{
			double denum = cutflow_3L_ttz->GetBinContent(2);
			//cout << "denum: " << denum << endl;
			double num = cutflow_3L_ttz->GetBinContent(i);
			//cout << "num: " << num << endl; 
			double eff = num/denum; 
			//cout << "Efficiency for bin " << i << " is: " << eff << endl; 
			eff_3L_ttz->SetBinContent(i,eff);

		}
		eff_3L_ttz->GetXaxis()->SetBinLabel(2,"initial");
		eff_3L_ttz->GetYaxis()->SetTitle("eff");
		eff_3L_ttz->GetXaxis()->SetBinLabel(3, "3L");
		c3L_ttz->cd(); 
		eff_3L_ttz->Draw("h");
		c3L_ttz->SaveAs("FCNC_plot_eff_3L_ttz.png");
	}
	
	
	
	//////////////////// SSdilepton
	//////////////////////////////////////:
	char plotTitle_eff_SSdilepton_ttw[30];
	sprintf(plotTitle_eff_SSdilepton_ttw,"The efficiency for %s channel: ttw",channelchar); 
	TH1F* eff_SSdilepton_ttw = new TH1F("eff_SSdilepton_ttw", plotTitle_eff_total, 6, -0.5,5.5); 
	eff_SSdilepton_ttw->Sumw2();
	
	TH1F* cutflow_SSdilepton_ttw;
	
	char plotTitle_eff_SSdilepton_ttt[30];
	sprintf(plotTitle_eff_SSdilepton_ttt,"The efficiency for %s channel: ttt",channelchar); 
	TH1F* eff_SSdilepton_ttt = new TH1F("eff_SSdilepton_ttt", plotTitle_eff_total, 6, -0.5,5.5); 
	eff_SSdilepton_ttt->Sumw2();
	
	TH1F* cutflow_SSdilepton_ttt;
	
	char plotTitle_eff_SSdilepton_wz[30];
	sprintf(plotTitle_eff_SSdilepton_wz,"The efficiency for %s channel: wz",channelchar); 
	TH1F* eff_SSdilepton_wz = new TH1F("eff_SSdilepton_wz", plotTitle_eff_total, 6, -0.5,5.5); 
	eff_SSdilepton_wz->Sumw2();
	
	TH1F* cutflow_SSdilepton_wz;
	
	
	
	char plotTitle_eff_SSdilepton_ttz[30];
	sprintf(plotTitle_eff_SSdilepton_ttz,"The efficiency for %s channel: ttz",channelchar); 
	TH1F* eff_SSdilepton_ttz = new TH1F("eff_SSdilepton_ttz", plotTitle_eff_total, 6, -0.5,5.5); 
	eff_SSdilepton_ttz->Sumw2();
	
	TH1F* cutflow_SSdilepton_ttz;
	
	TCanvas *cSSdilepton_ttw = new TCanvas();
	TCanvas *cSSdilepton_ttt = new TCanvas();
	TCanvas *cSSdilepton_wz = new TCanvas();
	
	TCanvas *cSSdilepton_ttz = new TCanvas();
	
	if(channel.find("SSdilepton")!=string::npos)
	{
	
	     
	     	cutflow_SSdilepton_ttw= (TH1F*) fin->Get("cutflow_ttw");
	
		for(int i=2;i<cutflow_SSdilepton_ttw->GetSize()-2; i++)
		{
			double denum = cutflow_SSdilepton_ttw->GetBinContent(2);
			//cout << "denum: " << denum << endl;
			double num = cutflow_SSdilepton_ttw->GetBinContent(i);
			//cout << "num: " << num << endl; 
			double eff = num/denum; 
			//cout << "Efficiency for bin " << i << " is: " << eff << endl; 
			eff_SSdilepton_ttw->SetBinContent(i,eff);

		}
		eff_SSdilepton_ttw->GetXaxis()->SetBinLabel(2,"initial");
		eff_SSdilepton_ttw->GetYaxis()->SetTitle("eff");
		eff_SSdilepton_ttw->GetXaxis()->SetBinLabel(3, "SSdilepton");
		cSSdilepton_ttw->cd(); 
		eff_SSdilepton_ttw->Draw("h");
		cSSdilepton_ttw->SaveAs("FCNC_plot_eff_SSdilepton_ttw.png");
	
		cutflow_SSdilepton_ttt= (TH1F*) fin->Get("cutflow_ttt");
	
		for(int i=2;i<cutflow_SSdilepton_ttt->GetSize()-2; i++)
		{
			double denum = cutflow_SSdilepton_ttt->GetBinContent(2);
			//cout << "denum: " << denum << endl;
			double num = cutflow_SSdilepton_ttt->GetBinContent(i);
			//cout << "num: " << num << endl; 
			double eff = num/denum; 
			//cout << "Efficiency for bin " << i << " is: " << eff << endl; 
			eff_SSdilepton_ttt->SetBinContent(i,eff);

		}
		eff_SSdilepton_ttt->GetXaxis()->SetBinLabel(2,"initial");
		eff_SSdilepton_ttt->GetYaxis()->SetTitle("eff");
		eff_SSdilepton_ttt->GetXaxis()->SetBinLabel(3, "SSdilepton");
		cSSdilepton_ttt->cd(); 
		eff_SSdilepton_ttt->Draw("h");
		cSSdilepton_ttt->SaveAs("FCNC_plot_eff_SSdilepton_ttt.png");
		
		cutflow_SSdilepton_wz= (TH1F*) fin->Get("cutflow_wz");
	
		for(int i=2;i<cutflow_SSdilepton_wz->GetSize()-2; i++)
		{
			double denum = cutflow_SSdilepton_wz->GetBinContent(2);
			//cout << "denum: " << denum << endl;
			double num = cutflow_SSdilepton_wz->GetBinContent(i);
			//cout << "num: " << num << endl; 
			double eff = num/denum; 
			//cout << "Efficiency for bin " << i << " is: " << eff << endl; 
			eff_SSdilepton_wz->SetBinContent(i,eff);

		}
		eff_SSdilepton_wz->GetXaxis()->SetBinLabel(2,"initial");
		eff_SSdilepton_wz->GetYaxis()->SetTitle("eff");
		eff_SSdilepton_wz->GetXaxis()->SetBinLabel(3, "SSdilepton");
		cSSdilepton_wz->cd(); 
		eff_SSdilepton_wz->Draw("h");
		cSSdilepton_wz->SaveAs("FCNC_plot_eff_SSdilepton_wz.png");
		
		
		
		cutflow_SSdilepton_ttz= (TH1F*) fin->Get("cutflow_ttz");
	
		for(int i=2;i<cutflow_SSdilepton_ttz->GetSize()-2; i++)
		{
			double denum = cutflow_SSdilepton_ttz->GetBinContent(2);
			//cout << "denum: " << denum << endl;
			double num = cutflow_SSdilepton_ttz->GetBinContent(i);
			//cout << "num: " << num << endl; 
			double eff = num/denum; 
			//cout << "Efficiency for bin " << i << " is: " << eff << endl; 
			eff_SSdilepton_ttz->SetBinContent(i,eff);

		}
		eff_SSdilepton_ttz->GetXaxis()->SetBinLabel(2,"initial");
		eff_SSdilepton_ttz->GetYaxis()->SetTitle("eff");
		eff_SSdilepton_ttz->GetXaxis()->SetBinLabel(3, "SSdilepton");
		cSSdilepton_ttz->cd(); 
		eff_SSdilepton_ttz->Draw("h");
		cSSdilepton_ttz->SaveAs("FCNC_plot_eff_SSdilepton_ttz.png");
	}
	
	
	
	/////// 1L3B
	////////////////////
	char plotTitle_eff_1L3B_wjets[30];
	sprintf(plotTitle_eff_1L3B_wjets,"The efficiency for %s channel: wjets",channelchar); 
	TH1F* eff_1L3B_wjets = new TH1F("eff_1L3B_wjets", plotTitle_eff_total, 6, -0.5,5.5); 
	eff_1L3B_wjets->Sumw2();
	
	TH1F* cutflow_1L3B_wjets;
	
	char plotTitle_eff_1L3B_ttbar[30];
	sprintf(plotTitle_eff_1L3B_ttbar,"The efficiency for %s channel: ttbar",channelchar); 
	TH1F* eff_1L3B_ttbar = new TH1F("eff_1L3B_ttbar", plotTitle_eff_total, 6, -0.5,5.5); 
	eff_1L3B_ttbar->Sumw2();
	
	TH1F* cutflow_1L3B_ttbar;
	
	
	
	TCanvas *c1L3B_wjets = new TCanvas();
	TCanvas *c1L3B_ttbar = new TCanvas();
	
	
	if(channel.find("1L3B")!=string::npos)
	{
	
	     
	     	cutflow_1L3B_wjets= (TH1F*) fin->Get("cutflow_wjets");
	
		for(int i=2;i<cutflow_1L3B_wjets->GetSize()-2; i++)
		{
			double denum = cutflow_1L3B_wjets->GetBinContent(2);
			//cout << "denum: " << denum << endl;
			double num = cutflow_1L3B_wjets->GetBinContent(i);
			//cout << "num: " << num << endl; 
			double eff = num/denum; 
			//cout << "Efficiency for bin " << i << " is: " << eff << endl; 
			eff_1L3B_wjets->SetBinContent(i,eff);

		}
		eff_1L3B_wjets->GetXaxis()->SetBinLabel(2,"initial");
		eff_1L3B_wjets->GetYaxis()->SetTitle("eff");
		eff_1L3B_wjets->GetXaxis()->SetBinLabel(3, "1L3B");
		c1L3B_wjets->cd(); 
		eff_1L3B_wjets->Draw("h");
		c1L3B_wjets->SaveAs("FCNC_plot_eff_1L3B_wjets.png");
	
		cutflow_1L3B_ttbar= (TH1F*) fin->Get("cutflow_ttbar");
	
		for(int i=2;i<cutflow_1L3B_ttbar->GetSize()-2; i++)
		{
			double denum = cutflow_1L3B_ttbar->GetBinContent(2);
			//cout << "denum: " << denum << endl;
			double num = cutflow_1L3B_ttbar->GetBinContent(i);
			//cout << "num: " << num << endl; 
			double eff = num/denum; 
			//cout << "Efficiency for bin " << i << " is: " << eff << endl; 
			eff_1L3B_ttbar->SetBinContent(i,eff);

		}
		eff_1L3B_ttbar->GetXaxis()->SetBinLabel(2,"initial");
		eff_1L3B_ttbar->GetYaxis()->SetTitle("eff");
		eff_1L3B_ttbar->GetXaxis()->SetBinLabel(3, "1L3B");
		c1L3B_ttbar->cd(); 
		eff_1L3B_ttbar->Draw("h");
		c1L3B_ttbar->SaveAs("FCNC_plot_eff_1L3B_ttbar.png");
		

	}
	

}

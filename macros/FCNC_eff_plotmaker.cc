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
	TH1F* cutflow;
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

}

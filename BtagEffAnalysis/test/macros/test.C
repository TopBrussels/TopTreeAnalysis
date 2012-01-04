{
   TCanvas *cst1 = new TCanvas("cst1","cst1",700,400);
   cst1->Divide(2,1);

   TH1F * hst11 = new TH1F("hst11", "", 20, -10, 10);
   hst11->Sumw2();
   hst11->FillRandom("gaus", 1000);
   hst11->SetFillColor(kViolet);
   hst11->SetLineColor(kViolet);

   TH1F * hst12 = new TH1F("hst12", "", 20, -10, 10);
   hst12->FillRandom("gaus", 500);
   hst12->SetFillColor(kBlue);
   hst12->SetLineColor(kBlue);

   THStack st1("st1", "st1");
   st1.Add(hst11);
   st1.Add(hst12);

   cst1->cd(1); st1.Draw();
   cst1->cd(2); st1.Draw("hist");

   return cst1;
}

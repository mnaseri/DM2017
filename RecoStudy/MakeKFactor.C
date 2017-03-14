#include "../interface/CMS_lumi.h"
#include <algorithm>
#include "TROOT.h"
#include <vector>
#include "TFile.h"
#include "TString.h"
#include "TTreeReader.h"
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>


vector<float> bosonPt {50, 70, 90, 110, 130, 150, 170, 200, 230, 260, 290, 320, 350, 390, 430, 470, 510, 550, 590, 640, 690, 740, 790, 840, 900, 960, 1020, 1090, 1160, 1250};

enum class Sample {zll, wjet};
float lumi_ = 1;


float GetXsec(const TString & file){
//LO
if(file.Contains("WJetsToLNu.root")) return 50690.0;
//NLO
if(file.Contains("DYJetsToLL_M-50_FXFX_ext2-v1")) return 5670.;
if(file.Contains("DYJetsToLL_Pt-100To250_FXFX-v2 1 pb")) return 1;
if(file.Contains("DYJetsToLL_Pt-100To250_FXFX_ext1-v1")) return 1;
if(file.Contains("DYJetsToLL_Pt-100To250_FXFX_ext2-v1")) return 1;
if(file.Contains("DYJetsToLL_Pt-100To250_FXFX_ext5-v1")) return 1;
if(file.Contains("DYJetsToLL_Pt-250To400_FXFX-v1")) return 3.047;
if(file.Contains("DYJetsToLL_Pt-250To400_FXFX_ext1-v1")) return 3.047;
if(file.Contains("DYJetsToLL_Pt-250To400_FXFX_ext2-v1")) return 3.047;
if(file.Contains("DYJetsToLL_Pt-250To400_FXFX_ext5-v1")) return 3.047;
if(file.Contains("DYJetsToLL_Pt-400To650_FXFX-v1")) return 0.3921;
if(file.Contains("DYJetsToLL_Pt-400To650_FXFX_ext1-v1")) return 0.3921;
if(file.Contains("DYJetsToLL_Pt-400To650_FXFX_ext2-v1")) return 0.3921;
if(file.Contains("DYJetsToLL_Pt-650ToInf_FXFX-v1")) return 0.03636;
if(file.Contains("DYJetsToLL_Pt-650ToInf_FXFX_ext1-v1")) return 0.03636;
if(file.Contains("DYJetsToLL_Pt-650ToInf_FXFX_ext2-v1")) return 0.03636;
if(file.Contains("WJetsToLNu_FXFX")) return 60290.0;
if(file.Contains("WJetsToLNu_Pt-100To250_FXFX")) return 676.3;
//if(file.Contains("WJetsToLNu_Pt-100To250_FXFX_ext1-v1")) return 676.3;
//if(file.Contains("WJetsToLNu_Pt-100To250_FXFX_ext4-v1")) return 676.3;
if(file.Contains("WJetsToLNu_Pt-250To400_FXFX")) return 23.94;
//if(file.Contains("WJetsToLNu_Pt-250To400_FXFX_ext1-v1")) return 23.94;
//if(file.Contains("WJetsToLNu_Pt-250To400_FXFX_ext4-v1")) return 23.94;
if(file.Contains("WJetsToLNu_Pt-400To600_FXFX")) return 3.031;
//if(file.Contains("WJetsToLNu_Pt-400To600_FXFX_ext1-v1")) return 3.031;
if(file.Contains("WJetsToLNu_Pt-600ToInf_FXFX")) return 0.4524;
//if(file.Contains("WJetsToLNu_Pt-600ToInf_FXFX_ext1-v1")) return 0.4524;

}


template <typename T1, typename T2> const bool Contains( TTreeReaderValue<T1> & Vec, const T2& Element )
{
    if (std::find(Vec->begin(), Vec->end(), Element) != Vec->end())
        return true;
    return false;
}
                                                                                                                       
void makeKFactor(string inputDIR_LO, string inputDIR_NLO, string outputDIR, Sample sample){
//gROOT->SetBatch(kTRUE);
setTDRStyle();
float scale_lo = 1;
float scale_nlo = 1;

system(("mkdir -p "+outputDIR).c_str());
vector<TTree*> tree_LO;
vector<TTree*> tree_NLO;
vector<TFile*> file_LO;
vector<TFile*> file_NLO;
vector<float>  eff_LO;
vector<float>  eff_NLO;

system(("ls "+inputDIR_LO+" | grep root > list.txt").c_str());
ifstream infileLO ("list.txt");
string line;
if(infileLO.is_open()){
  while(!infileLO.eof()){
    getline(infileLO,line);
    if(!TString(line).Contains("root") || line == "") continue;
    file_LO.push_back(TFile::Open((inputDIR_LO+"/"+line).c_str(),"READ"));
    tree_LO.push_back((TTree*) file_LO.back()->FindObjectAny("EventTree"));
  }
}

infileLO.close();

system("rm list.txt");


system(("ls "+inputDIR_NLO+" | grep root > list.txt").c_str());
ifstream infileNLO ("list.txt");
if(infileNLO.is_open()){
  while(!infileNLO.eof()){
    getline(infileNLO,line);
    if(not TString(line).Contains("root") or line == "") continue;
    file_NLO.push_back(TFile::Open((inputDIR_NLO+"/"+line).c_str(),"READ"));
    tree_NLO.push_back((TTree*) file_NLO.back()->FindObjectAny("EventTree"));
  }
}

infileNLO.close();
system("rm list.txt");

//calculate sumwgt


vector<double> sumwgt_lo;
int ifile = 0;
for(auto tree : tree_LO){
  float xsec= GetXsec(file_LO.at(ifile)->GetName());

  TTreeReader reader (tree);
  TTreeReaderValue<float> wgt (reader,"genWeight");
  TTreeReaderValue<std::vector<int>> mcMomPID(reader,"mcMomPID");
  cout<<"Calculate sumwgt for LO file "<<file_LO.at(ifile)->GetName()<<endl;
  double sumwgt = 0;
  float efficiency=1;
  while(reader.Next()){
    // filter away bad events with no matching
    //if(!Contains(mcMomPID,24)&&!Contains(mcMomPID,-24)) continue;
    sumwgt += *wgt;
  }
  cout<<"Tree LO with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<endl;
  TH1F * Histo_WLO = (TH1F*) file_LO.at(ifile)->Get("hcount");
  efficiency=sumwgt/Histo_WLO->GetBinContent(2);
  sumwgt = Histo_WLO->GetBinContent(2); 
  std::cout<<"efficiency::"<<efficiency<<std::endl;
  sumwgt_lo.push_back(sumwgt);
  eff_LO.push_back(efficiency);
  ifile++;
}


vector<double> sumwgt_nlo;
ifile = 0;

for(auto tree : tree_NLO){
  float xsec= GetXsec(file_NLO.at(ifile)->GetName());
  TTreeReader reader (tree);
  TTreeReaderValue<float> wgt (reader,"genWeight");
  TTreeReaderValue<std::vector<int>> mcMomPID(reader,"mcMomPID");
  cout<<"Calculate sumwgt for NLO file "<<file_NLO.at(ifile)->GetName()<<endl;
  double sumwgt = 0;
  float efficiency=1;
  while(reader.Next()){
    // filter away bad events with no matching
 //   if(!Contains(mcMomPID,24)&&!Contains(mcMomPID,-24)) continue;
    sumwgt += *wgt;
  }
  cout<<"Tree NLO with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<endl;
  TH1F * Histo_WNLO = (TH1F*) file_NLO.at(ifile)->Get("hcount");
  efficiency=sumwgt/Histo_WNLO->GetBinContent(2);
  sumwgt = Histo_WNLO->GetBinContent(2);
  std::cout<<"efficiency::"<<efficiency<<std::endl;
  sumwgt_nlo.push_back(sumwgt);
  eff_NLO.push_back(efficiency);
  ifile++;
}

TH1F* bosonPt_LO = new TH1F("bosonPt_LO","",bosonPt.size()-1,&bosonPt[0]);
TH1F* bosonPt_NLO = new TH1F("bosonPt_NLO","",bosonPt.size()-1,&bosonPt[0]);

bosonPt_LO->Sumw2();
bosonPt_NLO->Sumw2();

ifile = 0;

for(auto tree: tree_LO){
  float xsec= GetXsec(file_LO.at(ifile)->GetName());
  TTreeReader reader (tree);
  TTreeReaderValue<float> wgt (reader,"genWeight");
  TTreeReaderValue<std::vector<int>> mcPID (reader,"mcPID");
  TTreeReaderValue<std::vector<int>> mcMomPID (reader,"mcMomPID");
  TTreeReaderValue<std::vector<float>> mcMomPt (reader,"mcMomPt");
  TTreeReaderValue<std::vector<float>> mcMomEta (reader,"mcMomEta");
  TTreeReaderValue<std::vector<float>> mcMomPhi (reader,"mcMomPhi");
  TTreeReaderValue<std::vector<float>> mcMomMass (reader,"mcMomMass");
  TTreeReaderValue<float> genMET (reader,"genMET");
  TTreeReaderValue<float> genMETPhi (reader,"genMETPhi");
  TTreeReaderValue<std::vector<int>> jetGenPartonID (reader,"jetGenPartonID");
  TTreeReaderValue<std::vector<float>> jetGenPt (reader,"jetGenPt");
  TTreeReaderValue<std::vector<float>> jetGenEta (reader,"jetGenEta");
  TTreeReaderValue<std::vector<float>> jetGenPhi (reader,"jetGenPhi");
  TTreeReaderValue<std::vector<float>> jetGenEn (reader,"jetGenEn");
  TTreeReaderValue<std::vector<int>> jetGenPartonMomID  (reader,"jetGenPartonMomID");
  TTreeReaderValue<std::vector<float>> mcMass (reader,"mcMass");
  TTreeReaderValue<std::vector<float>> mcPt (reader,"mcPt");
  TTreeReaderValue<std::vector<float>> mcEta (reader,"mcEta");
  TTreeReaderValue<std::vector<float>> mcPhi (reader,"mcPhi");
  TTreeReaderValue<std::vector<int>> mcStatus (reader,"mcStatus");
  TTreeReaderValue<int> nMC (reader,"nMC");

  cout<<"Loop on LO file "<<file_LO.at(ifile)->GetName()<<endl;
  while(reader.Next()){
    float mindphi = 0.5;
    bool isclosejet=false;
    bool isgoodevt=false;
    float WBosonPt=-100.0;
    for (int igen=0;igen < *nMC; igen++){
                if (fabs(mcPID->at(igen)) ==24   && (mcStatus->at(igen) ==22 || mcStatus->at(igen) ==44 || mcStatus->at(igen) ==62)){ 
                   WBosonPt= mcPt->at(igen); // In inclusive we have status 62||22||44 while in HTbins we have just 22
                   isgoodevt=true; 
                }
    }

 
    if (sample == Sample::wjet && isgoodevt ){   
      for(size_t ijet = 0; ijet < jetGenPhi->size(); ijet++){
        if(jetGenPhi->at(ijet)<-10) continue;
        float dphi = fabs(*genMETPhi-jetGenPhi->at(ijet));
        if(dphi > TMath::Pi())
          dphi = 2*TMath::Pi()-dphi;
        if(dphi < mindphi)
          {isclosejet=true;break;}
      }



      //if(sample == Sample::wjet and (fabs(*l1id) == 15 or fabs(*l1id) == 16 ) continue; // skip taus
      if(jetGenPt->at(0) > 100){
        bosonPt_LO->Fill(WBosonPt,lumi_*(*wgt)*(xsec)*scale_lo/(sumwgt_lo.at(ifile)));
      //bosonPt_LO->Fill(wzpt1,(*wgt));
      }

    }
  }
  ifile++;
}

ifile = 0;
for(auto tree: tree_NLO){
  float xsec= GetXsec(file_NLO.at(ifile)->GetName());
  
  std::cout<<"xsec::"<<xsec<<std::endl;
  TTreeReader reader (tree);
  TTreeReaderValue<float> wgt (reader,"genWeight");
  TTreeReaderValue<std::vector<int>> mcPID (reader,"mcPID");
  TTreeReaderValue<std::vector<int>> mcMomPID (reader,"mcMomPID");
  TTreeReaderValue<std::vector<float>> mcMomPt (reader,"mcMomPt");
  TTreeReaderValue<std::vector<float>> mcMomEta (reader,"mcMomEta");
  TTreeReaderValue<std::vector<float>> mcMomPhi (reader,"mcMomPhi");
  TTreeReaderValue<std::vector<float>> mcMomMass (reader,"mcMomMass");
  TTreeReaderValue<float> genMET (reader,"genMET");
  TTreeReaderValue<float> genMETPhi (reader,"genMETPhi");
  TTreeReaderValue<std::vector<int>> jetGenPartonID (reader,"jetGenPartonID");
  TTreeReaderValue<std::vector<float>> jetGenPt (reader,"jetGenPt");
  TTreeReaderValue<std::vector<float>> jetGenEta (reader,"jetGenEta");
  TTreeReaderValue<std::vector<float>> jetGenPhi (reader,"jetGenPhi");
  TTreeReaderValue<std::vector<float>> jetGenEn (reader,"jetGenEn");
  TTreeReaderValue<std::vector<int>> jetGenPartonMomID  (reader,"jetGenPartonMomID");
  TTreeReaderValue<std::vector<float>> mcMass (reader,"mcMass");
  TTreeReaderValue<std::vector<float>> mcPt (reader,"mcPt");
  TTreeReaderValue<std::vector<float>> mcEta (reader,"mcEta");
  TTreeReaderValue<std::vector<float>> mcPhi (reader,"mcPhi");
  TTreeReaderValue<std::vector<int>> mcStatus (reader,"mcStatus");
  TTreeReaderValue<int> nMC (reader,"nMC");
  cout<<"Loop on NLO file "<<file_NLO.at(ifile)->GetName()<<endl;
  while(reader.Next()){
    float mindphi = 0.5;         
    bool isclosejet=false;   
    bool isgoodevt=false;
    float WBosonPt=-100.0; 
    for (int igen=0;igen < *nMC; igen++){
                if (fabs(mcPID->at(igen))==24   && mcStatus->at(igen)==22){
                   WBosonPt= mcPt->at(igen); // In inclusive we have status 62||22||44 while in HTbins we have just 22
                   isgoodevt=true;
                }
    }

    if (sample == Sample::wjet && isgoodevt){
      for(size_t ijet = 0; ijet < jetGenPhi->size(); ijet++){
        if(jetGenPhi->at(ijet)<-10) continue;
        float dphi = fabs(*genMETPhi-jetGenPhi->at(ijet));
        if(dphi > TMath::Pi())
          dphi = 2*TMath::Pi()-dphi;
        if(dphi < mindphi)
          {isclosejet=true;break;}
      }



      //if(sample == Sample::wjet and (fabs(*l1id) == 15 or fabs(*l1id) == 16 ) continue; // skip taus
     
      if(jetGenPt->at(0) > 100){
         if (WBosonPt <= 100.00) bosonPt_NLO->Fill(WBosonPt, lumi_*(*wgt)*scale_nlo/ ( sumwgt_nlo.at(0)/ GetXsec("WJetsToLNu_FXFX.root")));
         else if(WBosonPt > 100.00 && WBosonPt<=250) bosonPt_NLO->Fill(WBosonPt, lumi_*(*wgt)*scale_nlo/ (sumwgt_nlo.at(1)/ GetXsec("WJetsToLNu_Pt-100To250_FXFX") +  sumwgt_nlo.at(0)/ GetXsec("WJetsToLNu_FXFX.root")));
         else if(WBosonPt > 250.0 && WBosonPt<=400) bosonPt_NLO->Fill(WBosonPt, lumi_*(*wgt)*scale_nlo/ (sumwgt_nlo.at(2)/ GetXsec("WJetsToLNu_Pt-250To400_FXFX") +  sumwgt_nlo.at(0)/ GetXsec("WJetsToLNu_FXFX.root")));
         else if(WBosonPt > 400.0 && WBosonPt<=600) bosonPt_NLO->Fill(WBosonPt, lumi_*(*wgt)*scale_nlo/ (sumwgt_nlo.at(3)/ GetXsec("WJetsToLNu_Pt-400To600_FXFX") +  sumwgt_nlo.at(0)/ GetXsec("WJetsToLNu_FXFX.root")));
         else if(WBosonPt > 600.0 ) bosonPt_NLO->Fill(WBosonPt, lumi_*(*wgt)*scale_nlo/ (sumwgt_nlo.at(4)/ GetXsec("WJetsToLNu_Pt-600ToInf_FXFX") +  sumwgt_nlo.at(0)/ GetXsec("WJetsToLNu_FXFX.root")));

     //   bosonPt_NLO->Fill(WBosonPt,lumi_*(*wgt)*(xsec)*scale_nlo/(sumwgt_nlo.at(ifile)));
     //   bosonPt_NLO->Fill(wzpt1,(*wgt));
      }

    }
  }
  ifile++;
}
      

TFile* output = new TFile((outputDIR+"/output.root").c_str(),"RECREATE");
output->cd();
bosonPt_NLO->Write();
bosonPt_LO->Write();
TCanvas* canvas = new TCanvas("canvas","canvas",600,625);
TPad* pad1 = new TPad("pad1","",0.,0.3,1,1);
pad1->Draw();
canvas->cd();
TPad* pad2 = new TPad("pad2","",0.,0.,1,0.3);
pad2->Draw();
canvas->cd();
canvas->cd();
pad1->cd();

TH1F* Pt_NLO = (TH1F*) bosonPt_NLO->Clone("Pt_NLO");
TH1F* Pt_LO = (TH1F*) bosonPt_LO->Clone("Pt_LO");

std::cout<<"integral::"<<Pt_NLO->Integral()<<std::endl;
std::cout<<"sumwgt_nlo::"<<sumwgt_nlo.at(0)<<std::endl;

//Pt_NLO->Scale((1.0*60290)/(sumwgt_nlo.at(0)),"width");
//Pt_LO->Scale((1.0*50690)/(sumwgt_lo.at(0)),"width");
Pt_NLO->Scale(1,"width");
Pt_LO->Scale(1,"width");

Pt_NLO->GetXaxis()->SetTitle("Boson p_{T} [GeV]");
Pt_NLO->GetYaxis()->SetTitle("d#sigma/dp_{T}[pb/GeV]");
Pt_NLO->GetYaxis()->SetTitleOffset(1.1);
Pt_NLO->GetXaxis()->SetTitleOffset(1.1);
Pt_NLO->GetYaxis()->SetRangeUser((Pt_NLO->GetMinimum())*0.7,(Pt_NLO->GetMaximum())*1.3);

TH1F* Pt_NLO_band = (TH1F*) Pt_NLO->Clone("Pt_NLO_band");
Pt_NLO_band->SetFillColor(kBlack);
Pt_NLO_band->SetFillStyle(3002);

TH1F* Pt_LO_band = (TH1F*) Pt_LO->Clone("Pt_LO_band");
Pt_LO_band->SetFillColor(kRed);
Pt_LO_band->SetFillStyle(3002);



Pt_NLO->SetLineColor(kBlack);
Pt_NLO->SetMarkerColor(kBlack);
Pt_NLO->SetLineWidth(2);
Pt_NLO->SetMarkerSize(1);
Pt_NLO->Draw("hist");
Pt_NLO_band->Draw("E2 same");
Pt_NLO->Draw("hist same");



Pt_LO->SetLineColor(kRed);
Pt_LO->SetMarkerColor(kRed);
Pt_LO->SetLineWidth(2);
Pt_LO->SetMarkerSize(1);
Pt_LO->Draw("hist same");
Pt_LO_band->Draw("E2 same");





TLegend leg (0.6,0.6,0.9,0.9);
leg.SetFillColor(0);
leg.SetFillStyle(0);
leg.SetBorderSize(0);
leg.AddEntry(Pt_NLO,"W+jets NLO","FL");
leg.AddEntry(Pt_LO,"W+jets LO","FL");
leg.Draw("same");

TLatex* latex = new TLatex ();
latex->SetNDC();
latex->SetTextSize(0.6*gPad->GetTopMargin());
latex->SetTextFont(42);
latex->SetTextAlign(31);
if(sample == Sample::zll)
  latex->DrawLatex(0.25,0.95,"DY+jets");
else if(sample == Sample::wjet)
  latex->DrawLatex(0.25,0.95,"W+jets");
pad2->cd();

TH1F* ratio_1 = (TH1F*) Pt_NLO->Clone("kfactor");
ratio_1->Divide(Pt_LO);
TH1F* ratio_2 = (TH1F*) Pt_LO->Clone("kfactor");
ratio_2->Divide(Pt_LO);

ratio_1->GetYaxis()->SetRangeUser(0,4);


ratio_1->SetLineColor(kBlack);
ratio_1->SetMarkerColor(kBlack);
ratio_1->SetMarkerStyle(20);
ratio_1->GetYaxis()->SetTitle("Ratio");
ratio_1->Draw("EPsame");

ratio_2->SetLineColor(kRed);
ratio_2->SetMarkerColor(kRed);
ratio_2->SetMarkerStyle(20);
ratio_2->GetYaxis()->SetTitle("Ratio");
ratio_2->Draw("hist same");




if(sample == Sample::zll){
  canvas->SaveAs((outputDIR+"/kfactor_zll.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/kfactor_zll.pdf").c_str(),"pdf");
}
else if(sample == Sample::wjet){
  canvas->SaveAs((outputDIR+"/kfactor_wjet.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/kfactor_wjet.pdf").c_str(),"pdf");
}
}


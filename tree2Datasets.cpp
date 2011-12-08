#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include "TH2.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"

#include "RooFit.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"

#include <TCanvas.h>
#include "TH2.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include <TVector3.h>
#include <math.h>

using namespace RooFit;
using namespace std;


bool isAccept(const TLorentzVector* aMuon) {
   // *USE* muon kinematical cuts (eta dependent momentum / pT cuts )
/*   return (fabs(aMuon->Eta()) < 2.4 &&
           ((fabs(aMuon->Eta()) < 1.3 && aMuon->Pt() > 3.3) ||
           (fabs(aMuon->Eta()) > 1.3 && fabs(aMuon->Eta()) < 2.2 && aMuon->P() > 2.9) ||
           (fabs(aMuon->Eta()) > 2.2 && aMuon->Pt() > 0.8)));
*/
   // *REMOVE* muon kinematical cuts (eta dependent momentum / pT cuts )
   // by just returning TRUE
     return true;
}

double CorrectMass(const TLorentzVector* mu1,const TLorentzVector* mu2, int mode){  
  double CMass=0;
  const double mumass=0.105658;
  double k1,k2;
  double pt1=mu1->Pt();
  double pt2=mu2->Pt();
  double eta1=mu1->Eta();
  double eta2=mu2->Eta();
  if (mode==1){
    k1=1.0009;//constant scale correction
    k2=1.0009;
  }
  if (mode==2){
    k1=1.0019-0.0004*pt1;
    k2=1.0019-0.0004*pt2; // pt dependent correction
  }
  if (mode==3){
      double a0=0.00038; //3.8 * pow(10,-4);
      double a1=0.0;
      double a2=0.0003; //3.0 * pow(10,-4);
      double a3=0.0;

      k1=1+a0+a1*fabs(eta1)+a2*eta1*eta1+a3*pt1;
      k2=1+a0+a1*fabs(eta2)+a2*eta2*eta2+a3*pt2;// pt and eta dependent
  }

  if (mode == 4){
      double a0=1.002;
      double a1=-0.002;
      double a2=0.001;
      double a3=-0.0001;

      k1=a0+a1*fabs(eta1)+a2*eta1*eta1+a3*pt1;
      k2=a0+a1*fabs(eta2)+a2*eta2*eta2+a3*pt2;// pt and eta dependent
  }

  TVector3 mom1=mu1->Vect();
  TVector3 mom2=mu2->Vect();
  mom1=k1*mom1; 
  mom2=k2*mom2;
  double E1=sqrt(mom1.Mag2()+(mumass*mumass));
  double E2=sqrt(mom2.Mag2()+(mumass*mumass));
  TVector3 momtot=mom1+mom2;
  CMass=sqrt((E1+E2)*(E1+E2)-momtot.Mag2());
  return CMass;
}



int main(int argc, char* argv[]) {
  
  const double Jpsi_MassMin=2.6;
  const double Jpsi_MassMax=3.5;
  const double Jpsi_PtMin=0;
  const double Jpsi_PtMax=100;
  const double Jpsi_YMin=0;
  const double Jpsi_YMax=2.4;
  const double Jpsi_CtMin = -3.0;
  const double Jpsi_CtMax = 3.5;
  
  char fileName[100];
  
  if ( argc != 5 && argc != 2 ){
    char msg[300];
    sprintf(msg,"Usage1: %s [input file] [directory name for output files] [start event #] [end event # (-1 for total)]",argv[0]);
    cout << msg << endl; 
    sprintf(msg,"Usage2: %s [input file] [directory name for output files]]",argv[0]);
    cout << msg << endl; 
    return 1;
  }
  strcpy(fileName,argv[1]);


  TFile *file= TFile::Open(fileName);
  TTree * Tree=(TTree*)file->Get("myTree");

  Int_t           Centrality;
  Int_t           Reco_QQ_size;
  Int_t           Reco_QQ_type[20];   //[Reco_QQ_size]
  Int_t           Reco_QQ_sign[20];   //[Reco_QQ_size]
  Int_t           Reco_QQ_trig[20];   //[Reco_QQ_size]
  Int_t           Gen_QQ_type[20];
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_QQ_mupl_4mom;
  TClonesArray    *Reco_QQ_mumi_4mom;
  Float_t         Reco_QQ_ctau[20];   //[Reco_QQ_size]
  Float_t         Reco_QQ_ctauErr[20];   //[Reco_QQ_size]
  Float_t         Reco_QQ_ctauTrue[20];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[20];   //[Reco_QQ_size]

  TBranch        *b_Centrality;   //!
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_QQ_type;   //!
  TBranch        *b_Reco_QQ_sign;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Gen_QQ_type;
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_QQ_mupl_4mom;   //!
  TBranch        *b_Reco_QQ_mumi_4mom;   //!
  TBranch        *b_Reco_QQ_ctau;   //!
  TBranch        *b_Reco_QQ_ctauErr;   //!
  TBranch        *b_Reco_QQ_ctauTrue;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!

  TLorentzVector* JP= new TLorentzVector;
  TLorentzVector* m1P= new TLorentzVector;
  TLorentzVector* m2P= new TLorentzVector;
  
  double vprob, theCt, theCtErr, theCtTrue;
  int trig,trig2,theCat,Jq,genType;

  static const unsigned int centRegions = 6;
  int centLimits[centRegions+1] = {0, 4, 8, 12, 16, 20, 40};
//  static const unsigned int centRegions = 2;
//  int centLimits[centRegions+1] = {0, 8, 40};

//  static const unsigned int ptRegions = 3;
//  float ptLimits[ptRegions+1] = {0.0, 6.5, 10.0, 30.0};
  
  static const unsigned int rapRegions = 3;
  float rapLimits[rapRegions+1] = {Jpsi_YMin,1.2,1.6,Jpsi_YMax};
//   static const unsigned int rapRegions = 1;
//   float rapLimits[rapRegions+1] = {Jpsi_YMin,Jpsi_YMax};
/*  RooDataSet* dataJpsi[rapRegions+1][centRegions+1];
  RooDataSet* dataJpsiSame[rapRegions+1][centRegions+1];
  RooDataSet* dataPsip[rapRegions+1][centRegions+1];*/
  RooDataSet* dataJpsi[centRegions+1];
  RooDataSet* dataJpsiSame[centRegions+1];
  RooDataSet* dataPsip[centRegions+1];
  RooRealVar* Jpsi_Mass;
  RooRealVar* Psip_Mass;      
  RooRealVar* Jpsi_Pt;
  RooRealVar* Jpsi_Ct;
  RooRealVar* Jpsi_CtErr;
  RooRealVar* Jpsi_CtTrue;
  RooRealVar* Jpsi_Y;
  RooCategory* Jpsi_Type;
  RooCategory* Jpsi_Sign;
  RooCategory* MCType;

  Jpsi_Mass = new RooRealVar("Jpsi_Mass","J/#psi mass",Jpsi_MassMin,Jpsi_MassMax,"GeV/c^{2}");
  Psip_Mass = new RooRealVar("Psip_Mass","#psi' mass",3.3,Jpsi_MassMax,"GeV/c^{2}");
  Jpsi_Pt = new RooRealVar("Jpsi_Pt","J/#psi pt",Jpsi_PtMin,Jpsi_PtMax,"GeV/c");
  Jpsi_Y = new RooRealVar("Jpsi_Y","J/#psi y",-Jpsi_YMax,Jpsi_YMax);
  Jpsi_Type = new RooCategory("Jpsi_Type","Category of Jpsi_");
  Jpsi_Sign = new RooCategory("Jpsi_Sign","Charge combination of Jpsi_");
  Jpsi_Ct = new RooRealVar("Jpsi_Ct","J/#psi c#tau",Jpsi_CtMin,Jpsi_CtMax,"mm");
  Jpsi_CtErr = new RooRealVar("Jpsi_CtErr","J/#psi c#tau error",-1.,1.,"mm");
//  MCType = new RooCategory("MCType","Type of generated Jpsi_");
//  Jpsi_CtTrue = new RooRealVar("Jpsi_CtTrue","J/#psi c#tau true",Jpsi_CtMin,Jpsi_CtMax,"mm");

  Jpsi_Type->defineType("GG",0);
  Jpsi_Type->defineType("GT",1);
  Jpsi_Type->defineType("TT",2);

  Jpsi_Sign->defineType("OS",0);
  Jpsi_Sign->defineType("PP",1);
  Jpsi_Sign->defineType("MM",2);

//  MCType->defineType("PR",0);
//  MCType->defineType("NP",1);
  
  Reco_QQ_4mom = 0;
  Reco_QQ_mupl_4mom = 0;
  Reco_QQ_mumi_4mom = 0;

  Tree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  Tree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  Tree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  Tree->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
  Tree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  Tree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  Tree->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
  Tree->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
  Tree->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
  Tree->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
//  Tree->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
//  Tree->SetBranchAddress("Reco_QQ_ctauTrue", Reco_QQ_ctauTrue, &b_Reco_QQ_ctauTrue);
  Tree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

  RooArgList varlist(*Jpsi_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*MCType,*Jpsi_Ct,*Jpsi_CtErr,*Jpsi_CtTrue);
  RooArgList varlistSame(*Jpsi_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*MCType,*Jpsi_Ct,*Jpsi_CtErr,*Jpsi_CtTrue);
  RooArgList varlist2(*Psip_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*MCType,*Jpsi_Ct,*Jpsi_CtErr,*Jpsi_CtTrue);
    
  for (unsigned int j = 0; j <= centRegions; j++) {
    dataJpsi[j] = new RooDataSet("dataJpsi","A sample",varlist);
    dataJpsiSame[j] = new RooDataSet("dataJpsiSame","A sample",varlistSame);
    dataPsip[j] = new RooDataSet("dataPsip","A sample",varlist2);
  }
/*    for (unsigned int i = 0; i <= rapRegions; i++) {
      dataJpsi[i][j] = new RooDataSet("dataJpsi","A sample",varlist);
      dataJpsiSame[i][j] = new RooDataSet("dataJpsiSame","A sample",varlistSame);
      dataPsip[i][j] = new RooDataSet("dataPsip","A sample",varlist2);
    }*/

  int initev = 0;
  int nevt = Tree->GetEntries();
  if (argc == 5) {
    initev = atoi(argv[3]);
    if (argv[4] != -1)
      nevt = atoi(argv[4]);
  }

  for (int ev=initev; ev<nevt; ++ev) {
    if (ev%5000==0) cout << ">>>>> EVENT " << ev << " / " << Tree->GetEntries()<<  endl;
    
    Tree->GetEntry(ev);

    int theCentrality=Centrality;
    for (int i=0; i<Reco_QQ_size; ++i) {
      JP = (TLorentzVector*) Reco_QQ_4mom->At(i);
      m1P = (TLorentzVector*) Reco_QQ_mupl_4mom->At(i);
      m2P = (TLorentzVector*) Reco_QQ_mumi_4mom->At(i);
      vprob = Reco_QQ_VtxProb[i];
      trig = ((Reco_QQ_trig[i] & 1) == 1) ? 1 : 0; // HLT_HIL1DoubleMu0_HighQ
      trig2 = ((Reco_QQ_trig[i] & 4) == 4) ? 1 : 0; // HLT_HIL2DoubleMu3
      theCat = Reco_QQ_type[i];
      Jq = Reco_QQ_sign[i];
      theCt = Reco_QQ_ctau[i];
      theCtErr = Reco_QQ_ctauErr[i];
//      genType = Gen_QQ_type[i];
//      theCtTrue = Reco_QQ_ctauTrue[i];
      
      double theMass =JP->M();
//      double CMass=CorrectMass(m1P,m2P,4);
      //      //cout << " Mass " << theMass << "   Corrected " << CMass << endl; 
//      if (CMass!=0) theMass=CMass;
      
      double theRapidity=JP->Rapidity();
      double thePt=JP->Pt();
      
//       cout << Jq << endl;
//       cout << thePt<< endl;
//       cout << theMass<< endl;
//       cout << theRapidity << endl;
//       cout << theCt<< endl;
//       cout << vprob << endl;
//       cout << "Triggered " << trig << "\t" << trig2 << endl;
      bool ok1=isAccept(m1P);
      bool ok2=isAccept(m2P);
      

      if (theMass > Jpsi_MassMin && theMass < Jpsi_MassMax && 
//      	  theCt > Jpsi_CtMin && theCt < Jpsi_CtMax && 
//	        thePt > Jpsi_PtMin && thePt < Jpsi_PtMax && 
//	        fabs(theRapidity) > Jpsi_YMin && fabs(theRapidity) < Jpsi_YMax &&
          ok2 && ok1 &&
      	  (trig == 1 || trig2 == 1) &&
      	  vprob > 0.001) {// &&
      
	Jpsi_Pt->setVal(thePt); 
	Jpsi_Y->setVal(theRapidity); 
	Jpsi_Mass->setVal(theMass);
	Psip_Mass->setVal(theMass);
	Jpsi_Ct->setVal(theCt);
	Jpsi_CtErr->setVal(theCtErr);
	Jpsi_Type->setIndex(theCat,kTRUE);
  if (Jq == 0){
    Jpsi_Sign->setIndex(Jq,kTRUE);
  } else {
    Jpsi_Sign->setIndex(Jq,kTRUE);
  }
	RooArgList varlist_tmp(*Jpsi_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*Jpsi_Ct,*Jpsi_CtErr);
	RooArgList varlist2_tmp(*Psip_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*Jpsi_Ct,*Jpsi_CtErr);

//	Jpsi_CtTrue->setVal(theCtTrue);
//	MCType->setIndex(genType,kTRUE);
//	RooArgList varlist_tmp(*Jpsi_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*MCType,*Jpsi_Ct,*Jpsi_CtErr,*Jpsi_CtTrue);
//	RooArgList varlist2_tmp(*Psip_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*MCType,*Jpsi_Ct,*Jpsi_CtErr,*Jpsi_CtTrue);


	for (unsigned int j = 0; j <= centRegions; j++) {
	  if ( (j==centRegions && theCentrality < centLimits[j] && theCentrality >= centLimits[0]) ||
	       (theCentrality < centLimits[j+1] && theCentrality >= centLimits[j]) ) {
      if (Jq == 0) {
        if (theMass < 3.5) dataJpsi[j]->add(varlist_tmp);
        if (theMass > 3.3) dataPsip[j]->add(varlist2_tmp);
      } else {
        if (theMass < 3.5) dataJpsiSame[j]->add(varlist_tmp);
      }
	  }
	}
	    /*for (unsigned int k = 0; k <= rapRegions; k++) {
	      if ( (k==rapRegions && fabs(theRapidity) < rapLimits[k] && fabs(theRapidity) >= rapLimits[0]) ||
		   (fabs(theRapidity) < rapLimits[k+1] && fabs(theRapidity) >= rapLimits[k]) ) {
          //if (theMass < 3.5) dataJpsi[k][j]->add(varlist_tmp);
          if (Jq == 0) {
            if (theMass < 3.5) dataJpsi[k][j]->add(varlist_tmp);
            if (theMass > 3.3) dataPsip[k][j]->add(varlist2_tmp);
          } else {
            if (theMass < 3.5) dataJpsiSame[k][j]->add(varlist_tmp);
          }
	      }
	    }*/
      }
    }
  }

  TFile *Out[centRegions+1];
  for (unsigned int j = 0; j <= centRegions; j++) {
    char namefile[300];
    std::cout << centLimits[j] << std::endl;
    if (j==centRegions) {
      sprintf(namefile,"%s/Data2011_cent%d-%d.root",argv[2],
        int(centLimits[0]*2.5),int(centLimits[j]*2.5));
    } else {
      sprintf(namefile,"%s/Data2011_cent%d-%d.root",argv[2],
        int(centLimits[j]*2.5),int(centLimits[j+1]*2.5));
    }
    Out[j] = new TFile(namefile,"RECREATE");
    Out[j]->cd();
    dataJpsi[j]->Write();
    dataJpsiSame[j]->Write();
    dataPsip[j]->Write();
    Out[j]->Close();
  }

  
  
/*  TFile* Out[rapRegions+1][centRegions+1];
  for (unsigned int j = 0; j <= centRegions; j++) {
    std::cout << centLimits[j] << std::endl;
    for (unsigned int i = 0; i <= rapRegions; i++) {
      if (j==centRegions && i==rapRegions) {
        sprintf(namefile,"datasets_np/Data2011_rap%.1f-%.1f_cent%d-%d.root",
          rapLimits[0],rapLimits[i],
          int(centLimits[0]*2.5),int(centLimits[j]*2.5));
      }
      else if (j==centRegions) {
        sprintf(namefile,"datasets_np/Data2011_rap%.1f-%.1f_cent%d-%d.root",
          rapLimits[i],rapLimits[i+1],
          int(centLimits[0]*2.5),int(centLimits[j]*2.5));
      }
      else if (i==rapRegions) {
        sprintf(namefile,"datasets_np/Data2011_rap%.1f-%.1f_cent%d-%d.root",
          rapLimits[0],rapLimits[i],
          int(centLimits[j]*2.5),int(centLimits[j+1]*2.5));
      }
      else {
        sprintf(namefile,"datasets_np/Data2011_rap%.1f-%.1f_cent%d-%d.root",
          rapLimits[i],rapLimits[i+1],
          int(centLimits[j]*2.5),int(centLimits[j+1]*2.5));
      }
      Out[i][j] = new TFile(namefile,"RECREATE");
      Out[i][j]->cd();
      dataJpsi[i][j]->Write();
      dataJpsiSame[i][j]->Write();
      dataPsip[i][j]->Write();
      Out[i][j]->Close();
    }
  }*/

  return 0;

}



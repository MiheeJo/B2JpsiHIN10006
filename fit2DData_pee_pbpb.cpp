#include <iostream>
#include <sstream>

#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TRandom.h>

#include <RooFit.h>
#include <RooGlobalFunc.h>
#include <RooCategory.h>
#include <RooHistPdfConv.h>
#include <RooWorkspace.h>
#include <RooBinning.h>
#include <RooHistPdf.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>

using namespace RooFit;
using namespace std;

bool superImpose = false;
bool analyticBlifetime = true;
bool narrowSideband = false;
bool oneGaussianResol = false;

void getOptRange(string &ran,float *min,float *max);
void setWSRange(RooWorkspace *ws, float lmin, float lmax, float errmin, float errmax);
RooBinning setCtBinning(float lmin,float lmax);
void defineMassBkg(RooWorkspace *ws);
void defineMassSig(RooWorkspace *ws);
void getMCTrueLifetime(RooWorkspace *ws, RooDataSet *redMCCutNP, float *bgmcVal, float *bctauVal, string titlestr);
void defineCTResol(RooWorkspace *ws);
void defineCTBkg(RooWorkspace *ws);
void defineCTSig(RooWorkspace *ws, RooDataSet *redMCCutNP, string titlestr);
RooDataHist* subtractSidebands(RooWorkspace* ws, RooDataHist* all, RooDataHist* side, float scalefactor, string varName);

int main (int argc, char* argv[]) {

  gROOT->Macro("/afs/cern.ch/user/m/miheejo/JpsiStyle.C");
  string FileName, FileNameMC1, FileNameMC2;
  int isGG = 0;
  bool prefitMass = false;
  bool prefitSignalCTau = false;
  bool prefitBkg = false;
  bool fracfix = false;
  bool isMB = false;
  bool isMBPT = false;
  bool isPT = false;
  string prange, yrange, crange, lrange, errrange;

  // *** Check options
  for (int i=1; i<argc; i++) {
    char *tmpargv = argv[i];
    switch (tmpargv[0]) {
      case '-':{
        switch (tmpargv[1]) {
          case 'f':
            FileName = argv[i+1];
            cout << "Fitted data file: " << FileName << endl;
            break;
          case 'm':
            FileNameMC1 = argv[i+1];
            cout << "MC data file 1: " << FileNameMC1 << endl;
            FileNameMC2 = argv[i+2];
            cout << "MC data file 2: " << FileNameMC2 << endl;
            break;
          case 'p':
            prange = argv[i+1];
            cout << "pT range: " << prange << " GeV/c" << endl;
            break;
          case 'y':
            yrange = argv[i+1];
            cout << "Rapidity range: " << yrange << " rad" << endl;
            break;
          case 't':
            crange = argv[i+1];
            cout << "Centrality range: " << crange << " %" << endl;
            break;
          case 'l':
            lrange = argv[i+1];
            cout << "l(J/psi) range: " << lrange << " mm" << endl; 
            break;
          case 'e':
            errrange = argv[i+1];
            cout << "Range for sigma_l(J/psi) is " << errrange << " mm" << endl;
            break;    
          case 'u':
            prefitMass = true;
            cout << "Turn on: signal(=data, depends on muon pair type) mass pre-fitting" << endl;
            break;
          case 'c':
            prefitSignalCTau = true;
            cout << "Turn on: prompt signal ctau distribution pre-fitting" << endl;
            break;
          case 'b':
            prefitBkg = true;
            cout << "Turn on: Background(=sideband data) ctau distribution pre-fitting" << endl;
            break;
          case 'v':
            isPT = true;
            cout << "Fit for pT6.5-30.0" << endl;
            break;
          case 'x':
            if (0 == atoi(argv[i+1])) {
              isMB = true;
              cout << "Fit for minbias" << endl;
            } else if (1 == atoi(argv[i+1])) {
              isMBPT = true;
              cout << "Fit for _mb_pt" << endl;
            } else if (2 == atoi(argv[i+1])) {
              isMBPT = false;
              cout << "Fit for _mb" << endl;
            }
            break;
          case 'z':
            if (0 == atoi(argv[i+1])) {
              fracfix = true;
              cout << "Fix frac2,3 BEFORE ctau bkg fitting" << endl;
            } else {
              fracfix = false;
              cout << "Fix frac2,3 AFTER ctau bkg fitting" << endl;
            }
            break;

        }
      }
    }
  }// End check options
 
  float pmin=0, pmax=0, ymin=0, ymax=0, lmin=0, lmax=0, cmin=0, cmax=0, errmin=0, errmax=0;
  getOptRange(prange,&pmin,&pmax);
  getOptRange(lrange,&lmin,&lmax);
  getOptRange(errrange,&errmin,&errmax);
  getOptRange(crange,&cmin,&cmax);
  getOptRange(yrange,&ymin,&ymax);

  if (isMB && isPT) {
    cout << "Please use isMB or isPT, not both!" << endl;
    return 1;
  }

  string dirPre;   //Name of directory to store results
  if (fracfix) {
    dirPre = "./fracfix";
  } else {
    dirPre = "./fracfree";
  }


  // *** TFile for saving fitting results
/*  string resultFN;
  resultFN = dirPre + "/rap" + yrange + "_cent" + crange + "_pT" + prange + "/fitResult.root";
  TFile *resultF = new TFile(resultFN.c_str(),"RECREATE");*/

  // *** Read MC and Data files
  TFile *fInMC = new TFile(FileNameMC1.c_str());   //Non-prompt J/psi MC
  cout << FileNameMC1.c_str() << endl;
  RooDataSet *dataMC;
  if (fInMC->IsZombie()) { cout << "CANNOT open MC1 root file\n"; return 1; }
  fInMC->cd();
  dataMC = (RooDataSet*)fInMC->Get("dataJpsi");
  dataMC->SetName("dataMC");

  TFile *fInMC2 = new TFile(FileNameMC2.c_str());  //Prompt J/psi MC
  cout << FileNameMC2.c_str() << endl;
  if (fInMC2->IsZombie()) { cout << "CANNOT open MC2 root file\n"; return 1; }
  fInMC2->cd();
  RooDataSet *dataMC2 = (RooDataSet*)fInMC2->Get("dataJpsi");
  dataMC2->SetName("dataMC2");

  TFile fInData(FileName.c_str());
  cout << FileName.c_str() << endl;
  if (fInData.IsZombie()) { cout << "CANNOT open data root file\n"; return 1; }
  fInData.cd();
  RooDataSet *data = (RooDataSet*)fInData.Get("dataJpsi");
  data->SetName("data");

  // Create workspace to play with
  RooWorkspace *ws = new RooWorkspace("workspace");

  // Reduce "dataMC" with given ranges/cuts
  char reduceDS[300];
  sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_CtErr > %.2f && Jpsi_CtErr < %.2f",pmin,pmax,ymin,ymax,errmin,errmax);
  cout << "reduceDS for MC and data: " << reduceDS << endl;

  RooDataSet *redMC = (RooDataSet*)dataMC->reduce(reduceDS);
  ws->import(*redMC);
  RooDataSet *redMC2 = (RooDataSet*)dataMC2->reduce(reduceDS);
  ws->import(*redMC2);

  RooDataSet *redData = (RooDataSet*)data->reduce(reduceDS);
  ws->import(*redData);
 
  setWSRange(ws,lmin,lmax,errmin,errmax);

  // Draw data
  ws->var("Jpsi_Mass")->SetTitle("J/#psi mass");
  ws->var("Jpsi_Ct")->SetTitle("#font[12]{l}_{J/#psi}");
  ws->var("Jpsi_CtTrue")->setBins(2000);

  // Test true lifetimes for non-prompt J/psi
  RooPlot *trueFrame = ws->var("Jpsi_CtTrue")->frame();
  ws->data("dataMC")->plotOn(trueFrame,DataError(RooAbsData::SumW2),Cut("MCType==MCType::NP"));

  TCanvas c0;
  c0.cd(); trueFrame->Draw();
  string titlestr;
  titlestr = dirPre + "/rap" + yrange + "_cent" + crange + "_pT" + prange + "/testTrueLife_Lin.pdf";
  c0.SaveAs(titlestr.c_str());

  ws->var("Jpsi_Mass")->setBins(60);
//  ws->var("Jpsi_Mass")->setBins(133);
  ws->var("Jpsi_CtErr")->setBins(25);
  if (pmin > 40.) ws->var("Jpsi_CtErr")->setBins(8);


  // Define binning for true lifetime
  RooBinning rb(-0.1,4.0);
  rb.addUniform(5,-0.1,0.0);
  rb.addUniform(100,0.0,0.5);
  rb.addUniform(15,0.5,1.0);
  rb.addUniform(20,1.0,2.5);
  rb.addUniform(5,2.5,4.0);
  if (analyticBlifetime) {
    cout << "analytic B lifetime option is turned on!\n";
    ws->var("Jpsi_CtTrue")->setBins(200);
  } else {
    cout << "analytic B lifetime option is NOT turned on!\n";
    ws->var("Jpsi_CtTrue")->setBinning(rb);
  }

  // Define binning for lifetime
  RooBinning rb2 = setCtBinning(lmin,lmax);
  ws->var("Jpsi_Ct")->setBinning(rb2);
  
  // Define ctau binning for plotting (coarser bin)
  RooBinning rb3(-lmin,lmax);
  rb3.addBoundary(-1.0);
  rb3.addBoundary(-0.7);
  rb3.addBoundary(-0.6);
  rb3.addBoundary(-0.5);
  rb3.addUniform(5,-0.5,-0.2);
  rb3.addUniform(15,-0.2,0.2);
  rb3.addUniform(5,0.2,0.5);
  rb3.addUniform(5,0.5,1.0);

  // Additional cuts on data and get sub-datasets/histograms
  RooDataSet *redDataCut;
  string reduceDSstr;
  if (isGG == 0) {
    reduceDSstr = "Jpsi_Type == Jpsi_Type::GG &&\
                  (MCType != MCType::NP || Jpsi_CtTrue>0.0001) &&\
                  (MCType == MCType::PR || MCType == MCType::NP)";
    redDataCut = (RooDataSet*)redData->reduce("Jpsi_Type == Jpsi_Type::GG");
    cout << "Using Global-Global muon pairs in this processing!\n";
  } else if (isGG == 1) {
    reduceDSstr = "(Jpsi_Type == Jpsi_Type::GT || Jpsi_Type == Jpsi_Type::GG) &&\
                  (MCType != MCType::NP || Jpsi_CtTrue > 0.0001) &&\
                  (MCType == MCType::PR || MCType == MCType::NP)";
    redDataCut = (RooDataSet*)redData->reduce("Jpsi_Type == Jpsi_Type::GT || Jpsi_Type == Jpsi_Type::GG");
  } else {
    reduceDSstr = "(MCType != MCType::NP || Jpsi_CtTrue>0.0001) &&\
                   (MCType == MCType::PR || MCType == MCType::NP)";
    redDataCut = (RooDataSet*)redData->reduce("Jpsi_Ct < 600000.");
  }

  RooDataHist *binDataCut = new RooDataHist("binDataCut","binDataCut",RooArgSet( *(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct")),*(ws->var("Jpsi_CtErr")) ), *redDataCut);
  RooDataHist *binDataCutCtErr = new RooDataHist("binDataCutCtErr","binDataCutCtErr",RooArgSet( *(ws->var("Jpsi_CtErr")) ), *redDataCut);
  cout << "DATA :: N events to fit: " << binDataCut->sumEntries() << endl;

  // *** Get MC sub-datasets and its histograms corresponds to data
  RooDataSet *redMCCut = (RooDataSet*) redMC->reduce(reduceDSstr.c_str());
  RooDataSet *redMCCutNP = (RooDataSet*) redMCCut->reduce(RooArgSet(*(ws->var("Jpsi_CtTrue"))),"MCType == MCType::NP");
  RooDataSet *redMCCut2 = (RooDataSet*) redMC2->reduce(reduceDSstr.c_str());
  RooDataSet *redMCCutPR = (RooDataSet*) redMCCut2->reduce("MCType == MCType::PR");

  RooDataHist *binMCCut = new RooDataHist("binMCCut","MC distribution for signal",RooArgSet( *(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct")) ),*redMCCut);
  cout << "MC :: N events to fit: " << binMCCut->sumEntries() << endl;
  RooDataHist *binMCCutPR = new RooDataHist("binMCCutPR","MC distribution for PR signal",RooArgSet( *(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct")) ),*redMCCutPR);
  RooDataHist *binMCCutNP = new RooDataHist("binMCCutNP","MC distribution for NP signal",RooArgSet(*(ws->var("Jpsi_CtTrue"))),*redMCCutNP);

  // *** Define signal, sideband datasets
  RooDataSet *redDataCutSIG = (RooDataSet*) redDataCut->reduce("Jpsi_Mass > 2.9 && Jpsi_Mass < 3.3");
  RooDataSet *redDataCutSB;
  if (narrowSideband) {
    redDataCutSB = (RooDataSet*) redDataCut->reduce("Jpsi_Mass < 2.8 || Jpsi_Mass > 3.4");
  } else {
    redDataCutSB = (RooDataSet*) redDataCut->reduce("Jpsi_Mass < 2.9 || Jpsi_Mass > 3.3");
  }
  RooDataHist *binDataCutSB = new RooDataHist("binDataCutSB","Data distribution for background",RooArgSet( *(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct")) ),*redDataCutSB);

  RooDataHist *binDataCutCtErrSB = new RooDataHist("binDataCutCtErrSB","Data ct error distribution for bkg",RooArgSet(*(ws->var("Jpsi_CtErr"))),*redDataCutSB);
  RooDataHist *binDataCutCtErrSIG = new RooDataHist("binDataCutCtErrSIG","Data ct error distribution for sig",RooArgSet(*(ws->var("Jpsi_CtErr"))),*redDataCutSIG);
  RooHistPdf errPdfBkg("errPdfBkg","Error PDF bkg",RooArgSet(*(ws->var("Jpsi_CtErr"))),*binDataCutCtErrSB);  ws->import(errPdfBkg);

  // Draw ctau error background pdf
  RooPlot *errframe2 = ws->var("Jpsi_CtErr")->frame();
  binDataCutCtErrSB->plotOn(errframe2,DataError(RooAbsData::SumW2));
  ws->pdf("errPdfBkg")->plotOn(errframe2,LineColor(kBlue),Normalization(binDataCutCtErrSB->sumEntries(),RooAbsReal::NumEvent));

  c0.Clear(); c0.cd(); c0.SetLogy(0); errframe2->Draw();
  titlestr = dirPre + "/rap" + yrange + "_cent" + crange + "_pT" + prange + "/testErrPdfBkg_Lin.pdf";
  c0.SaveAs(titlestr.c_str());
  c0.SetLogy(1); errframe2->Draw();
  titlestr = dirPre + "/rap" + yrange + "_cent" + crange + "_pT" + prange + "/testErrPdfBkg_Log.pdf";
  c0.SaveAs(titlestr.c_str());


  // *** Define PDFs with parameters (mass and ctau)
  // J/psi mass parameterization
  defineMassBkg(ws);
  defineMassSig(ws);
  // J/psi CTau parameterization
  defineCTResol(ws);              // R(l) : resolution function
  defineCTBkg(ws);                // theta(l') convolution R(l')
  titlestr = dirPre + "/rap" + yrange + "_cent" + crange + "_pT" + prange + "/testTrueLifeFit_Log.pdf";
  defineCTSig(ws,redMCCutNP,titlestr); // F_B(l) : R(l') convolution X_mc(l')
  RooProdPdf bkgCtTot_PEE( "bkgCtTot_PEE","PDF with PEE", *(ws->pdf("errPdfBkg")), Conditional( *(ws->pdf("bkgCtTot")), RooArgList(*(ws->var("Jpsi_Ct"))) ) );
  ws->import(bkgCtTot_PEE);
//  ws->factory("PROD::totBKG(polFunct,bkgCtTot)"); // F_bkg(l) : exp*theta(l')
  ws->factory("PROD::totBKG(expFunct,bkgCtTot)"); // F_bkg(l) : exp*theta(l')

  // Binning for invariant mass distribution
  RooBinning rbm(2.6,3.5);
  rbm.addUniform(45,2.6,3.5);

  string partTit, partFile;
  if (isGG == 0) { partTit = "glb-glb"; partFile = "GG"; }
  else if (isGG = 1) { partTit = "glb-trk"; partFile = "GT"; }
  else { partTit = "all"; partFile = "ALL"; }


  // Set some fitting variables to constant. It depends on the prefitting options.
  if (prefitMass) {
//    ws->var("enne")->setConstant(kTRUE);
//    ws->var("alpha")->setConstant(kTRUE);

/*    if (!isPT) {
      double inputNS[2] = {0};
      string inputFN;
      //_pt
      inputFN = "./mb/rap" + yrange + "_cent" + crange + "_pT6.5-30.0/2D_GG.txt";
      ifstream input(inputFN.c_str());
      if (!input.good()) {cout << "Fail to open MB mass fit results file." << endl; return 1;}
      string tmp;
      double inputP[6]={0}, inputTmp[2]={0};
      input >> tmp >> inputNS[0] >> inputNS[1]; //NSig
      input >> tmp >> inputTmp[0] >> inputTmp[1]; //NBkg
      for (int p=0; p<3; p++) {
        input >> tmp >> inputP[p]; cout << tmp << " " << inputP[p] << endl;
      }
      input >> tmp >> inputTmp[0] >> inputTmp[1]; //sigmaSig1
      for (int p=3; p<6; p++) {
        input >> tmp >> inputP[p]; cout << tmp << " " << inputP[p] << endl;
      }
//      ws->var("coefExp")->setVal(inputP[0]);
      ws->var("coeffGaus")->setVal(inputP[1]);
      ws->var("meanSig1")->setVal(inputP[2]);
      ws->var("sigmaSig1")->setVal(inputTmp[0]);
      ws->var("alpha")->setVal(inputP[4]);
      ws->var("enne")->setVal(inputP[5]);
      ws->var("alpha")->setConstant(kTRUE);
      ws->var("enne")->setConstant(kTRUE);
//      ws->var("coefExp")->setConstant(kTRUE);
      ws->var("coeffGaus")->setConstant(kTRUE);
      ws->var("meanSig1")->setConstant(kTRUE);
      ws->var("sigmaSig1")->setConstant(kTRUE);
      input.close();
    }

*/
/*    if (!isMB) {
      double inputNS[2] = {0};
      string inputFN;
      if (isMBPT) {
        //_mb_pt
        cout << "_mb_pt\n";
        inputFN =  dirPre + "_rap" + yrange + "_cent0-100_pT" + prange + "_2D_GG.txt";
      } else if (!isMBPT) {
        //_mb
        cout << "_mb\n";
        inputFN = dirPre + "_rap" + yrange + "_cent0-100_pT6.5-30.0_2D_GG.txt";
      }
      ifstream input(inputFN.c_str());
      if (!input.good()) {cout << "Fail to open MB mass fit results file." << endl; return 1;}
      string tmp;
      double inputP[8][2]={0}, inputTmp[2]={0};
      input >> tmp >> inputNS[0] >> inputNS[1]; //NSig
      input >> tmp >> inputTmp[0] >> inputTmp[1]; //NBkg
      for (int p=0; p<8; p++) {   //Mass signal parameters
        input >> tmp >> inputP[p][0] >> inputP[p][1];
        cout << tmp << " " << inputP[p][0] << endl;
      }

      ws->var("coeffGaus")->setVal(inputP[1][0]);
      ws->var("meanSig1")->setVal(inputP[2][0]);
      ws->var("sigmaSig1")->setVal(inputP[3][0]);
      ws->var("sigmaSig2")->setVal(inputP[4][0]);
      ws->var("alpha")->setVal(inputP[5][0]);
      ws->var("enne")->setVal(inputP[6][0]);
      ws->var("enneW")->setVal(inputP[7][0]);

      ws->var("alpha")->setConstant(kTRUE);
      ws->var("enne")->setConstant(kTRUE);
      ws->var("enneW")->setConstant(kTRUE);
      ws->var("coeffGaus")->setConstant(kTRUE);
      ws->var("meanSig1")->setConstant(kTRUE);
      ws->var("sigmaSig1")->setConstant(kTRUE);
      ws->var("sigmaSig2")->setConstant(kTRUE);
      input.close();
    }*/
    
    ws->factory("SUM::sigMassPDF(NSig[500.0,10.0,10000000.0]*sigCBWNG1,NBkg[2000.,10.,100000000.0]*expFunct)");
//    ws->factory("SUM::sigMassPDF(NSig[500.0,10.0,10000.0]*sigCBWNG1,NBkg[2000.,10.,100000.0]*polFunct)");
//    ws->factory("SUM::sigMassPDF(NSig[500.0,10.0,10000.0]*signalG1,NBkg[2000.,10.,100000.0]*expFunct)");
    RooFitResult *fitM = ws->pdf("sigMassPDF")->fitTo(*redDataCut,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(4));
//    resultF->cd();
//    fitM->Write("sigMassPDF");

    RooPlot *mframe_wob = ws->var("Jpsi_Mass")->frame();
    redDataCut->plotOn(mframe_wob,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(1),Binning(rbm));
    titlestr = "2D fit for" + partTit + "muons (mass projection), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
    mframe_wob->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
    mframe_wob->GetXaxis()->CenterTitle(1);
    const double max = mframe_wob->GetMaximum() * 1.3;
    mframe_wob->SetMaximum(max);

    ws->pdf("sigMassPDF")->plotOn(mframe_wob,DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
//    ws->pdf("sigMassPDF")->plotOn(mframe_wob,Components("polFunct"),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("sigMassPDF")->plotOn(mframe_wob,Components("expFunct"),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
//    ws->pdf("sigMassPDF")->plotOn(mframe_wob,Components("polFunct"),LineColor(kBlue),LineStyle(7),LineWidth(5),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("sigMassPDF")->plotOn(mframe_wob,Components("expFunct"),LineColor(kBlue),LineStyle(7),LineWidth(5),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("sigMassPDF")->plotOn(mframe_wob,LineColor(kBlack),LineWidth(2),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));

    redDataCut->plotOn(mframe_wob,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(1),Binning(rbm));
    TCanvas c00;
    c00.cd(); mframe_wob->Draw();
    titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_2D_" + partFile + "massfit_wob.pdf";
    c00.SaveAs(titlestr.c_str());

  } else {
    RooRealVar NSig("NSig","dummy total signal events",0.);
    ws->import(NSig);
  }

  Double_t NSig_fin = ws->var("NSig")->getVal();
  Double_t ErrNSig_fin = ws->var("NSig")->getError();

  float bc = ws->var("coefExp")->getVal();
  float scaleF = (exp(2.9*bc)-exp(3.3*bc))/(exp(2.6*bc)-exp(2.9*bc)+exp(3.3*bc)-exp(3.6*bc));
  RooDataHist* binSubtrData = subtractSidebands(ws,binDataCutCtErrSIG,binDataCutCtErrSB,scaleF,"Jpsi_CtErr");
  binSubtrData->SetName("binSubtrData");
  RooHistPdf errPdfSig("errPdfSig","Error PDF signal",RooArgSet(*(ws->var("Jpsi_CtErr"))),*binSubtrData);  ws->import(errPdfSig);

  if (prefitMass) {
    ws->var("alpha")->setConstant(kTRUE);
    ws->var("enne")->setConstant(kTRUE);
    ws->var("enneW")->setConstant(kTRUE);
    ws->var("coeffGaus")->setConstant(kTRUE);
    ws->var("sigmaSig1")->setConstant(kTRUE);
    ws->var("sigmaSig2")->setConstant(kTRUE);
    ws->var("meanSig1")->setConstant(kTRUE);
    ws->var("NSig")->setConstant(kTRUE);
    ws->var("NBkg")->setConstant(kTRUE);

    RooFormulaVar fBkg("fBkg","@0/(@0+@1)",RooArgList(*(ws->var("NBkg")),*(ws->var("NSig"))));   ws->import(fBkg);
    ws->factory("PROD::totSIGPR(sigCBWNG1,sigPR)"); 
    ws->factory("PROD::totSIGNP(sigCBWNG1,sigNP)");
    RooProdPdf totSIGPR_PEE( "totSIGPR_PEE","PDF with PEE", *(ws->pdf("errPdfSig")), Conditional( *(ws->pdf("totSIGPR")), RooArgList(*(ws->var("Jpsi_Ct")), *(ws->var("Jpsi_Mass"))) ) );    ws->import(totSIGPR_PEE);
    RooProdPdf totSIGNP_PEE( "totSIGNP_PEE","PDF with PEE", *(ws->pdf("errPdfSig")), Conditional( *(ws->pdf("totSIGNP")), RooArgList(*(ws->var("Jpsi_Ct")), *(ws->var("Jpsi_Mass"))) ) );    ws->import(totSIGNP_PEE);
    RooProdPdf totBKG_PEE( "totBKG_PEE","PDF with PEE", *(ws->pdf("errPdfBkg")), Conditional( *(ws->pdf("totBKG")), RooArgList(*(ws->var("Jpsi_Ct")), *(ws->var("Jpsi_Mass"))) ) );    ws->import(totBKG_PEE);
    ws->factory("RSUM::totPDF_PEE(fBkg*totBKG_PEE,Bfrac[0.2,0.,1.]*totSIGNP_PEE,totSIGPR_PEE)");  //Final F(l,m)

    RooPlot *errframe3 = ws->var("Jpsi_CtErr")->frame();
    binDataCutCtErr->plotOn(errframe3,DataError(RooAbsData::SumW2));
    ws->pdf("errPdfSig")->plotOn(errframe3,LineColor(kBlue),Normalization(binDataCutCtErr->sumEntries(),RooAbsReal::NumEvent));
    
    c0.Clear(); c0.cd(); c0.SetLogy(0); errframe3->Draw();
    titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_testErrPdfSig_Lin.pdf";
    c0.SaveAs(titlestr.c_str());
    c0.SetLogy(1); errframe3->Draw();
    titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_testErrPdfSig_Log.pdf";
    c0.SaveAs(titlestr.c_str());

  } else {
    ws->factory("PROD::totSigPR(sigCBG1,sigPR)");
    ws->factory("PROD::totSigNP(sigCBG1,sigNP)");
    ws->factory("SUM::totPDF(NSigPR[4000.0,10.0,1000000.0]*totSigPR,NSigNP[900.0,10.,1000000.]*totSigNP,NBkg[1400.,10.,1000000.]*totBKG)");   //Final F(l,m)
  }

  // *** Prefit ctau distribution of signal
  if (prefitSignalCTau) {
    RooProdPdf sigPR_PEE("sigPR_PEE","PDF with PEE", *(ws->pdf("errPdfSig")),Conditional(*(ws->pdf("sigPR")), RooArgList(*(ws->var("Jpsi_Ct")))));  ws->import(sigPR_PEE);
    RooFitResult *fitSigPrPee = ws->pdf("sigPR_PEE")->fitTo(*redMCCutPR,Range("promptfit"),SumW2Error(kTRUE),ConditionalObservables(RooArgSet(*(ws->var("Jpsi_CtErr")))),NumCPU(4));
//    resultF->cd();
//    fitSigPrPee->Write("sigPR_PEE");
    if (ws->var("sigmaResSigW")) ws->var("sigmaResSigW")->setConstant(kTRUE);
    ws->var("meanResSigW")->setConstant(kTRUE);

    RooPlot *errframe3 = ws->var("Jpsi_CtErr")->frame();
    redMCCutPR->plotOn(errframe3,DataError(RooAbsData::SumW2));
    ws->pdf("sigPR_PEE")->plotOn(errframe3,LineColor(kBlue),Normalization(redMCCutPR->sumEntries(),RooAbsReal::NumEvent));
    
    c0.Clear(); c0.cd(); c0.SetLogy(0); errframe3->Draw();
    titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_testErrPdfSigPRPee_Lin.pdf";
    c0.SaveAs(titlestr.c_str());
    c0.SetLogy(1); errframe3->Draw();
    titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_testErrPdfSigPRPee_Log.pdf";
    c0.SaveAs(titlestr.c_str());


    RooRealVar *CtWeighted = new RooRealVar("CtWeighted","#font[12]{l}_{J/#psi} / #sigma( #font[12]{l}_{J/#psi} )",-5.,5.);
    ws->import(*CtWeighted);

    const RooArgSet *thisRow = (RooArgSet*)redMCCutPR->get(0); 
    RooArgSet *newRow = new RooArgSet(*CtWeighted);
    RooDataSet *tempJpsi = new RooDataSet("tempJpsi","ctau weighted data",*newRow);

    for (Int_t iSamp = 0; iSamp < redMCCutPR->numEntries(); iSamp++) {
      thisRow = (RooArgSet*)redMCCutPR->get(iSamp);
      RooRealVar *myct = (RooRealVar*)thisRow->find("Jpsi_Ct");
      RooRealVar *mycterr = (RooRealVar*)thisRow->find("Jpsi_CtErr");
      CtWeighted->setVal(myct->getVal()/mycterr->getVal());
      RooArgSet *tempRow = new RooArgSet(*CtWeighted);
      tempJpsi->add(*tempRow);
    }

    if (oneGaussianResol) {
      ws->factory("Gaussian::tempsigPR(CtWeighted,meanResSigW,sigmaResSigN)");
    } else {
      ws->factory("Gaussian::tempresGW(CtWeighted,meanResSigW,sigmaResSigW)");
      ws->factory("Gaussian::tempresGN(CtWeighted,meanResSigW,sigmaResSigN)");
      ws->factory("SUM::tempsigPR(fracRes*tempresGW,tempresGN)");
    }  

    titlestr = "Prompt resolution fit for" + partTit + "muons (ctau projection), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
    RooPlot *tframePR = ws->var("CtWeighted")->frame();
    tframePR->SetTitle(titlestr.c_str());
    tframePR->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
    tframePR->GetXaxis()->SetRangeUser(-1.5,-3);
    tempJpsi->plotOn(tframePR,DataError(RooAbsData::SumW2),Binning(rb2));
    ws->pdf("tempsigPR")->plotOn(tframePR,LineColor(kBlue),Normalization(tempJpsi->sumEntries(),RooAbsReal::NumEvent));

    c0.Clear();   c0.cd();   c0.SetLogy(0);   tframePR->Draw();
    titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_2D_" + partFile + "resolfit_Lin.pdf";
    c0.SaveAs(titlestr.c_str());
    c0.cd();   c0.SetLogy(1);   tframePR->Draw();
    titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_2D_" + partFile + "resolfit_Log.pdf";
    c0.SaveAs(titlestr.c_str());

    /* // Get sum of squared residual to fisher's F-test
    RooHist *hresid;
    TH1 *hres = redMCCutPR->createHistogram("hres",*ws->var("CtWeighted"),Binning(rb2));
    unsigned int nBins = hres->GetNbinsX();
    unsigned int nFullBinsResid = 0;
    hresid = tframePR->residHist();
    hresid->SetName("hresid");
    double *yresid = hresid->GetY();
    for (unsigned int i=0; i < nBins; i++) {
      if (hres->GetBinContent(i) == 0) continue;
      nFullBinsResid++;
      cout << "Residual of bin " << i << " = " << yresid[i] << endl;
      RSS = RSS + pow(yresid[i],2);
    }
    cout << "Residual sum of squares : " << RSS << endl;
    cout << "nBins: " << nBins << endl;
    cout << "nFullBinsResid: " << nFullBinsResid << endl;
    // End of fisher's F-test */
  
    if (ws->var("sigmaResSigW")) ws->var("sigmaResSigW")->setConstant(kFALSE);
    ws->var("meanResSigW")->setConstant(kFALSE);

  }

  // *** Prefit ctau distribution of sideband
  if(prefitBkg){
    cout << "Prefitting background on " << redDataCutSB->sumEntries() << " events " << endl;

    ws->var("fpm")->setConstant(kTRUE);
    if(prefitSignalCTau){
      if (ws->var("fracRes")) ws->var("fracRes")->setConstant(kTRUE);
      ws->var("meanResSigW")->setConstant(kTRUE);
      if (ws->var("sigmaResBkgW")) ws->var("sigmaResBkgW")->setVal(ws->var("sigmaResSigW")->getVal());
      if (ws->var("sigmaResBkgN")) ws->var("sigmaResBkgN")->setVal(ws->var("sigmaResSigN")->getVal());
    }

    cout << "DATA :: N events to fit on the sidebands: " << redDataCutSB->sumEntries() << endl;
    RooFitResult *fitBkfCtTotPee= ws->pdf("bkgCtTot_PEE")->fitTo(*redDataCutSB,SumW2Error(kTRUE),Minos(0),NumCPU(4),ConditionalObservables(RooArgSet(*(ws->var("Jpsi_CtErr")))),NumCPU(4));
//    resultF->cd();
//    fitBkfCtTotPee->Write("bkgCtTot_PEE");
    ws->var("fpm")->setConstant(kTRUE);
    ws->var("fLiving")->setConstant(kTRUE);
    ws->var("fbkgCtTot")->setConstant(kTRUE);
    // ws->var("fracResBkg")->setConstant(kTRUE);
    if (ws->var("sigmaResBkgN")) ws->var("sigmaResBkgN")->setConstant(kTRUE);
    if (ws->var("sigmaResBkgW")) ws->var("sigmaResBkgW")->setConstant(kTRUE);
    if (ws->var("meanResBkgW")) ws->var("meanResBkgW")->setConstant(kTRUE);
    ws->var("lambdap")->setConstant(kTRUE);
    ws->var("lambdam")->setConstant(kTRUE);
    ws->var("lambdasym")->setConstant(kTRUE);

    if(prefitSignalCTau){
      if (ws->var("fracRes")) ws->var("fracRes")->setConstant(kFALSE);
      //ws->var("meanResSigN")->setConstant(kFALSE);
      ws->var("meanResSigW")->setConstant(kFALSE);
      //ws->var("sigmaResSigN")->setConstant(kFALSE);
      //ws->var("sigmaResSigW")->setConstant(kFALSE);
    }

    RooPlot *tframe1 = ws->var("Jpsi_Ct")->frame();
    titlestr = "2D fit for" + partTit + "muons (J/ #psi c  #tau projection, mass sidebands), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
    tframe1->SetTitle(titlestr.c_str());
    redDataCutSB->plotOn(tframe1,DataError(RooAbsData::SumW2),Binning(rb2));
    ws->pdf("bkgCtTot_PEE")->plotOn(tframe1,ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCutCtErr,kTRUE),NumCPU(4),Normalization(redDataCutSB->sumEntries(),RooAbsReal::NumEvent));
    
    c0.Clear();   c0.cd();    c0.SetLogy(0);  tframe1->Draw();
    titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_2D_" + partFile + "timeside_Lin.pdf";
    c0.SaveAs(titlestr.c_str());
    c0.cd();    c0.SetLogy(1);    tframe1->Draw();
    titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_2D_" + partFile + "timeside_Log.pdf";
    c0.SaveAs(titlestr.c_str()); 
    
  }


  // Fix below bkg variables in any case
  ws->var("fpm")->setConstant(kTRUE);
  ws->var("fLiving")->setConstant(kTRUE);

  // *** Get NSig, NBkg, Bfraction and their errors
  Double_t NBkg_fin = ws->var("NBkg")->getVal();
  Double_t ErrNBkg_fin = ws->var("NBkg")->getError();
  Double_t NSigPR_fin, ErrNSigPR_fin;
  Double_t NSigNP_fin, ErrNSigNP_fin;
  Double_t Bfrac_fin, ErrBfrac_fin;
  int nFitPar;
  Double_t theNLL;

  if (prefitMass) {
    RooFitResult *rfr;
    if (redDataCut->sumEntries() < 5000) {
      rfr = ws->pdf("totPDF_PEE")->fitTo(*redDataCut,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(4),ConditionalObservables(RooArgSet(*(ws->var("Jpsi_CtErr")))));
//      resultF->cd();
//      rfr->Write("totPDF_PEE");
    } else {
      rfr = ws->pdf("totPDF_PEE")->fitTo(*binDataCut,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(4),ConditionalObservables(RooArgSet(*(ws->var("Jpsi_CtErr")))));
//      resultF->cd();
//      rfr->Write("totPDF_PEE");
    }
    nFitPar = rfr->floatParsFinal().getSize();
    theNLL = rfr->minNll();
    Bfrac_fin = ws->var("Bfrac")->getVal();
    ErrBfrac_fin = ws->var("Bfrac")->getError();
    // ws->var("sigmaResSigN")->setConstant(kTRUE);
    // ws->var("sigmaResSigW")->setConstant(kTRUE);
     if (ws->var("fracRes")) ws->var("fracRes")->setConstant(kTRUE);
    ws->var("meanResSigW")->setConstant(kTRUE);
    // ws->var("Bfrac")->setConstant(kTRUE);

    NSigNP_fin = NSig_fin*Bfrac_fin;
    NSigPR_fin = NSig_fin*(1.0-Bfrac_fin);
    ErrNSigNP_fin = NSigNP_fin*sqrt(pow(ErrNSig_fin/NSig_fin,2) + pow(ErrBfrac_fin/Bfrac_fin,2));
    ErrNSigPR_fin = NSigPR_fin*sqrt(pow(ErrNSig_fin/NSig_fin,2) + pow(ErrBfrac_fin/(1.0-Bfrac_fin),2));

  } else {
    RooFitResult *rfr = ws->pdf("totSim")->fitTo(*redDataCut,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(4));
//    resultF->cd();
//    rfr->Write("totSim");
    nFitPar = rfr->floatParsFinal().getSize();   
    
    NSigNP_fin = ws->var("NSigNP")->getVal();
    NSigPR_fin = ws->var("NSigPR")->getVal();
    ErrNSigNP_fin = ws->var("NSigNP")->getError();
    ErrNSigPR_fin = ws->var("NSigPR")->getError();

    Bfrac_fin = NSigNP_fin/(NSigNP_fin + NSigPR_fin);
    ErrBfrac_fin = sqrt(pow(NSigNP_fin*ErrNSigPR_fin,2) + pow(NSigPR_fin*ErrNSigNP_fin,2))/pow(NSigNP_fin + NSigPR_fin,2);
  }

  // *** Plot various fit results and data points
  // Temporary variables for plotting
  RooRealVar tmpVar1("tmpVar1","tmpVar1",NSigNP_fin);
  RooRealVar tmpVar2("tmpVar2","tmpVar2",NBkg_fin);

  // Mass plot with bfraction
  RooPlot *mframe = ws->var("Jpsi_Mass")->frame();
  redDataCut->plotOn(mframe,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(1),Binning(rbm));
  titlestr = "2D fit for" + partTit + "muons (mass projection), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
  mframe->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  mframe->GetXaxis()->CenterTitle(1);
//  mframe->GetYaxis()->SetTitle("Events/(0.015 GeV/c^{2})");
  const double max = mframe->GetMaximum() * 1.3;
  mframe->SetMaximum(max);

  if (prefitMass) {
    // Fill color
    ws->pdf("totPDF_PEE")->plotOn(mframe,DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
//    RooAddPdf tmpPDF("tmpPDF","tmpPDF",RooArgList(*(ws->pdf("signalG1")),*(ws->pdf("expFunct"))),RooArgList(tmpVar1,tmpVar2));
//    RooAddPdf tmpPDF("tmpPDF","tmpPDF",RooArgList(*(ws->pdf("sigCBWNG1")),*(ws->pdf("polFunct"))),RooArgList(tmpVar1,tmpVar2));
    RooAddPdf tmpPDF("tmpPDF","tmpPDF",RooArgList(*(ws->pdf("sigCBWNG1")),*(ws->pdf("expFunct"))),RooArgList(tmpVar1,tmpVar2));
    tmpPDF.plotOn(mframe,LineColor(kRed),DrawOption("F"),FillColor(kWhite),FillStyle(1001),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    tmpPDF.plotOn(mframe,LineColor(kRed),DrawOption("F"),FillColor(kRed),FillStyle(3444),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    gStyle->SetHatchesLineWidth(2);
//    ws->pdf("totPDF_PEE")->plotOn(mframe,Components("polFunct"),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("totPDF_PEE")->plotOn(mframe,Components("expFunct"),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));

    //Line color
//    ws->pdf("totPDF_PEE")->plotOn(mframe,Components("polFunct"),LineColor(kBlue),LineStyle(7),LineWidth(5),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("totPDF_PEE")->plotOn(mframe,Components("expFunct"),LineColor(kBlue),LineStyle(7),LineWidth(5),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    tmpPDF.plotOn(mframe,LineColor(kRed),LineStyle(9),LineWidth(5),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    ws->pdf("totPDF_PEE")->plotOn(mframe,LineColor(kBlack),LineWidth(2),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
  } else {
    ws->pdf("totPDF_PEE")->plotOn(mframe,Components("totSigNP,totBKG"),LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF_PEE")->plotOn(mframe,Components("totBKG"),LineColor(kBlue),LineColor(7),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF_PEE")->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  }

  redDataCut->plotOn(mframe,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(1),Binning(rbm));

  TLatex *t = new TLatex();
  t->SetNDC(); t->SetTextAlign(12);
/*  t->SetTextSize(0.05);
  t->DrawLatex(0.17,0.90,"CMS PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  t->SetTextSize(0.04);
  t->DrawLatex(0.17,0.83,"L_{int} =  70 #mub^{-1}");
  sprintf(reduceDS,"%.0f-%.0f%, |y| < %.1f",cmin,cmax,ymax);
  t->DrawLatex(0.17,0.76,reduceDS);
  sprintf(reduceDS,"%.1f < p_{T} < %.1f GeV/c",pmin,pmax);
  t->DrawLatex(0.17,0.70,reduceDS);
  sprintf(reduceDS,"#sigma = %0.0f MeV/c^{2}",ws->var("sigmaSig1")->getVal()*1000);
  t->DrawLatex(0.65,0.715,reduceDS);
//  sprintf(reduceDS,"#sigma = %0.0f #pm %0.0f MeV/c^{2}",ws->var("sigmaSig1")->getVal()*1000,ws->var("sigmaSig1")->getError()*1000);
*/
  Double_t fx[2], fy[2], fex[2], fey[2];
  TGraphErrors *gfake1 = new TGraphErrors(2,fx,fy,fex,fey);
  gfake1->SetMarkerStyle(20);
  gfake1->SetMarkerSize(1);
  TH1F hfake11 = TH1F("hfake1","hfake1",100,200,300);
  hfake11.SetLineColor(kBlue); hfake11.SetLineWidth(4); hfake11.SetLineStyle(7); hfake11.SetFillColor(kAzure-9); hfake11.SetFillStyle(1001);
  TH1F hfake21 = TH1F("hfake2","hfake2",100,200,300);
  hfake21.SetLineColor(kBlack); hfake21.SetLineWidth(4); hfake21.SetFillColor(kBlack); hfake21.SetFillStyle(3354);
  TH1F hfake31 = TH1F("hfake3","hfake3",100,200,300);
  hfake31.SetLineColor(kRed); hfake31.SetMarkerStyle(kCircle); hfake31.SetLineWidth(4); hfake31.SetMarkerColor(kRed); hfake31.SetLineStyle(9); hfake31.SetFillColor(kRed-7); hfake31.SetFillStyle(3444);

/*  TLegend * leg1 = new TLegend(0.16,0.50,0.57,0.65,NULL,"brNDC");
  leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetShadowColor(0);
  leg1->SetMargin(0.2);
  leg1->AddEntry(gfake1,"data","p");
  leg1->AddEntry(&hfake21,"total fit","lf");
  leg1->AddEntry(&hfake11,"background","lf");
  leg1->Draw("same");
  titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_2D_" + partFile + "massfit_wob.pdf";
  c00.SaveAs(titlestr.c_str());*/

  TCanvas c1;
  c1.cd(); mframe->Draw();
  t->SetTextSize(0.05);

  t->DrawLatex(0.17,0.90,"CMS PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  t->SetTextSize(0.04);
  t->DrawLatex(0.17,0.83,"L_{int} =  70 #mub^{-1}");
  sprintf(reduceDS,"%.0f-%.0f%, |y| < %.1f",cmin,cmax,ymax);
  t->DrawLatex(0.17,0.76,reduceDS);
  sprintf(reduceDS,"%.1f < p_{T} < %.1f GeV/c",pmin,pmax);
  t->DrawLatex(0.17,0.70,reduceDS);
  sprintf(reduceDS,"#sigma = %0.0f MeV/c^{2}",ws->var("sigmaSig1")->getVal()*1000);
  t->DrawLatex(0.65,0.715,reduceDS);


  TLegend * leg11 = new TLegend(0.16,0.48,0.57,0.66,NULL,"brNDC");
  leg11->SetFillStyle(0); leg11->SetBorderSize(0); leg11->SetShadowColor(0);
  leg11->SetMargin(0.2);
  leg11->AddEntry(gfake1,"data","p");
  leg11->AddEntry(&hfake21,"total fit","lf");
  leg11->AddEntry(&hfake31,"bkgd + non-prompt","lf"); 
  leg11->AddEntry(&hfake11,"background","lf");
  leg11->Draw("same");
  titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_2D_" + partFile + "massfit.pdf";
  c1.SaveAs(titlestr.c_str());

  ws->var("Jpsi_Ct")->setBinning(rb3);
  RooPlot *tframe = ws->var("Jpsi_Ct")->frame();
  titlestr = "2D fit for" + partTit + "muons (c#tau projection), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
  tframe->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  tframe->GetXaxis()->CenterTitle(1);
  tframe->GetYaxis()->SetTitle("Events / (0.088 mm)");
//  tframe->SetTitleOffset(1.35,"Y");

  RooHist *hpull;
  redDataCut->plotOn(tframe,DataError(RooAbsData::SumW2),Binning(rb3),MarkerSize(1));
  if (prefitMass) {
    ws->pdf("totPDF_PEE")->plotOn(tframe,ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCutCtErr,kTRUE),DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
//    ws->pdf("totPDF_PEE")->plotOn(tframe,DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    hpull = tframe->pullHist();
    hpull->SetName("hpull");
    RooAddPdf tmpPDF2("tmpPDF2","tmpPDF2",RooArgList(*(ws->pdf("totSIGNP_PEE")),*(ws->pdf("bkgCtTot_PEE"))),RooArgList(tmpVar1,tmpVar2));
    tmpPDF2.plotOn(tframe,DrawOption("F"),FillColor(kWhite),FillStyle(1001),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    tmpPDF2.plotOn(tframe,DrawOption("F"),FillColor(kRed),FillStyle(3444),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    gStyle->SetHatchesLineWidth(2);
    ws->pdf("totPDF_PEE")->plotOn(tframe,Components("totBKG_PEE"),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCutCtErr,kTRUE),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("totPDF_PEE")->plotOn(tframe,LineColor(kBlack),LineWidth(2),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("totPDF_PEE")->plotOn(tframe,Components("totBKG_PEE"),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCutCtErr,kTRUE),LineColor(kBlue),LineWidth(5),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent),LineStyle(7));
    tmpPDF2.plotOn(tframe,LineColor(kRed),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent),LineWidth(5),LineStyle(9));
  } else {
    ws->pdf("totPDF")->plotOn(tframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF")->plotOn(tframe,Components("totSigNP,totBKG"),LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected),LineStyle(9));
    ws->pdf("totPDF")->plotOn(tframe,Components("totBKG"),LineColor(kBlue),Normalization(1.0,RooAbsReal::RelativeExpected),LineStyle(7));
    ws->pdf("totPDF")->plotOn(tframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  }
  redDataCut->plotOn(tframe,DataError(RooAbsData::SumW2),Binning(rb3),MarkerSize(1));

  tframe->SetMaximum(tframe->GetMaximum()*9); 
  tframe->SetMinimum(0.5); 
  tframe->Draw();

  ////////// Line only drawing is commented out to speed up
/*  redDataCut->plotOn(tframe,DataError(RooAbsData::SumW2),Binning(rb3),MarkerSize(1));
  if (prefitMass) {
   ws->pdf("totPDF_PEE")->plotOn(tframe,ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCutCtErr,kTRUE),LineColor(kBlack),LineWidth(2),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    hpull = tframe->pullHist();
    hpull->SetName("hpull");
    ws->pdf("totPDF_PEE")->plotOn(tframe,Components("totBKG_PEE"),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCutCtErr,kTRUE),LineColor(kBlue),LineWidth(5),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent),LineStyle(7));
    RooAddPdf tmpPDF2("tmpPDF2","tmpPDF2",RooArgList(*(ws->pdf("totSIGNP_PEE")),*(ws->pdf("bkgCtTot_PEE"))),RooArgList(tmpVar1,tmpVar2));
    tmpPDF2.plotOn(tframe,LineColor(kRed),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent),LineWidth(5),LineStyle(9));
//    ws->pdf("totPDF_PEE")->plotOn(tframe,LineColor(kBlack),LineWidth(2),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("totPDF_PEE")->plotOn(tframe,LineColor(kBlack),LineWidth(2),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCutCtErr,kTRUE),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));  
  } else {
    ws->pdf("totPDF")->plotOn(tframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
    hpull = tframe->pullHist();
    hpull->SetName("hpull");
    ws->pdf("totPDF")->plotOn(tframe,Components("totSigNP,totBKG"),LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected),LineStyle(9));
    ws->pdf("totPDF")->plotOn(tframe,Components("totBKG"),LineColor(kBlue),Normalization(1.0,RooAbsReal::RelativeExpected),LineStyle(7));
    ws->pdf("totPDF")->plotOn(tframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  }*/

  double chi2 = 0, unNormChi2 = 0;
  int dof = 0;
  double *ypulls = hpull->GetY();
  unsigned int nFullBins = 0;
  int nBins = ws->var("Jpsi_Ct")->getBinning().numBins();
  TH1 *hpullhist = redDataCut->createHistogram("hpullhist",*ws->var("Jpsi_Ct"),Binning(ws->var("Jpsi_Ct")->getBinning()));
  
  for (unsigned int i = 0; i < nBins; i++) {
    if (hpullhist->GetBinContent(i+1) == 0) continue;
    nFullBins++;
//    cout << "Pull of bin " << i << " = " << ypulls[i] << endl;
    chi2 += ypulls[i]*ypulls[i];
  }
  unNormChi2 = chi2;
  dof = nFullBins - nFitPar;
  chi2 /= (nFullBins - nFitPar);
  for (unsigned int i = 0; i < nBins; i++) {
    if (fabs(ypulls[i]) < 0.0001) ypulls[i] = 999.; 
    hpull->SetPointError(i,0.,0.,0.,0.);
  } 


  // WITH RESIDUALS
  TCanvas* c2 = new TCanvas("c2","The Canvas",200,10,600,880);
  c2->cd();
  TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.03,0.95,0.35);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.24);
  pad2->Draw();

  pad1->cd(); /* pad1->SetLogy(1); */ tframe->Draw();

  t->SetTextSize(0.05);
  t->DrawLatex(0.44,0.9,"CMS PbPb  #sqrt{s_{NN}} = 2.76 TeV"); 
  sprintf(reduceDS,"%.0f-%.0f%, |y| < %.1f",cmin,cmax,ymax);
  t->SetTextSize(0.04);
  t->DrawLatex(0.47,0.82,reduceDS);
  sprintf(reduceDS,"%.1f < p_{T} < %.1f GeV/c",pmin,pmax);
  t->DrawLatex(0.47,0.73,reduceDS);
  t->DrawLatex(0.47,0.66,"L_{int} =  70 #mub^{-1}"); 

  TLegend * leg = new TLegend(0.455,0.35,0.85,0.55,NULL,"brNDC");
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetShadowColor(0);
  leg->SetMargin(0.2);
  leg->AddEntry(gfake1,"data","p");
  leg->AddEntry(&hfake21,"total fit","lf");
  leg->AddEntry(&hfake31,"bkgd + non-prompt","lf"); 
  leg->AddEntry(&hfake11,"background","lf");
  leg->Draw("same"); 

  RooPlot* tframepull =  ws->var("Jpsi_Ct")->frame(Title("Pull Distribution")) ;
  tframepull->GetYaxis()->SetTitle("Pull");
  tframepull->SetLabelSize(0.08,"XYZ");
  tframepull->SetTitleSize(0.1,"XYZ");
  tframepull->SetTitleOffset(0.55,"Y");
//  tframepull->SetTitleOffset(0.95,"X");
  tframepull->addPlotable(hpull,"P") ;
  tframepull->SetMaximum(-(tframepull->GetMinimum())); 
  tframepull->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  tframepull->GetXaxis()->CenterTitle(1);

  pad2->cd(); tframepull->Draw();


  TLatex *t2 = new TLatex();
  t2->SetNDC(); t2->SetTextAlign(22);
  t2->SetTextSize(0.07);
  sprintf(reduceDS,"#chi^{2}/dof = %.2f/%d",unNormChi2,dof);
  t2->DrawLatex(0.76,0.90,reduceDS);
  
  c2->Update();

  titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_2D_" + partFile + "timefit_Lin.pdf";
  c2->SaveAs(titlestr.c_str());

  TCanvas* c2a = new TCanvas("c2a","The Canvas",200,10,600,880);
  c2a->cd();
  TPad *pad1a = new TPad("pad1a","This is pad1",0.05,0.35,0.95,0.97);
  pad1a->SetBottomMargin(0);
  pad1a->Draw();
  TPad *pad2a = new TPad("pad2a","This is pad2",0.05,0.03,0.95,0.35);
  pad2a->SetTopMargin(0);
  pad2a->SetBottomMargin(0.24);
  pad2a->Draw();

  pad1a->cd(); pad1a->SetLogy(1);
  tframe->SetMaximum(tframe->GetMaximum()*9); 
  tframe->SetMinimum(0.5); 
  tframe->Draw();
  t->SetTextSize(0.05);
  t->DrawLatex(0.17,0.9,"CMS PbPb  #sqrt{s_{NN}} = 2.76 TeV"); 
  t->SetTextSize(0.04);
  t->DrawLatex(0.17,0.819,"L_{int} =  70 #mub^{-1}"); 
  sprintf(reduceDS,"%.0f-%.0f%, |y| < %.1f",cmin,cmax,ymax);
  t->DrawLatex(0.5,0.819,reduceDS);
  sprintf(reduceDS,"%.1f < p_{T} < %.1f GeV/c",pmin,pmax);
  t->DrawLatex(0.5,0.759,reduceDS);

  leg->SetX1NDC(0.49);
  leg->SetY1NDC(0.55);
  leg->SetX2NDC(0.92);
  leg->SetY2NDC(0.72);
  leg->Draw("same");

  pad2a->cd(); tframepull->Draw();

  sprintf(reduceDS,"#chi^{2}/dof = %.2f/%d",unNormChi2,dof);
  t2->DrawLatex(0.76,0.90,reduceDS);
  
  c2a->Update();
  titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_2D_" + partFile + "timefit_Log.pdf";
  c2a->SaveAs(titlestr.c_str());

  TCanvas* c2b = new TCanvas("c2b","The Canvas",200,10,540,546);
  c2b->cd(); c2b->Draw(); c2b->SetLogy(1);

/*
  ws->var("Jpsi_Ct")->setBinning(rb3);
  RooPlot *tframefill = ws->var("Jpsi_Ct")->frame();
  tframefill->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  tframefill->GetXaxis()->CenterTitle(1);
  tframefill->GetYaxis()->SetTitle("Events / (0.088 mm)");

  redDataCut->plotOn(tframefill,DataError(RooAbsData::SumW2),Binning(rb3),MarkerSize(1));
  if (prefitMass) {
    ws->pdf("totPDF_PEE")->plotOn(tframefill,ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCutCtErr,kTRUE),DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
//    ws->pdf("totPDF_PEE")->plotOn(tframefill,DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    RooAddPdf tmpPDF2("tmpPDF2","tmpPDF2",RooArgList(*(ws->pdf("totSIGNP_PEE")),*(ws->pdf("bkgCtTot_PEE"))),RooArgList(tmpVar1,tmpVar2));
    tmpPDF2.plotOn(tframefill,DrawOption("F"),FillColor(kWhite),FillStyle(1001),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    tmpPDF2.plotOn(tframefill,DrawOption("F"),FillColor(kRed),FillStyle(3444),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    gStyle->SetHatchesLineWidth(2);
    ws->pdf("totPDF_PEE")->plotOn(tframefill,Components("totBKG_PEE"),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCutCtErr,kTRUE),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("totPDF_PEE")->plotOn(tframefill,LineColor(kBlack),LineWidth(2),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("totPDF_PEE")->plotOn(tframefill,Components("totBKG_PEE"),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCutCtErr,kTRUE),LineColor(kBlue),LineWidth(5),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent),LineStyle(7));
    tmpPDF2.plotOn(tframefill,LineColor(kRed),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent),LineWidth(5),LineStyle(9));
  } else {
    ws->pdf("totPDF")->plotOn(tframefill,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF")->plotOn(tframefill,Components("totSigNP,totBKG"),LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected),LineStyle(9));
    ws->pdf("totPDF")->plotOn(tframefill,Components("totBKG"),LineColor(kBlue),Normalization(1.0,RooAbsReal::RelativeExpected),LineStyle(7));
    ws->pdf("totPDF")->plotOn(tframefill,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  }
  redDataCut->plotOn(tframefill,DataError(RooAbsData::SumW2),Binning(rb3),MarkerSize(1));

  tframefill->SetMaximum(tframefill->GetMaximum()*9); 
  tframefill->SetMinimum(0.5); 
  tframefill->Draw();*/

  tframe->Draw();
  t->SetTextSize(0.05);
  t->DrawLatex(0.17,0.9,"CMS PbPb  #sqrt{s_{NN}} = 2.76 TeV"); 
  t->SetTextSize(0.04);
  t->DrawLatex(0.17,0.819,"L_{int} =  70 #mub^{-1}"); 
  sprintf(reduceDS,"%.0f-%.0f%, |y| < %.1f",cmin,cmax,ymax);
  t->DrawLatex(0.5,0.819,reduceDS);
  sprintf(reduceDS,"%.1f < p_{T} < %.1f GeV/c",pmin,pmax);
  t->DrawLatex(0.5,0.759,reduceDS);
  leg->SetX1NDC(0.49);
  leg->SetY1NDC(0.55);
  leg->SetX2NDC(0.92);
  leg->SetY2NDC(0.72);
  leg->Draw("same");

  titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_2D_" + partFile + "timefit_Log_wopull.pdf";
  c2b->SaveAs(titlestr.c_str());

/*  RooRealVar jpsiM = ws->var("Jpsi_Mass");
  RooRealVar jpsiCt = ws->var("Jpsi_Ct");
  RooDataHist hist2D("2DtotPDF","2D totPDF with mass and ctau",RooArgSet(jpsiM,jpsiCt),redDataCut);
  TH2F *htotPDF = (TH2F*)hist2D.createHistogram();
  TCanvas* c2c = new TCanvas("c2c","The Canvas",200,10,540,546);
  c2c->cd();
  htotPDF->Draw("COLZ");
  titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_2D_" + partFile + "2D_totFit.pdf";
  c2c->SaveAs(titlestr.c_str());
*/


  // To check values of fit parameters
  cout << endl << "J/psi yields:" << endl;
  cout << "PROMPT :     Fit : " << NSigPR_fin << " +/- " << ErrNSigPR_fin << endl;
  cout << "NON-PROMPT : Fit : " << NSigNP_fin << " +/- " << ErrNSigNP_fin << endl;
  cout << "Bfraction : Fit : " << Bfrac_fin << " +/- " << ErrBfrac_fin << endl;
//  cout << "Resolution : Fit : " << resol*1000. << " +/- " << Errresol*1000. << " mum" << endl;

  titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_2D_" + partFile +".txt";

  ofstream outputFile(titlestr.c_str());
  if (!outputFile.good()) {cout << "Fail to open result file." << endl; return 1;}
  outputFile
  << "NSig "         << NSig_fin                          << " " << ErrNSig_fin << "\n"
  << "NBkg "         << NBkg_fin                          << " " << ErrNBkg_fin << "\n"
  << "coefExp "      << ws->var("coefExp")->getVal()      << " " << ws->var("coefExp")->getError() << "\n"
  << "coeffGaus "    << ws->var("coeffGaus")->getVal()    << " " << ws->var("coeffGaus")->getError() << "\n"
  << "meanSig1 "     << ws->var("meanSig1")->getVal()     << " " << ws->var("meanSig1")->getError() << "\n"
  << "sigmaSig1 "    << ws->var("sigmaSig1")->getVal()    << " " << ws->var("sigmaSig1")->getError() << "\n"
  << "sigmaSig2 "    << ws->var("sigmaSig2")->getVal()    << " " << ws->var("sigmaSig2")->getError()<< "\n"
  << "alpha "        << ws->var("alpha")->getVal()        << " " << ws->var("alpha")->getError() << "\n"
  << "enne "         << ws->var("enne")->getVal()         << " " << ws->var("enne")->getError() << "\n"
  << "enneW "        << ws->var("enneW")->getVal()        << " " << ws->var("enneW")->getError() << "\n"
  << "fracRes "      << ws->var("fracRes")->getVal()      << " " << ws->var("fracRes")->getError() <<  "\n"
//  << "fracRes2 "     << ws->var("fracRes2")->getVal()     << " " << ws->var("fracRes2")->getError() << "\n"
//  << "fracRes3 "     << ws->var("fracRes3")->getVal()     << " " << ws->var("fracRes3")->getError() << "\n"
  << "meanResSigW "  << ws->var("meanResSigW")->getVal()  << " " << ws->var("meanResSigW")->getError() << "\n"
  << "sigmaResSigW " << ws->var("sigmaResSigW")->getVal() << " " << ws->var("sigmaResSigW")->getError() <<  "\n"
  << "sigmaResSigN " << ws->var("sigmaResSigN")->getVal() << " " << ws->var("sigmaResSigN")->getError() << "\n"
//  << "sigmaResSigM " << ws->var("sigmaResSigM")->getVal() << " " << ws->var("sigmaResSigM")->getError() << "\n"
//  << "sigmaResSigO " << ws->var("sigmaResSigO")->getVal() << " " << ws->var("sigmaResSigO")->getError() << "\n"
  << "fLiving "      << ws->var("fLiving")->getVal()      << " " << ws->var("fLiving")->getError() << "\n"
  << "fpm "          << ws->var("fpm")->getVal()          << " " << ws->var("fpm")->getError() << "\n"
  << "fbkgCtTot "    << ws->var("fbkgCtTot")->getVal()    << " " << ws->var("fbkgCtTot")->getError() << "\n"
  << "lambdam "      << ws->var("lambdam")->getVal()      << " " << ws->var("lambdam")->getError() << "\n"
  << "lambdap "      << ws->var("lambdap")->getVal()      << " " << ws->var("lambdap")->getError() << "\n"
  << "lambdasym "    << ws->var("lambdasym")->getVal()    << " " << ws->var("lambdasym")->getError() << "\n"
  << "NLL "          << theNLL                            << endl
  << "PROMPT "       << NSigPR_fin                        << " " << ErrNSigPR_fin << endl
  << "NON-PROMPT "   << NSigNP_fin                        << " " << ErrNSigNP_fin << endl
  << "Bfraction "    << Bfrac_fin                         << " " << ErrBfrac_fin << endl;
//  << "Resolution "   << resol*1000.                       << " " << Errresol*1000. << endl
//  << "nFullBinsResid "<< nFullBinsResid                   << endl
//  << "RSS "          << RSS                               << endl;

  outputFile.close();
  
/*  titlestr = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_workspace.root";
  if (ws->importClassCode("RooHistPdfConv")) cout << "importClassCode() succeed" << endl;
  else cout << "importClassCode() failed" << endl;
  if (ws->writeToFile(titlestr.c_str(),kTRUE)) cout << "Workspace saved" << endl;
  else cout << "Workspace NOT saved" << endl;*/

  fInMC->Close();
  fInMC2->Close();
  fInData.Close();
//  resultF->Close();

/*  string resultRoot;
  resultRoot = dirPre + "_rap" + yrange + "_cent" + crange + "_pT" + prange + "_workspace.root";
  Bool_t ok = ws->writeToFile(resultRoot.c_str());
  if (!ok) { cout << "CANNOT write on workspace.root file\n"; }*/

  return 0;
}




/////////////////////////////////////////////////////////
//////////////////// Sub-routines ///////////////////////
/////////////////////////////////////////////////////////
void getOptRange(string &ran, float *min, float *max) {
  if (sscanf(ran.c_str(), "%f-%f", min, max) == 0) {
    cout << ran.c_str() << ": not valid!" << endl;
    assert(0);
  }
  return ;
}

void setWSRange(RooWorkspace *ws, float lmin, float lmax, float errmin, float errmax) {

  float minRangeForPF = -4*errmax;
  if (minRangeForPF < -lmin) minRangeForPF = -lmin;

  ws->var("Jpsi_Ct")->setRange("promptfit",minRangeForPF,4*errmax);
  // ws->var("Jpsi_Ct")->setRange("psipfit",-lmin-0.3,lmax-0.3);
  ws->var("Jpsi_CtTrue")->setRange(-0.1,4.0);
  ws->var("Jpsi_CtErr")->setRange(errmin,errmax);
  ws->var("Jpsi_Ct")->setRange(-lmin,lmax);

  ws->cat("MCType")->setRange("prompt","PR");
  ws->cat("MCType")->setRange("nonprompt","NP");

  ws->var("Jpsi_Mass")->setRange("all",2.6,3.5);
  ws->var("Jpsi_Mass")->setRange("left",2.6,2.9);
  ws->var("Jpsi_Mass")->setRange("right",3.3,3.5);

  ws->cat("Jpsi_Type")->setRange("glbglb","GG");
  ws->cat("Jpsi_Type")->setRange("glbtrk","GT");

  return;
}

RooBinning setCtBinning(float lmin,float lmax) {
  RooBinning rb2(-lmin,lmax);

  if (lmax+lmin>4.9) {
    rb2.addBoundary(-1.5);
    rb2.addBoundary(-1.0);
    rb2.addBoundary(-0.8);
    rb2.addBoundary(-0.6);
    rb2.addBoundary(-0.5);
    rb2.addUniform(6,-0.5,-0.2);
    rb2.addUniform(18,-0.2,0.2);
    rb2.addUniform(9,0.2,0.5);
    rb2.addUniform(5,0.5,1.0);
    rb2.addUniform(15,1.0,lmax);
  } else if (lmax+lmin > 4.4) {
    rb2.addBoundary(-1.5);
    rb2.addBoundary(-1.0);
    rb2.addBoundary(-0.8);
    rb2.addBoundary(-0.6);
    rb2.addBoundary(-0.5);
    rb2.addUniform(9,-0.5,-0.2);
    rb2.addUniform(36,-0.2,0.2);
    rb2.addUniform(12,0.2,0.5);
    rb2.addUniform(5,0.5,1.0);
    rb2.addUniform(5,1.0,lmax);
  } else {
    rb2.addBoundary(-1.0);
    rb2.addBoundary(-0.7);
    rb2.addBoundary(-0.6);
    rb2.addBoundary(-0.5);
    rb2.addUniform(9,-0.5,-0.2);
    rb2.addUniform(40,-0.2,0.2);
    rb2.addUniform(12,0.2,0.5);
    rb2.addUniform(10,0.5,1.0);
    rb2.addUniform(4,1.0,lmax);
  }
  return rb2;
}

void defineMassBkg(RooWorkspace *ws) {
  // 1st order polynomial
  ws->factory("Polynomial::polFunct(Jpsi_Mass,{coefPol1[-0.05,-150.,150.]})");
  // Exponential
  ws->factory("Exponential::expFunct(Jpsi_Mass,coefExp[-1.,-3.,1.])");

  return;
}

void defineMassSig(RooWorkspace *ws) {
  //////// Candidates for signal
  // Normal gaussians
  ws->factory("Gaussian::signalG1(Jpsi_Mass,meanSig1[3.0975,3.05,3.15],sigmaSig1[0.02,0.008,0.2])");
//  ws->factory("Gaussian::signalG2(Jpsi_Mass,meanSig2[3.0975,3.05,3.15],sigmaSig2[0.03,0.008,0.2])");
//  ws->factory("Gaussian::signalGm1s2(Jpsi_Mass,meanSig1,sigmaSig2)");

  // Crystall Ball
  ws->factory("CBShape::signalCB(Jpsi_Mass,meanSig1,sigmaSig1,alpha[0.5,0.,3.],enne[5.,1.,30.])");
  ws->factory("CBShape::signalCB2(Jpsi_Mass,meanSig1,sigmaSig2[0.03,0.008,0.04],alpha,enne)");
  ws->factory("CBShape::signalCBWN(Jpsi_Mass,meanSig1,sigmaSig1,alpha,enneW[5.,1.,50.])");
  ws->factory("CBShape::signalCB2WN(Jpsi_Mass,meanSig1,sigmaSig2,alpha,enneW)");

  //////// Sum of signal functions
  // Sum of gaussian with different mean
//  ws->factory("SUM::sigPDF(coeffGaus[0.05,0.,1.]*signalG1,signalG2)");

  // Sum of gaussian with same mean
//  ws->factory("SUM::sigPDFm1(coeffGaus*signalG1,signalGm1s2)");

  // Sum of gaussian 1 and a crystall ball
  ws->factory("SUM::sigCBG1(coeffGaus[0.1,0.05,1.]*signalG1,signalCB)");
  // Sum of gaussian 1 and crystall ball 2
  ws->factory("SUM::sigCB2G1(coeffGaus*signalG1,signalCB2)");
  // Sum of gaussian 1 and crystall ball with wide n
  ws->factory("SUM::sigCBWNG1(coeffGaus*signalG1,signalCBWN)");
  // Sum of gaussian 1 and crystall ball 2 with wide n
  ws->factory("SUM::sigCB2WNG1(coeffGaus*signalG1,signalCB2WN)");
  return;
}


void defineCTResol(RooWorkspace *ws) {
  if (oneGaussianResol) {
    ws->factory("GaussModel::sigPR(Jpsi_Ct,meanResSigW[0.,-0.01,0.01],sigmaResSigN[0.8,0.6,2.0],one[1.0],Jpsi_CtErr)");
  } else {
    ws->factory("GaussModel::resGW(Jpsi_Ct,meanResSigW[0.,-0.01,0.01],sigmaResSigW[2.3,1.3,3.5],one[1.0],Jpsi_CtErr)");
    ws->factory("GaussModel::resGN(Jpsi_Ct,meanResSigW,sigmaResSigN[0.8,0.6,1.1],one,Jpsi_CtErr)");
    ws->factory("AddModel::sigPR({resGW,resGN},{fracRes[0.05,0.001,0.3]})");
  }

  return;
}

void defineCTBkg(RooWorkspace *ws) {
  ws->factory("Decay::bkg2(Jpsi_Ct,lambdap[0.42,0.05,1.5],sigPR,RooDecay::SingleSided)");
  ws->factory("Decay::bkg3(Jpsi_Ct,lambdam[0.79,0.01,1.5],sigPR,RooDecay::Flipped)");
  ws->factory("Decay::bkg4(Jpsi_Ct,lambdasym[0.69,0.02,5.0],sigPR,RooDecay::DoubleSided)");

  ws->factory("SUM::bkgPart1(fpm[1.0,0.0,1.0]*bkg2,bkg3)");
  ws->factory("SUM::bkgPart2(fLiving[0.9,0.0,1.0]*bkgPart1,bkg4)");
  ws->factory("SUM::bkgCtTot(fbkgCtTot[0.29,0.0,1.0]*sigPR,bkgPart2)");

  return;
}

void getMCTrueLifetime(RooWorkspace *ws, RooDataSet *redMCCutNP, float *bgmcVal, float *bctauVal, string titlestr) {

  RooFitResult *fitBMcTrue = ws->pdf("bMCTrue")->fitTo(*redMCCutNP,Minos(0),SumW2Error(kTRUE),NumCPU(4));
//  resultF->cd();
//  fitBMcTrue->Write("bMCTrue");
//  bgmcVal[0] = ws->var("Gmc1")->getVal();
//  bgmcVal[1] = ws->var("Gmc2")->getVal();
  *bgmcVal = ws->var("Gmc")->getVal();
  *bctauVal = ws->var("bTau")->getVal();

  // *** test True Lifetime fit
  RooPlot *trueframef = ws->var("Jpsi_CtTrue")->frame();
  redMCCutNP->plotOn(trueframef);
  ws->pdf("bMCTrue")->plotOn(trueframef,LineColor(kBlue),Normalization(redMCCutNP->sumEntries(),RooAbsReal::NumEvent));

  TCanvas c0f;
  c0f.cd(); c0f.SetLogy(1); trueframef->Draw();
  c0f.SaveAs(titlestr.c_str());
  // *** end test True Lifetimes
  return ;
}

void defineCTSig(RooWorkspace *ws, RooDataSet *redMCCutNP, string titlestr) {
  if (analyticBlifetime) {
    ws->factory("GaussModel::bresGTrue(Jpsi_CtTrue,mean[0.0],Gmc[0.002,0.00001,0.02])");
//    ws->factory("GaussModel::bresGTrue2(Jpsi_CtTrue,mean,Gmc2[0.002,0.00001,0.02])");
//    ws->factory("AddModel::bresGTrue({bresGTrue1,bresGTrue2},{fracCT[0.5,0,1]})");
    ws->factory("Decay::bMCTrue(Jpsi_CtTrue,bTau[0.01,0.01,1.0],bresGTrue,RooDecay::SingleSided)");
//    float GmcVal[2], bTauVal;
    float GmcVal, bTauVal;

//    float GmcValF = sqrt(pow(GmcVal[0],2)+pow(GmcVal[1]*ws->var("fracCT")->getVal(),2));
    getMCTrueLifetime(ws, redMCCutNP, &GmcVal, &bTauVal, titlestr);
    RooRealVar gmc("gmc","Sigma of MC Gaussian",GmcVal);
    RooRealVar btauFix("btauFix","Slope of MC exponential",bTauVal);   ws->import(btauFix);
    RooFormulaVar bResSigN("bResSigN", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws->var("sigmaResSigN")), *(ws->var("Jpsi_CtErr")),gmc));  ws->import(bResSigN);
    if (oneGaussianResol) {
      ws->factory("GaussModel::bresG(Jpsi_Ct,meanResSigW,bResSigN)");
    } else {
      RooFormulaVar bResSigW("bResSigW", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws->var("sigmaResSigW")), *(ws->var("Jpsi_CtErr")),gmc));  ws->import(bResSigW);
      
      ws->factory("GaussModel::bresGN(Jpsi_Ct,meanResSigW,bResSigN)");
      ws->factory("GaussModel::bresGW(Jpsi_Ct,meanResSigW,bResSigW)");
      ws->factory("AddModel::bresG({bresGW,bresGN},{fracRes})");
    }
    // fix tau_B
    // ws->factory("Decay::sigNP(Jpsi_Ct,btauFix,bresG,RooDecay::SingleSided)");
    // float tau_B
    ws->factory("Decay::sigNP(Jpsi_Ct,bTau,bresG,RooDecay::SingleSided)");
    
  } else {
    RooDataHist *binMCCutNP = new RooDataHist("redMCCutNP","MC distribution for NP signal",RooArgSet(*(ws->var("Jpsi_CtTrue"))),*redMCCutNP);
    if (oneGaussianResol) {
      RooHistPdfConv sigNP("sigNP","Non-prompt signal with narrow gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigN")),*(ws->var("one")),*(ws->var("Jpsi_CtErr")),*binMCCutNP);  ws->import(sigNP);
    } else {
      RooHistPdfConv sigNPW("sigNPW","Non-prompt signal with wide gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigW")),*(ws->var("one")),*(ws->var("Jpsi_CtErr")),*binMCCutNP);  ws->import(sigNPW);
      RooHistPdfConv sigNPN("sigNPN","Non-prompt signal with narrow gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigN")),*(ws->var("one")),*(ws->var("Jpsi_CtErr")),*binMCCutNP);  ws->import(sigNPN);
      RooAddPdf sigNP("sigNP","Non-prompt signal",RooArgSet(*(ws->pdf("sigNPW")),*(ws->pdf("sigNPN"))),RooArgSet(*(ws->var("fracRes"))));  ws->import(sigNP); 
    }
  }
  
  return;
}

RooDataHist* subtractSidebands(RooWorkspace* ws, RooDataHist* all, RooDataHist* side, float scalefactor, string varName = "Jpsi_CtErr") {
  const RooArgSet* aRow;
  const RooArgSet* aRowS;
 
  if (all->numEntries() != side->numEntries()) {
    cout << "ERROR subtractSidebands : different binning!" << endl;
    return 0;
  }

  RooDataHist* subtrData = new RooDataHist("subtrData","Subtracted data",RooArgSet(*(ws->var(varName.c_str())))); 
  for (Int_t i=0; i<all->numEntries(); i++) {
    aRow = all->get(i);
    aRowS = side->get(i);
    RooRealVar* thisVar = (RooRealVar*)aRow->find(varName.c_str());
    ws->var(varName.c_str())->setVal(thisVar->getVal());
    float newWeight = all->weight(*aRow,0,false) - scalefactor*side->weight(*aRowS,0,false);
    if (newWeight <= 2.0) newWeight = 2.0;
    subtrData->add(RooArgSet(*(ws->var(varName.c_str()))),newWeight);
  }
  return subtrData;

}

#define EDBR2PKUTree_cxx
#include "EDBR2PKUTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
//#include "BTagCalibrationStandalone.h"
//
vector<Double_t> generate_weights(TH1* data_npu_estimated, Int_t isForSynch){
	// see SimGeneral/MixingModule/python/mix_2015_25ns_Startup_PoissonOOTPU_cfi.pyy; copy and paste from there:
	const Double_t npu_probs[50] = {
/*
		4.8551E-07,
		1.74806E-06,
		3.30868E-06,
		1.62972E-05,
		4.95667E-05,
		0.000606966,
		0.003307249,
		0.010340741,
		0.022852296,
		0.041948781,
		0.058609363,
		0.067475755,
		0.072817826,
		0.075931405,
		0.076782504,
		0.076202319,
		0.074502547,
		0.072355135,
		0.069642102,
		0.064920999,
		0.05725576,
		0.047289348,
		0.036528446,
		0.026376131,
		0.017806872,
		0.011249422,
		0.006643385,
		0.003662904,
		0.001899681,
		0.00095614,
		0.00050028,
		0.000297353,
		0.000208717,
		0.000165856,
		0.000139974,
		0.000120481,
		0.000103826,
		8.88868E-05,
		7.53323E-05,
		6.30863E-05,
		5.21356E-05,
		4.24754E-05,
		3.40876E-05,
		2.69282E-05,
		2.09267E-05,
		1.5989E-05,
		4.8551E-06,
		2.42755E-06,
		4.8551E-07,
		2.42755E-07,
		1.21378E-07,
		4.8551E-08
*/
/*
                0.000108643,
                0.000388957,
                0.000332882,
                0.00038397,
                0.000549167,
                0.00105412,
                0.00459007,
                0.0210314,
                0.0573688,
                0.103986,
                0.142369,
                0.157729,
                0.147685,
                0.121027,
                0.08855,
                0.0582866,
                0.0348526,
                0.019457,
                0.0107907,
                0.00654313,
                0.00463195,
                0.00370927,
                0.0031137,
                0.00261141,
                0.00215499,
                0.00174491,
                0.00138268,
                0.00106731,
                0.000798828,
                0.00057785,
                0.00040336,
                0.00027161,
                0.000176535,
                0.00011092,
                6.75502e-05,
                4.00323e-05,
                2.32123e-05,
                1.32585e-05,
                7.51611e-06,
                4.25902e-06,
                2.42513e-06,
                1.39077e-06,
                8.02452e-07,
                4.64159e-07,
                2.67845e-07,
                1.5344e-07,
                8.68966e-08,
                4.84931e-08,
                2.6606e-08,
                1.433e-08	};
*/
		0.000829312873542,
 		0.00124276120498,
 		0.00339329181587,
 		0.00408224735376,
 		0.00383036590008,
		0.00659159288946,
 		0.00816022734493,
 		0.00943640833116,
 		0.0137777376066,
 		0.017059392038,
 		0.0213193035468,
 		0.0247343174676,
 		0.0280848773878,
 		0.0323308476564,
 		0.0370394341409,
 		0.0456917721191,
 		0.0558762890594,
 		0.0576956187107,
 		0.0625325287017,
 		0.0591603758776,
 		0.0656650815128,
 		0.0678329011676,
 		0.0625142146389,
 		0.0548068448797,
 		0.0503893295063,
 		0.040209818868,
 		0.0374446988111,
 		0.0299661572042,
 		0.0272024759921,
 		0.0219328403791,
 		0.0179586571619,
 		0.0142926728247,
 		0.00839941654725,
 		0.00522366397213,
 		0.00224457976761,
 		0.000779274977993,
 		0.000197066585944,
 		7.16031761328e-05,
 		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
 		0.0,
 		0.0,
		0.0 };
	if (isForSynch==0) { //OFFICIAL RECIPE
		vector<Double_t> result(50);
		Double_t s = 0.0;
		for(Int_t npu=0; npu<50; ++npu){
			Double_t npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));              
			result[npu] = npu_estimated / npu_probs[npu];
			s += npu_estimated;
		}
		// normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
		for(Int_t npu=0; npu<50; ++npu){
			result[npu] /= s;
		}
		return result;
	}
	else { //THIS IS FOR THE SYNCH ONLY. THIS IS NOT THE OFFICIAL RECIPE!
		vector<Double_t> result(60);
		for(Int_t npu=0; npu<60; ++npu){
			if (data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu))==NULL)
			  result[npu] = 0.;
			else {
				Double_t npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));            
				result[npu] = npu_estimated;
			}
		}
		return result;
	}

}


Double_t bsv (Int_t cud, Double_t x ) // cud=1,2,3 for central,up,down; x for pt
{
  double result=1.0;

  if (cud==1) {  //central
   result=0.892452;
  }
  else if (cud==2) { //up
   if(x<30) {result=1;}
   else if(x<50)  {result=0.892452+0.017849041149020195;}
   else if(x<70)  {result=0.892452+0.017849041149020195;}
   else if(x<100) {result=0.892452+0.017849041149020195;}
   else if(x<140) {result=0.892452+0.020885121077299118;}
   else if(x<200) {result=0.892452+0.025080939754843712;}
   else if(x<300) {result=0.892452+0.10671335458755493;}
   else if(x<670) {result=0.892452+0.16398745775222778;}
   else {result=1;}
  }
  else if (cud==3) {//down
   if(x<30) {result=1;}
   else if(x<50)  {result=0.892452-0.017849041149020195;}
   else if(x<70)  {result=0.892452-0.017849041149020195;}
   else if(x<100) {result=0.892452-0.017849041149020195;}
   else if(x<140) {result=0.892452-0.020885121077299118;}
   else if(x<200) {result=0.892452-0.025080939754843712;}
   else if(x<300) {result=0.892452-0.10671335458755493;}
   else if(x<670) {result=0.892452-0.16398745775222778;}
   else {result=1;}
  }
  return result;  
}



Double_t csv (Int_t cud, Double_t x ) // c-jet; cud=1,2,3 for central,up,down; x for pt
{
  double result=1.0;

  if (cud==1) {  //central
   result=0.892452;
  }
  else if (cud==2) { //up 
   if(x<30) {result=1;}
   else if(x<50)  {result=0.892452+0.03569808229804039;}
   else if(x<70)  {result=0.892452+0.03569808229804039;}
   else if(x<100) {result=0.892452+0.03569808229804039;}
   else if(x<140) {result=0.892452+0.041770242154598236;}
   else if(x<200) {result=0.892452+0.050161879509687424;}
   else if(x<300) {result=0.892452+0.21342670917510986;}
   else if(x<670) {result=0.892452+0.32797491550445557;}
   else {result=1;}
  }
  else if (cud==3) {//down
   if(x<30) {result=1;}
   else if(x<50)  {result=0.892452-0.03569808229804039 ;}
   else if(x<70)  {result=0.892452-0.03569808229804039 ;}
   else if(x<100) {result=0.892452-0.03569808229804039 ;}
   else if(x<140) {result=0.892452-0.041770242154598236 ;}
   else if(x<200) {result=0.892452-0.050161879509687424 ;}
   else if(x<300) {result=0.892452-0.21342670917510986 ;}
   else if(x<670) {result=0.892452-0.32797491550445557 ;}
   else {result=1;}
  }
  return result;
}

Double_t lsv (Int_t cud, Double_t x ) // light flavor; cud=1,2,3 for central,up,down; x for pt
{
  double result=1.0;

  if (cud==1) {  //central
   result=0.924144-0.000861952*x+3.46078e-06*x*x-2.4028e-09*x*x*x;
  }
  else if (cud==2) { //up 
   result=1.02632+-0.000806342*x+3.54292e-06*x*x+-2.50001e-09*x*x*x;
  }
  else if (cud==3) {//down
   result=0.821939-0.00091531*x+3.37305e-06*x*x-2.30372e-09*x*x*x;
  }
  return result;
}


void EDBR2PKUTree::Loop(TString channelname, Double_t XS, Double_t totaleventnumber, Int_t IsData) {

	std::vector<Double_t> weights_pu1; //these are made with our recipe
	std::vector<Double_t> weights_pu2; //these are made with the official recipe
	TFile* pileupFile1 = TFile::Open("pileupDataRun2016B_63mb_80X.root");  
	TH1F* pileupHisto1 = (TH1F*)pileupFile1->Get("pileup");  
	weights_pu1 = generate_weights(pileupHisto1,0);
	pileupFile1->Close();

	//  TFile* pileupFile2 = TFile::Open("puweights.root");  
	TFile* pileupFile2 = TFile::Open("PUxSynch.root");  
	TH1F *pileupHisto2 = (TH1F*)pileupFile2->Get("puweights");
	weights_pu2 = generate_weights(pileupHisto2,1);
	pileupFile2->Close();

	//TFile * input1 = new TFile ("puweights.root");
	//TH1F* hR1= (TH1F*)input1->Get("puweights");
	//zixu
	TFile * input1 = new TFile ("puweight.root");	
	TH1F* hR1= (TH1F*)input1->Get("h2");
	//TFile * input1 = new TFile ("test_mu.root");
	//TH1F* hR1= (TH1F*)input1->Get("hRatio"); //"pileup");//hRatio");


	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();

	Double_t n_deltaRlepjet = 0; 
	Double_t n_delPhijetlep = 0; 
	Double_t ntau = 0;
	Double_t number_qq = 0; 
	Double_t nmassVhad = 0; 
	Double_t nptVlepJEC = 0;
	Double_t nID_e = 0;
	Double_t npt_e = 0;
	Double_t nmet_e = 0; 
	Double_t nnum_bJet_e = 0; 
	Double_t n_delPhijetmet = 0; 

	Double_t nID_mu = 0;
	Double_t npt_mu = 0;
	Double_t nmet_mu = 0; 
	Double_t nnum_bJet_mu = 0; 
	//Double_t nbtb_mu = 0; 

	Double_t nptVhad = 0;
	Double_t yields = 0;
	//TLorentzVector jetV, genjetV;
	//some constants inside this analysis
	Double_t pi_2=1.57079632679;
	Long64_t npp = fChain->GetEntries("theWeight>0.");
	Long64_t nmm = fChain->GetEntries("theWeight<0.");
	cout<<"npp="<<npp<<" nmm="<<nmm<<" totaleventnumber="<<totaleventnumber<<endl;

	Double_t nn;
	Double_t eff_and_pu_Weight;
	Double_t eff_and_pu_Weight1;
	Float_t Identical_lumiWeight = XS;//All the events inside a sample are same lumiweight
	//Float_t Identical_lumiWeight = XS/totaleventnumber;//All the events inside a sample are same lumiweight

	Long64_t nbytes = 0, nb = 0;
	//for (Long64_t jentry=0; jentry<10;jentry++)
	for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
      //             if (jentry%5000==0)
      //                {std::cout<<jentry<<std::endl;}
		nb = fChain->GetEntry(jentry); 
		nbytes += nb;
		pfMET             = Float_t(met);
		pfMETPhi          = Float_t(metPhi);
		l_pt              = Float_t(ptlep1);
		l_eta             = Float_t(etalep1);
		l_phi             = Float_t(philep1);
		ptVhad            = Float_t(ptVhad);
		jet_eta           = Float_t(yVhad);
		jet_phi           = Float_t(phiVhad);
		jet_mass_pruned   = Float_t(massVhadJEC);
                jet_mass_puppi    = Float_t(jetAK8puppi_sdJEC);
                jet_mass_puppi_un = Float_t(jetAK8puppi_sd);
                jet_tau2tau1_puppi    = Float_t(jetAK8puppi_tau21);
                jet_pt_puppi      = Float_t(jetAK8puppi_ptJEC);
		jetAK8_mass       = Float_t(jetAK8_mass);
		jet_mass_softdrop = Float_t(sdropJEC);
		jet_tau2tau1      = Float_t(tau21);
		W_pt              = Float_t(ptVlepJEC);
		W_eta             = Float_t(yVlep);
		W_phi             = Float_t(phiVlep);
		m_lvj             = Float_t(candMassJEC);
		fjet2_pt          = Float_t(jet2_pt);
		fjet2_btag        = Float_t(jet2_btag);
		fjet3_pt          = Float_t(jet3_pt);
		fjet3_btag        = Float_t(jet3_btag);
                mtVlepnew         = Float_t(sqrt(2*ptlep1*met*(1.0-cos(philep1-metPhi))));
                
                //puppi+softdrop recorrected by Thea's "JEC"
                Double_t gencorrect=1.0;
                Double_t recocorrect_0eta1p3=1.0;
                Double_t recocorrect_1p3eta2p5=1.0;

                gencorrect=1.0-0.321*pow(jet_pt_puppi*0.0354,-1.1);
                recocorrect_0eta1p3=1.09-1.69e-04*jet_pt_puppi+3.34e-07*pow(jet_pt_puppi,2)-2.47e-10*pow(jet_pt_puppi,3)+7.8e-14*pow(jet_pt_puppi,4)-8.83e-18*pow(jet_pt_puppi,5);
                recocorrect_1p3eta2p5=1.3-7.76e-04*jet_pt_puppi+1.11e-06*pow(jet_pt_puppi,2)-6.79e-10*pow(jet_pt_puppi,3)+1.87e-13*pow(jet_pt_puppi,4)-1.9e-17*pow(jet_pt_puppi,5);
                if (fabs(jet_eta)<=1.3){jet_mass_puppi_corr=jet_mass_puppi_un*gencorrect*recocorrect_0eta1p3;}
                else if (fabs(jet_eta)<2.5 && fabs(jet_eta)>1.3){jet_mass_puppi_corr=jet_mass_puppi_un*gencorrect*recocorrect_1p3eta2p5;}


		// GEN-RECO match
		//deltaRleplep = deltaR(etalep1,philep1,etalep2,philep2);
		//deltaRWlepGen = deltaR(etaGenVlep, phiGenVlep, yVlep, phiVlep);
		//jetV.SetPtEtaPhiM(ptVhad, yVhad, phiVhad, massVhad);
		//genjetV.SetPtEtaPhiM(ptGenVhad, etaGenVhad, phiGenVhad, massGenVhad);
		//deltaRWhadGen = deltaR(etaGenVhad, phiGenVhad, yVhad, phiVhad);
		Double_t deltaRWhadGen = sqrt(pow(etaGenVhad-yVhad,2) + pow(phiGenVhad-phiVhad,2));

		//Weight Calculation
		Int_t bin = hR1->FindBin(npT);
		pileupWeight = hR1->GetBinContent(bin);		

		eff_and_pu_Weight = 0;
		eff_and_pu_Weight1 = 0;
		if(IsData>0) {
			if(npT < weights_pu1.size()){
				eff_and_pu_Weight = weights_pu1[npT];
			}
			if(npT < weights_pu2.size()){
				eff_and_pu_Weight1 = weights_pu2[npT];
			}
		}

                trigger_eff=1.0;
/*
                if (ptlep1>55. && ptlep1<70.)
                   {if (etalep1<-2.0){trigger_eff=0.812;}
                    else if (etalep1<-1.56){trigger_eff=0.814;}
                    else if (etalep1<-1.44){trigger_eff=0.781;}
                    else if (etalep1<-0.8){trigger_eff=0.802;}
                    else if (etalep1<0.0){trigger_eff=0.869;}
                    else if (etalep1<0.8){trigger_eff=0.862;}
                    else if (etalep1<1.44){trigger_eff=0.879;}
                    else if (etalep1<1.56){trigger_eff=0.642;}
                    else if (etalep1<2.0){trigger_eff=0.843;}
                    else if (etalep1>=2.0){trigger_eff=0.845;}
                   }
                else if (ptlep1>=70. && ptlep1<100.)
                   {if (etalep1<-2.0){trigger_eff=0.788;}
                    else if (etalep1<-1.56){trigger_eff=0.809;}
                    else if (etalep1<-1.44){trigger_eff=0.833;}
                    else if (etalep1<-0.8){trigger_eff=0.888;}
                    else if (etalep1<0.0){trigger_eff=0.883;}
                    else if (etalep1<0.8){trigger_eff=0.882;}
                    else if (etalep1<1.44){trigger_eff=0.884;}
                    else if (etalep1<1.56){trigger_eff=0.806;}
                    else if (etalep1<2.0){trigger_eff=0.833;}
                    else if (etalep1>=2.0){trigger_eff=0.853;}
                   }
                else if (ptlep1>=100. && ptlep1<120.)
                   {if (etalep1<-2.0){trigger_eff=0.919;}
                    else if (etalep1<-1.56){trigger_eff=0.871;}
                    else if (etalep1<-1.44){trigger_eff=0.894;}
                    else if (etalep1<-0.8){trigger_eff=0.908;}
                    else if (etalep1<0.0){trigger_eff=0.902;}
                    else if (etalep1<0.8){trigger_eff=0.924;}
                    else if (etalep1<1.44){trigger_eff=0.927;}
                    else if (etalep1<1.56){trigger_eff=0.728;}
                    else if (etalep1<2.0){trigger_eff=0.87;}
                    else if (etalep1>=2.0){trigger_eff=0.867;}
                   }
                else if (ptlep1>=120. && ptlep1<180.)
                   {if (etalep1<-2.0){trigger_eff=0.955;}
                    else if (etalep1<-1.56){trigger_eff=0.986;}
                    else if (etalep1<-1.44){trigger_eff=0.887;}
                    else if (etalep1<-0.8){trigger_eff=0.924;}
                    else if (etalep1<0.0){trigger_eff=0.941;}
                    else if (etalep1<0.8){trigger_eff=0.934;}
                    else if (etalep1<1.44){trigger_eff=0.956;}
                    else if (etalep1<1.56){trigger_eff=0.81;}
                    else if (etalep1<2.0){trigger_eff=0.977;}
                    else if (etalep1>=2.0){trigger_eff=0.948;}
                   }
                else if (ptlep1>=180. && ptlep1<250.)
                   {if (etalep1<-2.0){trigger_eff=1.0;}
                    else if (etalep1<-1.56){trigger_eff=0.984;}
                    else if (etalep1<-1.44){trigger_eff=0.918;}
                    else if (etalep1<-0.8){trigger_eff=0.934;}
                    else if (etalep1<0.0){trigger_eff=0.981;}
                    else if (etalep1<0.8){trigger_eff=0.937;}
                    else if (etalep1<1.44){trigger_eff=0.931;}
                    else if (etalep1<1.56){trigger_eff=0.939;}
                    else if (etalep1<2.0){trigger_eff=0.968;}
                    else if (etalep1>=2.0){trigger_eff=0.954;}
                   }
                else if (ptlep1>=250.)
                   {if (etalep1<-2.0){trigger_eff=1.0;}
                    else if (etalep1<-1.56){trigger_eff=1.0;}
                    else if (etalep1<-1.44){trigger_eff=1.0;}
                    else if (etalep1<-0.8){trigger_eff=0.927;}
                    else if (etalep1<0.0){trigger_eff=0.953;}
                    else if (etalep1<0.8){trigger_eff=0.959;}
                    else if (etalep1<1.44){trigger_eff=0.972;}
                    else if (etalep1<1.56){trigger_eff=1.0;}
                    else if (etalep1<2.0){trigger_eff=1.0;}
                    else if (etalep1>=2.0){trigger_eff=1.0;}
                   }
*/


/*
      TFile * input_trigger = new TFile ("EleEff-TH2D.root");
      input_trigger->cd("");
      TH2D* HLTeff= (TH2D*) input_trigger->Get("2Dh");

      Double_t ele_pt=ptlep1;
      if (ele_pt>1000) {ele_pt=999.0;}
      int elebin=HLTeff->FindBin(etalep1,ele_pt);
      trigger_eff=HLTeff->GetBinContent(elebin); 
      input_trigger->Close();
*/
                IDweight=1.0;
/*
                if (ptlep1<20.)
                   {if (fabs(etalep1)<0.8){IDweight=0.992611;}
                    else if (fabs(etalep1)<1.44){IDweight=1.11702;}
                    else if (fabs(etalep1)<1.56){IDweight=0.825893;}
                    else if (fabs(etalep1)<2.0){IDweight=0.486322;}
                    else if (fabs(etalep1)>=2.0){IDweight=0.795031;}
                   }
                else if (ptlep1<30.)
                   {if (fabs(etalep1)<0.8){IDweight=0.970909;}
                    else if (fabs(etalep1)<1.44){IDweight=0.949343;}
                    else if (fabs(etalep1)<1.56){IDweight=0.812865;}
                    else if (fabs(etalep1)<2.0){IDweight=0.604607;}
                    else if (fabs(etalep1)>=2.0){IDweight=0.794971;}
                   }
                else if (ptlep1<40.)
                   {if (fabs(etalep1)<0.8){IDweight=0.975867;}
                    else if (fabs(etalep1)<1.44){IDweight=0.975232;}
                    else if (fabs(etalep1)<1.56){IDweight=0.819502;}
                    else if (fabs(etalep1)<2.0){IDweight=0.694099;}
                    else if (fabs(etalep1)>=2.0){IDweight=0.853312;}
                   }
                else if (ptlep1<50.)
                   {if (fabs(etalep1)<0.8){IDweight=0.974359;}
                    else if (fabs(etalep1)<1.44){IDweight=0.972603;}
                    else if (fabs(etalep1)<1.56){IDweight=0.819805;}
                    else if (fabs(etalep1)<2.0){IDweight=0.740027;}
                    else if (fabs(etalep1)>=2.0){IDweight=0.861777;}
                   }
                else if (ptlep1>=50.)
                   {if (fabs(etalep1)<0.8){IDweight=0.984752;}
                    else if (fabs(etalep1)<1.44){IDweight=0.978093;}
                    else if (fabs(etalep1)<1.56){IDweight=0.82263;}
                    else if (fabs(etalep1)<2.0){IDweight=0.80829;}
                    else if (fabs(etalep1)>=2.0){IDweight=0.889186;}
                   }

*/
/*
      TFile * input_ID = new TFile ("elesf.root");
      input_ID->cd("");
      TH2D* hele= (TH2D*) input_ID->Get("elesf");

      Double_t elept=ptlep1;
      if (elept>200) {elept=199.0;}
      int elebin0=hele->FindBin(etalep1,elept);
      IDweight=hele->GetBinContent(elebin0);   
      input_ID->Close();
*/
		//cout << "pileupWeight:"<<pileupWeight<< " eff_and_pu_Weight:" << eff_and_pu_Weight << " eff_and_pu_Weight1:" << eff_and_pu_Weight1 << endl;
		if(theWeight>0) nn=1;
		else nn= -1;
		if(npp>0) lumiWeight=Identical_lumiWeight/(npp-nmm);
		else lumiWeight=Identical_lumiWeight/nentries;
		weight_nobtag=lumiWeight*triggerWeight*eff_and_pu_Weight*nn*trigger_eff;
		//weight=lumiWeight*triggerWeight*pileupWeight*nn;
		if (IsData>1 ) weight_nobtag = weight_nobtag*1.21;

		//lumiWeight=Identical_lumiWeight;
		//if(npp>0) weight=lumiWeight*triggerWeight*pileupWeight/(npp-nmm)*nn*0.04024;//0.00559;
		//else weight=lumiWeight*triggerWeight*pileupWeight/nentries*0.04024;//0.00559;
		if ( IsData==0 ) weight_nobtag=1;
		//Weight Calculation Done



//--BSF----------------------------
      btagweight_center=1.0;
      btagweight_up=1.0;
      btagweight_down=1.0;
      
/*
      double  bweight=1.0, dweight=1.0;
      double  bweightup=1.0;
      double  bweightdown=1.0;

      double  beff=1.0, ceff=1.0, leff=1.0;

      for(Int_t i=0; i<8; i++)  {
       if(ak4jet_pt[i]>30 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8[i]>=0.8 ) {
         if(abs(ak4jet_hf[i])==5 && ak4jet_pt[i]<670.) 
          { 
            double jet_sf    = bsv(1,ak4jet_pt[i]); 
            double jet_sfu    = bsv(2,ak4jet_pt[i]);
            double jet_sfd    = bsv(3,ak4jet_pt[i]);
//            int bbin=hb->FindBin(ak4jet_pt[i],ak4jet_eta[i]);
//            beff=hb->GetBinContent(bbin);
           if (ak4jet_pt[i]>30. && ak4jet_pt[i]<50.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.702463;}
               else if (ak4jet_eta[i]<-0.8){beff=0.628582;}
               else if (ak4jet_eta[i]>0.8){beff=0.628031;}
              }
           else if (ak4jet_pt[i]>=50. && ak4jet_pt[i]<70.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.737214;}
               else if (ak4jet_eta[i]<-0.8){beff=0.665792;}
              else if (ak4jet_eta[i]>0.8){beff=0.667027;}
              }
           else if (ak4jet_pt[i]>=70. && ak4jet_pt[i]<100.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.75137;}
               else if (ak4jet_eta[i]<-0.8){beff=0.67818;}
               else if (ak4jet_eta[i]>0.8){beff=0.682461;}
              }
           else if (ak4jet_pt[i]>=100. && ak4jet_pt[i]<140.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.749506;}
               else if (ak4jet_eta[i]<-0.8){beff=0.674191;}
               else if (ak4jet_eta[i]>0.8){beff=0.678941;}
              }
           else if (ak4jet_pt[i]>=140. && ak4jet_pt[i]<200.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.73442;}
               else if (ak4jet_eta[i]<-0.8){beff=0.662146;}
               else if (ak4jet_eta[i]>0.8){beff=0.665269;}
              }
           else if (ak4jet_pt[i]>=200. && ak4jet_pt[i]<300.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.69194;}
               else if (ak4jet_eta[i]<-0.8){beff=0.614117;}
               else if (ak4jet_eta[i]>0.8){beff=0.616712;}
              }
           else if (ak4jet_pt[i]>=300. && ak4jet_pt[i]<670.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.620711;}
               else if (ak4jet_eta[i]<-0.8){beff=0.579713;}
               else if (ak4jet_eta[i]>0.8){beff=0.567196;}
              }

             if(ak4jet_icsv[i]>0.800) {
                   bweight=bweight*beff*jet_sf; dweight=dweight*beff;
                   bweightup=bweightup*beff*jet_sfu; 
                   bweightdown=bweightdown*beff*jet_sfd; 
                   }
             else {
                   bweight=bweight*(1-beff*jet_sf); dweight=dweight*(1-beff);
                   bweightup=bweightup*(1-beff*jet_sfu); 
                   bweightdown=bweightdown*(1-beff*jet_sfd); 
                  }
          }
        else if(abs(ak4jet_hf[i])==4 && ak4jet_pt[i]<670.)
          {   
            double jet_sf    = csv(1,ak4jet_pt[i]);
            double jet_sfu    = csv(2,ak4jet_pt[i]);
            double jet_sfd    = csv(3,ak4jet_pt[i]);
//            int cbin=hc->FindBin(ak4jet_pt[i],ak4jet_eta[i]);
//            ceff=hc->GetBinContent(cbin);
           if (ak4jet_pt[i]>30. && ak4jet_pt[i]<50.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.198287;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.172776;}
               else if (ak4jet_eta[i]>0.8){ceff=0.169512;}
              }
           else if (ak4jet_pt[i]>=50. && ak4jet_pt[i]<70.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.200463;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.174732;}
              else if (ak4jet_eta[i]>0.8){ceff=0.173456;}
              }
           else if (ak4jet_pt[i]>=70. && ak4jet_pt[i]<100.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.203827;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.179102;}
               else if (ak4jet_eta[i]>0.8){ceff=0.181135;}
              }
           else if (ak4jet_pt[i]>=100. && ak4jet_pt[i]<140.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.201279;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.179472;}
               else if (ak4jet_eta[i]>0.8){ceff=0.18073;}
              }
           else if (ak4jet_pt[i]>=140. && ak4jet_pt[i]<200.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.195573;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.183093;}
               else if (ak4jet_eta[i]>0.8){ceff=0.182443;}
              }
           else if (ak4jet_pt[i]>=200. && ak4jet_pt[i]<300.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.180287;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.166444;}
               else if (ak4jet_eta[i]>0.8){ceff=0.166148;}
              }
           else if (ak4jet_pt[i]>=300. && ak4jet_pt[i]<670.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.161175;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.161221;}
               else if (ak4jet_eta[i]>0.8){ceff=0.162861;}
              }

             if(ak4jet_icsv[i]>0.800) {
                   bweight=bweight*ceff*jet_sf; dweight=dweight*ceff;
                   bweightup=bweightup*ceff*jet_sfu; 
                   bweightdown=bweightdown*ceff*jet_sfd; 
                   }
             else {
                   bweight=bweight*(1-ceff*jet_sf); dweight=dweight*(1-ceff);
                   bweightup=bweightup*(1-ceff*jet_sfu); 
                   bweightdown=bweightdown*(1-ceff*jet_sfd); 
                  }
          }
        else if (abs(ak4jet_hf[i])==0 && ak4jet_pt[i]<1000.)
          {  
            double jet_sf    = lsv(1,ak4jet_pt[i]);
            double jet_sfu    = lsv(2,ak4jet_pt[i]);
            double jet_sfd    = lsv(3,ak4jet_pt[i]);
//            int lbin=hl->FindBin(ak4jet_pt[i],ak4jet_eta[i]);
//            leff=hl->GetBinContent(lbin);
           if (ak4jet_pt[i]>=20. && ak4jet_pt[i]<200.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.018909;}
               else if (ak4jet_eta[i]<-0.8){leff=0.021948;}
               else if (ak4jet_eta[i]>0.8){leff=0.020182;}
              }
           else if (ak4jet_pt[i]>=200. && ak4jet_pt[i]<500.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.015463;}
               else if (ak4jet_eta[i]<-0.8){leff=0.021675;}
               else if (ak4jet_eta[i]>0.8){leff=0.021937;}
              }
           else if (ak4jet_pt[i]>=500.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.025439;}
               else if (ak4jet_eta[i]<-0.8){leff=0.030005;}
               else if (ak4jet_eta[i]>0.8){leff=0.029121;}
              }

             if(ak4jet_icsv[i]>0.800) {
                   bweight=bweight*leff*jet_sf; dweight=dweight*leff;
                   bweightup=bweightup*leff*jet_sfu; 
                   bweightdown=bweightdown*leff*jet_sfd; 
                   }

             else {
                   bweight=bweight*(1-leff*jet_sf); dweight=dweight*(1-leff);
                   bweightup=bweightup*(1-leff*jet_sfu); 
                   bweightdown=bweightdown*(1-leff*jet_sfd); 
                  }
          }
      //std::cout<<" yy1 "<<abs(ak4jet_hf[i])<<" "<<dweight<<" "<<bweight<<std::endl;
        } 
      //std::cout<<" yy2 "<<abs(ak4jet_hf[i])<<" "<<ak4jet_pt[i]<<" "<<fabs(ak4jet_eta[i])<<" "<<ak4jet_IDLoose[i]<<" "<<deltaRAK4AK8[i]<<std::endl;
     }
//     input->Close();
//     std::cout<<" xx "<<bweight/dweight<<" "<<bweightup/dweight<<" "<<bweightdown/dweight<<std::endl;

     btagweight_center=bweight/dweight;
     btagweight_up=bweightup/dweight;
     btagweight_down=bweightdown/dweight;
*/
//     std::cout<<btagweight_center<<std::endl;

     
//--BSF----------------------------

                if(theWeight>0) nn=1;
                else nn= -1;
                if(npp>0) lumiWeight=Identical_lumiWeight/(npp-nmm);
                else lumiWeight=Identical_lumiWeight/nentries;
                weight=lumiWeight*triggerWeight*eff_and_pu_Weight*nn*trigger_eff*btagweight_center*IDweight;
      //          std::cout<<IDweight<<"        "<<btagweight_center<<"	"<<trigger_eff<<std::endl;
                if (IsData>1 ) weight = weight*1.21;
                if ( IsData==0 ) weight=1;

		//number of bjet calculation
		num_bJet=0.;
		num_bJet_loose=0.;
		num_bJet_tight=0.;
		for(Int_t i=0; i<8; i++)  {
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.800 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8[i]>=0.8 ) {num_bJet=num_bJet+1;}
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.46 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8[i]>=0.8 ) {num_bJet_loose=num_bJet_loose+1;}
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.935 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8[i]>=0.8 ) {num_bJet_tight=num_bJet_tight+1;}
		}
		nbtag=num_bJet;
		//number of bjet calculation Done

		Int_t nLooseLep=nLooseEle+nLooseMu;//the tight Lep included

		Double_t isAnaHP=1.;
		Double_t isAnaLP=1.;
		Double_t isAnaNP=1.;
		Double_t isTTBarControl=1.;
		Int_t tmp_categoryID_channel=0;
		if( channelname=="el" ){
		tmp_categoryID_channel=-1;// -1 for el; 1 for mu

			//HP: 0<tau21<=0.5;
			if (isAnaHP>0 && lep==11 && (HLT_Ele1>0 || HLT_Ele3>0) && nLooseLep==1){ nID_e = nID_e+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptlep1>55 && fabs(etalep1)<2.5){ npt_e = npt_e+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && MET_et>80) { nmet_e = nmet_e+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptVlepJEC > 200.) { nptVlepJEC = nptVlepJEC +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptVhad>200 && fabs(yVhad)<2.4 && IDLoose>0 ){ nptVhad = nptVhad+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && num_bJet<1){ nnum_bJet_e = nnum_bJet_e +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && deltaRlepjet>pi_2) {n_deltaRlepjet = n_deltaRlepjet+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && fabs(delPhijetmet) >2.0)  {n_delPhijetmet = n_delPhijetmet+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && fabs(delPhijetlep)>2.0) {n_delPhijetlep = n_delPhijetlep +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && tau21>0. && tau21<=0.45) {ntau = ntau+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && massVhadJEC>40 && massVhadJEC <150)// && m_lvj>200 && m_lvj<5000)
			{
				nmassVhad = nmassVhad +1;
				yields = yields + weight;
				(*file_cutflow)<<"event:"<<event<<endl;
			} else{ isAnaHP=-1; }

			//LP: 0.5<tau21<=0.75;
			if ( lep==11 && (HLT_Ele1>0 || HLT_Ele3>0) && nLooseLep==1 && ptlep1>55 && fabs(etalep1)<2.5 && MET_et>80 && ptVlepJEC > 200.  && ptVhad>200 && fabs(yVhad)<2.4 && IDLoose>0 && num_bJet<1 && deltaRlepjet>pi_2 && fabs(delPhijetmet) >2.0 && fabs(delPhijetlep)>2.0 && tau21>0.45 && tau21<=0.75 && (( massVhadJEC >40 &&  massVhadJEC< 150 )) )
			{ isAnaLP=1.; } 
			else{ isAnaLP=-1.; }
			//NP: 0.75<tau21
			if ( lep==11 && (HLT_Ele1>0 || HLT_Ele3>0) && nLooseLep==1 && ptlep1>55 && fabs(etalep1)<2.5 && MET_et>80 && ptVlepJEC > 200.  && ptVhad>200 && fabs(yVhad)<2.4 && IDLoose>0 && num_bJet<1 && deltaRlepjet>pi_2 && fabs(delPhijetmet) >2.0 && fabs(delPhijetlep)>2.0 && tau21>0.75 && (( massVhadJEC >40 &&  massVhadJEC< 150 )) )
			{ isAnaNP=1.; } 
			else{ isAnaNP=-1.; }


			//TTbar control
			if ( lep==11 && (HLT_Ele1>0 || HLT_Ele3>0) && nLooseLep==1 && ptlep1>55 && fabs(etalep1)<2.5 && MET_et>80 && ptVlepJEC > 200. &&ptVhad>200 && fabs(yVhad)<2.4 && IDLoose>0  && num_bJet>0 &&  massVhadJEC>40 && massVhadJEC <150)
			{ isTTBarControl=1.; } 
			else{ isTTBarControl=-1.; }
		}
		else if( channelname=="mu" ){
		tmp_categoryID_channel=1;// -1 for el; 1 for mu
			//HP: 0<tau21<=0.5;
			if (isAnaHP>0 && lep==13 && HLT_Mu1>0 && trackIso/ptlep1<0.1 && fabs(etalep1)<2.1 && nLooseLep==1 ) { nID_mu = nID_mu+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptlep1>50){ npt_mu = npt_mu+1; }else{ isAnaHP=-1; }
			if (isAnaHP>0 && MET_et>40) { nmet_mu = nmet_mu+1; }else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptVlepJEC>200) { nptVlepJEC = nptVlepJEC +1;} else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptVhad>200 && fabs(yVhad)<2.4 && IDLoose>0){ nptVhad = nptVhad+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && num_bJet<1){ nnum_bJet_mu = nnum_bJet_mu +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && deltaRlepjet>pi_2) {n_deltaRlepjet = n_deltaRlepjet+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && fabs(delPhijetmet) >2.0)  {n_delPhijetmet = n_delPhijetmet+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && fabs(delPhijetlep)>2.0) {n_delPhijetlep = n_delPhijetlep +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && tau21>0. && tau21<=0.45) {ntau = ntau+1;} else{ isAnaHP=-1; }
			if (isAnaHP>0 && (( massVhadJEC >40&& massVhadJEC< 150 )))// && m_lvj>100 && m_lvj<5000 )
			{ 
				nmassVhad = nmassVhad +1; 
				(*file_cutflow)<<"event:"<<event<<endl;
			} else{ isAnaHP=-1; }

			//LP: 0.5<tau21<=0.75;
			if (lep==13 && HLT_Mu1>0 && trackIso/ptlep1<0.1 && fabs(etalep1)<2.1 && nLooseLep==1 && ptlep1>50 && MET_et>40 && ptVlepJEC>200 && ptVhad>200 && fabs(yVhad)<2.4 && IDLoose>0 && num_bJet<1 && deltaRlepjet>pi_2 && fabs(delPhijetmet) >2.0 && fabs(delPhijetlep)>2.0 && tau21>0.45 && tau21<=0.75 && (( massVhadJEC >40&& massVhadJEC< 150 )))
			{ isAnaLP=1.; } 
			else{ isAnaLP=-1.; }

			//NP: 0.75<tau21;
			if (lep==13 && HLT_Mu1>0 && trackIso/ptlep1<0.1 && fabs(etalep1)<2.1 && nLooseLep==1 && ptlep1>50 && MET_et>40 && ptVlepJEC>200 && ptVhad>200 && fabs(yVhad)<2.4 && IDLoose>0 && num_bJet<1 && deltaRlepjet>pi_2 && fabs(delPhijetmet) >2.0 && fabs(delPhijetlep)>2.0 && tau21>0.75  && (( massVhadJEC >40&& massVhadJEC< 150 )))
			{ isAnaNP=1.; } 
			else{ isAnaNP=-1.; }

			//TTbar control
			if (lep==13 && HLT_Mu1>0 && trackIso/ptlep1<0.1 && fabs(etalep1)<2.1 && nLooseLep==1 && ptlep1>50 && MET_et>40 && ptVlepJEC>200 && ptVhad>200 && fabs(yVhad)<2.4 && IDLoose>0 && num_bJet>0 && massVhadJEC>40 && massVhadJEC <150)
			{ isTTBarControl=1.; } 
			else{ isTTBarControl=-1.; }

		}else{
			cout<<"We don't know channelname:"<<channelname<<endl;
		}

		Int_t tmp_categoryID_eventselection=0;
		if(isAnaHP>0)tmp_categoryID_eventselection=1;
		else if(isAnaLP>0)tmp_categoryID_eventselection=2;
		else if(isAnaNP>0)tmp_categoryID_eventselection=4;
		else if(isTTBarControl>0)tmp_categoryID_eventselection=3;
		else tmp_categoryID_eventselection=100;

		CategoryID=tmp_categoryID_channel* tmp_categoryID_eventselection;

		isMatch=1.;
		if(deltaRWhadGen >= 0.3) isMatch=-1;
		//cout << "massVhad" << massVhad << "jet_mass_pruned " << jet_mass_pruned << endl;
		if(tau21<=0){vTagID=2;}
		else if(tau21>0.45 && tau21<=0.60){vTagID=1;}
		else if(tau21>0.60 && tau21<=0.75){vTagID=0;}
		else if(tau21>0.75 && tau21<=1){vTagID=-1;}
		else {vTagID=-2;}

		if(TMath::Abs(CategoryID)<10) ExTree->Fill();
	}

	if(channelname=="el"){ 
		std::cout << "nID_e" << nID_e << "; npt_e" << npt_e << "; nmet_e" << nmet_e << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<<"; nnum_bJet_e" << nnum_bJet_e <<"; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad << "; number_qq" << number_qq << "; yields " << yields << std::endl;
		(*file_cutflow) << "nID_e" << nID_e << "; npt_e" << npt_e << "; nmet_e" << nmet_e << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<<"; nnum_bJet_e" << nnum_bJet_e <<"; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad << "; number_qq" << number_qq << std::endl;
	}
	if(channelname=="mu"){
		std::cout << "nID_mu" << nID_mu << "; npt_mu" << npt_mu << "; nmet_mu" << nmet_mu << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<< "; nnum_bJet_mu" << nnum_bJet_mu << "; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad<< "; number_qq" << number_qq << std::endl;
		(*file_cutflow) << "nID_mu" << nID_mu << "; npt_mu" << npt_mu << "; nmet_mu" << nmet_mu << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<< "; nnum_bJet_mu" << nnum_bJet_mu << "; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad<< "; number_qq" << number_qq << std::endl;
	}

}

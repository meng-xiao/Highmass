
#include "external_cConstants.h"
//float getDVBF2jetsConstant(float ZZMass){
//  float par[9]={
//    1.876,
//    -55.488,
//    403.32,
//    0.3906,
//    80.8,
//    27.7,
//    -0.06,
//    54.97,
//    309.96
//  };
//  float kappa =
//    pow(1.-atan((ZZMass-par[1])/par[2])*2./TMath::Pi(), par[0])
//    + par[3]*exp(-pow((ZZMass-par[4])/par[5], 2))
//    + par[6]*exp(-pow((ZZMass-par[7])/par[8], 2));
//  float constant = kappa/(1.-kappa);
//  return constant;
//}

void data(){


	TString treename [3]={"ZZTree/candTree","ZZTreelooseEle/candTree","ZZTreetle/candTree"};
	TString newtreename[3]={"","_rse","_tle"};
	TFile* fnew = new TFile("data.root","recreate");
	int count=0;
	for(int t =0;t<2;t++){
		TTree *tnew = new TTree("SelectedTree"+newtreename[t],"SelectedTree"+newtreename[t]);
		TChain *tqqzz= new TChain(treename[t]);
		tqqzz->Add("root://lxcms03//data3/Higgs/170222/AllData/ZZ4lAnalysis.root");
		float ZZPt,ZZMass;
		vector<float> *LepPt=new vector<float>;
		short Z1Flav,Z2Flav;
		short nCleanedJetsPt30;
		float pvbf_VAJHU_old;
		float phjj_VAJHU_old;
		float bkg_VAMCFM,p0plus_VAJHU;
		short ZZsel;
		float TLE_dR_Z;
		short nExtraLep;
		short nCleanedJetsPt30BTagged_bTagSF;
		vector<short> *LepLepId=0;
		tqqzz->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal",&pvbf_VAJHU_old);
		tqqzz->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal",&phjj_VAJHU_old);
		tqqzz->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_JHUGen",&p0plus_VAJHU);
		tqqzz->SetBranchAddress("p_QQB_BKG_MCFM",&bkg_VAMCFM);
		tqqzz->SetBranchAddress("ZZPt",&ZZPt);
		tqqzz->SetBranchAddress("ZZMass",&ZZMass);
		tqqzz->SetBranchAddress("Z1Flav",&Z1Flav);
		tqqzz->SetBranchAddress("Z2Flav",&Z2Flav);
		tqqzz->SetBranchAddress("ZZsel",&ZZsel);
		tqqzz->SetBranchAddress("TLE_dR_Z",&TLE_dR_Z);
		tqqzz->SetBranchAddress("LepLepId",&LepLepId);
		tqqzz->SetBranchAddress("LepPt",&LepPt);
		tqqzz->SetBranchAddress("nCleanedJetsPt30",&nCleanedJetsPt30);
		tqqzz->SetBranchAddress("nExtraLep",&nExtraLep);
		tqqzz->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF",&nCleanedJetsPt30BTagged_bTagSF);
		int chan;
		int vbfcate;
		float dbkg_kin;
		if(t==0)
		tnew->Branch("mreco",&ZZMass,"mreco/F");
		else
		tnew->Branch("mreco_rse",&ZZMass,"mreco_rse/F");

		tnew->Branch("dbkg_kin",&dbkg_kin,"dbkg_kin/F");
		tnew->Branch("chan",&chan,"chan/I");
		tnew->Branch("vbfcate",&vbfcate,"vbfcate/I");
		for(int i=0;i<tqqzz->GetEntries();i++){
			tqqzz->GetEntry(i);
			if(ZZsel!=120 && t>0)
				continue;
			if(t>0)
				if(ZZMass<300)
					continue;

			if(t!=2){
				if(abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)==121){
					chan=2;
				}
				else if (abs(Z1Flav)==abs(Z2Flav) && abs(Z1Flav)!=121){
					chan=1;
				}
				else{
					chan=3;
				}
			}
			bool patle = true;
			if(t==2){
				for (int k=0 ; k<LepLepId->size() ; k++){ 
					if (abs(LepLepId->at(k))==22 && LepPt->at(k)<=30) 
						patle = false;
				}
				if(TLE_dR_Z<=1.6)
					patle = false;
				if(abs(Z1Flav*Z2Flav)==29282)
					chan =2;
				else if(abs(Z1Flav*Z2Flav)==40898)
					chan =3;
			}
			if(!patle)
				continue;	
			float c_Mela2j = getDVBF2jetsConstant(ZZMass);
			float vbfMela= 1./(1.+ c_Mela2j*phjj_VAJHU_old/pvbf_VAJHU_old);
			float WP_VBF2j = getDVBF2jetsWP(ZZMass, 0);
			if(t==1)
				vbfcate =2;
			else{
				//cout<< vbfMela<<"\t"<<WP_VBF2j<< "\t"<<nCleanedJetsPt30<<"\t"<<nExtraLep<<"\t"<<nCleanedJetsPt30BTagged_bTagSF<<endl;
				//		 if(vbfMela> (1.043-460./(ZZMass+634.)) && nCleanedJetsPt30>=2)
				if( nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && vbfMela>WP_VBF2j ){
					vbfcate=1;
					if(ZZMass>118&&ZZMass<130)
						count++;
				}

				else
					vbfcate=0;
			}

			short ZZFlav = Z1Flav*Z2Flav;
			dbkg_kin = p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFM*getDbkgkinConstant(ZZFlav, ZZMass));
			tnew->Fill();
		}
		tqqzz->SetLineColor(2);
		tqqzz->SetMarkerColor(2);
	fnew->cd();
	tnew->Write();
	}
	fnew->Close();
	cout<<count<<endl;
}

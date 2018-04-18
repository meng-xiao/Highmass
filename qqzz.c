#include "external_cConstants.h"
void qqzz(){

	TString treename [3]={"ZZTree/candTree","ZZTreelooseEle/candTree","ZZTreetle/candTree"};
	TString newtreename[3]={"","_rse","_tle"};
	TFile* fnew = new TFile("qqzz_80_Moriond.root","recreate");
	        int nbins=32;
	        double xbin[33]={
			110,120,130,140,150,160,170,180,190,200,
			210,220,230,240,250,270,290,310,330,350,
			370,390,420,460,500,540,580,620,660,700,
			750,800,3500
	                };
	
		TFile *input_file= TFile::Open("root://lxcms03://data3/Higgs/170222/ZZTo4l/ZZ4lAnalysis.root");
		TH1F *hCounters= (TH1F*)input_file->Get("ZZTree/Counters");
		float gen_sum_weights = hCounters->GetBinContent(40);
		input_file->Close();

	for(int t =0;t<2;t++){
		TChain *tqqzz= new TChain(treename[t]);
		tqqzz->Add("root://lxcms03://data3/Higgs/170222/ZZTo4l/ZZ4lAnalysis.root");
		

		//TH2F *temp_zz = new TH2F("temp_zz"+newtreename[t],"",389,110.,4000.,30,0.,1.); 
		TH2F *temp_zz_4e = new TH2F("vartemp_zz_4e"+newtreename[t],"",nbins,xbin,30,0.,1.); 
		TH2F *temp_zz_4mu = new TH2F("vartemp_zz_4mu"+newtreename[t],"",nbins,xbin,30,0.,1.); 
		TH2F *temp_zz_2e2mu = new TH2F("vartemp_zz_2e2mu"+newtreename[t],"",nbins,xbin,30,0.,1.); 
		//		TH2F *temp_zz = new TH2F("temp_zz"+newtreename[t],"",nbins,xbin,30,0.,1.); 
		//tqqzz->Add("root://lxcms03://data3/Higgs/160624/ZZTo4l/ZZ4lAnalysis.root");
		float ZZPt,ZZMass;
		float xsec,KFactorEWKqqZZ,overallEventWeight,KFactorQCDqqZZ_M;
		vector<float> *LepPt=new vector<float>;
		short Z1Flav,Z2Flav;
		short nCleanedJetsPt30;
		float pvbf_VAJHU_old;
		float phjj_VAJHU_old;
		float bkg_VAMCFM,p0plus_VAJHU;
		short ZZsel;
		float TLE_dR_Z;
		vector<short> *LepLepId=0;
		short nExtraLep;
		short nCleanedJetsPt30BTagged_bTagSF;
		tqqzz->SetBranchAddress("nExtraLep",&nExtraLep);
		tqqzz->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF",&nCleanedJetsPt30BTagged_bTagSF);

		tqqzz->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal",&pvbf_VAJHU_old);
		tqqzz->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal",&phjj_VAJHU_old);
		tqqzz->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_JHUGen",&p0plus_VAJHU);
		tqqzz->SetBranchAddress("p_QQB_BKG_MCFM",&bkg_VAMCFM);
		tqqzz->SetBranchAddress("ZZPt",&ZZPt);
		tqqzz->SetBranchAddress("ZZMass",&ZZMass);
		tqqzz->SetBranchAddress("Z1Flav",&Z1Flav);
		tqqzz->SetBranchAddress("Z2Flav",&Z2Flav);
		tqqzz->SetBranchAddress("KFactor_EW_qqZZ",&KFactorEWKqqZZ);
		//tqqzz->SetBranchAddress("KFactorEWKqqZZ",&KFactorEWKqqZZ);
		tqqzz->SetBranchAddress("xsec",&xsec);
		tqqzz->SetBranchAddress("ZZsel",&ZZsel);
		tqqzz->SetBranchAddress("TLE_dR_Z",&TLE_dR_Z);
		tqqzz->SetBranchAddress("LepLepId",&LepLepId);
		tqqzz->SetBranchAddress("LepPt",&LepPt);
		tqqzz->SetBranchAddress("overallEventWeight",&overallEventWeight);
		tqqzz->SetBranchAddress("KFactor_QCD_qqZZ_M",&KFactorQCDqqZZ_M);
		//tqqzz->SetBranchAddress("KFactorQCDqqZZ_M",&KFactorQCDqqZZ_M);
		tqqzz->SetBranchAddress("nCleanedJetsPt30",&nCleanedJetsPt30);
		float weight, weight_up, weight_dn;
		float weight_vbf, weight_vbf_up, weight_vbf_dn;
		int chan;
		int vbfcate;
		float dbkg_kin;
		TTree *tnew = new TTree("SelectedTree"+newtreename[t],"SelectedTree"+newtreename[t]);
		tnew->Branch("mreco",&ZZMass,"mreco/F");
		tnew->Branch("dbkg_kin",&dbkg_kin,"dbkg_kin/F");
		tnew->Branch("weight",&weight,"weight/F");
		tnew->Branch("weight_up",&weight_up,"weight_up/F");
		tnew->Branch("weight_dn",&weight_dn,"weight_dn/F");
		tnew->Branch("weight_vbf",&weight_vbf,"weight_vbf/F");
		tnew->Branch("weight_vbf_up",&weight_vbf_up,"weight_vbf_up/F");
		tnew->Branch("weight_vbf_dn",&weight_vbf_dn,"weight_vbf_dn/F");
		tnew->Branch("chan",&chan,"chan/I");
		tnew->Branch("vbfcate",&vbfcate,"vbfcate/I");
		for(int i=0;i<tqqzz->GetEntries();i++){
			tqqzz->GetEntry(i);
			if(ZZsel!=120 && t>0)
				continue;
			if(t>0)
				if(ZZMass<300)
					continue;
			weight= xsec*KFactorEWKqqZZ * overallEventWeight*KFactorQCDqqZZ_M/gen_sum_weights;
			//      double eff = 6.548111e-03 - 5.866523e-06*ZZMass*TMath::Gaus((ZZMass-2.432632e+02)/2.272477e+01);
//			double eff = 0.0139011;
			double eff = 0.0155;
			if(t>0)
				eff=0;
			double rho = ZZPt/(LepPt->at(0)+LepPt->at(1)+LepPt->at(2)+LepPt->at(3)); 

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
			short ZZFlav = Z1Flav*Z2Flav;
			dbkg_kin = p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFM*getDbkgkinConstant(ZZFlav, ZZMass));
			if(chan==1)
			temp_zz_4mu->Fill(ZZMass,dbkg_kin,weight);
			else if(chan==2)
			temp_zz_4e->Fill(ZZMass,dbkg_kin,weight);
			else
			temp_zz_2e2mu->Fill(ZZMass,dbkg_kin,weight);

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
			//float vbfMela = pvbf_VAJHU_old / ( phjj_VAJHU_old*0.06 + pvbf_VAJHU_old );
			//float vbfMela = pvbf_VAJHU_old / ( phjj_VAJHU_old + pvbf_VAJHU_old );
		float c_Mela2j = getDVBF2jetsConstant(ZZMass);
		float vbfMela= 1./(1.+ c_Mela2j*phjj_VAJHU_old/pvbf_VAJHU_old);
		float WP_VBF2j = getDVBF2jetsWP(ZZMass, 0);
		   if( nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && vbfMela>WP_VBF2j )
//			if(vbfMela>WP_VBF2j  && nCleanedJetsPt30>=2)
				vbfcate=1;
			else
				vbfcate=0;

			if(rho<0.3){
				weight_dn = weight*(1-abs((KFactorQCDqqZZ_M-1)*(KFactorEWKqqZZ-1))); 
				weight_up = weight*(1+abs((KFactorQCDqqZZ_M-1)*(KFactorEWKqqZZ-1))); 
			}
			else{
				weight_dn = weight*(1-abs((KFactorEWKqqZZ-1))); 
				weight_up = weight*(1+abs((KFactorEWKqqZZ-1))); 
			}
			
			weight*= (1-eff);
			weight_dn*= (1-eff);
			weight_up*= (1-eff);

			weight_vbf= weight*eff;
			weight_vbf_dn= weight_dn*eff;
			weight_vbf_up= weight_up*eff;
			

			//weight_vbf = weight;
			//weight_vbf_up = weight_up;
			//weight_vbf_dn = weight_dn;

			tnew->Fill();
		}
		//		tnew->Draw("ZZMass","weight_up","same");
		//		tnew->Draw("ZZMass","weight_dn","same");
		tqqzz->SetLineColor(2);
		tqqzz->SetMarkerColor(2);
		//		tqqzz->Draw("ZZMass","xsec*KFactorEWKqqZZ * overallEventWeight*KFactorQCDqqZZ_M","same");
		fnew->cd();
		for(int binx=0;binx<temp_zz_4e->GetXaxis()->GetNbins();binx++){
			double inttmp1 = temp_zz_4mu->Integral(binx+1,binx+1);
			double inttmp2 = temp_zz_4e->Integral(binx+1,binx+1);
			double inttmp3 = temp_zz_2e2mu->Integral(binx+1,binx+1);
			for(int biny=0;biny<temp_zz_4e->GetNbinsY();biny++){
			if(inttmp1!=0 )
				temp_zz_4mu->SetBinContent(binx+1, biny+1, temp_zz_4mu->GetBinContent(binx+1,biny+1)/inttmp1);
			if(inttmp2!=0 )
				temp_zz_4e->SetBinContent(binx+1, biny+1, temp_zz_4e->GetBinContent(binx+1,biny+1)/inttmp2);
			if(inttmp3!=0 )
				temp_zz_2e2mu->SetBinContent(binx+1, biny+1, temp_zz_2e2mu->GetBinContent(binx+1,biny+1)/inttmp3);
			}

		}
		temp_zz_4e->Smooth();
		temp_zz_4mu->Smooth();
		temp_zz_2e2mu->Smooth();
	        int nbinsfinal=455;
	        double xbinfinal[456]={
			110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,192,194,196,198,200,202,204,206,208,210,212,214,216,218,220,222,224,226,228,230,232,234,236,238,240,242,244,246,248,250,252,254,256,258,260,262,264,266,268,270,272,274,276,278,280,282,284,286,288,290,292,294,296,298,300,302,304,306,308,310,312,314,316,318,320,322,324,326,328,330,332,334,336,338,340,342,344,346,348,350,352,354,356,358,360,362,364,366,368,370,372,374,376,378,380,382,384,386,388,390,392,394,396,398,400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550,560,570,580,590,600,610,620,630,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800,810,820,830,840,850,860,870,880,890,900,910,920,930,940,950,960,970,980,990,1000,1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,1110,1120,1130,1140,1150,1160,1170,1180,1190,1200,1210,1220,1230,1240,1250,1260,1270,1280,1290,1300,1310,1320,1330,1340,1350,1360,1370,1380,1390,1400,1410,1420,1430,1440,1450,1460,1470,1480,1490,1500,1510,1520,1530,1540,1550,1560,1570,1580,1590,1600,1610,1620,1630,1640,1650,1660,1670,1680,1690,1700,1710,1720,1730,1740,1750,1760,1770,1780,1790,1800,1810,1820,1830,1840,1850,1860,1870,1880,1890,1900,1910,1920,1930,1940,1950,1960,1970,1980,1990,2000,2010,2020,2030,2040,2050,2060,2070,2080,2090,2100,2110,2120,2130,2140,2150,2160,2170,2180,2190,2200,2210,2220,2230,2240,2250,2260,2270,2280,2290,2300,2310,2320,2330,2340,2350,2360,2370,2380,2390,2400,2410,2420,2430,2440,2450,2460,2470,2480,2490,2500,2510,2520,2530,2540,2550,2560,2570,2580,2590,2600,2610,2620,2630,2640,2650,2660,2670,2680,2690,2700,2710,2720,2730,2740,2750,2760,2770,2780,2790,2800,2810,2820,2830,2840,2850,2860,2870,2880,2890,2900,2910,2920,2930,2940,2950,2960,2970,2980,2990,3000,3010,3020,3030,3040,3050,3060,3070,3080,3090,3100,3110,3120,3130,3140,3150,3160,3170,3180,3190,3200,3210,3220,3230,3240,3250,3260,3270,3280,3290,3300,3310,3320,3330,3340,3350,3360,3370,3380,3390,3400,3410,3420,3430,3440,3450,3460,3470,3480,3490,3500
	                };
		//TH2F *temp_zz_rebin= new TH2F("temp_zz"+newtreename[t],"",nbinsfinal,xbinfinal,30,0.,1.); 
		TH2F *temp_zz_rebin_4e= new TH2F("temp_zz_4e"+newtreename[t],"",1695,110,3500,30,0.,1.); 
		TH2F *temp_zz_rebin_4mu= new TH2F("temp_zz_4mu"+newtreename[t],"",1695,110,3500,30,0.,1.); 
		TH2F *temp_zz_rebin_2e2mu= new TH2F("temp_zz_2e2mu"+newtreename[t],"",1695,110,3500,30,0.,1.); 

		for(int binx=0;binx<temp_zz_rebin_4e->GetXaxis()->GetNbins();binx++){
			int b = temp_zz_4e->GetXaxis()->FindBin(temp_zz_rebin_4e->GetXaxis()->GetBinCenter(binx+1));
			for(int biny=0;biny<30;biny++){
				float bc_4e = temp_zz_4e->GetBinContent(b,biny+1);
				float bc_4mu = temp_zz_4mu->GetBinContent(b,biny+1);
				float bc_2e2mu = temp_zz_2e2mu->GetBinContent(b,biny+1);
				temp_zz_rebin_4e->SetBinContent(binx+1,biny+1,bc_4e);
				temp_zz_rebin_4mu->SetBinContent(binx+1,biny+1,bc_4mu);
				temp_zz_rebin_2e2mu->SetBinContent(binx+1,biny+1,bc_2e2mu);
			}
		}
		temp_zz_rebin_4e->Draw("colz");
		gPad->Print("zztemplate_4e.png");
		temp_zz_rebin_4mu->Draw("colz");
		gPad->Print("zztemplate_4mu.png");
		temp_zz_rebin_2e2mu->Draw("colz");
		gPad->Print("zztemplate_2e2mu.png");

		temp_zz_rebin_4e->Write();	
		temp_zz_rebin_4mu->Write();	
		temp_zz_rebin_2e2mu->Write();	
		tnew->Write();
	}
	fnew->Close();
}

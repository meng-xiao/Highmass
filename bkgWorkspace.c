#include "ProcessNormalization.cc+"
//#include "VerticalInterpPdf.cc+"

using namespace RooFit;

// highmass = 0 onshell
// highmass = 1 offshell
// highmass = 2 search 
//
//
// cate_vbf =0 ggH
// cate_vbf =1 VBF
// cate_vbf =2 RSE
// cate_vbf =3 TLE
//

void dosomething(TString chan="2e2mu", int cate_vbf=1, int highmass=0,int is2D=0){
	// double lumi = 9.235;
	double lumi = 35.8;
	// double lumi = 10;

	double parzx_all[5][6]={
		//		4e
		1,141.9,21.3,0,0,0,
		//4mu
		1,130.4,15.6,0,0,0,
		//2e2mu
		0.45,131.1,18.1,0.55,133.8,18.9,
		//4e RSE
		0.298623315141,238.804039513,38.5525989967,4.83114145422,-0.0097489713697,0,
		//2e2mu RSE
		0.000171310428786,209.221006175,26.5346636174,11.193766044,-0.00296426709129,0
	};
	double parzx_rse_all[2][8]{
		//4e
		2.18216e+01  ,2.99737e+02   ,-2.48222e+02   ,1.38710e+01   ,2.76729e-03   ,2.02956e+02   ,2.54751e+01   ,4.84386e-05 , 
			//2e2mu RSE
			2.19408e+00   ,-1.25292e+02  ,-8.72399e+01   ,1.50523e+01   ,8.52961e-03   ,2.40365e+02   ,4.26689e+01   ,6.19292e-05   
	};
	double parzx[6]={0.};
	double parzx_rse[8]={0.};

	if(cate_vbf!=2){
		if (chan=="4e") 	{
			for (int i=0;i<6;i++){parzx[i]=parzx_all[0][i];}

		} 
		if (chan=="4mu")  {
			for (int i=0;i<6;i++){parzx[i]=parzx_all[1][i];}
		}
		if (chan=="2e2mu") {
			for (int i=0;i<6;i++){parzx[i]=parzx_all[2][i];}
		} 
	}
	else{
		if (chan=="2e2mu") {
			for (int i=0;i<8;i++){parzx_rse[i]=parzx_rse_all[1][i];}
		} 
		if (chan=="4e") {
			for (int i=0;i<8;i++){parzx_rse[i]=parzx_rse_all[0][i];}
		} 
		}
		ofstream yields("yields.txt",std::fstream::app);



		RooWorkspace w("w");
	//	double recolowarr[3]={105.,104,110};
	//	double recohigharr[3]={140.,1604.,3500};
	//	const int reconbinsarr[3]={35,750,1695};
		double recolowarr[3]={110,110,300};
		double recohigharr[3]={3500.,3500.,3500.};
		const int reconbinsarr[3]={1695,1695,1600};
	//	const int reconbinsarr[3]={3390,3390,3200};


		const double low_reco=recolowarr[cate_vbf];
		const double high_reco=recohigharr[cate_vbf];
		const int nbins_reco=reconbinsarr[cate_vbf];



		RooRealVar* mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",125,low_reco,high_reco);
		if(cate_vbf==2){
			mreco->SetName("mreco_rse");
			mreco->SetTitle("mreco_rse");
		}

		RooRealVar* dbkg= new RooRealVar("dbkg_kin","Dbkg_{kin} ",0.5,0.,1.);
		RooPlot* frame= mreco->frame(low_reco,high_reco) ;
		if(highmass==2){
			mreco->setBins(nbins_reco);
			dbkg->setBins(30);
		}

		//qqzz pdf
		//TChain *tqqzz= new TChain("ZZTree/candTree");

		// simplify: for ggF, this is qqZZ bkg, for VBF, this should be ggZZ bkg
		// simplify: if no systematics wanted, comment out _up and _dn related stuff
		TString treename[4]={"","","_rse","_tle"};
		TChain *tqqzz= new TChain("SelectedTree"+treename[cate_vbf]);

		TH1F *hqqzz= new TH1F ("hqqzz","",nbins_reco,low_reco,high_reco);
		TH1F *hqqzz_up= new TH1F ("hqqzz_up","",nbins_reco,low_reco,high_reco);
		TH1F *hqqzz_dn= new TH1F ("hqqzz_dn","",nbins_reco,low_reco,high_reco);

		//tqqzz->Add("root://lxcms03://data3/Higgs/160225/ZZTo4l/ZZ4lAnalysis.root");
		//	if(cate_vbf!=2)
		tqqzz->Add("qqzz_80_Moriond.root");
		//	else
		//	tqqzz->Add("qqzz_80.root");

		int channum=1;
		if(chan=="2e2mu")
			channum=3;
		else if(chan=="4e")
			channum=2;
		else
			channum=1;

		//	tqqzz->Draw("mreco>>hqqzz",Form("weight*(chan==%d&&vbfcate==%d)",channum,cate_vbf));
		//	tqqzz->Draw("mreco>>hqqzz_up",Form("weight_up*(chan==%d&&vbfcate==%d)",channum,cate_vbf));
		//	tqqzz->Draw("mreco>>hqqzz_dn",Form("weight_dn*(chan==%d&&vbfcate==%d)",channum,cate_vbf));


		RooRealVar* wt= new RooRealVar("weight","wt",125,0.000,1000.);
		RooRealVar* wt_up= new RooRealVar("weight_up","wt_up",125,0.000,1000.);
		RooRealVar* wt_dn= new RooRealVar("weight_dn","wt_dn",125,0.000,1000.);


		float ZZMass;
		float weight, weight_up, weight_dn;
		int channel;
		tqqzz->SetBranchAddress("mreco",&ZZMass);
		

		if(cate_vbf!=1){
			tqqzz->SetBranchAddress("weight",&weight);
			tqqzz->SetBranchAddress("weight_up",&weight_up);
			tqqzz->SetBranchAddress("weight_dn",&weight_dn);
			tqqzz->Draw("mreco>>hqqzz",Form("weight*(chan==%d)",channum));
			tqqzz->Draw("mreco>>hqqzz_dn",Form("weight_dn*(chan==%d)",channum));
			tqqzz->Draw("mreco>>hqqzz_up",Form("weight_up*(chan==%d)",channum));
		}
		else{
			tqqzz->SetBranchAddress("weight_vbf",&weight);
			tqqzz->SetBranchAddress("weight_vbf_up",&weight_up);
			tqqzz->SetBranchAddress("weight_vbf_dn",&weight_dn);
			tqqzz->Draw("mreco>>hqqzz",Form("weight_vbf*(chan==%d)",channum));
			tqqzz->Draw("mreco>>hqqzz_dn",Form("weight_vbf_dn*(chan==%d)",channum));
			tqqzz->Draw("mreco>>hqqzz_up",Form("weight_vbf_up*(chan==%d)",channum));
		}

		tqqzz->SetBranchAddress("chan",&channel);


		if(cate_vbf!=2){
			cout<<chan<<"\t"<<cate_vbf<<"\t"<<hqqzz->Integral()*lumi*1000.<<endl;
			cout<<chan<<"\t"<<cate_vbf<<"\t"<<hqqzz_up->Integral()*lumi*1000.<<endl;
			cout<<chan<<"\t"<<cate_vbf<<"\t"<<hqqzz_dn->Integral()*lumi*1000.<<endl;
		}
		else{
			cout<<chan<<"\t"<<cate_vbf<<"\t"<<hqqzz->Integral()*lumi*1000.<<endl;  /// Remember why 1.6?
			cout<<chan<<"\t"<<cate_vbf<<"\t"<<hqqzz_up->Integral()*lumi*1000<<endl;
			cout<<chan<<"\t"<<cate_vbf<<"\t"<<hqqzz_dn->Integral()*lumi*1000<<endl;
		}

		TTree *cuttree =new TTree("SelectedTree","SelectedTree");
		TString mrecob = Form("%s/F",mreco->GetName());
		cuttree->Branch(mreco->GetName(),&ZZMass,mrecob);
		cuttree->Branch("weight",&weight,"weight/F");
		cuttree->Branch("weight_up",&weight_up,"weight_up/F");
		cuttree->Branch("weight_dn",&weight_dn,"weight_dn/F");

		for(int i =0;i<tqqzz->GetEntries();i++){
			tqqzz->GetEntry(i);
			if(channel==channum){
				cuttree->Fill();
			}
		}


		//Remember to revert
		double rho =1;
		if(highmass==0)
			rho=4;
		RooDataSet bkgdata ("bkgdata"+chan+Form("_%d",cate_vbf),"",cuttree,RooArgSet(*mreco,*wt),"weight");
		RooKeysPdf qqzzpdf_1d("bkg_qqzz_1d","",*mreco,bkgdata,RooKeysPdf::MirrorBoth,rho);

		TFile *ff = new TFile("qqzz_80_Moriond.root");
		TH2F *temp_zz=(TH2F*)ff->Get("temp_zz_"+chan);

		for(int bx=0;bx<temp_zz->GetNbinsX();bx++){
			for(int by=0;by<temp_zz->GetNbinsY();by++){
				if(temp_zz->GetBinContent(bx+1,by+1)==0){
					temp_zz->SetBinContent(bx+1,by+1,1.0e-10);
				}
			}
		}
		RooDataHist* template_sig= new RooDataHist("temp_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*mreco,*dbkg),temp_zz);
		RooHistPdf* pdf_2d_sig = new RooHistPdf("pdf_2d_sig_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*mreco,*dbkg),*template_sig);

		RooProdPdf *qqzzpdf_2d= new RooProdPdf("bkg_qqzz_2d_"+chan+Form("_%d",cate_vbf),"",qqzzpdf_1d,Conditional(*pdf_2d_sig,*dbkg));
		//RooProdPdf *qqzzpdf_2d= new RooProdPdf("bkg_qqzz_2d","",qqzzpdf_1d,*pdf_2d_sig);


		//    RooDataSet bkgdata_up("bkgdata_up"+chan+Form("_%d",cate_vbf),"",cuttree,RooArgSet(*mreco,*wt_up),"weight_up");
		//    RooKeysPdf qqzzpdf_up("bkg_qqzz_EWKscale_VVUp","",*mreco,bkgdata_up,RooKeysPdf::MirrorLeft,rho);
		//
		//    RooDataSet bkgdata_dn("bkgdata_dn"+chan+Form("_%d",cate_vbf),"",cuttree,RooArgSet(*mreco,*wt_dn),"weight_dn");
		//    RooKeysPdf qqzzpdf_dn("bkg_qqzz_EWKscale_VVDown","",*mreco,bkgdata_dn,RooKeysPdf::MirrorLeft);



		//    bkgdata.plotOn(frame,Binning(35));
		//    qqzzpdf.plotOn(frame);
		//    frame->Draw();
		

		/*	
			RooDataHist* qqzzhist= new RooDataHist("qqzzhist"+chan+Form("_%d",cate_vbf),"qqzzhist"+chan+Form("_%d",cate_vbf),RooArgSet(*mreco),hqqzz);
		//simplify: important, if VBF, change name "bkg_qqzz" to "bkg_ggzz"
		RooHistPdf* qqzzpdf= new RooHistPdf("bkg_qqzz","bkg_qqzz",RooArgSet(*mreco),*qqzzhist);

		RooDataHist* qqzzhist_up= new RooDataHist("qqzzhist_up"+chan+Form("_%d",cate_vbf),"qqzzhist_up"+chan+Form("_%d",cate_vbf),RooArgSet(*mreco),hqqzz_up);
		RooHistPdf* qqzzpdf_up= new RooHistPdf("bkg_qqzz_EWKscale_VVUp","bkg_qqzz_EWKscale_VVUp",RooArgSet(*mreco),*qqzzhist_up);
		RooDataHist* qqzzhist_dn= new RooDataHist("qqzzhist_dn"+chan+Form("_%d",cate_vbf),"qqzzhist_dn"+chan+Form("_%d",cate_vbf),RooArgSet(*mreco),hqqzz_dn);
		RooHistPdf* qqzzpdf_dn= new RooHistPdf("bkg_qqzz_EWKscale_VVDown","bkg_qqzz_EWKscale_VVDown",RooArgSet(*mreco),*qqzzhist_dn);
		*/	

		//  RooRealVar *ewkvar=new RooRealVar("EWKscale_VV","EWKscale_VV",-7,7);
		//  ProcessNormalization *int_qqzz;
		//	int_qq= new ProcessNormalization("bkg_qqzz_norm","bkg_qqzz_norm",1.);
		//	double hqqzz_dn_int = hqqzz_dn->Integral();
		//	double hqqzz_up_int = hqqzz_up->Integral();
		//	double hqqzz_norm_int = hqqzz->Integral();
		//  int_qq->addAsymmLogNormal(hqqzz_dn_int/hqqzz_norm_int, hqqzz_up_int/hqqzz_norm_int, *ewkvar);
		//	cout << hqqzz_dn_int/hqqzz_norm_int<<" "<< hqqzz_up_int/hqqzz_norm_int<<endl;

		//zjet pdf

		TF1 *fsum ;
		int parsize =6;
		double exp;
		fsum=new TF1("fsum"+chan,"landau( 0 )+ landau(3)",50,3500);
		if(cate_vbf==2){
			//fsum=new TF1("fsum"+chan,"landau( 0 )*(1 + exp( pol1(3)))",50,3500);
			parsize=8;
			fsum= new TF1("fsum"+chan, "[0]*TMath::Landau(x, [1], [2]) + [7]*(1 + TMath::Exp([3] - [4]*x))*TMath::Landau(x, [5], [6])",50,3500);
		}

		if(cate_vbf!=2)
			fsum->SetParameters(parzx);
		else
			fsum->SetParameters(parzx_rse);

		if(cate_vbf!=2){
			if(chan=="4e")
				exp=9.8*lumi/12.9;
			else if(chan=="2e2mu")
				exp = 20.4*lumi/12.9;
			else
				exp=10.2*lumi/12.9;
		}
		else{
			if(chan=="4e")
				exp =19; 
			else if(chan=="2e2mu")
				exp= 25; 
		}

		//cout <<fsum->GetExpFormula("p")<<endl;
		TString formuO= fsum->GetExpFormula("p");
		cout<<formuO<<endl;
		formuO.ReplaceAll("x",mreco->GetName());
		cout<<formuO<<endl;

		if(cate_vbf==1){
			yields<< chan<< " "<< cate_vbf<<" ZX "<< 0.0762466*exp*fsum->Integral(low_reco,high_reco)/fsum->Integral(50,3500)<<endl;
			cout<< chan<< " "<< cate_vbf<<" ZX "<< 0.0762466*exp*fsum->Integral(low_reco,high_reco)/fsum->Integral(50,3500)<<endl;
		}
		else if(cate_vbf==0){
			yields<< chan<< " "<< cate_vbf<<" ZX "<< (1-0.0762466)*exp*fsum->Integral(low_reco,high_reco)/fsum->Integral(50,3500)<<endl;
			cout<< chan<< " "<< cate_vbf<<" ZX "<< (1-0.0762466)*exp*fsum->Integral(low_reco,high_reco)/fsum->Integral(50,3500)<<endl;
		}
		else{
			cout<< chan<< " "<< cate_vbf<<" ZX "<< exp*fsum->Integral(300,high_reco)/fsum->Integral(50,3500)<<endl;
			yields<< chan<< " "<< cate_vbf<<" ZX "<< exp*fsum->Integral(300,high_reco)/fsum->Integral(50,3500)<<endl;
		}
		fsum->Draw();

		RooArgList* parlist_zx=new RooArgList("parlist_zx");
		const int parsize_const = parsize;
		RooConstVar *pars_zx[parsize_const];

		parlist_zx->add(*mreco);
		for (int pa=0; pa<parsize; pa++){ 
			if(cate_vbf!=2)
			pars_zx[pa]=new RooConstVar(Form("par_zx_w%d_%s_%d",pa,chan.Data(),cate_vbf),"",parzx[pa]);	
			else
			pars_zx[pa]=new RooConstVar(Form("par_zx_w%d_%s_%d",pa,chan.Data(),cate_vbf),"",parzx_rse[pa]);	
			cout << pars_zx[pa]->getVal()<<endl;

			parlist_zx->add(*pars_zx[pa]);
		}


		//RooGenericPdf 
		TString form="@1*TMath::Landau(@0,@2,@3)+@4*TMath::Landau(@0,@5,@6)";
		if(cate_vbf==2)
			form= "@1*TMath::Landau(@0, @2, @3) + @8*(1 + TMath::Exp(@4 - @5*@0))*TMath::Landau(@0, @6, @7)"; 
//		RooGenericPdf *zjetpdf = new RooGenericPdf("bkg_zjet","bkg_zjet","@1*TMath::Landau(@0,@2,@3)+@4*TMath::Landau(@0,@5,@6)",*parlist_zx);

		RooGenericPdf * zjetpdf_1d;
		zjetpdf_1d=new RooGenericPdf("bkg_zjet_1d","bkg_zjet_1d",form,*parlist_zx);

		TFile *fzjet = new TFile("zx.root");
		TH2F *temp_zjet=(TH2F*)fzjet->Get("zx_"+chan);

		RooDataHist* template_zx= new RooDataHist("temp_zx_"+chan+Form("_%d",cate_vbf),"",RooArgSet(*mreco,*dbkg),temp_zjet);
		RooHistPdf* pdf_2d_zx= new RooHistPdf("pdf_2d_zx"+chan+Form("_%d",cate_vbf),"",RooArgSet(*mreco,*dbkg),*template_zx);

		RooProdPdf *zjetpdf_2d= new RooProdPdf("bkg_zjet_2d_"+chan+Form("_%d",cate_vbf),"",*zjetpdf_1d,Conditional(*pdf_2d_zx,*dbkg));

		//zjetpdf_1d->plotOn(frame);
		//frame->Draw();
		//return;
		//	TCanvas *c2=new TCanvas("c2","",800,600);
		//	fsum->Draw("");
		//
		//	return;
		//data
		TChain *tdata = new TChain("SelectedTree"+treename[cate_vbf]);
		tdata->Add("data.root");
		TTree* reducetree= tdata->CopyTree(Form("chan==%d&&vbfcate==%d",channum,cate_vbf));
		RooDataSet* data_obs_1d= new RooDataSet("data_obs_1d","data_obs_1d",reducetree,*mreco);
		RooDataSet* data_obs_2d = new RooDataSet("data_obs_2d","data_obs_2d",reducetree,RooArgSet(*mreco,*dbkg));

// * Remember to revert
		if(is2D){
			data_obs_2d->SetNameTitle("data_obs","data_obs");
			qqzzpdf_2d->SetNameTitle("bkg_qqzz","bkg_qqzz");
			zjetpdf_2d->SetNameTitle("bkg_zjet","bkg_zjet");
			w.import(*qqzzpdf_2d);
			w.import(*zjetpdf_2d);
			w.import(*data_obs_2d);
		}
		else{
			data_obs_1d->SetNameTitle("data_obs","data_obs");
			qqzzpdf_1d.SetNameTitle("bkg_qqzz","bkg_qqzz");
			zjetpdf_1d->SetNameTitle("bkg_zjet","bkg_zjet");
			w.import(qqzzpdf_1d);
			w.import(*zjetpdf_1d);
			w.import(*data_obs_1d);
		}
		//		cout<<data_obs->sumEntries()<<endl;

		//	}

		//	qqzzpdf_2d.setNameTitle("bkg_qqzz","bkg_qqzz");
		//	w.import(qqzzpdf,RecycleConflictNodes());
		//	w.import(*qqzzpdf_2d,RecycleConflictNodes());
		//simplif no systematics, comment out begin
		//	w.import(qqzzpdf_up,RecycleConflictNodes());
		//	w.import(qqzzpdf_dn,RecycleConflictNodes());
		//	w.import(*int_qq,RecycleConflictNodes());
		//w.import(*zjetpdf,RecycleConflictNodes());


TFile *fwork ;
if(highmass==0)
	fwork= new TFile("workspace_2d_bkg_onshell/hzz4l_"+chan+Form("%dS_13TeV.input_func.root",cate_vbf),"recreate");
if(highmass==1)
	fwork= new TFile("workspace_2d_bkg_offshell/hzz4l_"+chan+Form("%dS_13TeV.input_func.root",cate_vbf),"recreate");
	if(highmass==2){
		if(is2D)
			fwork= new TFile("workspace_template2D_bkg_Moriond/hzz4l_"+chan+Form("%dS_13TeV.input_func.root",cate_vbf),"recreate");
		else
			fwork= new TFile("workspace_template1D_bkg_Moriond/hzz4l_"+chan+Form("%dS_13TeV.input_func.root",cate_vbf),"recreate");
	}
fwork->cd();
w.Write();
fwork->Close();

}
void bkgWorkspace(TString chan, int vbfcate, int highmass, int is2D){
	dosomething(chan,vbfcate,highmass,is2D);
}

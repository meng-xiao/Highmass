
#include "HZZ4L_RooHighmass_1D.cc+"
#include "HZZ4L_RooHighmass.cc+"
#include "HZZ2L2QRooPdfs.cc+"
#include "ProcessNormalization.cc+"
#include "VerticalInterpPdf.cc+"
#include "SplinePdf.cc+"
using namespace RooFit;

double dcbPara_resolved_2nd[6][11]={
	1000,1500,-1.10999, 0.00749822, -4.5463e-06, 4.04181, -0.00280538, 6.05505e-07, 3.98471, -0.00272925, 5.80127e-07, 
	1000,1500,0.0665556, 0.00364624, -2.39996e-06, 3.3315, -0.00288365, 8.64985e-07, 2.01928, -0.00113402, 2.81776e-07, 
	1000,1500,2.6701, 1.00153, -1.70739e-05, 44.0866, 0.918697, 2.43426e-05, -32.3866, 1.02066, -9.64553e-06, 
	1000,1500,15.9333, -0.0311407, 1.6586e-05, 3.3495, -0.0059731, 4.00217e-06, -6.32149, 0.00692155, -2.96044e-07, 
	1000,1500,1.62105, 0.00202981, 5.77829e-07, -2.19538, 0.00966267, -3.2386e-06, 5.18966, -0.00018405, 4.36387e-08, 
	1000,1500,-2.2513, 0.0303664, -2.88022e-07, -4.51134, 0.0348865, -2.54806e-06, 32.5499, -0.0145285, 1.39236e-05

//	1000,2000,0.120035, 0.00387941, -2.39538e-06, 2.84887, -0.00157826, 3.33457e-07, 13.9058, -0.0126352, 3.09769e-06, 
//	1000,2000,0.706763, 0.001371, -8.93629e-07, 2.36546, -0.00194639, 7.65068e-07, 1.25632, -0.000837254, 4.87783e-07, 
//	1000,2000,5.17318, -0.00807666, -2.12297e-06, 28.0614, -0.0538531, 2.07652e-05, -224.648, 0.198856, -4.24121e-05, 
//	1000,2000,13.0375, -0.0239911, 1.26205e-05, 6.93773, -0.0117916, 6.52068e-06, -84.7658, 0.079912, -1.64052e-05, 
//	1000,2000,-0.233604, 0.00931405, -4.26494e-06, 3.38117, 0.0020845, -6.50163e-07, -4.14694, 0.00961261, -2.53219e-06, 
//	1000,2000,5.2775, -0.000556645, 2.99758e-05, 3.00995, 0.00397846, 2.77082e-05, -240.441, 0.247429, -3.31545e-05

};


double dcbPara_merged_2nd[6][11]={
//	1000,2000,0.136929, 0.00258175, -1.40827e-06, 1.70736, -0.000559112, 1.62166e-07, 0.55962, 0.000588628, -1.24769e-07, 
//	1000,2000,-0.843714, 0.00835354, -4.7006e-06, 4.26184, -0.00185757, 4.04956e-07, 2.69182, -0.000287548, 1.24513e-08, 
//	1000,2000,12.9218, 0.00156746, 1.0281e-06, 13.8266, -0.00024214, 1.9329e-06, -2.04953, 0.015634, -2.03613e-06, 
//	1000,2000,5.92329, 0.0117041, -7.92576e-06, 14.7481, -0.00594552, 8.9905e-07, 16.5705, -0.00776792, 1.35465e-06, 
//	1000,2000,1.94528, -0.00449054, 4.75637e-06, -5.1446, 0.00968922, -2.33351e-06, 3.25791, 0.00128671, -2.32887e-07, 
//	1000,2000,10.8679, 0.0115416, 9.45775e-06, -2.51232, 0.038302, -3.92247e-06, 12.3401, 0.0234496, -2.09366e-07

	1000,1500,1.79404, -0.0011454, 6.50693e-07, 1.01854, 0.0004056, -1.24807e-07, 1.22566, 0.00012944, -3.27534e-08, 
	1000,1500,-1.40284, 0.00996887, -5.86499e-06, 5.91995, -0.00467671, 1.4578e-06, 2.76244, -0.000466697, 5.44607e-08, 
	1000,1500,-26.1949, 1.07972, -4.16358e-05, 24.2941, 0.978742, 8.8532e-06, -4.35604, 1.01694, -3.8802e-06, 
	1000,1500,10.7366, -0.00214211, -1.24418e-06, 14.652, -0.00997291, 2.67122e-06, 10.5603, -0.00451731, 8.52683e-07, 
	1000,1500,5.55773, -0.0146569, 1.11044e-05, -11.6884, 0.0198354, -6.14177e-06, 0.79069, 0.00319657, -5.95509e-07, 
	1000,1500,4.59191, 0.0348467, -5.74807e-06, 10.704, 0.0226225, 3.64019e-07, 6.08797, 0.0287772, -1.68755e-06
};

TString chanName [2]={"eeqq","mumuqq"};
TString jetName [2]={"Merged","Resolved"};
TString jetNameLower [2]={"merged","resolved"};
TString cateName [3]={"vbf-tagged","b-tagged","untagged"};
TString sampleName[2]={"ggH","VBF"};
TString dataName[2]={"ggH","qqH"};
TString workName[2]={"1D","2D"};
float btagU[2]={0.2,0.05};
float jesU[2]={0.1,0.05};

double lowEdge_reco[2]={500,500};
double lowEdge[2]={450,450};
//TString mzz_name[2]= {"zz2lJ_mass","zz2l2q_mass"};
TString mzz_name[2]= {"ZZMass","ZZMass"};
//float hxsec[2]={0.855*3.04,0.248*1.26844};
float bkgxsec[2]={3.19,0.642};
//float bkgxsec[2]={3.19,0.};
//float bkgxsec[2]={3.19,0.1};
float hxsec[2]={0.855*3.04,0.248*1.4568*1.15};
//float hxsec[2]={0.855*3.04,0.248};
//float hxsec[2]={0.855*3.04,0.};

//chan="eeqq , mumuqq",
//jetType="Merge, resolve"
//cate="vbf-tagged","b-tagged","untagged"
//sample="ggH, VBF"
TString hName[2]={
	"root://eoscms.cern.ch//store/user/xiaomeng/lxplusBackUp/rootfiles/rootfiles/125_sig.root",
	"root://eoscms.cern.ch//store/user/xiaomeng/lxplusBackUp/rootfiles/rootfiles/vbf/phantom_125_2e2mu.root"
};
void dosomething(int chan=0,int jetType=0, int cate=0, int sample=0,int is2D=1,double mass_d=450, double width_d=46.8){


	gStyle->SetOptStat(0);
	double sigfull, bkgfull,hfull,lumi,h_bkgfull;
	lumi = 35.8;
	sigfull = 2*0.0673*0.699*lumi;
	hfull = hxsec[sample]*lumi*0.699/0.0336*2;//3.04 LHC YR 
	bkgfull = bkgxsec[sample]*lumi*0.699/0.0336*2;


	double dcbPara_2nd[6][11];
	//	double effsig[11]={0.};
	double parzx[5]={0.};
	if(jetType==1){
		for (int i=0;i<6;i++){for(int j=0;j<11;j++){dcbPara_2nd[i][j]= dcbPara_resolved_2nd[i][j];}}
	}
	else{
		for (int i=0;i<6;i++){for(int j=0;j<11;j++){dcbPara_2nd[i][j]= dcbPara_merged_2nd[i][j];}}
	}
	//	for (int i=0;i<11;i++){effsig[i]=eff[2][i];}


	RooWorkspace w("w");

	// range of gen and reco
	const double low=lowEdge[jetType];
	const double high=3600;
	const int nbins=(high-low)/2;

	const double low_reco=lowEdge_reco[jetType];
	const double high_reco=3500;
	const int nbins_reco=(high_reco-low_reco);

	// define gen mass and reco mass
	RooRealVar* mzz = new RooRealVar("genmass","M_{ZZ} (GeV)",400,low,high);
	//RooRealVar* dbkg= new RooRealVar("Dspin0","Dbkg_{kin} ",0.5,0.,1.);
	RooRealVar* dbkg= new RooRealVar("Dbkg_0plus","Dbkg_{kin} ",0.5,0.,1.);
	RooRealVar* mreco= new RooRealVar(mzz_name[jetType],"M_{ZZ} (GeV)",400,low_reco,high_reco);

	RooPlot* frame= mreco->frame(Title("Z mass")) ;
	RooPlot* frame_dbkg= dbkg->frame(Title("")) ;
	RooRealVar* mean = new RooRealVar("mean","mean",450,100,3000);
	RooRealVar* sigma= new RooRealVar("sigma","sigma",50,0.,2000);
	RooRealVar *r= new RooRealVar("r","",1,0,1000);
	RooRealVar *vbfrac= new RooRealVar("fvbf","",0,0,1.);
	RooFormulaVar* signorm;
	if(sample==0)
		signorm =new RooFormulaVar("signorm_"+sampleName[sample],"@0*(1-@1)",RooArgSet(*r,*vbfrac));
	else{
		vbfrac->setVal(1);
		signorm =new RooFormulaVar("signorm_"+sampleName[sample],"@0*@1",RooArgSet(*r,*vbfrac));
	}

	//signal

	//TString floname = "ggh_input_spline.root";
	TString floname = "gghKfac_input_spline.root";
	TString spname ="sp_xsec_statnom";
	if(sample==1){
		//floname = "width_new_spline.root";
		floname = "width_new_spline_3500.root";
		spname = "br_2e2mu";
	}

	//	 TFile *flo=new TFile("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/Reduce/CMSSW_7_1_5/src/prepareInputs/"+floname,"read");
	TFile *flo=new TFile(floname,"read");
	TSpline3 *lo=(TSpline3*) flo->Get(spname);

	SplinePdf *pdf;

	TChain *da= new TChain("selectedTree");
	da->Add(Form("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_8_0_6/src/ZZMatrixElement/MELA/test/2l2q/%smH%.0f_gen.root",sampleName[sample].Data(),mass_d)); 

	TH1F *dag= new TH1F ("dag","",nbins_reco,low_reco,high_reco);
	TH1F *dag_noit= new TH1F ("dag_noit","",nbins_reco,low_reco,high_reco);

	TH1F *dag_dbkg= new TH1F ("dag_dbkg","",30,0,1.);
	TH1F *dag_dbkg_noit= new TH1F ("dag_dbkg_noit","",30,0.,1.);

	TString inta = "weightBSI";
	/*
	   da->Draw("ZZMass>>dag",Form("(jetType==%d)*cpsweight*weight*(%s)",jetType,inta.Data()));
	   da->Draw("ZZMass>>dag_noit",Form("(jetType==%d)*cpsweight*weight*(weightHBI+weightBS-weight2)",jetType));

	   da->Draw("Dbkg_0plus>>dag_dbkg",Form("(jetType==%d&&ZZMass>%f&&ZZMass<%f)*cpsweight*weight*(%s)",jetType,low_reco,high_reco,inta.Data()));
	   da->Draw("Dbkg_0plus>>dag_dbkg_noit",Form("(jetType==%d&&ZZMass>%f&&ZZMass<%f)*cpsweight*weight*(weightHBI+weightBS-weight2)",jetType,low_reco,high_reco));

	//da->Draw("ZZMass>>dag",Form("(jetType==%d)*cpsweight*weight*(weightBS-weight2)",jetType));
	//	da->Draw("ZZMass>>dag",Form("(jetType==%d&&lepFlav==%d&&tag==%d)*cpsweight*weight*weightBSI",jetType,chan,cate));
	//	da->Draw("ZZMass>>dag_noit",Form("(jetType==%d&&lepFlav==%d&&tag==%d)*cpsweight*weight*weightBS",jetType,chan,cate));
	//da->Draw("ZZMass>>dag","cpsweight*weight");
	RooDataHist* h160hist= new RooDataHist("h160hist","h160hist",RooArgSet(*mreco),dag);
	RooDataHist* h160hist_noit= new RooDataHist("h160hist_noit","h160hist_noit",RooArgSet(*mreco),dag_noit);

	RooDataHist* h160hist_dbkg= new RooDataHist("h160hist_dbkg","h160hist_dbkg",RooArgSet(*dbkg),dag_dbkg);
	RooDataHist* h160hist_dbkg_noit= new RooDataHist("h160hist_dbkg_noit","h160hist_dbkg_noit",RooArgSet(*dbkg),dag_dbkg_noit);

	TH1F *dag_gen= new TH1F ("dag_gen","",nbins,low,high);
	//da->Draw("GenHMass>>dag_gen",Form("(jetType==%d&&lepFlav==%d&&tag==%d)*cpsweight*weight*weightBSI",jetType,chan,cate));
	da->Draw("GenHMass>>dag_gen",Form("(jetType==%d)*cpsweight*weight*(%s)",jetType,inta.Data()));
	//	da->Draw("GenHMass>>dag_gen",Form("(jetType==%d)*cpsweight*weight*(weightBS-weight2)",jetType));
	RooDataHist* h160hist_gen= new RooDataHist("h160hist_gen","h160hist_gen",RooArgSet(*mzz),dag_gen);
	*/

	TFile *fkfactor = new TFile("Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
	//TSpline3  *ggZZ_kf = (TSpline3*)fkfactor->Get("sp_kfactor_Nominal");
	TGraph *ggZZ_kf = (TGraph*)fkfactor->Get("kfactor_Nominal");
	pdf=new SplinePdf("pdf_2e2mu","",*mzz,*mean,*sigma,*lo);

	RooPlot* frame_mzz= mzz->frame(Title("Z mass")) ;
	//h160hist_gen->plotOn(frame_mzz,Binning(nbins/2));
	//h160hist_noit->plotOn(frame,Binning(nbins_reco/10),MarkerColor(2));
	//h160hist->plotOn(frame,Binning(nbins_reco/10));

	//h160hist_dbkg_noit->plotOn(frame_dbkg,Binning(30),MarkerColor(2));
	//h160hist_dbkg->plotOn(frame_dbkg,Binning(30));

	//h125
	mean->setVal(125.);
	sigma->setVal(0.00407);

	TChain *t125 = new TChain ("SelectedTree");
	//t125->Add("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/highmass/rootfiles/125_sig.root");
	t125->Add(hName[sample]);
	double hall = t125->GetEntries("ZZMass!=0");
	double hsel = t125->GetEntries(Form("ZZMass>%f&&ZZMass<%f",low,high));
	double hnum = hsel/hall*hfull;

	TH1F *pdf125_hist=new TH1F("pdf125_hist","",nbins,low,high);
	double bwid = (high-low)/double(nbins);
	for(int k=0;k<nbins;k++){
		double x = pdf125_hist->GetBinCenter(k+1);
		mzz->setVal(x);
		double bc_125= pdf->getVal(*mzz)*hnum*bwid;
		if(x==125){	
			bc_125=0;
			for(int j =-50;j<50;j++){
				mzz->setVal(x+j*0.2*0.00407);
				bc_125+=pdf->getVal(*mzz)*hnum*0.2*0.00407;
			}
		}
		pdf125_hist->SetBinContent(k+1,bc_125);		
	}
	pdf125_hist->Draw();
	//return;


	// signal
	mean->setVal(mass_d);
	sigma->setVal(width_d);
	//pdf->plotOn(frame_mzz);
	//frame_mzz->Draw();
	//return;


	//bkg
	TChain *ggzz = new TChain("SelectedTree");
	ggzz->Add("/eos/cms/store/user/xiaomeng/lxplusBackUp/rootfiles/rootfiles/ggZZ_Bkg_xcheck.root");
	ggzz->Add("/eos/cms/store/user/xiaomeng/lxplusBackUp/rootfiles/rootfiles/ggZZ_ulascan_new.root");
	ggzz->Add("/eos/cms/store/user/xiaomeng/lxplusBackUp/rootfiles/rootfiles/ggZZ_ulascan.root");

	TChain *vbfbkg= new TChain("ZZTree/candTree");
	TChain *vbfbkg_failed= new TChain("ZZTree/candTree_failed");
	vbfbkg->Add("root://lxcms03//data3/Higgs/170222/VBFTo2e2muJJ_Contin_phantom128/ZZ4lAnalysis.root");
	vbfbkg_failed->Add("root://lxcms03//data3/Higgs/170222/VBFTo2e2muJJ_Contin_phantom128/ZZ4lAnalysis.root");

	//TChain *ggzz2= new TChain("selectedTree");
	//ggzz2->Add("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_8_0_6/src/ZZMatrixElement/MELA/test/2l2q/ggHmH1000_all.root");

	double bkgall = ggzz->GetEntries("");
	double bkgsel = ggzz->GetEntries(Form("ZZMass>%f&&ZZMass<%f",low,high));

	if(sample){
		bkgall = vbfbkg->GetEntries("")+vbfbkg_failed->GetEntries("");
		bkgsel = vbfbkg->GetEntries(Form("GenHMass>%f&&GenHMass<%f",low,high))+vbfbkg_failed->GetEntries(Form("GenHMass>%f&&GenHMass<%f",low,high));
	}
	double bkgnum = bkgsel/bkgall*bkgfull; 

	TH1F *genm= new TH1F ("genm","",nbins,low,high);
	if(!sample)
	ggzz->Draw("ZZMass>>genm",Form("ZZMass>%f",low));
	else{
	vbfbkg->Draw("GenHMass>>genm",Form("GenHMass>%f",low));
	vbfbkg_failed->Draw("GenHMass>>+genm",Form("GenHMass>%f",low));
	}
	//	ggzz2->Draw("GenHMass>>genm",Form("(GenHMass>%f)*weight2",low));


	genm->Scale(bkgnum/genm->Integral());

	//// resolution
	//TString formu_2nd=" (@0<@1)*(@3+@0*@4+@0*@0*@5 ) + ( @0>=@1 && @0<@2)*(@6+@0*@7+@0*@0*@8) + (@0>=@2)*(@9+@0*@10+@0*@0*@11)";	

	//RooArgList formuList_a1;
	//RooArgList formuList_a2;
	//RooArgList formuList_mean;
	//RooArgList formuList_n1;
	//RooArgList formuList_n2;
	//RooArgList formuList_sigma;
	//formuList_a1.add(*mzz);
	//formuList_a2.add(*mzz);
	//formuList_mean.add(*mzz);
	//formuList_n1.add(*mzz);
	//formuList_n2.add(*mzz);
	//formuList_sigma.add(*mzz);

	//RooConstVar* a1_p0_0_2nd[11] ;
	//RooConstVar* a2_p0_0_2nd[11] ;
	//RooConstVar* mean_p0_0_2nd[11] ;
	//RooConstVar* n1_p0_0_2nd[11] ;
	//RooConstVar* n2_p0_0_2nd[11] ;
	//RooConstVar* sigma_p0_0_2nd[11] ;
	//for(int i =0; i<11;i++){
	//	a1_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_a1_p0_0_2nd",jetName[jetType].Data(),i),Form("%s_%d_a1_p0_0_2nd",jetName[jetType].Data(),i),dcbPara_2nd[0][i]);
	//	a2_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_a2_p0_0_2nd",jetName[jetType].Data(),i),Form("%s_%d_a2_p0_0_2nd",jetName[jetType].Data(),i),dcbPara_2nd[1][i]);
	//	mean_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_mean_p0_0_2nd",jetName[jetType].Data(),i),Form("%s_%d_mean_p0_0_2nd",jetName[jetType].Data(),i),dcbPara_2nd[2][i]);
	//	n1_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_n1_p0_0_2nd",jetName[jetType].Data(),i),Form("%s_%d_n1_p0_0_2nd",jetName[jetType].Data(),i),dcbPara_2nd[3][i]);
	//	n2_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_n2_p0_0_2nd",jetName[jetType].Data(),i),Form("%s_%d_n2_p0_0_2nd",jetName[jetType].Data(),i),dcbPara_2nd[4][i]);
	//	sigma_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_sigma_p0_0_2nd",jetName[jetType].Data(),i),Form("%s_%d_sigma_p0_0_2nd",jetName[jetType].Data(),i),dcbPara_2nd[5][i]);

	//	formuList_a1.add(*a1_p0_0_2nd[i]);
	//	formuList_a2.add(*a2_p0_0_2nd[i]);
	//	formuList_mean.add(*mean_p0_0_2nd[i]);
	//	formuList_n1.add(*n1_p0_0_2nd[i]);
	//	formuList_n2.add(*n2_p0_0_2nd[i]);
	//	formuList_sigma.add(*sigma_p0_0_2nd[i]);
	//}

	//RooFormulaVar* a1_p0_2nd= new RooFormulaVar("a1_p0_2nd"+jetName[jetType],"a1_p0_2nd"+jetName[jetType],formu_2nd,formuList_a1);
	//RooFormulaVar* a2_p0_2nd= new RooFormulaVar("a2_p0_2nd"+jetName[jetType],"a2_p0_2nd"+jetName[jetType],formu_2nd,formuList_a2);
	//RooFormulaVar* mean_p0_2nd= new RooFormulaVar("mean_p0_2nd"+jetName[jetType],"mean_p0_2nd"+jetName[jetType],"("+formu_2nd+")-@0",formuList_mean);
//	//RooFormulaVar* mean_p0_2nd= new RooFormulaVar("mean_p0_2nd"+jetName[jetType],"mean_p0_2nd"+jetName[jetType],formu_2nd,formuList_mean);
	//RooFormulaVar* n1_p0_2nd= new RooFormulaVar("n1_p0_2nd"+jetName[jetType],"n1_p0_2nd"+jetName[jetType],formu_2nd,formuList_n1);
	//RooFormulaVar* n2_p0_2nd= new RooFormulaVar("n2_p0_2nd"+jetName[jetType],"n2_p0_2nd"+jetName[jetType],formu_2nd,formuList_n2);
	//RooFormulaVar* sigma_p0_2nd= new RooFormulaVar("sigma_p0_2nd"+jetName[jetType],"sigma_p0_2nd"+jetName[jetType],formu_2nd,formuList_sigma);

	//RooDoubleCB dcrReso_piece("dcrReso"+jetName[jetType],"Double Crystal ball ",*mreco,*mzz,*mean_p0_2nd,*sigma_p0_2nd,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);
	//RooDoubleCB dcrReso_piece_up("dcrReso"+jetName[jetType]+"_up","dcb up",*mreco,*mzz,*mean_p0_2nd,*sigma_p0_up,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);
	//RooDoubleCB dcrReso_piece_dn("dcrReso"+jetName[jetType]+"_dn","dcb dn",*mreco,*mzz,*mean_p0_2nd,*sigma_p0_dn,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);

	RooRealVar *mean_p0_2nd = new RooRealVar("mean_p0_2nd","",0,-500,500);
	RooRealVar *sigma_p0_2nd = new RooRealVar("sigma_p0_2nd","",0,-500,500);
	RooRealVar *n1_p0_2nd = new RooRealVar("n1_p0_2nd","",0,-500,500);
	RooRealVar *n2_p0_2nd = new RooRealVar("n2_p0_2nd","",0,-500,500);
	RooRealVar *a1_p0_2nd = new RooRealVar("a1_p0_2nd","",0,-500,500);
	RooRealVar *a2_p0_2nd = new RooRealVar("a2_p0_2nd","",0,-500,500);

	RooFormulaVar *sigma_p0_up = new RooFormulaVar("sigma_p0_up"+jetName[jetType],"","@0+0.1*@0",*sigma_p0_2nd);
	RooFormulaVar *sigma_p0_dn = new RooFormulaVar("sigma_p0_dn"+jetName[jetType],"","@0-0.1*@0",*sigma_p0_2nd);

	RooDoubleCB dcrReso_piece("dcrReso"+jetName[jetType],"Double Crystal ball ",*mreco,*mzz,*mean_p0_2nd,*sigma_p0_2nd,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);
	RooDoubleCB dcrReso_piece_up("dcrReso"+jetName[jetType]+"_up","dcb up",*mreco,*mzz,*mean_p0_2nd,*sigma_p0_up,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);
	RooDoubleCB dcrReso_piece_dn("dcrReso"+jetName[jetType]+"_dn","dcb dn",*mreco,*mzz,*mean_p0_2nd,*sigma_p0_dn,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);

	TFile *resoFile = new TFile("2l2qRes/2l2q_resolution_"+jetNameLower[jetType]+".root","read");
	TGraph *a1_gr = (TGraph*)resoFile->Get("a1"); 
	TGraph *a2_gr = (TGraph*)resoFile->Get("a2"); 
	TGraph *n1_gr = (TGraph*)resoFile->Get("n1"); 
	TGraph *n2_gr = (TGraph*)resoFile->Get("n2"); 
	TGraph *mean_gr = (TGraph*)resoFile->Get("mean"); 
	TGraph *sigma_gr = (TGraph*)resoFile->Get("sigma"); 

	mzz->setVal(mass_d);
	//cout<<mean_p0_2nd->getVal()<<endl;
	//cout<<sigma_p0_2nd->getVal()<<endl;
	//dcrReso_piece.plotOn(frame);
	//frame->Draw();
	//return;

	//TFile *ftemplate=new TFile("template_cond.root");

	TFile *ftemplate=new TFile(Form("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_8_0_6/src/ZZMatrixElement/MELA/test/template_2l2q_%s_%s_Moriond.root",jetName[jetType].Data(),sampleName[sample].Data()));
	TH2F *template_1= (TH2F*)ftemplate->Get("T_2D_1c");
	TH2F *template_2= (TH2F*)ftemplate->Get("T_2D_2c");
	TH2F *template_124= (TH2F*)ftemplate->Get("T_2D_124c");

	TH2F *template_1r= new TH2F("template1","",nbins,low,high,30,0.,1.); 
	TH2F *template_2r= new TH2F("template2","",nbins,low,high,30,0.,1.); 
	TH2F *template_124r= new TH2F("template124","",nbins,low,high,30,0.,1.); 

	TFile *ftemplate_dbkgup=new TFile(Form("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_8_0_6/src/ZZMatrixElement/MELA/test/template_2l2q_%s_up_%s_Moriond.root",jetName[jetType].Data(),sampleName[sample].Data()));
	TH2F *template_1_dbkgup= (TH2F*)ftemplate->Get("T_2D_1c");
	TH2F *template_2_dbkgup= (TH2F*)ftemplate->Get("T_2D_2c");
	TH2F *template_124_dbkgup= (TH2F*)ftemplate->Get("T_2D_124c");

	TH2F *template_1r_dbkgup= new TH2F("template1_dbkgup","",nbins,low,high,30,0.,1.); 
	TH2F *template_2r_dbkgup= new TH2F("template2_dbkgup","",nbins,low,high,30,0.,1.); 
	TH2F *template_124r_dbkgup= new TH2F("template124_dbkgup","",nbins,low,high,30,0.,1.); 

	TFile *ftemplate_dbkgdn=new TFile(Form("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_8_0_6/src/ZZMatrixElement/MELA/test/template_2l2q_%s_dn_%s_Moriond.root",jetName[jetType].Data(),sampleName[sample].Data()));
	TH2F *template_1_dbkgdn= (TH2F*)ftemplate->Get("T_2D_1c");
	TH2F *template_2_dbkgdn= (TH2F*)ftemplate->Get("T_2D_2c");
	TH2F *template_124_dbkgdn= (TH2F*)ftemplate->Get("T_2D_124c");

	TH2F *template_1r_dbkgdn= new TH2F("template1_dbkgdn","",nbins,low,high,30,0.,1.); 
	TH2F *template_2r_dbkgdn= new TH2F("template2_dbkgdn","",nbins,low,high,30,0.,1.); 
	TH2F *template_124r_dbkgdn= new TH2F("template124_dbkgdn","",nbins,low,high,30,0.,1.); 

	TH2F *conv_template_1= new TH2F("conv_template1","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_2= new TH2F("conv_template2","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_124= new TH2F("conv_template124","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 

	TH2F *conv_template_1_dbkgup= new TH2F("conv_template1_dbkgup","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_2_dbkgup= new TH2F("conv_template2_dbkgup","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_124_dbkgup= new TH2F("conv_template124_dbkgup","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 

	TH2F *conv_template_1_dbkgdn= new TH2F("conv_template1_dbkgdn","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_2_dbkgdn= new TH2F("conv_template2_dbkgdn","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_124_dbkgdn= new TH2F("conv_template124_dbkgdn","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 


	TH2F *conv_template_1_up= new TH2F("conv_template1_up","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_2_up= new TH2F("conv_template2_up","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_124_up= new TH2F("conv_template124_up","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 

	TH2F *conv_template_1_dn= new TH2F("conv_template1_dn","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_2_dn= new TH2F("conv_template2_dn","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_124_dn= new TH2F("conv_template124_dn","",nbins_reco/10,low_reco,high_reco,30,0.,1.); 


	TH2F *template_merge= new TH2F("template_merge","",nbins_reco,low_reco,high_reco,1,0.,1.); 

	//TFile *fphase_noweight=new TFile("fphase_ggH.root");
	//TGraph *cosfunc = (TGraph*)fphase_noweight->Get("cosspline");
	//TGraph *sinfunc = (TGraph*)fphase_noweight->Get("sinspline");
	TString phaseName[2]={"fphase_ggH.root","fgraph_vbf_phase.root"};
	TString cName[2]={"cosspline","cos"};
	TString sName[2]={"sinspline","sin"};
	TFile *fphase_noweight=new TFile(phaseName[sample]);
	TGraph *cosfunc = (TGraph*)fphase_noweight->Get(cName[sample]);
	TGraph *sinfunc = (TGraph*)fphase_noweight->Get(sName[sample]);

//	TFile *sigeff = new TFile ("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CreateDatacards_ZZ2l2q_ICHEP_spin0_v4/SigEff/2l2q_Efficiency_spin0_"+sampleName[sample]+".root");
	TFile *sigeff = new TFile ("V3/efficiency_2l2q_spin0_"+sampleName[sample]+".root");
	//	cout<<"spin0_ggH_"+chan+"_"+jetType+"_untagged"<<endl;
	TGraph *effsig=((TGraph*)(sigeff->Get("spin0_"+sampleName[sample]+"_"+chanName[chan]+"_"+jetName[jetType]+"_"+cateName[cate])));
	TGraph *effsig_btag=((TGraph*)(sigeff->Get("spin0_"+sampleName[sample]+"_"+chanName[chan]+"_"+jetName[jetType]+"_b-tagged")));
	TGraph *effsig_vbftag=((TGraph*)(sigeff->Get("spin0_"+sampleName[sample]+"_"+chanName[chan]+"_"+jetName[jetType]+"_vbf-tagged")));
	TGraph *effsig_untag=((TGraph*)(sigeff->Get("spin0_"+sampleName[sample]+"_"+chanName[chan]+"_"+jetName[jetType]+"_untagged")));

	TFile *bkgeff = new TFile ("eff_vbf_bkg.root");
	TGraph *effbkg = (TGraph*)(bkgeff->Get("eff_"+chanName[chan]+"_"+jetName[jetType]+"_"+cateName[cate])); 
	TGraph *effbkg_btag=((TGraph*)(bkgeff->Get("eff_"+chanName[chan]+"_"+jetName[jetType]+"_b-tagged")));
	TGraph *effbkg_vbftag=((TGraph*)(bkgeff->Get("eff_"+chanName[chan]+"_"+jetName[jetType]+"_vbf-tagged")));
	TGraph *effbkg_untag=((TGraph*)(bkgeff->Get("eff_"+chanName[chan]+"_"+jetName[jetType]+"_untagged")));

	//TF1 *effsig=((TF1*)((TGraph*)(sigeff->Get("spin0_"+sampleName[sample]+"_"+chanName[chan]+"_"+jetName[jetType]+"_"+cateName[cate])))->GetListOfFunctions()->First());

	float ratio_btag=0 ;
	float ratio_vbftag = 0;
	if(effsig_untag->Eval(mass_d)!=0)
		ratio_btag=-1*(effsig_btag->Eval(mass_d)/effsig_untag->Eval(mass_d));
	if((effsig_untag->Eval(mass_d)+effsig_btag->Eval(mass_d))!=0)
		ratio_vbftag=-1*(effsig_vbftag->Eval(mass_d)/(effsig_untag->Eval(mass_d)+effsig_btag->Eval(mass_d)));
	if(cate==1)
		ratio_btag = 1; 
	else if(cate==0)
		ratio_vbftag = 1; 




	//	TFile *bkgeff = new TFile ("eff_2l2q_bkg.root");
	//	TGraph *effbkg=(TGraph*)(bkgeff->Get("eff_"+jetName[jetType]));

	cout << "Gen "<<sigfull<<endl;
	TAxis *xax = template_1->GetXaxis();


	//	for (int i =0; i< template_1->GetNbinsX();i++){
	for (int i =0; i< nbins;i++){

		//double x = template_1->GetXaxis()->GetBinCenter(i+1);
		double x = low+ (high-low)/double(nbins)*i; 

		mzz->setVal(x);
		double bc =sigfull*(fabs((mass_d-x))<0.1);
		if(width_d>0.5)
			bc = pdf->getVal(*mzz)*sigfull*bwid;
		double bc_bkg= genm->GetBinContent(i+1);
		double bc_125= pdf125_hist->GetBinContent(i+1);

		//	double effval = effsig1->Eval(x)+ effsig2->Eval(x)+effsig3->Eval(x)+effsig4->Eval(x)+effsig5->Eval(x)+effsig6->Eval(x);
		double effval = effsig->Eval(x); 
		double effbkg_val = effbkg->Eval(x);
		double kfac = ggZZ_kf->Eval(x);

		//important
		double fa_sig = bc*effval; // effval effcate
		//		double fa_sig = dag_gen->GetBinContent(i); 
		double fa_bkg = bc_bkg*effval; // effcate*kfac

		//if(sample)
		//	fa_bkg =bc_bkg*effbkg_val; // effcate*kfac
		
		//		double fa_bkg = bc_bkg; 
		//		double fa_sig = bc; // effval effcate
		//		double fa_125 = bc_125*effval; // effcate
		double fa_125 = bc_125*effval; // effcate
		if(sample==1){
			//			fa_bkg = fa_bkg/kfac;
			//fa_sig = fa_sig/x;
		}
		else{
			//			fa_sig = fa_sig/kfac;
			fa_bkg = fa_bkg*kfac;
			//			fa_125 = fa_125/kfac;
		}

		double a= (x*x-mean->getVal()*mean->getVal());
		double b = mean->getVal()*sigma->getVal();

		double cossig = a/TMath::Sqrt(a*a+b*b);
		double sinsig = b/TMath::Sqrt(a*a+b*b);

		double a125= (x*x-125.*125.);
		double b125 = 125*0.00407;

		double cossig125 = a125/TMath::Sqrt(a125*a125+b125*b125);
		double sinsig125 = b125/TMath::Sqrt(a125*a125+b125*b125);

		double maphase = x;
		if(x>1600){
			if(sample==0)
			maphase=1600.;
			else if(sample==1&&x>3000)
			maphase=3000.;

		}
		double cosfuncv=cosfunc->Eval(maphase);
		double sinfuncv=sinfunc->Eval(maphase);
		if(sample){
			cosfuncv = cosfuncv/1.76*2; 
			sinfuncv = sinfuncv/1.76*2; 
		}

		double sigbkgsqrt = TMath::Sqrt(fa_sig*fa_bkg); 
		double hbkgsqrt = TMath::Sqrt(fa_125*fa_bkg); 
		double sighsqrt = TMath::Sqrt(fa_sig*fa_125); 

		double inter_sig_bkg = 1.76*sigbkgsqrt*(cosfuncv*cossig+sinfuncv*sinsig);
		double inter_125_bkg = 1.76*hbkgsqrt* (cosfuncv*cossig125+sinfuncv*sinsig125); 
		double inter_sig_125 = 2*sighsqrt * (cossig125*cossig-sinsig125*sinsig); 


		int binN = xax->FindBin(x);
		for(int j = 0;j<template_1->GetNbinsY();j++){
			double cont1, cont2, cont124;
			double cont1_dbkgup, cont2_dbkgup, cont124_dbkgup;
			double cont1_dbkgdn, cont2_dbkgdn, cont124_dbkgdn;

			cont1 = template_1->GetBinContent(binN,j+1);
			cont2 = template_2->GetBinContent(binN,j+1);
			cont124 = template_124->GetBinContent(binN,j+1);

			cont1_dbkgup = template_1_dbkgup->GetBinContent(binN,j+1);
			cont2_dbkgup = template_2_dbkgup->GetBinContent(binN,j+1);
			cont124_dbkgup= template_124_dbkgup->GetBinContent(binN,j+1);

			cont1_dbkgdn = template_1_dbkgdn->GetBinContent(binN,j+1);
			cont2_dbkgdn = template_2_dbkgdn->GetBinContent(binN,j+1);
			cont124_dbkgdn= template_124_dbkgdn->GetBinContent(binN,j+1);

			//cont1 = 1.; 
			//cont2 = 1.; 
			//cont124 = 1.; 

			template_1r->SetBinContent(i+1,j+1,cont1*fa_sig);
			template_2r->SetBinContent(i+1,j+1,cont2*fa_bkg+cont1*fa_125+cont124*inter_125_bkg);

			if(width_d>0.5)
			template_124r->SetBinContent(i+1,j+1,cont124*inter_sig_bkg+cont1*inter_sig_125);
			else
			template_124r->SetBinContent(i+1,j+1,0);

			template_1r_dbkgup->SetBinContent(i+1,j+1,cont1_dbkgup*fa_sig);
			template_2r_dbkgup->SetBinContent(i+1,j+1,cont2_dbkgup*fa_bkg+cont1_dbkgup*fa_125+cont124_dbkgup*inter_125_bkg);
			if(width_d>0.5)
			template_124r_dbkgup->SetBinContent(i+1,j+1,cont124_dbkgup*inter_sig_bkg+cont1_dbkgup*inter_sig_125);
else
			template_124r_dbkgup->SetBinContent(i+1,j+1,0);

			template_1r_dbkgdn->SetBinContent(i+1,j+1,cont1_dbkgdn*fa_sig);
			template_2r_dbkgdn->SetBinContent(i+1,j+1,cont2_dbkgdn*fa_bkg+cont1_dbkgdn*fa_125+cont124_dbkgdn*inter_125_bkg);
			if(width_d>0.5)
			template_124r_dbkgdn->SetBinContent(i+1,j+1,cont124_dbkgdn*inter_sig_bkg+cont1_dbkgdn*inter_sig_125);
else
			template_124r_dbkgdn->SetBinContent(i+1,j+1,0);
		}
	}
	cout << "Gen sig"<<template_1r->Integral()<<endl;
	cout << "Gen bkg "<<template_2r->Integral()<<endl;
	cout << "Gen int"<<template_124r->Integral()<<endl;
	//return;

	TH1F *template_1_proj = (TH1F*)template_1r->ProjectionX();
	TH1F *template_2_proj = (TH1F*)template_2r->ProjectionX();
	TH1F *template_124_proj= (TH1F*)template_124r->ProjectionX();

	//template_2_proj->Draw();
	//genm->Draw("same");	
	//return;


	RooDataHist *pdfsig_hist_hist_gen= new RooDataHist("pdfsig_hist_hist_gen"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mzz),template_1_proj);
	RooHistFunc *pdfsig_hist_func_gen= new RooHistFunc("pdfsig_hist_func_gen"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mzz),*pdfsig_hist_hist_gen);

	RooDataHist* hgenm_gen= new RooDataHist("hgenm_gen"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mzz),template_2_proj);
	RooHistFunc* hgenm_func_gen= new RooHistFunc("hgenm_func_gen"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mzz),*hgenm_gen);

	RooDataHist *inter_hist_gen= new RooDataHist("inter_hist_gen"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mzz),template_124_proj);
	RooHistFunc *inter_func_gen= new RooHistFunc("inter_func_gen"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mzz),*inter_hist_gen);

	HZZ4L_RooHighmass_1D *bsi_hist_gen=new HZZ4L_RooHighmass_1D("bsi_hist_gen","bsi_hist_gen",*mzz,*signorm,RooArgList(*pdfsig_hist_func_gen,*hgenm_func_gen,*inter_func_gen)); 
	bsi_hist_gen->plotOn(frame_mzz);
	//	frame_mzz->SetMaximum(6000);
	frame_mzz->Draw();

	//	gPad->Print(Form("2l2q_closure/m%.0f_%s_%s_%s_%s_gen.png",mass_d,chanName[chan].Data(),jetName[jetType].Data(),cateName[cate].Data(),sampleName[sample].Data()));
	//	gPad->Print(Form("2l2q_closure/m%.0f_%s_%s_%s_%s_gen.pdf",mass_d,chanName[chan].Data(),jetName[jetType].Data(),cateName[cate].Data(),sampleName[sample].Data()));
	//	gPad->Print(Form("2l2q_closure/m%.0f_%s_BSI_%s_gen.png",mass_d,jetName[jetType].Data(),sampleName[sample].Data()));
	//	gPad->Print(Form("2l2q_closure/m%.0f_%s_BSI_%s_gen.pdf",mass_d,jetName[jetType].Data(),sampleName[sample].Data()));
	//	gPad->Print(Form("2l2q_closure/m%.0f_%s_S_%s_gen.png",mass_d,jetName[jetType].Data(),sampleName[sample].Data()));
	//	gPad->Print(Form("2l2q_closure/m%.0f_%s_S_%s_gen.pdf",mass_d,jetName[jetType].Data(),sampleName[sample].Data()));

	//	for(int i =0;i<nbins;i++){
	//		double ratio = dag_gen->GetBinContent(i+1)/template_2_proj->GetBinContent(i+1);
	//		cout<< dag_gen->GetBinCenter(i+1)<<"\t"<< ratio<<"\t"<<template_2_proj->GetBinContent(i+1)<<endl;
	//	}
	//	return;
	//template_1->Draw("colz");
	//return;
	//conv_template_1->Draw("colz");
	//		return;

	for(int m =0;m< nbins; m++){
		//			double genval = template_1->GetXaxis()->GetBinCenter(m+1);
		double genval = low+(high-low)/double(nbins)*m; 
		mzz->setVal(genval);
		a1_p0_2nd->setVal(a1_gr->Eval(genval));
		a2_p0_2nd->setVal(a2_gr->Eval(genval));
		n1_p0_2nd->setVal(n1_gr->Eval(genval));
		n2_p0_2nd->setVal(n2_gr->Eval(genval));
		mean_p0_2nd->setVal(mean_gr->Eval(genval));
		sigma_p0_2nd->setVal(sigma_gr->Eval(genval));

		int low_bound = genval + mean_p0_2nd->getVal()-sigma_p0_2nd->getVal()*30 ; 
		int high_bound = genval + mean_p0_2nd->getVal()+sigma_p0_2nd->getVal()*30 ; 
		if(low_bound<low_reco)
			low_bound=low_reco;
		if(high_bound>high_reco)
			high_bound=high_reco;
		//	for(int k = 0;k<nbins_reco;k++){
		//		double recoval = low_reco+(high_reco-low_reco)/double(nbins_reco)*k; 
		

		for(int k = low_bound;k<high_bound;k++){
			double recoval = k; 
			mreco->setVal(recoval);
			double reso = dcrReso_piece.getVal(*mreco); 
			//cout<<genval<<"\t"<< recoval<<"\t"<<reso<< "\t"<< a1_p0_2nd->getVal()<< "\t"<<a2_p0_2nd->getVal()<<"\t"<<(recoval-genval-mean_p0_2nd->getVal())/sigma_p0_2nd->getVal()<< endl;
			if(reso<1.E-5)
				continue;

			double reso_up = dcrReso_piece_up.getVal(*mreco); 
			double reso_dn = dcrReso_piece_dn.getVal(*mreco); 
			for(int j =0;j< template_1r->GetNbinsY(); j++){
				double dbkg_val = template_1r->GetYaxis()->GetBinCenter(j+1);
//				cout<< recoval<<"\t"<<genval<<"\t"<<template_1->GetBinContent(m+1,j+1)<<"\t"<<reso<<endl;
				conv_template_1->Fill(recoval, dbkg_val,template_1r->GetBinContent(m+1,j+1)*reso);
				conv_template_2->Fill(recoval, dbkg_val,template_2r->GetBinContent(m+1,j+1)*reso);
				conv_template_124->Fill(recoval, dbkg_val,template_124r->GetBinContent(m+1,j+1)*reso);

				conv_template_1_up->Fill(recoval, dbkg_val,template_1r->GetBinContent(m+1,j+1)*reso_up);
				conv_template_2_up->Fill(recoval, dbkg_val,template_2r->GetBinContent(m+1,j+1)*reso_up);
				conv_template_124_up->Fill(recoval, dbkg_val,template_124r->GetBinContent(m+1,j+1)*reso_up);

				conv_template_1_dn->Fill(recoval, dbkg_val,template_1r->GetBinContent(m+1,j+1)*reso_dn);
				conv_template_2_dn->Fill(recoval, dbkg_val,template_2r->GetBinContent(m+1,j+1)*reso_dn);
				conv_template_124_dn->Fill(recoval, dbkg_val,template_124r->GetBinContent(m+1,j+1)*reso_dn);

				conv_template_1_dbkgup->Fill(recoval, dbkg_val,template_1r_dbkgup->GetBinContent(m+1,j+1)*reso);
				conv_template_2_dbkgup->Fill(recoval, dbkg_val,template_2r_dbkgup->GetBinContent(m+1,j+1)*reso);
				conv_template_124_dbkgup->Fill(recoval, dbkg_val,template_124r_dbkgup->GetBinContent(m+1,j+1)*reso);

				conv_template_1_dbkgdn->Fill(recoval, dbkg_val,template_1r_dbkgdn->GetBinContent(m+1,j+1)*reso);
				conv_template_2_dbkgdn->Fill(recoval, dbkg_val,template_2r_dbkgdn->GetBinContent(m+1,j+1)*reso);
				conv_template_124_dbkgdn->Fill(recoval, dbkg_val,template_124r_dbkgdn->GetBinContent(m+1,j+1)*reso);

			}
		}
	}

	//cout << "reco sig "<< conv_template_1->Integral()<<endl;
	//cout << "reco bkg "<< conv_template_2->Integral()<<endl;
	TH1F *conv_template_1_proj = (TH1F*)conv_template_1->ProjectionX();
	TH1F *conv_template_2_proj = (TH1F*)conv_template_2->ProjectionX();
	TH1F *conv_template_124_proj= (TH1F*)conv_template_124->ProjectionX();

	TH1F *conv_template_1_proj_up = (TH1F*)conv_template_1_up->ProjectionX();
	TH1F *conv_template_2_proj_up = (TH1F*)conv_template_2_up->ProjectionX();
	TH1F *conv_template_124_proj_up= (TH1F*)conv_template_124_up->ProjectionX();

	TH1F *conv_template_1_proj_dn = (TH1F*)conv_template_1_dn->ProjectionX();
	TH1F *conv_template_2_proj_dn = (TH1F*)conv_template_2_dn->ProjectionX();
	TH1F *conv_template_124_proj_dn= (TH1F*)conv_template_124_dn->ProjectionX();
	//	conv_template_1->Draw("colz");
	if(is2D){
		RooDataHist *pdfsig_hist_hist= new RooDataHist("pdfsig_hist_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mreco,*dbkg),conv_template_1);
		RooHistFunc *pdfsig_hist_func= new RooHistFunc("pdfsig_hist_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mreco,*dbkg),*pdfsig_hist_hist);

		RooDataHist* hgenm= new RooDataHist("hgenm"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mreco,*dbkg),conv_template_2);
		RooHistFunc* hgenm_func= new RooHistFunc("hgenm_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mreco,*dbkg),*hgenm);

		RooDataHist *inter_hist= new RooDataHist("inter_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mreco,*dbkg),conv_template_124);
		RooHistFunc *inter_func= new RooHistFunc("inter_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mreco,*dbkg),*inter_hist);

		HZZ4L_RooHighmass *bsi_hist=new HZZ4L_RooHighmass("bsi_hist","bsi_hist",*mreco,*dbkg,*signorm,RooArgList(*pdfsig_hist_func,*hgenm_func,*inter_func)); 
		bsi_hist->SetNameTitle(dataName[sample],dataName[sample]);

		// Resolution up
		RooDataHist *pdfsig_hist_hist_up= new RooDataHist("pdfsig_hist_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResUp","",RooArgSet(*mreco,*dbkg),conv_template_1_up);
		RooHistFunc *pdfsig_hist_func_up= new RooHistFunc("pdfsig_hist_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResUp","",RooArgSet(*mreco,*dbkg),*pdfsig_hist_hist_up);

		RooDataHist* hgenm_up= new RooDataHist("hgenm"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResUp","",RooArgSet(*mreco,*dbkg),conv_template_2_up);
		RooHistFunc* hgenm_func_up= new RooHistFunc("hgenm_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResUp","",RooArgSet(*mreco,*dbkg),*hgenm_up);

		RooDataHist *inter_hist_up= new RooDataHist("inter_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResUp","",RooArgSet(*mreco,*dbkg),conv_template_124_up);
		RooHistFunc *inter_func_up= new RooHistFunc("inter_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResUp","",RooArgSet(*mreco,*dbkg),*inter_hist_up);

		HZZ4L_RooHighmass *bsi_hist_up=new HZZ4L_RooHighmass("bsi_hist_2l2qJetResUp","bsi_hist_2l2qJetResUp",*mreco,*dbkg,*signorm,RooArgList(*pdfsig_hist_func_up,*hgenm_func_up,*inter_func_up)); 
		bsi_hist_up->SetNameTitle(dataName[sample]+"_CMS_2l2q"+jetName[jetType]+"ResUp",dataName[sample]+"_CMS_2l2q"+jetName[jetType]+"ResUp");

		// Resolution dn 
		RooDataHist *pdfsig_hist_hist_dn= new RooDataHist("pdfsig_hist_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResDown","",RooArgSet(*mreco,*dbkg),conv_template_1_dn);
		RooHistFunc *pdfsig_hist_func_dn= new RooHistFunc("pdfsig_hist_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResDown","",RooArgSet(*mreco,*dbkg),*pdfsig_hist_hist_dn);

		RooDataHist* hgenm_dn= new RooDataHist("hgenm"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResDown","",RooArgSet(*mreco,*dbkg),conv_template_2_dn);
		RooHistFunc* hgenm_func_dn= new RooHistFunc("hgenm_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResDown","",RooArgSet(*mreco,*dbkg),*hgenm_dn);

		RooDataHist *inter_hist_dn= new RooDataHist("inter_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResDown","",RooArgSet(*mreco,*dbkg),conv_template_124_dn);
		RooHistFunc *inter_func_dn= new RooHistFunc("inter_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResDown","",RooArgSet(*mreco,*dbkg),*inter_hist_dn);

		HZZ4L_RooHighmass *bsi_hist_dn=new HZZ4L_RooHighmass("bsi_hist_2l2qJetResDown","bsi_hist_2l2qJetResDown",*mreco,*dbkg,*signorm,RooArgList(*pdfsig_hist_func_dn,*hgenm_func_dn,*inter_func_dn)); 
		bsi_hist_dn->SetNameTitle(dataName[sample]+"_CMS_2l2q"+jetName[jetType]+"ResDown",dataName[sample]+"_CMS_2l2q"+jetName[jetType]+"ResDown");

		// dbkg up
		RooDataHist *pdfsig_hist_hist_dbkgup= new RooDataHist("pdfsig_hist_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qsigMELAUp","",RooArgSet(*mreco,*dbkg),conv_template_1_dbkgup);
		RooHistFunc *pdfsig_hist_func_dbkgup= new RooHistFunc("pdfsig_hist_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qsigMELAUp","",RooArgSet(*mreco,*dbkg),*pdfsig_hist_hist_dbkgup);

		RooDataHist* hgenm_dbkgup= new RooDataHist("hgenm"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qsigMELAUp","",RooArgSet(*mreco,*dbkg),conv_template_2_dbkgup);
		RooHistFunc* hgenm_func_dbkgup= new RooHistFunc("hgenm_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qsigMELAUp","",RooArgSet(*mreco,*dbkg),*hgenm_dbkgup);

		RooDataHist *inter_hist_dbkgup= new RooDataHist("inter_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qsigMELAUp","",RooArgSet(*mreco,*dbkg),conv_template_124_dbkgup);
		RooHistFunc *inter_func_dbkgup= new RooHistFunc("inter_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qsigMELAUp","",RooArgSet(*mreco,*dbkg),*inter_hist_dbkgup);

		HZZ4L_RooHighmass *bsi_hist_dbkgup=new HZZ4L_RooHighmass("bsi_hist_2l2qsigMELAUp","bsi_hist_2l2qsigMELAUp",*mreco,*dbkg,*signorm,RooArgList(*pdfsig_hist_func_dbkgup,*hgenm_func_dbkgup,*inter_func_dbkgup)); 
		bsi_hist_dbkgup->SetNameTitle(dataName[sample]+"_CMS_2l2q_"+jetName[jetType]+"_sigMELAUp",dataName[sample]+"_CMS_2l2q_"+jetName[jetType]+"_sigMELAUp");

		//dbkg dn
		RooDataHist *pdfsig_hist_hist_dbkgdn= new RooDataHist("pdfsig_hist_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qsigMELADn","",RooArgSet(*mreco,*dbkg),conv_template_1_dbkgdn);
		RooHistFunc *pdfsig_hist_func_dbkgdn= new RooHistFunc("pdfsig_hist_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qsigMELADn","",RooArgSet(*mreco,*dbkg),*pdfsig_hist_hist_dbkgdn);

		RooDataHist* hgenm_dbkgdn= new RooDataHist("hgenm"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qsigMELADn","",RooArgSet(*mreco,*dbkg),conv_template_2_dbkgdn);
		RooHistFunc* hgenm_func_dbkgdn= new RooHistFunc("hgenm_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qsigMELADn","",RooArgSet(*mreco,*dbkg),*hgenm_dbkgdn);

		RooDataHist *inter_hist_dbkgdn= new RooDataHist("inter_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qsigMELADn","",RooArgSet(*mreco,*dbkg),conv_template_124_dbkgdn);
		RooHistFunc *inter_func_dbkgdn= new RooHistFunc("inter_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qsigMELADn","",RooArgSet(*mreco,*dbkg),*inter_hist_dbkgdn);

		HZZ4L_RooHighmass *bsi_hist_dbkgdn=new HZZ4L_RooHighmass("bsi_hist_2l2qsigMELADown","bsi_hist_2l2qsigMELADown",*mreco,*dbkg,*signorm,RooArgList(*pdfsig_hist_func_dbkgdn,*hgenm_func_dbkgdn,*inter_func_dbkgdn)); 
		bsi_hist_dbkgdn->SetNameTitle(dataName[sample]+"_CMS_2l2q_"+jetName[jetType]+"_sigMELADown",dataName[sample]+"_CMS_2l2q_"+jetName[jetType]+"_sigMELADown");

		w.import(*bsi_hist,RecycleConflictNodes());
		w.import(*bsi_hist_up,RecycleConflictNodes());
		w.import(*bsi_hist_dn,RecycleConflictNodes());
		w.import(*bsi_hist_dbkgup,RecycleConflictNodes());
		w.import(*bsi_hist_dbkgdn,RecycleConflictNodes());
	}
	else{
		RooDataHist *pdfsig_hist_hist= new RooDataHist("pdfsig_hist_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mreco),conv_template_1_proj);
		RooHistFunc *pdfsig_hist_func= new RooHistFunc("pdfsig_hist_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mreco),*pdfsig_hist_hist);

		RooDataHist* hgenm= new RooDataHist("hgenm"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mreco),conv_template_2_proj);
		RooHistFunc* hgenm_func= new RooHistFunc("hgenm_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mreco),*hgenm);

		RooDataHist *inter_hist= new RooDataHist("inter_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mreco),conv_template_124_proj);
		RooHistFunc *inter_func= new RooHistFunc("inter_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"",RooArgSet(*mreco),*inter_hist);

		HZZ4L_RooHighmass_1D *bsi_hist=new HZZ4L_RooHighmass_1D("bsi_hist","bsi_hist",*mreco,*signorm,RooArgList(*pdfsig_hist_func,*hgenm_func,*inter_func)); 
		bsi_hist->SetNameTitle(dataName[sample],dataName[sample]);

		// Resolution up 
		RooDataHist *pdfsig_hist_hist_up= new RooDataHist("pdfsig_hist_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResUp","",RooArgSet(*mreco),conv_template_1_proj_up);
		RooHistFunc *pdfsig_hist_func_up= new RooHistFunc("pdfsig_hist_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResUp","",RooArgSet(*mreco),*pdfsig_hist_hist_up);

		RooDataHist* hgenm_up= new RooDataHist("hgenm"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResUp","",RooArgSet(*mreco),conv_template_2_proj_up);
		RooHistFunc* hgenm_func_up= new RooHistFunc("hgenm_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResUp","",RooArgSet(*mreco),*hgenm_up);

		RooDataHist *inter_hist_up= new RooDataHist("inter_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResUp","",RooArgSet(*mreco),conv_template_124_proj_up);
		RooHistFunc *inter_func_up= new RooHistFunc("inter_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResUp","",RooArgSet(*mreco),*inter_hist_up);

		HZZ4L_RooHighmass_1D *bsi_hist_up=new HZZ4L_RooHighmass_1D("bsi_hist","bsi_hist",*mreco,*signorm,RooArgList(*pdfsig_hist_func_up,*hgenm_func_up,*inter_func_up)); 
		bsi_hist_up->SetNameTitle(dataName[sample]+"_CMS_2l2q"+jetName[jetType]+"ResUp",dataName[sample]+"_CMS_2l2q"+jetName[jetType]+"ResUp");

		// Resolution dn 
		RooDataHist *pdfsig_hist_hist_dn= new RooDataHist("pdfsig_hist_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResDown","",RooArgSet(*mreco),conv_template_1_proj_dn);
		RooHistFunc *pdfsig_hist_func_dn= new RooHistFunc("pdfsig_hist_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResDown","",RooArgSet(*mreco),*pdfsig_hist_hist_dn);

		RooDataHist* hgenm_dn= new RooDataHist("hgenm"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResDown","",RooArgSet(*mreco),conv_template_2_proj_dn);
		RooHistFunc* hgenm_func_dn= new RooHistFunc("hgenm_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResDown","",RooArgSet(*mreco),*hgenm_dn);

		RooDataHist *inter_hist_dn= new RooDataHist("inter_hist"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResDown","",RooArgSet(*mreco),conv_template_124_proj_dn);
		RooHistFunc *inter_func_dn= new RooHistFunc("inter_func"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample]+"_2l2qJetResDown","",RooArgSet(*mreco),*inter_hist_dn);

		HZZ4L_RooHighmass_1D *bsi_hist_dn=new HZZ4L_RooHighmass_1D("bsi_hist","bsi_hist",*mreco,*signorm,RooArgList(*pdfsig_hist_func_dn,*hgenm_func_dn,*inter_func_dn)); 
		bsi_hist_dn->SetNameTitle(dataName[sample]+"_CMS_2l2q"+jetName[jetType]+"ResDown",dataName[sample]+"_CMS_2l2q"+jetName[jetType]+"ResDown");
		w.import(*bsi_hist,RecycleConflictNodes());
		w.import(*bsi_hist_up,RecycleConflictNodes());
		w.import(*bsi_hist_dn,RecycleConflictNodes());
	}
	//		RooAbsPdf* bsi_histy= bsi_hist->createProjection(*mreco) ;
	//		RooAbsPdf* bsi_histx= bsi_hist->createProjection(*dbkg) ;
	////	return;
	cout << "Reco Sig "<<conv_template_1_proj->Integral()<<endl;
	cout << "Reco bkg "<<conv_template_2_proj->Integral()<<endl;
	cout << "Reco int"<<conv_template_124_proj->Integral()<<endl;
	//return;



	//hgenm_func->plotOn(frame_dbkg);
	//pdfsig_hist_func->plotOn(frame_dbkg);
	//inter_func->plotOn(frame_dbkg);
	//bsi_hist->plotOn(frame,Normalization(h160hist->sum(false),RooAbsReal::NumEvent));
	//bsi_histx->plotOn(frame,Normalization(h160hist->sum(false),RooAbsReal::NumEvent));
	//bsi_histy->plotOn(frame_dbkg,Normalization(h160hist->sum(false),RooAbsReal::NumEvent));
	////	h160hist_noit->plotOn(frame,Binning(nbins_reco/20),MarkerColor(2));
	////bsi_hist->plotOn(frame_dbkg);
	//frame->Draw();
	////frame_dbkg->Draw();

	//gPad->Print(Form("2l2q_closure/m%.0f_%s_%s_%s_%s.png",mass_d,chanName[chan].Data(),jetName[jetType].Data(),cateName[cate].Data(),sampleName[sample].Data()));
	//gPad->Print(Form("2l2q_closure/m%.0f_%s_%s_%s_%s.pdf",mass_d,chanName[chan].Data(),jetName[jetType].Data(),cateName[cate].Data(),sampleName[sample].Data()));
	//	gPad->Print(Form("2l2q_closure/m%.0f_%s_BSI_%s.png",mass_d,jetName[jetType].Data(),sampleName[sample].Data()));
	//	gPad->Print(Form("2l2q_closure/m%.0f_%s_BSI_%s.pdf",mass_d,jetName[jetType].Data(),sampleName[sample].Data()));
	//	frame_dbkg->Draw();
	//	gPad->Print(Form("2l2q_closure/m%.0f_%s_dbkg_BSI_%s.png",mass_d,jetName[jetType].Data(),sampleName[sample].Data()));
	//	gPad->Print(Form("2l2q_closure/m%.0f_%s_dbkg_BSI_%s.pdf",mass_d,jetName[jetType].Data(),sampleName[sample].Data()));
	//	gPad->Print(Form("2l2q_closure/m%.0f_%s_S_%s.png",mass_d,jetName[jetType].Data(),sampleName[sample].Data()));
	//	gPad->Print(Form("2l2q_closure/m%.0f_%s_S_%s.pdf",mass_d,jetName[jetType].Data(),sampleName[sample].Data()));
	//	return;



	ProcessNormalization* int_sig;
	ProcessNormalization* int_bkg;
	ProcessNormalization* int_int;

	float int_sig_norm = conv_template_1->Integral();
	float int_bkg_norm = conv_template_2->Integral();
	float int_int_norm = conv_template_124->Integral();

	float int_sig_up = conv_template_1_up->Integral();
	float int_bkg_up = conv_template_2_up->Integral();
	float int_int_up = conv_template_124_up->Integral();

	float int_sig_dn = conv_template_1_dn->Integral();
	float int_bkg_dn = conv_template_2_dn->Integral();
	float int_int_dn = conv_template_124_dn->Integral();

	cout<< int_int_norm <<"\t"<<int_int_dn<<"\t"<<int_int_up<<endl;
	int_sig = new ProcessNormalization("int_sig"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"int_sig"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],int_sig_norm);
	int_bkg = new ProcessNormalization("int_bkg"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"int_bkg"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],int_bkg_norm); 
	int_int = new ProcessNormalization("int_int"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],"int_int"+jetName[jetType]+chanName[chan]+cateName[cate]+sampleName[sample],int_int_norm); 

	RooRealVar *sysrr=new RooRealVar("CMS_2l2q"+jetName[jetType]+"Res","2l2q"+jetName[jetType]+"Res",-7,7);
	RooRealVar *kfacrr=new RooRealVar("kfactor_ggZZ","kfactor_ggZZ",-7,7);
	if(int_sig_norm!=0)
		int_sig->addAsymmLogNormal(int_sig_dn/int_sig_norm, int_sig_up/int_sig_norm, *sysrr);
	if(int_bkg_norm!=0){
		int_bkg->addAsymmLogNormal(int_bkg_dn/int_bkg_norm, int_bkg_up/int_bkg_norm, *sysrr);
		//if(sample==0)
		int_bkg->addAsymmLogNormal(0.9, 1.1, *kfacrr);
	}
	if(int_int_norm!=0 && fabs(int_int_norm/int_sig_norm)>1.E-4){
		int_int->addAsymmLogNormal(int_int_dn/int_int_norm, int_int_up/int_int_norm, *sysrr);
		//if(sample==0)
		int_int->addAsymmLogNormal(0.95, 1.05, *kfacrr);
	}

	RooRealVar *btag= new RooRealVar("BTAG_"+jetNameLower[jetType],"BTAG_"+jetNameLower[jetType],0,-3,3); 
	RooRealVar *jes= new RooRealVar("JES","JES",0,-3,3); 

	float errbtag = btagU[jetType]*ratio_btag;
	float errjes = jesU[sample]*ratio_vbftag;
	cout<< errbtag<<"\t"<<errjes<<endl;

	string formula = Form("(@0*@3+@1+@2*sqrt(@3))*(@4*(%.3f)+1)*(1+@5*(%.3f))",errbtag, errjes);

	RooFormulaVar *ggH_norm=new RooFormulaVar(dataName[sample]+"_norm",dataName[sample]+"_norm",formula.c_str(),RooArgList(*int_sig,*int_bkg,*int_int,*signorm,*btag,*jes));
	w.import(*ggH_norm,RecycleConflictNodes());
	cout << "final "<<ggH_norm->getVal()<<endl;

	TFile *fwork ;
	if(width_d<0.1)
		fwork= new TFile(Form("workspace_template%s_2l2q_Moriond/hzz4l_%d_%d_%d_%dS_13TeV.input_func_%2.0f_%2.2f.root",workName[is2D].Data(),chan,jetType,cate,sample,mass_d,width_d),"recreate");
	else
		fwork= new TFile(Form("workspace_template%s_2l2q_Moriond/hzz4l_%d_%d_%d_%dS_13TeV.input_func_%2.0f_%2.1f.root",workName[is2D].Data(),chan,jetType,cate,sample,mass_d,width_d),"recreate");
	fwork->cd();
	w.Write();
	fwork->Close();

	}
	void construct_2l2q(double m=600,double w=123.0,int cate=0, int is2D=1,double inteval=2,int loo=1){
		//gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
		// gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
		// gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");
		gROOT->ProcessLine(".x tdrstyle.cc");
		for(int i=0;i<loo;i++){
			dosomething(0,0,0,cate,is2D,m,w);
			dosomething(1,0,0,cate,is2D,m,w);
			dosomething(0,1,0,cate,is2D,m,w);
			dosomething(1,1,0,cate,is2D,m,w);
	      		dosomething(0,0,1,cate,is2D,m,w);
			dosomething(1,0,1,cate,is2D,m,w);
			dosomething(0,1,1,cate,is2D,m,w);
			dosomething(1,1,1,cate,is2D,m,w);
			dosomething(0,0,2,cate,is2D,m,w);
			dosomething(1,0,2,cate,is2D,m,w);
			dosomething(1,1,2,cate,is2D,m,w);
			dosomething(0,1,2,cate,is2D,m,w);
			m +=inteval; 
		}
	}

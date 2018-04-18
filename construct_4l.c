
#include "HZZ4L_RooHighmass_1D.cc+"
#include "HZZ4L_RooHighmass.cc+"
#include "HZZ2L2QRooPdfs.cc+"
#include "ProcessNormalization.cc+"
#include "VerticalInterpPdf.cc+"
#include "SplinePdf.cc+"
using namespace RooFit;
double dcbPara_2e2mu_2nd[6][11]={
	/*
	   1000,2000,0.796001, 0.00121753, -4.04361e-07, 1.26493, 0.000279677, 6.45654e-08, -0.625592, 0.0021702, -4.08065e-07,
	   1000,2000,1.62422, 0.000557576, -3.24469e-07, 2.0048, -0.000203594, 5.61161e-08, 0.966682, 0.000834528, -2.03414e-07,
	   1000,2000,0.00417567, -0.00301974, 8.98306e-07, -1.65418, 0.000296969, -7.60047e-07, 1.98469, -0.0033419, 1.49671e-07,
	   1000,2000,3.27247, -0.00189252, 5.20826e-07, 2.80943, -0.000966426, 5.77799e-08, 4.8379, -0.0029949, 5.64897e-07,
	   1000,2000,3.23411, 0.00205592, -1.27997e-07, 1.73047, 0.0050632, -1.63164e-06, 14.1017, -0.00730808, 1.46118e-06,
	   1000,2000,0.318145, 0.0100607, 4.59762e-06, -0.3251, 0.0113472, 3.95437e-06, -18.4583, 0.0294803, -5.78928e-07
	   */

	1000,2000,0.779987, 0.00134559, -3.60197e-07, 1.10251, 0.000700543, -3.76727e-08, -0.825775, 0.00262883, -5.19744e-07,
	1000,2000,1.67285, 0.000100963, 1.04155e-07, 1.30519, 0.000836286, -2.63506e-07, 2.28023, -0.00013876, -1.97447e-08,
	1000,2000,-0.149739, -0.00279488, 5.71049e-07, 0.304491, -0.00370334, 1.02528e-06, -4.24097, 0.000842119, -1.11086e-07,
	1000,2000,2.99058, -0.00128841, 4.28005e-09, 3.34155, -0.00199035, 3.55253e-07, 4.15905, -0.00280785, 5.59629e-07,
	1000,2000,3.54188, 0.00393276, -3.24087e-06, 8.37189, -0.00572726, 1.58914e-06, 0.359371, 0.00228526, -4.1399e-07,
	1000,2000,0.310114, 0.0100951, 4.6391e-06, -1.01784, 0.012751, 3.31115e-06, -15.5681, 0.0273012, -3.26413e-07

};
double dcbPara_4mu_2nd[6][11]=
{
	/*
	   1000,2000,1.19829, 0.000371215, 2.99251e-07, 0.556051, 0.00165569, -3.42986e-07, 0.145479, 0.00206626, -4.45629e-07,
	   1000,2000,1.83825, -2.55154e-05, -2.65552e-08, 1.89271, -0.000134432, 2.79031e-08, 0.977273, 0.000781005, -2.00956e-07,
	   1000,2000,0.159507, -0.00284656, 4.13771e-07, 3.19728, -0.00892209, 3.45154e-06, -36.0397, 0.0303149, -6.3577e-06,
	   1000,2000,1.83379, 0.00162983, -1.76978e-06, 4.53299, -0.00376858, 9.29427e-07, 2.22756, -0.00146315, 3.53069e-07,
	   1000,2000,3.26367, -0.0013482, 3.61841e-06, -4.48684, 0.0141528, -4.1321e-06, 19.7713, -0.0101053, 1.93244e-06,
	   1000,2000,-0.19907, 0.00940217, 1.0038e-05, -7.33146, 0.0236669, 2.90565e-06, -14.5777, 0.0309132, 1.09408e-06
	   */

	1000,2000,1.23952, 8.20617e-05, 5.85123e-07, -0.311558, 0.00318423, -9.6596e-07, 5.08132, -0.00220866, 3.82261e-07,
	1000,2000,1.90846, -0.000474028, 4.07841e-07, 0.620535, 0.00210183, -8.80086e-07, 9.2037, -0.00648134, 1.26571e-06,
	1000,2000,0.16386, -0.00287822, 4.47869e-07, 2.80211, -0.00815472, 3.08612e-06, -35.0484, 0.0296958, -6.37651e-06,
	1000,2000,1.80762, 0.00183435, -1.95209e-06, 4.92739, -0.0044052, 1.16768e-06, 0.741145, -0.000218948, 1.2112e-07,
	1000,2000,2.61507, 0.00353264, -1.84442e-06, 13.0793, -0.0173958, 8.61979e-06, -61.9304, 0.0576139, -1.01326e-05,
	1000,2000,-0.183005, 0.00929162, 1.02552e-05, -8.32009, 0.0255658, 2.11815e-06, -12.1115, 0.0293572, 1.17029e-06


};
double dcbPara_4e_2nd[6][11]=
{
	/*	
		1000,2000,0.94329, 0.000246632, -3.71884e-08, 1.00242, 0.000128368, 2.19433e-08, 1.06071, 7.00765e-05, 3.65162e-08,
		1000,2000,1.86047, -0.000488889, 2.42265e-07, 1.80599, -0.000379925, 1.87783e-07, 0.239024, 0.00118704, -2.03959e-07,
		1000,2000,-0.58501, -0.0012942, 6.12243e-07, -3.62988, 0.00479555, -2.43263e-06, 8.63908, -0.00747342, 6.34614e-07,
		1000,2000,5.935, -0.00866883, 4.86395e-06, 0.0461084, 0.00310894, -1.02494e-06, 5.13685, -0.0019818, 2.4775e-07,
		1000,2000,2.96678, 0.00137339, 1.09836e-06, -2.24232, 0.0117916, -4.11074e-06, 23.1239, -0.0135746, 2.23081e-06,
		1000,2000,1.27682, 0.00873798, -1.09829e-06, 3.8105, 0.00367062, 1.4354e-06, -3.17133, 0.0106524, -3.10063e-07
		*/

	1000,2000,0.820461, 0.000974084, -4.39409e-07, 1.45941, -0.000303822, 1.99544e-07, 0.542795, 0.000612796, -2.96107e-08,
	1000,2000,1.93006, -0.00103167, 7.81493e-07, 1.13399, 0.000560478, -1.45809e-08, 0.123982, 0.00157048, -2.67082e-07,
	1000,2000,-0.853194, -0.00110413, 4.11286e-09, -0.333603, -0.00214331, 5.23704e-07, 0.39652, -0.00287344, 7.06234e-07,
	1000,2000,6.08942, -0.0102987, 6.10307e-06, -1.83225, 0.00554461, -1.8186e-06, 8.11201, -0.00439966, 6.67463e-07,
	1000,2000,2.98608, 0.00571114, -4.06742e-06, 7.25445, -0.00282559, 2.00946e-07, 9.46499, -0.00503613, 7.53581e-07,
	1000,2000,1.23564, 0.00884888, -1.03523e-06, 2.8199, 0.00568035, 5.49037e-07, 0.0497875, 0.00845046, -1.4349e-07

};
double dcbPara_rse_4e_2nd[6][11]={
	//	1000,2000,1.19262, -0.000667809, 5.36599e-07, 0.43979, 0.000837851, -2.16231e-07, 3.21441, -0.00193677, 4.77424e-07, 
	//	1000,2000,1.61631, -5.49678e-05, 7.0388e-08, 1.57029, 3.70722e-05, 2.4368e-08, 3.06399, -0.00145663, 3.97793e-07, 
	//	1000,2000,-2.05537, 0.00359008, -2.8143e-06, 0.285642, -0.00109194, -4.73284e-07, -5.86174, 0.00505544, -2.01013e-06, 
	//	1000,2000,2.63838, 0.000180588, -2.8094e-07, 2.92796, -0.000398572, 8.6395e-09, 1.24401, 0.00128538, -4.12348e-07, 
	//	1000,2000,3.85484, 0.0051618, -3.11686e-06, 5.98042, 0.00091064, -9.91285e-07, 14.3802, -0.00748914, 1.10866e-06, 
	//	1000,2000,2.59881, 0.00292441, 2.30583e-06, 0.0443096, 0.00803341, -2.48673e-07, 8.55604, -0.00047832, 1.87926e-06 
	//4e
	1000,2000,1.28198, -0.000800173, 8.13379e-07, 0.121824, 0.00152013, -3.46774e-07, 0.837972, 0.000803985, -1.67737e-07,
	1000,2000,2.70798, -0.0033039, 2.13291e-06, 0.256486, 0.00159908, -3.18585e-07, 2.42084, -0.000565271, 2.22504e-07,
	1000,2000,-1.41565, 0.000481046, -1.44466e-06, 1.73532, -0.0058209, 1.70631e-06, -7.28826, 0.00320269, -5.49588e-07,
	1000,2000,3.51468, -0.00217425, 5.86534e-07, 3.27893, -0.00170274, 3.50782e-07, 3.58753, -0.00201134, 4.27931e-07,
	1000,2000,-0.957121, 0.0120793, -7.63449e-06, 7.82858, -0.00549207, 1.15121e-06, 3.83404, -0.00149753, 1.52578e-07,
	1000,2000,2.96621, 0.00459584, 2.03718e-06, 0.150321, 0.0102276, -7.78714e-07, 8.93172, 0.00144622, 1.41663e-06
};
double dcbPara_rse_2e2mu_2nd[6][11]=
{
	//	1000,2000,1.15302, 0.000116426, 2.96574e-07, 0.64029, 0.00114189, -2.16156e-07, 0.699032, 0.00108314, -2.0147e-07, 
	//	1000,2000,2.06924, -0.00029104, 1.16066e-07, 1.86057, 0.0001263, -9.26045e-08, 3.08304, -0.00109617, 2.13013e-07, 
	//	1000,2000,-0.320145, -0.00206645, 1.3675e-07, -0.877443, -0.000951854, -4.20548e-07, 6.11059, -0.00793989, 1.32646e-06, 
	//	1000,2000,2.03043, 0.00162785, -1.59229e-06, 4.21386, -0.00273901, 5.91142e-07, 3.09987, -0.00162502, 3.12644e-07, 
	//	1000,2000,3.03098, -0.000334323, 2.01689e-06, -0.246812, 0.00622126, -1.2609e-06, -9.65545, 0.0156299, -3.61306e-06, 
	//	1000,2000,-0.0479635, 0.0112779, 3.21573e-06, 2.12083, 0.00694031, 5.38453e-06, -31.425, 0.0404861, -3.00193e-06
	//2e2mu
	1000,2000,1.05507, 0.0004398, 1.22621e-07, 0.424062, 0.00170181, -5.08384e-07, 3.51185, -0.00138598, 2.63564e-07,
	1000,2000,0.686783, 0.0028037, -1.5382e-06, 2.66987, -0.00116246, 4.44878e-07, -3.20192, 0.00470932, -1.02307e-06,
	1000,2000,-0.799399, -0.00109496, -9.43704e-07, 3.05267, -0.00879909, 2.90836e-06, -21.1501, 0.0154037, -3.14234e-06,
	1000,2000,3.05656, -0.00191973, 5.45926e-07, 3.30836, -0.00242333, 7.97722e-07, -1.32157, 0.0022066, -3.59759e-07,
	1000,2000,7.32362, -0.0143608, 8.87835e-06, -4.65324, 0.00959291, -3.09851e-06, 17.0498, -0.0121101, 2.32725e-06,
	1000,2000,0.128783, 0.0121983, 2.87379e-06, 1.73759, 0.00898074, 4.4826e-06, -22.7136, 0.0334319, -1.63019e-06
};
double eff[3][11]={

	/*
	// ggH  4e
	-4.446093E+00, 4.591943E+00, -3.587753E+02, 3.396792E+02, 2.255649E+00, 8.701775E-04, -2.591608E-07, 2.980908E-01, 1.108333E+02, 8.950199E+01, 1.697673E-11,
	// ggH  2e2mu
	-4.479593E+00, 4.563339E+00, -3.803250E+02, 3.375862E+02, 4.963176E+00, 2.306610E-03, -8.351650E-07, 6.989435E-01, 6.915814E+01, 1.040211E+02, 8.912866E-11,
	// ggH  4mu
	-4.458068E+00, 4.582843E+00, -7.488067E+03, 4.547506E+03, 1.123692E+01, -4.069213E-03, 1.134668E-06, -3.384899E-01, 7.094597E+01, 3.862078E+01, -1.147461E-10
	*/


	// ggH  4e
	-4.409092E+00, 4.629159E+00, -3.499506E+02, 3.458808E+02, 1.584749E+00, 5.822571E-04, -1.884481E-07, 2.117278E-01, 1.355489E+02, 7.628844E+01, 1.331265E-11,
	// ggH  4mu
	-4.460109E+00, 4.580740E+00, -7.560276E+03, 4.554498E+03, 1.170131E+01, -3.898090E-03, 1.026058E-06, -3.831164E-01, 6.646347E+01, 4.041114E+01, -9.867291E-11,
	// ggH  2e2mu
	-4.485994E+00, 4.557560E+00, -3.779531E+02, 3.229685E+02, 5.752549E+00, 3.435045E-03, -1.399594E-06, 6.674296E-01, 6.925980E+01, 1.011219E+02, 1.699300E-10


};
double eff_vbf[3][11]={
	/*
	// VBF  4e
	-4.439290E+00, 4.601427E+00, -3.197092E+02, 2.908201E+02, 1.981189E+00, 9.519061E-04, -3.058299E-07, 1.409187E-01, 1.313160E+02, 6.803395E+01, 2.495282E-11,
	// VBF  2e2mu
	-4.481897E+00, 4.560994E+00, -3.857275E+02, 3.260331E+02, 5.196698E+00, 2.741781E-03, -9.661342E-07, 5.555572E-01, 7.482259E+01, 1.007256E+02, 9.658647E-11,
	// VBF  4mu
	-4.446348E+00, 4.593324E+00, -8.046137E+03, 4.979096E+03, 9.585262E+00, -3.042051E-03, 7.532618E-07, -3.148329E-01, 7.269329E+01, 3.878612E+01, -6.804696E-11,
	*/

	// VBF  4e
	-4.457080E+00, 4.583149E+00, -3.220522E+02, 2.772364E+02, 2.590067E+00, 1.438727E-03, -5.122000E-07, 1.369010E-01, 1.296028E+02, 6.672205E+01, 4.702222E-11,
	// VBF  4mu
	-4.464713E+00, 4.575797E+00, -7.792695E+03, 4.494984E+03, 1.001845E+01, -1.132429E-03, -6.647936E-08, -4.162270E-01, 5.929335E+01, 4.422993E+01, 4.317671E-11,
	// VBF  2e2mu
	-4.487157E+00, 4.555744E+00, -3.852006E+02, 3.212818E+02, 6.250167E+00, 3.400357E-03, -1.287049E-06, 6.712138E-01, 6.316531E+01, 1.033601E+02, 1.402357E-10

};
TString sampleName[2]={"ggH","VBF"};
TString dataName[2]={"ggH","qqH"};
TString workName[2]={"1D","2D"};
TString bkgName[2][2]={
	"root://eoscms.cern.ch//store/user/xiaomeng/lxplusBackUp/rootfiles/rootfiles/ggZZ_Bkg_xcheck.root","root://eoscms.cern.ch//store/user/xiaomeng/lxplusBackUp/rootfiles/rootfiles/ggZZ_4e.root",
	"root://eoscms.cern.ch//store/user/xiaomeng/lxplusBackUp/rootfiles/rootfiles/vbf/phantom_bkg.root","root://eoscms.cern.ch//store/user/xiaomeng/lxplusBackUp/rootfiles/rootfiles/vbf/phantom_bkg_4e.root"
};
TString hName[2][2]={
	"root://eoscms.cern.ch//store/user/xiaomeng/lxplusBackUp/rootfiles/rootfiles/125_sig.root","root://eoscms.cern.ch//store/user/xiaomeng/lxplusBackUp/rootfiles/rootfiles/125_sig_4e.root",
	"root://eoscms.cern.ch//store/user/xiaomeng/lxplusBackUp/rootfiles/rootfiles/vbf/phantom_125_2e2mu.root","root://eoscms.cern.ch//store/user/xiaomeng/lxplusBackUp/rootfiles/rootfiles/vbf/phantom_125_4e.root"
};
float bkgxsec[2][2]={3.19,1.59,0.642,0.310}; // ggH/VBF 2e2mu/4e

//float bkgxsec[2][2]={0,0,0.1,0.05}; // ggH/VBF 2e2mu/4e
//float hxsec[2][2]={0.855*3.04,0.452*3.04,0.248*1.26844,0.133*1.26844};

//float hxsec[2][2]={0.855*3.425*1.13,0.452*3.425*1.13,0.248*1.26844,0.133*1.26844};// Remember to revert it back
float hxsec[2][2]={0.855*3.425*1.15,0.452*3.425*1.15,0.248*1.4568*1.15,0.133*1.4568*1.15};// Remember to revert it back
//float hxsec[2][2]={0.,0.,0,0};// Remember to revert it back

void dosomething(TString chan="2e2mu",double mass_d=450, double width_d=46.8,int sample =0, int cate_vbf=0,int is2D=1){
	int chanN = 0;
	if(chan!="2e2mu")
		chanN=1;

	double sigfull, bkgfull,hfull,lumi,signum;
	lumi = 35.8;
	//	sigfull = 0.0673*0.0673*lumi*1000*0.6138650761079254; //Remember to revert it back
	sigfull = 0.0673*0.0673*lumi; //Remember to revert it back

	TChain *t125 = new TChain ("SelectedTree");
	TChain *ggzz = new TChain("SelectedTree");

	t125->Add(hName[sample][chanN]);
	ggzz->Add(bkgName[sample][chanN]);
	hfull = hxsec[sample][chanN]*lumi;
	bkgfull = bkgxsec[sample][chanN]*lumi;


	double dcbPara_2nd[6][11];
	double effsig[11]={0.};

	if (chan=="4e") 	{
		if(cate_vbf!=2){
			for (int i=0;i<6;i++){for(int j=0;j<11;j++){dcbPara_2nd[i][j]= dcbPara_4e_2nd[i][j];}}
		}
		else{
			for (int i=0;i<6;i++){for(int j=0;j<11;j++){dcbPara_2nd[i][j]= dcbPara_rse_4e_2nd[i][j];}}
		}
		if(sample)
			for (int i=0;i<11;i++){effsig[i]=eff_vbf[0][i];}
		else
			for (int i=0;i<11;i++){effsig[i]=eff[0][i];}
	} 
	if (chan=="4mu")  {
		for (int i=0;i<6;i++){for(int j=0;j<11;j++){dcbPara_2nd[i][j]= dcbPara_4mu_2nd[i][j];}}
		if(sample)
			for (int i=0;i<11;i++){effsig[i]=eff_vbf[1][i];}
		else
			for (int i=0;i<11;i++){effsig[i]=eff[1][i];}
	}
	if (chan=="2e2mu") {
		if(cate_vbf!=2){
			for (int i=0;i<6;i++){for(int j=0;j<11;j++){dcbPara_2nd[i][j]= dcbPara_2e2mu_2nd[i][j];}}
		}
		else{
			for (int i=0;i<6;i++){for(int j=0;j<11;j++){dcbPara_2nd[i][j]= dcbPara_rse_2e2mu_2nd[i][j];}}
		}
		for (int i=0;i<6;i++){for(int j=0;j<11;j++){dcbPara_2nd[i][j]= dcbPara_2e2mu_2nd[i][j];}}
		if(sample)
			for (int i=0;i<11;i++){effsig[i]=eff_vbf[2][i];}
		else
			for (int i=0;i<11;i++){effsig[i]=eff[2][i];}
	}

	RooWorkspace w("w");
	// range of gen and reco
	const double low=100;
	const double high=3600;
	const int nbins=int(high-low);


	double lowvalue [3]={110.,110.,300.};

	const double low_reco=lowvalue[cate_vbf];
	const double high_reco=3500;
	const int nbins_reco=int(high_reco-low_reco);

	// define gen mass and reco mass
	RooRealVar* mzz = new RooRealVar("ZZMass_gen","M_{ZZ} (GeV)",400,low,high);
	RooRealVar* dbkg= new RooRealVar("dbkg_kin","Dbkg_{kin} ",0.5,0.,1.);
	RooRealVar* mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",400,low_reco,high_reco);
	if(cate_vbf==2){
		mreco->SetName("mreco_rse");
		mreco->SetTitle("mreco_rse");
	}


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
	TString floname = "gghKfac_input_spline.root";
	//	TString floname = "ggh_input_spline.root";
	TString spname ="sp_xsec_statnom";
	TString spname_4e ="sp_xsec_statnom_4e";
	if(sample==1){
		floname = "width_new_spline_3500.root";
		spname = "br_2e2mu";
		spname_4e= "br_4e";
	}
	TFile *flo=new TFile(floname,"read");
	TSpline3 *lo=(TSpline3*) flo->Get(spname);
	TSpline3 *lo_4e=(TSpline3*) flo->Get(spname_4e);

	double pole2e2mu ,pole4e;
	pole2e2mu = lo->Eval(mass_d);
	pole4e= lo_4e->Eval(mass_d);

	SplinePdf *pdf;
	SplinePdf *pdf_other;

	if(chan=="2e2mu"){
		pdf=new SplinePdf("pdf_"+chan,"",*mzz,*mean,*sigma,*lo);
		pdf_other=new SplinePdf("pdf_other","",*mzz,*mean,*sigma,*lo_4e);
		if(width_d<0.5)
			signum = sigfull*pole2e2mu/(pole2e2mu+pole4e*2);	
		else{
			double integralsum = (pdf->createIntegral(*mzz))->getVal();
			double integralsum_other= (pdf_other->createIntegral(*mzz))->getVal();
			signum = integralsum/(integralsum+integralsum_other*2.)*sigfull;
			cout<<integralsum/(integralsum+integralsum_other*2.)<<endl;
		}
	}
	else{
		pdf=new SplinePdf("pdf_"+chan,"",*mzz,*mean,*sigma,*lo_4e);
		pdf_other=new SplinePdf("pdf_other","",*mzz,*mean,*sigma,*lo);
		if(width_d<0.5)
			signum = sigfull*pole4e/(pole2e2mu+pole4e*2);	
		else{
			double integralsum = (pdf->createIntegral(*mzz))->getVal();
			double integralsum_other= (pdf_other->createIntegral(*mzz))->getVal();
			signum = integralsum/(integralsum*2+integralsum_other)*sigfull;
			cout<<integralsum/(integralsum*2+integralsum_other)<<endl;
		}
	}


	TFile *fkfactor = new TFile("Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
	//TGraph *ggZZ_kf = (TGraph*)fkfactor->Get("sp_kfactor_Nominal");
	TGraph *ggZZ_kf = (TGraph*)fkfactor->Get("kfactor_Nominal");

	RooPlot* frame_mzz= mzz->frame(Title("Z mass")) ;

	mean->setVal(125.);
	sigma->setVal(0.00407);

	double hall = t125->GetEntries("ZZMass!=0");
	double hsel = t125->GetEntries(Form("ZZMass>%f&&ZZMass<%f",low,high));
	double hnum = hsel/hall*hfull;

	TH1F *pdf125_hist=new TH1F("pdf125_hist","",nbins,low,high);
	double bwid = (high-low)/double(nbins);
	for(int k=0;k<nbins;k++){
		//		double x = pdf125_hist->GetBinCenter(k+1);
		double x = low+k*bwid; 
		//double kfac = ggZZ_kf->Eval(x);
		mzz->setVal(x);
		double bc_125= pdf->getVal(*mzz)*hnum*bwid; //Remember to revert it back
		if(x==125){	
			bc_125=0;
			for(int j =-50;j<50;j++){
				mzz->setVal(x+j*0.2*0.00407);
				bc_125+=pdf->getVal(*mzz)*hnum*0.2*0.00407; //Remember to revert it back
			}
		}
		pdf125_hist->SetBinContent(k+1,bc_125);		
	}
	//	cout<< pdf125_hist->Integral(1,31)<<endl;
	//	return;


	// signal
	mean->setVal(mass_d);
	sigma->setVal(width_d);



	//bkg
	double bkgall = ggzz->GetEntries("");
	double bkgsel = ggzz->GetEntries(Form("ZZMass>%f&&ZZMass<%f",low,high));
	double bkgnum = bkgsel/bkgall*bkgfull; 

	TH1F *genm= new TH1F ("genm","",nbins,low,high);
	ggzz->Draw("ZZMass>>genm",Form("ZZMass>%f",low));

	genm->Scale(bkgnum/genm->Integral());

	// resolution
	TString formu_2nd=" (@0<@1)*(@3+@0*@4+@0*@0*@5 ) + ( @0>=@1 && @0<@2)*(@6+@0*@7+@0*@0*@8) + (@0>=@2)*(@9+@0*@10+@0*@0*@11)";	

	RooArgList formuList_a1;
	RooArgList formuList_a2;
	RooArgList formuList_mean;
	RooArgList formuList_n1;
	RooArgList formuList_n2;
	RooArgList formuList_sigma;
	formuList_a1.add(*mzz);
	formuList_a2.add(*mzz);
	formuList_mean.add(*mzz);
	formuList_n1.add(*mzz);
	formuList_n2.add(*mzz);
	formuList_sigma.add(*mzz);

	RooConstVar* a1_p0_0_2nd[11] ;
	RooConstVar* a2_p0_0_2nd[11] ;
	RooConstVar* mean_p0_0_2nd[11] ;
	RooConstVar* n1_p0_0_2nd[11] ;
	RooConstVar* n2_p0_0_2nd[11] ;
	RooConstVar* sigma_p0_0_2nd[11] ;
	for(int i =0; i<11;i++){
		a1_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_a1_p0_0_2nd",chan.Data(),i),Form("%s_%d_a1_p0_0_2nd",chan.Data(),i),dcbPara_2nd[0][i]);
		a2_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_a2_p0_0_2nd",chan.Data(),i),Form("%s_%d_a2_p0_0_2nd",chan.Data(),i),dcbPara_2nd[1][i]);
		mean_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_mean_p0_0_2nd",chan.Data(),i),Form("%s_%d_mean_p0_0_2nd",chan.Data(),i),dcbPara_2nd[2][i]);
		n1_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_n1_p0_0_2nd",chan.Data(),i),Form("%s_%d_n1_p0_0_2nd",chan.Data(),i),dcbPara_2nd[3][i]);
		n2_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_n2_p0_0_2nd",chan.Data(),i),Form("%s_%d_n2_p0_0_2nd",chan.Data(),i),dcbPara_2nd[4][i]);
		sigma_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_sigma_p0_0_2nd",chan.Data(),i),Form("%s_%d_sigma_p0_0_2nd",chan.Data(),i),dcbPara_2nd[5][i]);

		formuList_a1.add(*a1_p0_0_2nd[i]);
		formuList_a2.add(*a2_p0_0_2nd[i]);
		formuList_mean.add(*mean_p0_0_2nd[i]);
		formuList_n1.add(*n1_p0_0_2nd[i]);
		formuList_n2.add(*n2_p0_0_2nd[i]);
		formuList_sigma.add(*sigma_p0_0_2nd[i]);
	}

	RooFormulaVar* a1_p0_2nd= new RooFormulaVar("a1_p0_2nd"+chan,"a1_p0_2nd"+chan,formu_2nd,formuList_a1);
	RooFormulaVar* a2_p0_2nd= new RooFormulaVar("a2_p0_2nd"+chan,"a2_p0_2nd"+chan,formu_2nd,formuList_a2);
	RooFormulaVar* mean_p0_2nd= new RooFormulaVar("mean_p0_2nd"+chan,"mean_p0_2nd"+chan,formu_2nd,formuList_mean);
	RooFormulaVar* n1_p0_2nd= new RooFormulaVar("n1_p0_2nd"+chan,"n1_p0_2nd"+chan,formu_2nd,formuList_n1);
	RooFormulaVar* n2_p0_2nd= new RooFormulaVar("n2_p0_2nd"+chan,"n2_p0_2nd"+chan,formu_2nd,formuList_n2);
	RooFormulaVar* sigma_p0_2nd= new RooFormulaVar("sigma_p0_2nd"+chan,"sigma_p0_2nd"+chan,formu_2nd,formuList_sigma);

	RooFormulaVar *sigma_p0_up = new RooFormulaVar("sigma_p0_up"+chan,"","@0+0.2*@0",*sigma_p0_2nd);
	RooFormulaVar *sigma_p0_dn = new RooFormulaVar("sigma_p0_dn"+chan,"","@0-0.2*@0",*sigma_p0_2nd);

	RooDoubleCB dcrReso_piece("dcrReso"+chan,"Double Crystal ball ",*mreco,*mzz,*mean_p0_2nd,*sigma_p0_2nd,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);
	RooDoubleCB dcrReso_piece_up("dcrReso"+chan+"_up","dcb up",*mreco,*mzz,*mean_p0_2nd,*sigma_p0_up,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);
	RooDoubleCB dcrReso_piece_dn("dcrReso"+chan+"_dn","dcb dn",*mreco,*mzz,*mean_p0_2nd,*sigma_p0_dn,*a1_p0_2nd,*n1_p0_2nd,*a2_p0_2nd,*n2_p0_2nd);

	//TFile *ftemplate=new TFile("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_8_0_6/src/ZZMatrixElement/MELA/test/template_com7_Rebin_var2.root");
	TFile *ftemplate=new TFile("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_8_0_6/src/ZZMatrixElement/MELA/test/template_"+chan+"_"+sampleName[sample]+"_new.root");
	//	TFile *ftemplate=new TFile("template_cond.root");

	TH2F *template_1= (TH2F*)ftemplate->Get("T_2D_1_s");
	TH2F *template_2= (TH2F*)ftemplate->Get("T_2D_2_s");
	TH2F *template_124= (TH2F*)ftemplate->Get("T_2D_124_s");
	TH2F *template_124Im= (TH2F*)ftemplate->Get("T_2D_124_Ims");

	TH2F *template_1r= new TH2F("template1","",nbins,low,high,30,0.,1.); 
	TH2F *template_2r= new TH2F("template2","",nbins,low,high,30,0.,1.); 
	TH2F *template_124r= new TH2F("template124","",nbins,low,high,30,0.,1.); 

	//	TH1F *template_1r= new TH1F("template1","",nbins,low,high); 
	//	TH1F *template_2r= new TH1F("template2","",nbins,low,high); 
	//	TH1F *template_124r= new TH1F("template124","",nbins,low,high); 

	TH2F *conv_template_1= new TH2F("conv_template1","",nbins_reco/2,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_2= new TH2F("conv_template2","",nbins_reco/2,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_124= new TH2F("conv_template124","",nbins_reco/2,low_reco,high_reco,30,0.,1.); 

	TH2F *conv_template_1_up= new TH2F("conv_template1_up","",nbins_reco/2,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_2_up= new TH2F("conv_template2_up","",nbins_reco/2,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_124_up= new TH2F("conv_template124_up","",nbins_reco/2,low_reco,high_reco,30,0.,1.); 

	TH2F *conv_template_1_dn= new TH2F("conv_template1_dn","",nbins_reco/2,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_2_dn= new TH2F("conv_template2_dn","",nbins_reco/2,low_reco,high_reco,30,0.,1.); 
	TH2F *conv_template_124_dn= new TH2F("conv_template124_dn","",nbins_reco/2,low_reco,high_reco,30,0.,1.); 

	TString phaseName[2]={"fphase_ggH.root","fgraph_vbf_phase.root"};
	TString cName[2]={"cosspline","cos"};
	TString sName[2]={"sinspline","sin"};
	TString gName[2]={"ggZZ","VBF"};

	TFile *fphase_noweight=new TFile(phaseName[sample]);
	TGraph *cosfunc = (TGraph*)fphase_noweight->Get(cName[sample]);
	TGraph *sinfunc = (TGraph*)fphase_noweight->Get(sName[sample]);

	TFile *fxsec=new TFile("xsec.root");
	TGraph *vbfxs= (TGraph*)fxsec->Get("vbf");
	TGraph *whxs= (TGraph*)fxsec->Get("wh");
	TGraph *zhxs= (TGraph*)fxsec->Get("zh");

	//	TFile *bkgeff = new TFile ("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/Reduce/CMSSW_7_1_5/src/160726/bkg_reg_eff.root");
	TFile *bkgeff = new TFile ("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/Efficiency/HighMassCombInputs/efficiency/bkg_reg_eff.root");
	TGraph *effwhat=  (TGraph*)bkgeff->Get(gName[sample]+"_reg_"+chan);
	if(cate_vbf==2){
		TFile *bkgeff_rse = new TFile("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/Efficiency/HighMassCombInputs/efficiency/bkg_rse_eff.root");
		//		if(!sample)
		effwhat = (TGraph*)bkgeff_rse->Get(gName[sample]+"_rse_"+chan);	
		//		else
		//			effwhat = (TGraph*)bkgeff_rse->Get(gName[sample]+"_rse_4e");	
	}
	TFile *sigeff_rse = new TFile("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/Efficiency/HighMassCombInputs/efficiency/eff_sig_rse.root");
	TGraph *effsig_rse = (TGraph*)sigeff_rse->Get(sampleName[sample]+"_"+chan); 

	TAxis *xax = template_1->GetXaxis();
	double inter_all;

	for (int i =0; i< nbins;i++){

		double x = low+ (high-low)/double(nbins)*i;

		mzz->setVal(x);

		double bc =signum*(fabs((mass_d-x))<0.1);
		if(width_d>0.5)
			bc = pdf->getVal(*mzz)*signum*bwid;
		double bc_bkg= genm->GetBinContent(i+1);
		double bc_125= pdf125_hist->GetBinContent(i+1);

		double effval = (effsig[0]+effsig[1]*TMath::Erf( (x-effsig[2])/effsig[3] ))*(effsig[4]+effsig[5]*x+effsig[6]*x*x+effsig[10]*x*x*x)+effsig[7]*TMath::Gaus(x,effsig[8],effsig[9]);
		double effval_bkg = effwhat->Eval(x);
		double m4l = x;
		double effcate_sig, effcate_sig_ori, effcate_bkg;
		if(!sample){
			effcate_sig= 0.0394; //Remember to revert it
			effcate_bkg=  (m4l<93 ? 0.086851-11.5824/(93+63.8731) : 0.086851-11.5824/(m4l+63.8731)); 
			//effcate_bkg= 0.0387026;
			//	effcate_bkg= 0;
			//	effcate_sig= 0;
		}
		else{
			//			effcate_sig = 0.476383 - 4.94232e-05*m4l;
			//			effcate_bkg = (m4l < 200 ? (3.849250e-01 + 2.677220e-03*(m4l-200)) : (m4l < 300 ? 3.849250e-01+1.520690e-04*(m4l-200) + -1.379680e-07*(m4l-200)*(m4l-200) : +1.426950e+02/(m4l+3.971120e+02) - 1.426950e+02/(300+3.971120e+02) + 3.849250e-01+1.520690e-04*100 + -1.379680e-07*100*100 )) ;
			effcate_sig = 0.489193 - 5.21995e-05*m4l;
			effcate_bkg = (m4l < 200 ? (0.392357 + 0.002486*(m4l-200)) : (m4l < 300 ? 0.392357+0.000229678*(m4l-200) + -1.18962e-06*(m4l-200)*(m4l-200) : 171.93/(m4l+437.284) - 171.93/(300+437.284) + 0.392357+0.000229678*100 + -1.18962e-06*100*100 ));
		}
		effcate_sig_ori= effcate_sig;

		double vbffrac=1;
		if(x<2000 && sample)
			vbffrac = vbfxs->Eval(x)*0.955/(whxs->Eval(x)*0.654+zhxs->Eval(x)*0.669+vbfxs->Eval(x)*0.955);

		if(cate_vbf==1){
			effcate_sig= effcate_sig*vbffrac;
			effcate_sig_ori = effcate_sig_ori*vbffrac;
		}
		else if(cate_vbf==0){
			effcate_sig= (1-effcate_sig)*vbffrac + (1-vbffrac)*1./0.7;
			effcate_bkg = 1-effcate_bkg;
			effcate_sig_ori = (1-effcate_sig_ori)*vbffrac+(1-vbffrac);
		}
		else{
			effval = effsig_rse->Eval(x);
			effval_bkg= effval_bkg;
			if(x<300){
				effval=0;
				effval_bkg=0;
			}
			effcate_sig=1;
			effcate_sig_ori=1;
			effcate_bkg=1;
		}

		double kfac = ggZZ_kf->Eval(x);
		if(sample)
			kfac=1;

		//important
		double fa_sig = bc*effval*effcate_sig_ori; // effcate
		double fa_bkg = bc_bkg*effval_bkg*effcate_bkg*kfac; // effcate*kfac  Remember to revert
		//double fa_bkg = bc_bkg*effval_bkg*effcate_bkg; // effcate*kfac
		double fa_125 = bc_125*effval*effcate_sig; // effcate
		//		if(x==125)
		//			cout<< effval<<"\t"<<bc_125<<"\t"<<effcate_sig<<endl;
		//if(x<1000)
		//cout<< x<<"\t"<<fa_bkg<<"\t"<<effval_bkg<<"\t"<<kfac<<"\t"<<effcate_bkg<<endl;


		double a= (x*x-mean->getVal()*mean->getVal());
		double b = mean->getVal()*sigma->getVal();

		double cossig = a/TMath::Sqrt(a*a+b*b);
		double sinsig = b/TMath::Sqrt(a*a+b*b);
		double a125= (x*x-125.*125.);
		double b125 = 125*0.00407;

		double cossig125 = a125/TMath::Sqrt(a125*a125+b125*b125);
		double sinsig125 = b125/TMath::Sqrt(a125*a125+b125*b125);

		double maphase = x;
		if(x>1600 ){
			if(sample&&x>3000)
				maphase=3000;
			else if(!sample)
				maphase=1600.;
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

		double inter_sig_bkgRe = 1.76*sigbkgsqrt*(cosfuncv*cossig);
		double inter_sig_bkgIm = 1.76*sigbkgsqrt*(sinfuncv*sinsig);
		double inter_125_bkgRe = 1.76*hbkgsqrt* (cosfuncv*cossig125);
		double inter_125_bkgIm = 1.76*hbkgsqrt* (sinfuncv*sinsig125); 
		double inter_sig_125 = 2*sighsqrt * (cossig125*cossig-sinsig125*sinsig); 

		//inter_all+= inter_sig_bkgRe+inter_sig_bkgIm+inter_sig_125;
		//cout<< x<<"\t"<<fa_sig<<"\t"<<fa_bkg<<"\t"<<fa_125<<"\t"<<inter_sig_125<<"\t"<<inter_sig_bkgRe<<"\t"<<inter_sig_bkgIm<< "\t"<<inter_all<<endl;



		int binN = xax->FindBin(x);

		double int_1=template_1->Integral(binN,binN);
		double int_2=template_2->Integral(binN,binN);
		double int_124=template_124->Integral(binN,binN);
		double int_124Im=template_124Im->Integral(binN,binN);
		//cout<< binN<<"\t"<<int_1<<"\t"<<int_2<<"\t"<<int_124<<"\t"<<int_124Im<<endl;

		for(int j = 0;j<template_1->GetNbinsY();j++){
			double cont1, cont2, cont124, cont124Im;
			cont1 = template_1->GetBinContent(binN,j+1);
			cont2 = template_2->GetBinContent(binN,j+1);
			cont124 = template_124->GetBinContent(binN,j+1);
			cont124Im = template_124Im->GetBinContent(binN,j+1);


			//			if (int_1!=0)
			template_1r->SetBinContent(i+1,j+1,cont1*fa_sig);
			//			if (int_2!=0&&int_124!=0&&int_124Im!=0)
			template_2r->SetBinContent(i+1,j+1,cont2*fa_bkg+cont1*fa_125+cont124*inter_125_bkgRe+cont124Im*inter_125_bkgIm);
			//			if (int_124!=0&&int_1!=0) 
			if(width_d>0.5)	
				template_124r->SetBinContent(i+1,j+1,cont124*inter_sig_bkgRe+cont124Im*inter_sig_bkgIm+cont1*inter_sig_125);
			else
				template_124r->SetBinContent(i+1,j+1,0);  
			//				if(x<1000)
			//					cout<<cont2*fa_bkg+cont1*fa_125<<"\t"<<cont124*inter_125_bkgRe<<"\t"<<cont124Im*inter_125_bkgIm<<endl;
		}

	}
	cout<< "sig "<<template_1r->Integral()<<endl;
	cout<< "bkg "<<template_2r->Integral(19,31)<<endl;
	cout<< "bkg "<<template_2r->Integral()<<endl;
	cout<< "int"<<template_124r->Integral()<<endl;
	//return;

	for(int m =0;m< nbins; m++){
		double genval = low +bwid*m; //template_1->GetXaxis()->GetBinCenter(m+1);
		mzz->setVal(genval);
		int low_bound = genval + mean_p0_2nd->getVal()-sigma_p0_2nd->getVal()*30 ; 
		int high_bound = genval + mean_p0_2nd->getVal()+sigma_p0_2nd->getVal()*30 ; 
		if(low_bound<low_reco)
			low_bound=low_reco;
		if(high_bound>high_reco)
			high_bound=high_reco;

		//	for(int k = 0;k<nbins_reco;k++){
		for(int k = low_bound;k<high_bound;k++){
			//double recoval = low_reco+(high_reco-low_reco)/double(nbins_reco)*k; 
			double recoval = k; 
			mreco->setVal(recoval);
			double reso = dcrReso_piece.getVal(*mreco); 
			if(reso<1.E-5)
				continue;
			//			cout<<genval<<"\t"<< recoval<<"\t"<<reso<< "\t"<< a1_p0_2nd->getVal()<< "\t"<<a2_p0_2nd->getVal()<<"\t"<<(recoval-genval-mean_p0_2nd->getVal())/sigma_p0_2nd->getVal()<< endl;

			double reso_up = dcrReso_piece_up.getVal(*mreco); 
			double reso_dn = dcrReso_piece_dn.getVal(*mreco); 

			for(int j =0;j< conv_template_1->GetNbinsY(); j++){
				double dbkg_val = conv_template_1->GetYaxis()->GetBinCenter(j+1);

				conv_template_1->Fill(recoval, dbkg_val,template_1r->GetBinContent(m+1,j+1)*reso);
				conv_template_2->Fill(recoval, dbkg_val,template_2r->GetBinContent(m+1,j+1)*reso);
				conv_template_124->Fill(recoval, dbkg_val,template_124r->GetBinContent(m+1,j+1)*reso);

				conv_template_1_up->Fill(recoval, dbkg_val,template_1r->GetBinContent(m+1,j+1)*reso_up);
				conv_template_2_up->Fill(recoval, dbkg_val,template_2r->GetBinContent(m+1,j+1)*reso_up);
				conv_template_124_up->Fill(recoval, dbkg_val,template_124r->GetBinContent(m+1,j+1)*reso_up);

				conv_template_1_dn->Fill(recoval, dbkg_val,template_1r->GetBinContent(m+1,j+1)*reso_dn);
				conv_template_2_dn->Fill(recoval, dbkg_val,template_2r->GetBinContent(m+1,j+1)*reso_dn);
				conv_template_124_dn->Fill(recoval, dbkg_val,template_124r->GetBinContent(m+1,j+1)*reso_dn);
			}
		}
	}
	cout<< "sig "<<conv_template_1->Integral()<<endl;
	cout<< "bkg "<<conv_template_2->Integral()<<endl;
	cout<< "bkg "<<conv_template_2->Integral(5,10)<<endl;
	cout<< "int"<<conv_template_124->Integral()<<endl;
	//	return;
	//for(int i=0;i<nbins;i++){
	//	cout<< i<<"\t"<<template_2r->Integral(i+1,i+1)<<"\t"<<conv_template_2->Integral(i+1,i+1)<<endl;
	//}

	TH1F *conv_template_1_proj = (TH1F*)conv_template_1->ProjectionX();
	TH1F *conv_template_2_proj = (TH1F*)conv_template_2->ProjectionX();
	TH1F *conv_template_124_proj= (TH1F*)conv_template_124->ProjectionX();

	TH1F *conv_template_1_proj_up = (TH1F*)conv_template_1_up->ProjectionX();
	TH1F *conv_template_2_proj_up = (TH1F*)conv_template_2_up->ProjectionX();
	TH1F *conv_template_124_proj_up= (TH1F*)conv_template_124_up->ProjectionX();

	TH1F *conv_template_1_proj_dn = (TH1F*)conv_template_1_dn->ProjectionX();
	TH1F *conv_template_2_proj_dn = (TH1F*)conv_template_2_dn->ProjectionX();
	TH1F *conv_template_124_proj_dn= (TH1F*)conv_template_124_dn->ProjectionX();



	if(is2D){
		RooDataHist *pdfsig_hist_hist= new RooDataHist(Form("pdfsig_hist_hist%d_%s_%s",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),conv_template_1);
		RooHistFunc *pdfsig_hist_func= new RooHistFunc(Form("pdfsig_hist_func%d_%s_%s",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),*pdfsig_hist_hist);

		RooDataHist* hgenm= new RooDataHist(Form("hgenm_%d_%s_%s",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),conv_template_2);
		RooHistFunc* hgenm_func= new RooHistFunc(Form("hgenm_func_%d_%s_%s",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),*hgenm);

		RooDataHist *inter_hist= new RooDataHist(Form("inter_hist_%d_%s_%s",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),conv_template_124);
		RooHistFunc *inter_func= new RooHistFunc(Form("inter_func_%d_%s_%s",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),*inter_hist);

		HZZ4L_RooHighmass *bsi_hist=new HZZ4L_RooHighmass("bsi_hist","bsi_hist",*mreco,*dbkg,*signorm,RooArgList(*pdfsig_hist_func,*hgenm_func,*inter_func)); 
		bsi_hist->SetNameTitle(dataName[sample],dataName[sample]);

		RooDataHist *pdfsig_hist_hist_up= new RooDataHist(Form("pdfsig_hist_hist%d_%s_%s_Resup",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),conv_template_1_up);
		RooHistFunc *pdfsig_hist_func_up= new RooHistFunc(Form("pdfsig_hist_func%d_%s_%s_Resup",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),*pdfsig_hist_hist_up);

		RooDataHist* hgenm_up= new RooDataHist(Form("hgenm_%d_%s_%s_Resup",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),conv_template_2_up);
		RooHistFunc* hgenm_func_up= new RooHistFunc(Form("hgenm_func_%d_%s_%s_Resup",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),*hgenm_up);

		RooDataHist *inter_hist_up= new RooDataHist(Form("inter_hist_%d_%s_%s_Resup",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),conv_template_124_up);
		RooHistFunc *inter_func_up= new RooHistFunc(Form("inter_func_%d_%s_%s_Resup",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),*inter_hist_up);

		HZZ4L_RooHighmass *bsi_hist_up=new HZZ4L_RooHighmass("bsi_hist_Resup","bsi_hist_Resup",*mreco,*dbkg,*signorm,RooArgList(*pdfsig_hist_func_up,*hgenm_func_up,*inter_func_up)); 
		bsi_hist_up->SetNameTitle(dataName[sample]+"_Res"+chan+"Up",dataName[sample]+"Res"+chan+"Up");

		RooDataHist *pdfsig_hist_hist_dn= new RooDataHist(Form("pdfsig_hist_hist%d_%s_%s_Resdn",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),conv_template_1_dn);
		RooHistFunc *pdfsig_hist_func_dn= new RooHistFunc(Form("pdfsig_hist_func%d_%s_%s_Resdn",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),*pdfsig_hist_hist_dn);

		RooDataHist* hgenm_dn= new RooDataHist(Form("hgenm_%d_%s_%s_Resdn",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),conv_template_2_dn);
		RooHistFunc* hgenm_func_dn= new RooHistFunc(Form("hgenm_func_%d_%s_%s_Resdn",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),*hgenm_dn);

		RooDataHist *inter_hist_dn= new RooDataHist(Form("inter_hist_%d_%s_%s_Resdn",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),conv_template_124_dn);
		RooHistFunc *inter_func_dn= new RooHistFunc(Form("inter_func_%d_%s_%s_Resdn",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco,*dbkg),*inter_hist_dn);

		HZZ4L_RooHighmass *bsi_hist_dn=new HZZ4L_RooHighmass("bsi_hist_Resdn","bsi_hist_Resdn",*mreco,*dbkg,*signorm,RooArgList(*pdfsig_hist_func_dn,*hgenm_func_dn,*inter_func_dn)); 
		bsi_hist_dn->SetNameTitle(dataName[sample]+"_Res"+chan+"Down",dataName[sample]+"Res"+chan+"Down");


		w.import(*bsi_hist,RecycleConflictNodes());
		w.import(*bsi_hist_dn,RecycleConflictNodes());
		w.import(*bsi_hist_up,RecycleConflictNodes());
		w.importClassCode("HZZ4L_RooHighmass");
	}
	else{
		RooDataHist *pdfsig_hist_hist= new RooDataHist(Form("pdfsig_hist_hist%d_%s_%s",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),conv_template_1_proj);
		RooHistFunc *pdfsig_hist_func= new RooHistFunc(Form("pdfsig_hist_func%d_%s_%s",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),*pdfsig_hist_hist);

		RooDataHist* hgenm= new RooDataHist(Form("hgenm_%d_%s_%s",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),conv_template_2_proj);
		RooHistFunc* hgenm_func= new RooHistFunc(Form("hgenm_func_%d_%s_%s",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),*hgenm);

		RooDataHist *inter_hist= new RooDataHist(Form("inter_hist_%d_%s_%s",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),conv_template_124_proj);
		RooHistFunc *inter_func= new RooHistFunc(Form("inter_func_%d_%s_%s",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),*inter_hist);

		HZZ4L_RooHighmass_1D *bsi_hist=new HZZ4L_RooHighmass_1D("bsi_hist","bsi_hist",*mreco,*signorm,RooArgList(*pdfsig_hist_func,*hgenm_func,*inter_func)); 
		bsi_hist->SetNameTitle(dataName[sample],dataName[sample]);

		RooDataHist *pdfsig_hist_hist_up= new RooDataHist(Form("pdfsig_hist_hist%d_%s_%s_up",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),conv_template_1_proj_up);
		RooHistFunc *pdfsig_hist_func_up= new RooHistFunc(Form("pdfsig_hist_func%d_%s_%s_up",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),*pdfsig_hist_hist_up);

		RooDataHist* hgenm_up	= new RooDataHist(Form("hgenm_%d_%s_%s_up",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),conv_template_2_proj_up);
		RooHistFunc* hgenm_func_up= new RooHistFunc(Form("hgenm_func_%d_%s_%s_up",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),*hgenm_up);

		RooDataHist *inter_hist_up= new RooDataHist(Form("inter_hist_%d_%s_%s_up",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),conv_template_124_proj_up);
		RooHistFunc *inter_func_up= new RooHistFunc(Form("inter_func_%d_%s_%s_up",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),*inter_hist_up);

		HZZ4L_RooHighmass_1D *bsi_hist_up=new HZZ4L_RooHighmass_1D("bsi_hist_up","bsi_hist_up",*mreco,*signorm,RooArgList(*pdfsig_hist_func_up,*hgenm_func_up,*inter_func_up)); 
		bsi_hist_up->SetNameTitle(dataName[sample]+"_Res"+chan+"Up",dataName[sample]+"_Res"+chan+"Up");

		RooDataHist *pdfsig_hist_hist_dn= new RooDataHist(Form("pdfsig_hist_hist%d_%s_%s_dn",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),conv_template_1_proj_dn);
		RooHistFunc *pdfsig_hist_func_dn= new RooHistFunc(Form("pdfsig_hist_func%d_%s_%s_dn",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),*pdfsig_hist_hist_dn);

		RooDataHist* hgenm_dn	= new RooDataHist(Form("hgenm_%d_%s_%s_dn",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),conv_template_2_proj_dn);
		RooHistFunc* hgenm_func_dn= new RooHistFunc(Form("hgenm_func_%d_%s_%s_dn",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),*hgenm_dn);

		RooDataHist *inter_hist_dn= new RooDataHist(Form("inter_hist_%d_%s_%s_dn",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),conv_template_124_proj_dn);
		RooHistFunc *inter_func_dn= new RooHistFunc(Form("inter_func_%d_%s_%s_dn",cate_vbf,chan.Data(),dataName[sample].Data()),"",RooArgSet(*mreco),*inter_hist_dn);

		HZZ4L_RooHighmass_1D *bsi_hist_dn=new HZZ4L_RooHighmass_1D("bsi_hist_dn","bsi_hist_dn",*mreco,*signorm,RooArgList(*pdfsig_hist_func_dn,*hgenm_func_dn,*inter_func_dn)); 
		bsi_hist_dn->SetNameTitle(dataName[sample]+"_Res"+chan+"Down",dataName[sample]+"_Res"+chan+"Down");


		w.import(*bsi_hist,RecycleConflictNodes());
		w.import(*bsi_hist_up,RecycleConflictNodes());
		w.import(*bsi_hist_dn,RecycleConflictNodes());
		w.importClassCode("HZZ4L_RooHighmass_1D");
	}



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

	int_sig = new ProcessNormalization("int_sig"+chan+"_"+dataName[sample]+Form("%d",cate_vbf),"int_sig"+chan+"_"+dataName[sample]+Form("%d",cate_vbf),int_sig_norm);
	int_bkg = new ProcessNormalization("int_bkg"+chan+"_"+dataName[sample]+Form("%d",cate_vbf),"int_bkg"+chan+"_"+dataName[sample]+Form("%d",cate_vbf),int_bkg_norm); 
	int_int = new ProcessNormalization("int_int"+chan+"_"+dataName[sample]+Form("%d",cate_vbf),"int_int"+chan+"_"+dataName[sample]+Form("%d",cate_vbf),int_int_norm); 

	RooRealVar *sysrr=new RooRealVar("Res"+chan,"Res"+chan,-7,7);
	RooRealVar *kfacrr=new RooRealVar("kfactor_ggZZ","kfactor_ggZZ",-7,7);
	if(int_sig_norm!=0)
		int_sig->addAsymmLogNormal(int_sig_dn/int_sig_norm, int_sig_up/int_sig_norm, *sysrr);
	if(int_bkg_norm!=0){
		int_bkg->addAsymmLogNormal(int_bkg_dn/int_bkg_norm, int_bkg_up/int_bkg_norm, *sysrr);
		//if(sample==0)
		int_bkg->addAsymmLogNormal(0.9, 1.1, *kfacrr);
	}
	if(int_int_norm!=0){
		int_int->addAsymmLogNormal(int_int_dn/int_int_norm, int_int_up/int_int_norm, *sysrr);
		//if(sample==0)
		int_int->addAsymmLogNormal(0.95, 1.05, *kfacrr);
	}

	RooFormulaVar *ggH_norm=new RooFormulaVar(dataName[sample]+"_norm",dataName[sample]+"_norm","@0*@3+@1+@2*sqrt(@3)",RooArgList(*int_sig,*int_bkg,*int_int,*signorm));
	w.import(*ggH_norm,RecycleConflictNodes());

	TFile *fwork ;
	if(width_d<0.1)
		fwork= new TFile(Form("workspace_template%s_4l_Moriond/hzz4l_%s_%d_%dS_13TeV.input_func_%2.0f_%2.2f.root",workName[is2D].Data(),chan.Data(),sample,cate_vbf,mass_d,width_d),"recreate");
	else
		fwork= new TFile(Form("workspace_template%s_4l_Moriond/hzz4l_%s_%d_%dS_13TeV.input_func_%2.0f_%2.1f.root",workName[is2D].Data(),chan.Data(),sample,cate_vbf,mass_d,width_d),"recreate");
	fwork->cd();
	w.Write();
	fwork->Close();

	}
	void construct_4l(double m=450,double w=2,int sample =0,int is2D=1,double inteval=2,int loo=1){
		//void construct_4l(double m=450,double w=2,int sample =0,int is2D=1,TString chan="2e2mu",int cate=1){
		//	gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
		//	gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
		//	gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");

		gROOT->ProcessLine(".x tdrstyle.cc");
		gStyle->SetOptStat(0);
		//		dosomething(chan,m,w,sample,cate,is2D);
		for(int i=0;i<loo;i++){
			dosomething("2e2mu",m,w,sample,0,is2D);
			dosomething("4e",m,w,sample,0,is2D);
			dosomething("4mu",m,w,sample,0,is2D);
			dosomething("2e2mu",m,w,sample,1,is2D);
			dosomething("4e",m,w,sample,1,is2D);
			dosomething("4mu",m,w,sample,1,is2D);
			dosomething("2e2mu",m,w,sample,2,is2D);
			dosomething("4e",m,w,sample,2,is2D);
			m +=inteval; 
		}
	}

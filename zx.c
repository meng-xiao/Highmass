void zx(){
//	TFile *f = new TFile("zx.root","recreate");
	TChain *t=new TChain("candTree");
	t->Add("/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/4lOff/CMSSW_8_0_24_patch1/src/HZZ4l-plotter/ZX.root");
//	int entry=0;
//	int count=0;
//        for(int j =100;j<4000;j+=2){
//                int bentry = t->GetEntries(Form("ZZMass>=%d&&ZZMass<%d",j,j+2));
//                entry +=bentry;
//                if(entry>1000){
//                        cout<< j<<",";
//                        entry=0;
//                        count++;
//                        if(count%10==0)
//                                cout<<endl;
//                }
//        }
//        return;

	int nbins=30;
	double xbin[31]={
		108,114,118,122,126,130,134,138,142,146,
		150,154,158,162,166,170,176,182,188,194,
		200,208,216,224,234,246,260,278,302,338,
		3500
	};

	TH2F *zx_var[3];
	TH2F *zx_Rebin[3];
	TString chanName[3]={"4mu","4e","2e2mu"};

	for(int chan =0;chan<3;chan++){
//		for(int cate=0;cate<2;cate++){
			zx_var[chan]= new TH2F(Form("zxVar_%s",chanName[chan].Data()),"",nbins,xbin,30,0,1.); 
			t->Draw(Form("dbkg_kin:ZZMass>>zxVar_%s",chanName[chan].Data()),Form("(chan==%d)*weight",chan+1));
			cate =0;
			//t->Draw(Form("dbkg_kin:ZZMass>>zxVar_%s",chanName[chan].Data()),Form("(chan==%d&&vbfcate==%d)*weight",chan+1,cate));
			zx_Rebin[chan]= new TH2F(Form("zx_%s",chanName[chan].Data()),"",1695,110.,3500.,30,0,1.); 
			cout<< chanName[chan]<<"\t"<<zx_var[chan]->Integral(1,nbins)<<endl;
			zx_var[chan]->Smooth();

			for(int i = 0; i< zx_Rebin[chan]->GetXaxis()->GetNbins();i++){
				int binx = zx_var[chan]->GetXaxis()->FindBin(zx_Rebin[chan]->GetXaxis()->GetBinCenter(i+1));
				float integral = zx_var[chan]->Integral(binx,binx);
				for(int j = 0; j< zx_Rebin[chan]->GetYaxis()->GetNbins();j++){
					if(integral!=0)
					zx_Rebin[chan]->SetBinContent(i+1,j+1, zx_var[chan]->GetBinContent(binx,j+1)/integral);
					if(zx_Rebin[chan]->GetBinContent(i+1,j+1)==0)
						zx_Rebin[chan]->SetBinContent(i+1,j+1,1.E-10);
				}
			}
			zx_Rebin[chan]->Draw("colz");
			//zx_Rebin[chan]->ProjectionY("y",50,50)->Draw();
			//zx_var[chan][cate]->Draw("colz");
			gPad->Print(Form("%s.png",zx_Rebin[chan]->GetName()));
			//zx_Rebin[chan]->Smooth();
			//zx_Rebin[chan]->ProjectionY("y",50,50)->Draw();
			//gPad->Print(Form("%s_smooth.png",zx_Rebin[chan]->GetName()));
//			f->cd();
//			zx_Rebin[chan]->Write();

		}
//	f->Close();
}

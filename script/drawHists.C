#include "range.h"

const int bins = 5000;

const double npe_sh = 3.0e5;

const char *suffix[3] = {
	"",
	"_tag",
	"_atag"
};

const char epCut[50] = "ep>0.7";
const char edCut[50] = "ed>1.5";
const char dtCut[50] = "dt<400";

const double mmv_cut = 0.2;

void drawHists(){
	
	TFile *f = new TFile("./data/proced.root");

	TTree *tr = f->Get("tr");

	char shCut[255] = "";
	for(int i=0;i<n_range;++i)
		if(range[i] >= npe_sh)
			sprintf(shCut,"%s&&(dtlSH_%d>400||dtlSH_%d<0)",shCut,i,i); 
	for(int i=0;i<3;++i){

		char siteCut[50];
		if(i==0 || i==1)
			sprintf(siteCut,"site==%d",i+1);
		else
			sprintf(siteCut,"site==4");

		for(int j=0;j<n_range;++j){
			for(int s=0;s<3;++s){
				TH1D *h = new TH1D("h","Time since last muon",bins,0,5000);
				char s_tmp[255];
				sprintf(s_tmp,"_%d%s",j,suffix[s]);
				char buf[255];
				sprintf(buf,"dtlSH%s>>h",s_tmp);
				char isoCut[500] = {0};
				//cout << isoCut << strlen(isoCut) << endl;
				for(int k=0;k<n_range;++k){
					if(j==k) continue;
					for(int s1=0;s1<3;++s1){
						char s_tmp1[255];
						sprintf(s_tmp1,"_%d%s",k,suffix[s1]);
						if(strlen(isoCut) == 0)
							sprintf(isoCut,"abs(dtlSH%s-dtlSH%s)>%f",s_tmp,s_tmp1,mmv_cut);
						else
							sprintf(isoCut,"%s&&abs(dtlSH%s-dtlSH%s)>%f",isoCut,s_tmp,s_tmp1,mmv_cut);
					}
				}
				char cut[500];
	//			sprintf(cut,"%s&&dtlSH%s>0&&dtlSH%s<5000&&%s&&%s&&%s",siteCut,s_tmp,s_tmp,epCut,edCut,isoCut);
				sprintf(cut,"%s&&dtlSH%s>0&&dtlSH%s<5000&&%s&&%s&&%s",siteCut,s_tmp,s_tmp,epCut,edCut,dtCut);
				if(range[j] < npe_sh)
					strcat(cut,shCut);
				cout << cut << endl;
				tr->Draw(buf,cut);
				sprintf(buf,"./hists/EH%d_dtlSH%s.root",i+1,s_tmp);
				h->SaveAs(buf);
				sprintf(buf,"./plots/EH%d_dtlSH%s.png",i+1,s_tmp);
				c1->SaveAs(buf);
				h->Delete();
			}
		}
	}

}

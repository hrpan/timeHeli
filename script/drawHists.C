#include "range.h"
#include "../hists/slice.h"

const int bins = 500;
const double histMin = 1.5;
const double histMax = 5000;

const double npe_sh = 3.0e5;

const char *suffix[3] = {
	"",
	"_tag",
	"_atag"
};

//const char *gdCut = "ed>6&&dt<200";
const char *defCut = "ed>4";
const char *nHCut = "ep>1.5&&ep<12&&ed>1.9&&ed<2.7";
const char *unifiedCut = "dist<500";

void drawHists(){
	string filename;
	cin >> filename;
	TFile *f = new TFile(filename.c_str());

	TTree *tr = f->Get("tr");

    TFile *dump = new TFile("dump.root","RECREATE");
    
    //TTree *tr_pre = tr;
    //TTree *tr_pre = tr->CopyTree(nHCut);
    TTree *tr_pre = tr->CopyTree(unifiedCut);
    cout << "pre selection complete" << endl;

	char shCut[255] = "";
	for(int i=0;i<n_range;++i)
		if(range[i] >= npe_sh)
			sprintf(shCut,"%s&&(dtlSH_%d>400||dtlSH_%d<0)",shCut,i,i); 
	for(int site=0;site<3;++site){

		char siteCut[50];
		if(site==0 || site==1)
			sprintf(siteCut,"site==%d",site+1);
		else
			sprintf(siteCut,"site==4");

        TTree *tr_site = tr_pre->CopyTree(siteCut);
        //TTree *tr_site = tr->CopyTree(siteCut);
        cout << "site selection complete" << endl;

		for(int r=0;r<n_range;++r){
			for(int t=0;t<1;++t){
				TH1D *h0 = new TH1D("h", "Time since last muon", bins, histMin, histMax);

				char s_tmp[255];
				sprintf(s_tmp,"_%d%s",r,suffix[t]);

				char buf[255];
				sprintf(buf,"dtlSH%s>>h",s_tmp);

				char cut[500];
				sprintf(cut,"%s&&dtlSH%s>0&&dtlSH%s<5000",siteCut,s_tmp,s_tmp);
				cout << cut << " " << buf << endl;
				tr_site->Draw(buf,cut);
				sprintf(buf,"./hists/EH%d_dtlSH%s.root",site+1,s_tmp);
				h0->SaveAs(buf);
				sprintf(buf,"./plots/EH%d_dtlSH%s.png",site+1,s_tmp);
				c1->SaveAs(buf);
				h0->Delete();

				for(int type=0;type<slice_types;++type){
					double _start = slice_range[type][0];
					double _end = slice_range[type][1];
					double _step = ( _end - _start ) / slices[type];
					char *var_name = slice_vars[type];
					for(int _s=0;_s<slices[type];++_s){
						TH1D *h = new TH1D("h","Time since last muon",bins, histMin, histMax);

						char sliceCut[255];
						sprintf(sliceCut,"%s>%f&&%s<%f", var_name, _start + _s * _step, var_name, _start + (_s + 1) * _step);

						sprintf(cut,"%s&&dtlSH%s>0&&dtlSH%s<5000&&%s",siteCut,s_tmp,s_tmp,sliceCut);

						sprintf(buf,"dtlSH%s>>h",s_tmp);
						cout << cut << " " << buf << endl;
						tr_site->Draw(buf,cut);
						sprintf(buf,"./hists/EH%d_dtlSH%s_%d_%d.root",site+1,s_tmp,type,_s);
						h->SaveAs(buf);
						sprintf(buf,"./plots/EH%d_dtlSH%s_%d_%d.png",site+1,s_tmp,type,_s);
						c1->SaveAs(buf);
						h->Delete();
					}//slice
				}//type
			}//tag
		}//r
	}//site

}

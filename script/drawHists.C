#include "range.h"
#include "slice.h"
const int bins = 100;

const double npe_sh = 3.0e5;

const char *suffix[3] = {
	"",
	"_tag",
	"_atag"
};

const char edCut[50] = "ed>1.5";
const char dtCut[50] = "dt<200";

void drawHists(){
	string filename;
	cin >> filename;
	TFile *f = new TFile(filename.c_str());

	TTree *tr = f->Get("tr");

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

		for(int r=0;r<n_range;++r){
			for(int t=0;t<3;++t){
				for(int _var=0;_var<slice_types;++_var){
					double _start = slice_range[_var][0];
					double _end = slice_range[_var][1];
					double _step = ( _end - _start ) / slices[_var];
					char *var_name = slice_vars[_var];
					for(int _s=0;_s<slices[_var];++_s){
						TH1D *h = new TH1D("h","Time since last muon",bins,1.5,5000);

						char s_tmp[255];
						sprintf(s_tmp,"_%d%s",r,suffix[t]);

						char buf[255];
						sprintf(buf,"dtlSH%s>>h",s_tmp);

						char sliceCut[255];
						sprintf(sliceCut,"%s>%f&&%s<%f", var_name, _start + _s * _step, var_name, _start + (_s + 1) * _step);

						char cut[500];
						sprintf(cut,"%s&&dtlSH%s>0&&dtlSH%s<5000&&%s",siteCut,s_tmp,s_tmp,sliceCut);

						cout << cut << endl;
						tr->Draw(buf,cut);
						sprintf(buf,"./hists/EH%d_dtlSH%s_%d.root",site+1,s_tmp,e);
						h->SaveAs(buf);
						sprintf(buf,"./plots/EH%d_dtlSH%s_%d.png",site+1,s_tmp,e);
						c1->SaveAs(buf);
						h->Delete();
					}
				}
			}//tag
		}//r
	}//site

}

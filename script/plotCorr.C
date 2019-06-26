/*
 * plot dtlSH correlations
 */
#include "range.h"

char *epCut = "ep>3.5";

const double axis_title_size = 0.1;

void plotCorr(){
	string filename;
	cin >> filename;
	TFile *f = new TFile(filename.c_str(),"READ");
	TTree *tr = f->Get("tr");
	char *dir = "./plots/corr/";
	char *suffix[3] = {
		"",
		"_tag",
		"_atag"
	};

	Color_t cls[3] = {
		kRed,
		kBlue,
		kGreen
	};
	char *range_title[3] = {
		"low",
		"mid",
		"high"
	};
	char *tag[3] = {
		" ",
		" tagged ",
		" anti-tagged "
	};
	for(int site=0;site<3;++site){
		char siteCut[255];
		sprintf(siteCut,"site==%d",site==2?4:site+1);
		for(int r0=0;r0<n_range;++r0){
			for(int r1=r0+1;r1<n_range;++r1){
				for(int s0=0;s0<3;++s0){
					for(int s1=0;s1<3;++s1){
						char cut[255];
						sprintf(cut,"%s",siteCut);						
						cout << r0 << " " << r1 << " " << s0 << " " << s1 << endl;
						char buf[512];
						char _sfx0[255];
						sprintf(_sfx0,"_%d%s",r0,suffix[s0]);
						char _sfx1[255];
						sprintf(_sfx1,"_%d%s",r1,suffix[s1]);

						sprintf(buf,"EH%d %s energy%smuon vs %s energy%smuon", site+1, range_title[r0], tag[s0], range_title[r1], tag[s1]);		
						TH2D *h_2d = new TH2D("h_2d", buf, 100, 0, 5000, 100, 0, 5000);
						sprintf(buf,"Time since last %s energy%smuon [ms]", range_title[r0], tag[s0]);
						h_2d->GetYaxis()->SetTitle(buf);
						h_2d->GetYaxis()->SetTitleOffset(1.5);
						sprintf(buf,"Time since last %s energy%smuon [ms]", range_title[r1], tag[s1]);
						h_2d->GetXaxis()->SetTitle(buf);
						sprintf(buf,"dtlSH%s:dtlSH%s>>h_2d",_sfx0,_sfx1);
						tr->Draw(buf,cut,"colz");
						gPad->SetLogz(1);
						gPad->SetLogy(0);
						sprintf(buf,"%sEH%d_dtlSH%s%s_corr.png",dir,site+1,_sfx0,_sfx1);
						c1->SaveAs(buf);
						h_2d->Delete();

						sprintf(buf,"EH%d %s energy%smuon vs %s energy%smuon (difference)", site+1, range_title[r0], tag[s0], range_title[r1], tag[s1]);	
						TH1D *h_1d_0 = new TH1D("h_1d_0", buf, 1000, -500, 500);	
						sprintf(buf,"t_%s%s-t%s%s[ms]", range_title[r0], tag[s0], range_title[r1], tag[s1]);
						h_1d_0->GetXaxis()->SetTitle(buf);
						sprintf(buf,"dtlSH%s-dtlSH%s>>h_1d_0",_sfx0,_sfx1);
						sprintf(cut,"dtlSH%s>0&&dtlSH%s>0&&%s",_sfx0,_sfx1,siteCut);
						tr->Draw(buf,cut);
						gPad->SetLogy(1);
						sprintf(buf,"%sEH%d_dtlSH%s%s_diff.png",dir,site+1,_sfx0,_sfx1);
						c1->SaveAs(buf);
						h_1d_0->Delete();

						sprintf(buf,"EH%d %s energy%smuon vs %s energy%smuon (difference)", site+1, range_title[r0], tag[s0], range_title[r1], tag[s1]);	
						TH1D *h_1d_1 = new TH1D("h_1d_1", buf, 100, -5, 5);	
						sprintf(buf,"t_%s%s-t_%s%s[ms]", range_title[r0], tag[s0], range_title[r1], tag[s1]);
						h_1d_1->GetXaxis()->SetTitle(buf);
						sprintf(buf,"dtlSH%s-dtlSH%s>>h_1d_1",_sfx0,_sfx1);
						tr->Draw(buf,cut);
						gPad->SetLogy(1);
						sprintf(buf,"%sEH%d_dtlSH%s%s_diff_short.png",dir,site+1,_sfx0,_sfx1);
						c1->SaveAs(buf);
						h_1d_1->Delete();	
					}
				}
			}
		}

				
	}		
}

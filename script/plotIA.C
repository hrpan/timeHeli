/*
 * plot inter-arrival times of muons
 */
#include "range.h"

void plotIA(){
	TFile *f = new TFile("./data/proced.root","READ");
	TTree *tr = f->Get("tr");

	char *to_plot[5] = {
		"dtlSH_sdt",
		"dtlSH_sdt",
		"dtlSH_dt",
		"dtlSH_dt",
		"dtlSH_nPESum"
	};

	int h_opts[5][3] = {
		{ 1000, 0, 5000	},
		{ 100, 0, 3 },
		{ 1000, 0, 5000	},
		{ 100, 0, 3	},
		{ 1000, 3e3, 8e5 }
	};

	char *output_names[5] = {
		"sdt_long",
		"sdt_short",
		"dt_long",
		"dt_short",
		"nPESum"
	};

	bool plot_overlap[5] = {
		false,
		true,
		true,
		true,
		true
	};

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

	for(int site=0;site<3;++site){

		char siteCut[255];
		sprintf(siteCut,"site==%d",site==2?4:site+1);

		for(int p=0;p<5;++p){
			for(int r=0;r<n_range;++r){
				if(plot_overlap[p]){
					char plot[255];
					for(int s=0;s<3;++s){
						sprintf(plot,"%s_%d%s>>h%d(%d,%d,%d)", 
								to_plot[p], r, suffix[s], s,
								h_opts[p][0], h_opts[p][1], h_opts[p][2]);
						
						tr->Draw(plot, siteCut);
						gPad->SetLogy();
					}	
					h0->SetLineColor(cls[0]);
					h1->SetLineColor(cls[1]);
					h2->SetLineColor(cls[2]);

					h0->DrawNormalized();
					h1->DrawNormalized("same");
					h2->DrawNormalized("same");

					char output[255];
					sprintf(output,"./plots/EH%d_%s_%d.png", site+1, output_names[p], r);
					c1->SaveAs(output);
				}else{
					for(int s=0;s<3;++s){
						char plot[255];
						sprintf(plot,"%s_%d%s>>h(%d,%d,%d)", 
								to_plot[p], r, suffix[s],
								h_opts[p][0], h_opts[p][1], h_opts[p][2]);

						tr->Draw(plot, siteCut);
						gPad->SetLogy();
						
						char output[255];
						sprintf(output,"./plots/EH%d_%s_%d%s.png", site+1, output_names[p], r, suffix[s]);
						c1->SaveAs(output);
					}	
				}
			}
		}				
	}		
}

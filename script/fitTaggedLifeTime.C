#include "range.h"

void printResults(double results[3][n_range][6]){

	printf("%4s %5s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n","site","range","n_sin","n_sin_err","n_gd","n_gd_err","n_h","n_h_err","n_gd+h","gd/h","gd/h_err");
	for(int i=0;i<3;++i){
		for(int j=0;j<n_range;++j){
			double ratio = results[i][j][2] / results[i][j][4];
			double _gd = results[i][j][3] / results[i][j][2];
			double _h = results[i][j][5] / results[i][j][4];
			double ratio_err = ratio * sqrt( _gd * _gd + _h * _h );

			printf("%4d %5d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
					i,
					j,
					results[i][j][0],
					results[i][j][1],
					results[i][j][2],
					results[i][j][3],
					results[i][j][4],
					results[i][j][5],
					results[i][j][2] + results[i][j][4],
					ratio,
					ratio_err);
		}
	}
}

void printLatexTable(double results[3][n_range][6]){
	for(int i=0;i<3;++i){
		for(int j=0;j<n_range;++j){
			double ratio = results[i][j][2] / results[i][j][4];
			double _gd = results[i][j][3] / results[i][j][2];
			double _h = results[i][j][5] / results[i][j][4];
			double ratio_err = ratio * sqrt( _gd * _gd + _h * _h );

			if(j==n_range-1)
				printf("%4d & [%10.1e,%10s] & $%10.2f\\pm%10.2f$ & $%10.2f\\pm%10.2f$ & $%10.2f\\pm%10.2f$ & $%10.2f\\pm%10.2f$\n",i+1,range[j],"inf",results[i][j][0],results[i][j][1],results[i][j][2],results[i][j][3],results[i][j][4],results[i][j][5],ratio,ratio_err);
			else
				printf("%4d & [%10.1e,%10.1e] & $%10.2f\\pm%10.2f$ & $%10.2f\\pm%10.2f$ & $%10.2f\\pm%10.2f$ & $%10.2f\\pm%10.2f$\n",i+1,range[j],range[j+1],results[i][j][0],results[i][j][1],results[i][j][2],results[i][j][3],results[i][j][4],results[i][j][5],ratio,ratio_err);
			fflush(stdout);
		}

	}

}


void fitTaggedLifeTime(){
	gStyle->SetOptFit(0111);
	//TFile *f = new TFile("./data/muons_split.root");
	TFile *f = new TFile("./data/tagged_dt.root");
	TList *l = f->GetListOfKeys();

	TF1 *func = new TF1("func","[0]+([1]/[2])*exp(-x/[2])+([3]/[4])*exp(-x/[4])",0,800);
	func->SetParNames("const.","n_GD","GD lifetime","n_H","H lifetime");
	func->SetParLimits(0,0,1e5);
	func->SetParLimits(1,1e3,1e7);
	func->SetParLimits(2,20,40);
	func->SetParLimits(3,1e3,1e7);
	func->SetParLimits(4,100,300);

	double results[3][n_range][6];

	char buf[255];

	for(size_t i=0;i<l->GetEntries();++i){
		char *key = l->At(i)->GetName();
		cout << key << endl;
		TH1 *h = f->Get(key);
		h->SetAxisRange(20,800,"X");
		h->Fit(func,"L","",20,800);

		results[i/n_range][i%n_range][0] = func->GetParameter(0);
		results[i/n_range][i%n_range][1] = func->GetParError(0);
		results[i/n_range][i%n_range][2] = func->GetParameter(1);
		results[i/n_range][i%n_range][3] = func->GetParError(1);
		results[i/n_range][i%n_range][4] = func->GetParameter(3);
		results[i/n_range][i%n_range][5] = func->GetParError(3);		
		
		sprintf(buf,"./plots/%s.png",key);



		h->Draw();
		func->Draw("same");

		gPad->SetLogy();
		c1->SaveAs(buf);	
	}

	/*
	for(int i=0;i<3;++i){
		cout << "Site: " << i+1 << endl;
		for(int j=0;j<n_range;++j){
			cout << "Range: " << j << endl;
			TH1 *h = new TH1D("h","h",800,0,800);
			if(j==n_range-1)
				sprintf(buf,"site==%d&&nPESum>%f",i==2?4:i+1,range[j]);
			else
				sprintf(buf,"site==%d&&nPESum<%f&&nPESum>%f",i==2?4:i+1,range[j+1],range[j]);
			tr->Draw("nTag_dt>>h",buf);
			h->Fit(func,"L","same",20,800);	

			results[i][j][0] = func->GetParameter(0);
			results[i][j][1] = func->GetParError(0);
			results[i][j][2] = func->GetParameter(1);
			results[i][j][3] = func->GetParError(1);
			results[i][j][4] = func->GetParameter(3);
			results[i][j][5] = func->GetParError(3);

			if(j==n_range-1)
				sprintf(buf,"EH%d tagged_dt fit %.2e<npe",i+1,range[j]);
			else
				sprintf(buf,"EH%d tagged_dt fit %.2e<npe<%.2e",i+1,range[j],range[j+1]);
			h->SetTitle(buf);
			sprintf(buf,"./plots/tagged_dt_fit_%d_%d.png",i,j);
			gPad->SetLogy();
			c1->SaveAs(buf);
			h->Delete();
		}
	}
*/
	printLatexTable(results);
	printResults(results);

}

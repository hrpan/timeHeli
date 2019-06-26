#include "range.h"

const char *hists_prefix = "./hists";

const double tau_li = 257.23;
const double tau_he = 171.68;
const double tau_b12 = 29.142;
const double tau_n12 = 15.9;

const double fitMin = 1.5;
const double fitMax = 5000;


void printResults(double results[3][n_range][2][4]){
	printf("%4s%20s%10s%10s%10s%10s%10s%10s%10s%10s\n","site","range","N_Li","E_Li","N_Li_Tag","E_Li_Tag","eff.","eff. err","R_mu","R_mu err");
	for(int i=0;i<3;++i){
		for(int j=0;j<n_range;++j){
			double _n = results[i][j][0][0];
			double _n_err = results[i][j][0][1];
			double _n_tag = results[i][j][1][0];
			double _n_tag_err = results[i][j][1][1];
			double eff = _n_tag / _n;
			double eff_err = eff * sqrt( (_n_err/_n) * (_n_err/_n) + (_n_tag_err/_n_tag) * (_n_tag_err/_n_tag));
			double r_mu = results[i][j][0][2];
			double r_mu_err = results[i][j][0][3];
			if(j==n_range-1)
				printf("%4d%10.2e%10.2e%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.3e%10.3e\n",i+1,range[j],-1,_n,_n_err,_n_tag,_n_tag_err,eff,eff_err,r_mu,r_mu_err);
			else
				printf("%4d%10.2e%10.2e%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.3e%10.3e\n",i+1,range[j],range[j+1],_n,_n_err,_n_tag,_n_tag_err,eff,eff_err,r_mu,r_mu_err);
		}
	}
}

void printLatexTable(double results[3][n_range][2][4]){
	for(int i=0;i<3;++i){
		for(int j=0;j<n_range;++j){
			double _n = results[i][j][0][0];
			double _n_err = results[i][j][0][1];
			double _n_tag = results[i][j][1][0];
			double _n_tag_err = results[i][j][1][1];
			double eff = _n_tag / _n;
			double eff_err = eff * sqrt( (_n_err/_n) * (_n_err/_n) + (_n_tag_err/_n_tag) * (_n_tag_err/_n_tag));
			if(j==n_range-1)
				printf("%4d & [%10.1e,%10s] & $%10.2f\\pm%10.2f$ & $%10.2f\\pm%10.2f$ & $%10.2f\\pm%10.2f$\n",i+1,range[j],"inf",_n,_n_err,_n_tag,_n_tag_err,eff,eff_err);
			else
				printf("%4d & [%10.2e,%10.2e] & $%10.2f\\pm%10.2f$ & $%10.2f\\pm%10.2f$ & $%10.2f\\pm%10.2f$\n",i+1,range[j],range[j+1],_n,_n_err,_n_tag,_n_tag_err,eff,eff_err);
		}

	}

}

char *fitOptions = "LM";

void fit(){

	gStyle->SetOptFit(1);
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
	//[0]:mu rate
	//[1]:n_ibd
	//[2]:n_li9
	//[3]:li9 lifetime
	//[4]:n_he8
	//[5]:he8 lifetime
	//[6]:
	char _f_dc[255] = "[1] * [0] * exp(-[0] * x )";
	char _f_li[255] = "[2] * ([0] + 1 / [3]) * exp(-([0] + 1/[3]) * x)";
	char _f_he[255] = "[4] * ([0] + 1 / [5]) * exp(-([0] + 1/[5]) * x)";
	char _f_bo[255] = "[6] * ([0] + 1 / [7]) * exp(-([0] + 1/[7]) * x)";
	char _f_ni[255] = "[8] * ([0] + 1 / [9]) * exp(-([0] + 1/[9]) * x)";
	char _f_sum[255];
	sprintf(_f_sum,"%s+%s+%s+%s+%s",_f_dc,_f_li,_f_he,_f_bo,_f_ni);
	TF1 *func = new TF1("func",_f_sum,0,5000);

	func->SetParName(0,"mu rate");
	func->SetParName(1,"N_uncorr");
	func->SetParName(2,"N_Li9");
	func->SetParName(3,"li9 lifetime");
	func->SetParName(4,"N_He8");
	func->SetParName(5,"He8 lifetime");
	func->SetParName(6,"N_B12");
	func->SetParName(7,"b12 lifetime");
	func->SetParName(8,"N_N12");
	func->SetParName(9,"N12 lifetime");

	func->SetParLimits(0,0,0.1);
	func->SetParLimits(2,0,1e6);
	func->SetParLimits(4,0,1e6);
	func->SetParLimits(6,0,1e6);
	func->SetParLimits(8,0,1e6);
	func->SetParameter(1, 1e5);
	func->SetParameter(2, 1e3);
	//func->SetParLimits(1,0,1e8);
	//func->SetParLimits(2,0,1e5);
	
	func->FixParameter(3,tau_li);
	func->FixParameter(4,0);
	func->FixParameter(5,tau_he);
	//func->FixParameter(6,0);
	//func->SetParLimits(6,0,1e5);
	func->FixParameter(7,tau_b12);
	func->FixParameter(9,tau_n12);
	func->SetLineColor(kRed);

	TF1 *funcDC = new TF1("funcDC",_f_dc,0,5000);
	funcDC->SetLineColor(kGreen);
	
	double results[3][n_range][2][4];
	for(int site=0;site<3;++site){
		for(int range=0;range<n_range;++range){
			char buf[255];
			sprintf(buf,"%s/EH%d_dtlSH_%d.root", hists_prefix, site+1, range);
			cout << "Fitting: " << buf << endl;
			TFile *f = new TFile(buf,"READ");
			TH1 *h = f->Get("h");
            h->SetTitle("Time since last muon;[ms]");
			h->Draw("E");

//			func->SetParameter(1,h->GetEntries());
//			func->FixParameter(2,0);
//			h->Fit(func,fitOptions,"E",fitMin,fitMax);
//			func->ReleaseParameter(2);
			func->SetParameter(0,1e-3);
			h->Fit(func,fitOptions,"E",fitMin,fitMax);
	//		c1->SetLogy();
			c1->SetLogx();
			funcDC->SetParameters(func->GetParameters());
			funcDC->Draw("same");
			sprintf(buf,"./plots/fits/EH%d_%d_fit.png",site+1,range);
			c1->SaveAs(buf);
			f->Close();
	
			results[site][range][0][0] = func->GetParameter(2);
			results[site][range][0][1] = func->GetParError(2);
			results[site][range][0][2] = func->GetParameter(0);
			results[site][range][0][3] = func->GetParError(0);
/*	
			sprintf(buf,"%s/EH%d_dtlSH_%d_tag.root", hists_prefix, site+1, range);
			cout << "Fitting: " << buf << endl;
			f = new TFile(buf,"READ");
			h =(TH1*) f->Get("h");
			h->Draw("E");

//			func->SetParameter(1,h->GetEntries());
//			func->FixParameter(2,0);
//			h->Fit(func,fitOptions,"E",fitMin,fitMax);
//			func->ReleaseParameter(2);
			func->SetParameter(0,1e-3);
			h->Fit(func,fitOptions,"E",fitMin,fitMax);
	//		c1->SetLogy();
			c1->SetLogx();
			funcDC->SetParameters(func->GetParameters());
			funcDC->Draw("same");
			sprintf(buf,"./plots/fits/EH%d_%d_tag_fit.png",site+1,range);
			c1->SaveAs(buf);
			f->Close();

			results[site][range][1][0] = func->GetParameter(2);
			results[site][range][1][1] = func->GetParError(2);
			results[site][range][1][2] = func->GetParameter(0);
			results[site][range][1][3] = func->GetParError(0);

			sprintf(buf,"%s/EH%d_dtlSH_%d_atag.root", hists_prefix, site+1, range);
			cout << "Fitting: " << buf << endl;
			f = new TFile(buf,"READ");
			h =(TH1*) f->Get("h");
			h->Draw("E");

			func->SetParameter(0,1e-3);
			h->Fit(func,fitOptions,"E",fitMin,fitMax);
		//	c1->SetLogy();
			c1->SetLogx();
			funcDC->SetParameters(func->GetParameters());
			funcDC->Draw("same");
			sprintf(buf,"./plots/fits/EH%d_%d_atag_fit.png",site+1,range);
			c1->SaveAs(buf);
*/
			f->Close();

		
		}
	}
	printLatexTable(results);
	printResults(results);
}



#include "range.h"
const double tau_li9 = 256.366;
const double tau_li9_err = 0.866;

const double tau_he8 = 171.68;

const double tau_b12 = 29.142;
const double tau_b12_err = 0.0288;

const double minimizer_args[2] = { 100000, 1e-6};
//const double fitMin = 10;
//const double fitMax = 5000;

const int npars = 12;

const char *par_names[npars] = {
	"r_mu_tag", "r_mu_atag", "N_DC",
	"N_Li", "EPS_Li", "tau_Li",
	"N_Bo", "EPS_Bo", "tau_Bo",
	"N_He", "EPS_He", "tau_He"
}; 

const double par_init[npars][3] = {
	{ 1e-4, 0, 0 },
	{ 1e-4, 0, 1 },
	{ 1e6, 0, 0 },
	{ 1e3, 0, 1e6 },
	{ 0.5, 0, 1 },
	{ tau_li9, 0, 0 },
	{ 1e3, 0, 1e6 },
	{ 0.5, 0, 1 },
	{ tau_b12, 0, 0 },
	{ 1e3, 0, 1e6 },
	{ 0.5, 0, 1 },
	{ tau_he8, 0, 0 },
};

const bool result_print[npars] = {
	true, true, true,
	true, true, false,
	true, true, false,
	true, true, false
};

const bool latex_print[npars] = {
	false, false, false,
	true, true, false,
	true, true, false,
	true, true, false
};

const int npars_single = 8;
//[0]:mu rate
//[1]:n_dc
//[2]:n_li9
//[3]:t_li9
//[4]:n_b12
//[5]:t_b12
//[6]:n_he8
//[7]:t_he8

const char *_pdfs[4] = {
	"[1] * [0] * exp(-[0] * x )",
	"[2] * ([0] + 1 / [3]) * exp(-([0] + 1 / [3]) * x )",
	"[4] * ([0] + 1 / [5]) * exp(-([0] + 2 / [5]) * x )",
	"[6] * ([0] + 1 / [7]) * exp(-([0] + 1 / [7]) * x )"
};
/*
const char _f_dc[255] = "[1] * [2] * exp(-[2] * x )";
const char _f_li[255] = "[0] * ([2] + 1 / [4]) * exp(-([2] + 1 / [4]) * x )";
const char _f_bo[255] = "[3] * ([2] + 2 / [5] ) * exp(-([2] + 2 / [5] ) * x )";
const char _f_he[255] = "[6] * ([2] + 1 / [7]) * exp(-([2] + 1 / [7]) * x )";
*/
double pull(double x, double mean, double err){

	double diff = x - mean;

	return diff * diff / ( err * err);

}

void printResults(double results[3][n_range][npars][2]){
	char buf[512];
	sprintf(buf, "%4s %20s", "site", "range");
	for(int i=0;i<npars;++i){
		if(result_print[i])
			sprintf(buf, "%s %17s", buf, par_names[i]);
	}
	printf("%s\n", buf);	
	for(int i=0;i<3;++i){
		for(int j=0;j<n_range;++j){

			sprintf(buf, "%4d %10.2e%10.2e", i+1, range[j], j==n_range-1? -1:range[j+1]);

			for(int k=0;k<npars;++k){
				if(result_print[k])
					sprintf(buf, "%s %8.2e/%8.2e", buf, results[i][j][k][0], results[i][j][k][1]);
			}
			printf("%s\n",buf);
		}
		printf("================\n");
	}
}

void printLatexTable(double results[3][n_range][npars][2]){
	for(int i=0;i<3;++i){
		for(int j=0;j<n_range;++j){
			if(j==n_range-1)
				printf("%4d & [%10.1e,%10s] ",i+1,range[j],"inf");
			else
				printf("%4d & [%10.1e,%10.1e] ",i+1,range[j],range[j+1]);
			for(int k=0;k<npars;++k)
				if(latex_print[k])
					printf("& $%10.2f\\pm%10.2f$ ",results[i][j][k][0],results[i][j][k][1]);
			
			printf("\\\\ \n");
		}
		printf("\\hline \n");
	}

}

char *fitOptions = "LM";

ROOT::Math::IMultiGenFunction *fPNLL_1;
ROOT::Math::IMultiGenFunction *fPNLL_2;
ROOT::Math::IMultiGenFunction *fPNLL_3;

TF1 *func1;
TF1 *func2;
TF1 *func3;

TH1 *h1;
TH1 *h2;
TH1 *h3;

void wrap(int &npar, double *g, double &result, double *par, int flag){
/*const char *par_names[npars] = {
	"r_mu_tag", "r_mu_atag", "N_DC",
	"N_Li", "EPS_Li", "tau_Li",
	"N_Bo", "EPS_Bo", "tau_Bo",
	"N_He", "EPS_He", "tau_He"
};*/ 
	double r_mu_tag = par[0];
	double r_mu_atag = par[1];
	double n_dc = par[2];
	double n_li = par[3];
	double eps_li = par[4];
	double tau_li = par[5];
	double n_bo = par[6];
	double eps_bo = par[7];
	double tau_bo = par[8];
	double n_he = par[9];
	double eps_he = par[10];
	double tau_he = par[11];

	//[0]:mu rate
	//[1]:n_dc
	//[2]:n_li9
	//[3]:t_li9
	//[4]:n_b12
	//[5]:t_b12
	//[6]:n_he8
	//[7]:t_he8


	double par_0[npars_single];
	par_0[0] = r_mu_tag + r_mu_atag;
	par_0[1] = n_dc;
	par_0[2] = n_li;
	par_0[3] = tau_li;	
	par_0[4] = n_bo;
	par_0[5] = tau_bo;
	par_0[6] = n_he;
	par_0[7] = tau_he;

	double par_1[npars_single];
	par_1[0] = r_mu_tag;
	par_1[1] = n_dc + (1-eps_li) * n_li + (1-eps_bo) * n_bo + (1-eps_he) * n_he;
	par_1[2] = eps_li * n_li;
	par_1[3] = tau_li;	
	par_1[4] = eps_bo * n_bo;
	par_1[5] = tau_bo;
	par_1[6] = eps_he * n_he;
	par_1[7] = tau_he;

	double par_2[npars_single];
	par_2[0] = r_mu_atag;
	par_2[1] = n_dc + eps_li * n_li + eps_bo * n_bo + eps_he * n_he;
	par_2[2] = (1-eps_li) * n_li;
	par_2[3] = tau_li;
	par_2[4] = (1-eps_bo) * n_bo;
	par_2[5] = tau_bo;
	par_2[6] = (1-eps_he) * n_he;
	par_2[7] = tau_he;

//	double pull_term = pull(par[7], tau_li, tau_li_err) + pull(par[8], tau_b12, tau_b12_err);
	double pull_term = 0;

	result = (*fPNLL_1)(par_0) + (*fPNLL_2)(par_1) + (*fPNLL_3)(par_2) + pull_term;
}


void fillPars(double *par, TF1 *f_0, TF1 *f_1, TF1 *f_2){
	double r_mu_tag = par[0];
	double r_mu_atag = par[1];
	double n_dc = par[2];
	double n_li = par[3];
	double eps_li = par[4];
	double tau_li = par[5];
	double n_bo = par[6];
	double eps_bo = par[7];
	double tau_bo = par[8];
	double n_he = par[9];
	double eps_he = par[10];
	double tau_he = par[11];


	double par_0[npars_single];
	par_0[0] = r_mu_tag + r_mu_atag;
	par_0[1] = n_dc;
	par_0[2] = n_li;
	par_0[3] = tau_li;	
	par_0[4] = n_bo;
	par_0[5] = tau_bo;
	par_0[6] = n_he;
	par_0[7] = tau_he;

	double par_1[npars_single];
	par_1[0] = r_mu_tag;
	par_1[1] = n_dc + (1-eps_li) * n_li + (1-eps_bo) * n_bo + (1-eps_he) * n_he;
	par_1[2] = eps_li * n_li;
	par_1[3] = tau_li;	
	par_1[4] = eps_bo * n_bo;
	par_1[5] = tau_bo;
	par_1[6] = eps_he * n_he;
	par_1[7] = tau_he;

	double par_2[npars_single];
	par_2[0] = r_mu_atag;
	par_2[1] = n_dc + eps_li * n_li + eps_bo * n_bo + eps_he * n_he;
	par_2[2] = (1-eps_li) * n_li;
	par_2[3] = tau_li;
	par_2[4] = (1-eps_bo) * n_bo;
	par_2[5] = tau_bo;
	par_2[6] = (1-eps_he) * n_he;
	par_2[7] = tau_he;

	f_0->SetParameters(par_0);

	f_1->SetParameters(par_1);

	f_2->SetParameters(par_2);	

}

void plotHists(int site, int range, TH1 *h1, TF1 *f1, TH1 *h2, TF1 *f2, TH1 *h3, TF1 *f3){
	
	char buf[255];

	sprintf(buf, "./plots/cfit_%d_%d_0.png", site, range);
	h1->Draw("E1");
	f1->Draw("same");
	gPad->SetLogy();
	gPad->SetLogx();
	c1->SaveAs(buf);
	sprintf(buf, "./plots/cfit_%d_%d_1.png", site, range);
	h2->Draw("E1");
	f2->Draw("same");
	c1->SaveAs(buf);
	sprintf(buf, "./plots/cfit_%d_%d_2.png", site, range);
	h3->Draw("E1");
	f3->Draw("same");
	c1->SaveAs(buf);
}

struct Config{
	double fitMin;
	double fitMax;

	int fix_lifetime;
	int fix_B12;
	int fix_He8;

	int bound_eps;
};

void parseInputs(Config &cfg){

	cin >> cfg.fitMin;
	cout << "fitMin: " << cfg.fitMin << endl;

	cin >> cfg.fitMax;
	cout << "fitMax: " << cfg.fitMax << endl;

	cin >> cfg.fix_B12;
	cout << "fix B12: " << cfg.fix_B12 << endl;

	cin >> cfg.fix_He8;
	cout << "fix He8: " << cfg.fix_He8 << endl; 

	cin >> cfg.bound_eps;
	cout << "Bound tagging efficiencies: " << cfg.bound_eps << endl;

	cin >> cfg.fix_lifetime;
	cout << "Fix isotope lifetimes: " << cfg.fix_lifetime << endl;
}

void fitterParInit(int site, TFitter &minuit, Config &cfg){

	double step_ratio = 0.01;

	for(int i=0;i<npars;++i){
		if(cfg.fix_lifetime == 1 && strstr(par_names[i],"tau") != NULL){
			minuit.SetParameter(i,par_names[i],par_init[i][0],0,0,0);
			minuit.FixParameter(i);
		}else if( (cfg.fix_B12 == 1 && strstr(par_names[i],"Bo") != NULL) ||
				  (cfg.fix_He8 == 1 && strstr(par_names[i],"He") != NULL) ){
			minuit.SetParameter(i,par_names[i],0,0,0,0);
			minuit.FixParameter(i);
		}else if(cfg.bound_eps == 0 && strstr(par_names[i],"EPS") != NULL){
			minuit.SetParameter(i,par_names[i],par_init[i][0], par_init[i][0]*step_ratio,0,0);
		}else{
			minuit.SetParameter(i,par_names[i],par_init[i][0], par_init[i][0]*step_ratio,par_init[i][1],par_init[i][2]);
		}

	}
}

void fix_and_release(TFitter &minuit, Config &cfg){

	double _tmp[2] = {3, 0};
	minuit.ExecuteCommand("SET PARAMETER", _tmp, 2);
	minuit.FixParameter(2);
	minuit.FixParameter(4);

	if(!cfg.fix_B12){
		_tmp[0] = 6;
		minuit.ExecuteCommand("SET PARAMETER", _tmp, 2);
		minuit.FixParameter(5);
		minuit.FixParameter(6);
	}	

	minuit.ExecuteCommand("MINIMIZE", minimizer_args, 2);

	minuit.ReleaseParameter(2);
	minuit.ReleaseParameter(4);
	minuit.ExecuteCommand("MINIMIZE", minimizer_args, 2);

	if(!cfg.fix_B12){
		minuit.ReleaseParameter(5);
		minuit.ReleaseParameter(6);
		minuit.ExecuteCommand("MINIMIZE", minimizer_args, 2);
	}
	
	minuit.ExecuteCommand("MINOS", minimizer_args, 1);


}

void combinedFit(){

	Config cfg;

	parseInputs(cfg);

	gStyle->SetOptFit(1);
	char _f_sum[512];
	for(int i=0;i<4;++i)
		if(i==0)
			sprintf(_f_sum,"%s",_pdfs[i]);
		else
			sprintf(_f_sum,"%s+%s",_f_sum,_pdfs[i]);
//	sprintf(_f_sum,"%s+%s+%s+%s",_f_dc,_f_li,_f_bo,_f_he);
	
	func1 = new TF1("func1",_f_sum,0,5000);
	func2 = new TF1("func2",_f_sum,0,5000);
	func3 = new TF1("func3",_f_sum,0,5000);

	ROOT::Math::WrappedMultiTF1 wf1(*func1,1);
	ROOT::Math::WrappedMultiTF1 wf3(*func2,1);
	ROOT::Math::WrappedMultiTF1 wf2(*func3,1);
	
	double results[3][n_range][npars][2];
	
	ROOT::Fit::DataOptions opt;
	opt.fUseEmpty = true;
	ROOT::Fit::DataRange rangeD;
	rangeD.SetRange(cfg.fitMin, cfg.fitMax);

	TFitter minuit(npars);
	double p0 = 1;
	minuit.ExecuteCommand("SET PRINTOUT",&p0,1);
	p0 = 2;
	minuit.ExecuteCommand("SET STRATEGY",&p0,1);
	p0 = 0.5;
	minuit.ExecuteCommand("SET ERRORDEF",&p0,1);
	p0 = 1e-16;
	minuit.ExecuteCommand("SET EPSMACHINE",&p0,1);
	for(int site=0;site<3;++site){
		for(int range=0;range<n_range;++range){


			char buf[255];
			sprintf(buf,"./hists/EH%d_dtlSH_%d.root",site+1,range);
			TFile *f1 = new TFile(buf,"READ");
			h1 = (TH1*)f1->Get("h");
			ROOT::Fit::BinData data1(opt,rangeD);
			ROOT::Fit::FillData(data1, h1);
			ROOT::Fit::PoissonLLFunction PNLL1(data1,wf1);
			fPNLL_1 = &PNLL1;

			sprintf(buf,"./hists/EH%d_dtlSH_%d_tag.root",site+1,range);
			TFile *f2 = new TFile(buf,"READ");
			h2 = (TH1*)f2->Get("h");
			ROOT::Fit::BinData data2(opt,rangeD);
			ROOT::Fit::FillData(data2, h2);
			ROOT::Fit::PoissonLLFunction PNLL2(data2,wf2);
			fPNLL_2 = &PNLL2;

			sprintf(buf,"./hists/EH%d_dtlSH_%d_atag.root",site+1,range);
			TFile *f3 = new TFile(buf,"READ");
			h3 = (TH1*)f3->Get("h");
			ROOT::Fit::BinData data3(opt,rangeD);
			ROOT::Fit::FillData(data3, h3);
			ROOT::Fit::PoissonLLFunction PNLL3(data3,wf3);
			fPNLL_3 = &PNLL3;

			minuit.SetFCN(wrap);
			fitterParInit(site, minuit, cfg);
		
			minuit.ExecuteCommand("SIMPLEX", minimizer_args, 2);
			minuit.ExecuteCommand("MINIMIZE", minimizer_args, 2);
			minuit.ExecuteCommand("MINOS", minimizer_args, 1);
			//fix_and_release(minuit, cfg);
	
			double _pars[npars];
			for(int _i=0;_i<npars;++_i){
				_pars[_i] = minuit.GetParameter(_i);
			}
			fillPars(_pars, func1, func2, func3);

			plotHists(site, range, h1, func1, h2, func2, h3, func3);

			for(int _i=0;_i<npars;++_i){
				results[site][range][_i][0] = minuit.GetParameter(_i);
				results[site][range][_i][1] = minuit.GetParError(_i);
			}
		
		}
	}
	printLatexTable(results);
	printResults(results);
}



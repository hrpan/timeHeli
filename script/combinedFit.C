#include "range.h"

const size_t max_bins = 100000;

const double tau_li9 = 256.366;
const double tau_li9_err = 0.866;

const double tau_he8 = 171.68;

const double tau_b12 = 29.142;
const double tau_b12_err = 0.0288;

const double simplex_args[2] = {1000000, 0.1};
const double minimizer_args[2] = {1000000, 0.1};
const double minos_args[1] = {1000000};
const double improve_args[1] = {1000000};
const double seek_args[2] = {1000, 1.0};
double fitMin = 10;
double fitMax = 5000;

const int npars = 11;

const char *par_names[npars] = {
	"r_mu_tag", "r_mu_atag", "N_DC",
	"N_Li/He", "EPS_Li/He", "ratio_Li/He", "tau_Li", "tau_He",
	"N_Bo", "EPS_Bo", "tau_Bo",
}; 

const double par_init[npars][3] = {
	{ 1e-4, 0, 1 },
	{ 1e-4, 0, 1 },
	{ 1e5, 0, 0 },
	{ 1e3, 0, 0 },
	{ 0.5, 0, 1 },
	{ 0.5, 0, 1 },
	{ tau_li9, 0, 0 },
	{ tau_he8, 0, 0 },
	{ 1e3, 0, 0 },
	{ 0.5, 0, 1 },
	{ tau_b12, 0, 0 },
};

const bool result_print[npars] = {
	true, true, true,
	true, true, true,
	false, false, true,
	true, false
};

const bool latex_print[npars] = {
	false, false, false,
	true, true, true,
	false, false,
	true, true, false
};

const int npars_single = 8;
//[0]:mu rate
//[1]:n_dc
//[2]:n_li9he8
//[3]:t_li9
//[4]:t_he9
//[5]:ratio_lihe
//[6]:n_b12
//[7]:t_b12
const int n_pdfs = 3;
const char *_pdfs[n_pdfs] = {
	"[1] * [0] * exp(-[0] * x )",
	"[2] * ([5] * ([0] + 1 / [3]) * exp(-([0] + 1 / [3]) * x ) + (1 - [5]) * ([0] + 1 / [4]) * exp(-([0] + 1 / [4]) * x) )",
	"[6] * ([0] + 2 / [7]) * exp(-([0] + 2 / [7]) * x )",
};
/*
const char _f_dc[255] = "[1] * [2] * exp(-[2] * x )";
const char _f_li[255] = "[0] * ([2] + 1 / [4]) * exp(-([2] + 1 / [4]) * x )";
const char _f_bo[255] = "[3] * ([2] + 2 / [5] ) * exp(-([2] + 2 / [5] ) * x )";
const char _f_he[255] = "[6] * ([2] + 1 / [7]) * exp(-([2] + 1 / [7]) * x )";
*/

double r_mu_measured[3][n_range][3] = {0};

double eps_pull[3][2] = {0};

bool enable_eps_pull = false;

double pull(double x, double mean, double err){

	double diff = x - mean;

	return diff * diff / (2 * err * err);

}

double min(double a, double b){
	if(a > b)
		return b;
	else 
		return a;
}

double max(double a, double b){
	if(a > b)
		return a;
	else
		return b;
}

double sigmoid(double x){
	if(x > 0)
		return 1 / ( 1 + exp(-x));
	else{
		double _tmp = exp(x);
		return _tmp / ( 1 + _tmp );
	}
}

double neumaierSum( double *vec, size_t size ){
	double sum = vec[0];
	double c = 0;
	for(size_t i=1;i<size;++i){
		double t = sum + vec[i];
		if( fabs(sum) >= fabs(vec[i]) )
			c += (sum - t) + vec[i];
		else
			c += (vec[i] - t) + sum;
		sum = t;
	}
	return sum + c;
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
					sprintf(buf, "%s %8.2e/%8.2e", buf, fabs(results[i][j][k][0]), results[i][j][k][1]);
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
					printf("& $%10.2f\\pm%10.2f$ ",fabs(results[i][j][k][0]),results[i][j][k][1]);
			
			printf("\\\\ \n");
		}
		printf("\\hline \n");
	}

}

char *fitOptions = "LM";

ROOT::Math::IMultiGenFunction *fPNLL[3];

TF1 *func[3];

TH1 *h[3];

double h_cache[3][max_bins];

void wrap_test(int &npar, double *g, double &result, double *par, int flag){
/*const char *par_names[npars] = {
	"r_mu_tag", "r_mu_atag", "N_DC",
	"N_Li", "EPS_Li", "tau_Li",
	"N_Bo", "EPS_Bo", "tau_Bo",
	"N_He", "EPS_He", "tau_He"
};*/ 
	double r_mu_tag = par[0];
	double r_mu_atag = par[1];
	double n_dc = fabs(par[2]);
	double n_li = fabs(par[3]);
	double eps_li = par[4];
	double tau_li = par[5];
	double n_bo = fabs(par[6]);
	double eps_bo = par[7];
	double tau_bo = par[8];
	double n_he = fabs(par[9]);
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

	func1->SetParameters(par_0);

	func2->SetParameters(par_1);

	func3->SetParameters(par_2);	

	result = 0;
	static double cache[max_bins];
	for(int i=0;i<h1->GetNbinsX();++i){
		cache[i] = 0;
		double x = h1->GetBinCenter(i);
		if(x<fitMin) continue;
		else if(x>fitMax) break;
		double _f0 = func1->Eval(x);
		double _f1 = func2->Eval(x);
		double _f2 = func3->Eval(x);

		//double _cache0 = h_cache[0][i] == 0? _f0 : _f0 - h_cache[0][i] * log(_f0) + TMath::LnGamma(h_cache[0][i]+1);
		//double _cache1 = h_cache[1][i] == 0? _f1 : _f1 - h_cache[1][i] * log(_f1) + TMath::LnGamma(h_cache[1][i]+1);
		//double _cache2 = h_cache[2][i] == 0? _f2 : _f2 - h_cache[2][i] * log(_f2) + TMath::LnGamma(h_cache[2][i]+1);
		double _cache0 = h_cache[0][i] == 0? _f0 : _f0 - h_cache[0][i] * log(_f0);
		double _cache1 = h_cache[1][i] == 0? _f1 : _f1 - h_cache[1][i] * log(_f1);
		double _cache2 = h_cache[2][i] == 0? _f2 : _f2 - h_cache[2][i] * log(_f2);


		cache[i] = _cache0 + _cache1 + _cache2;
	}
	result = neumaierSum(cache,h1->GetNbinsX());
}


void wrap(int &npar, double *g, double &result, double *par, int flag){

/*const char *par_names[npars] = {
	"r_mu_tag", "r_mu_atag", "N_DC",
	"N_Li/He", "EPS_Li/He", "ratio_Li/He", "tau_Li", "tau_He",
	"N_Bo", "EPS_Bo", "tau_Bo",
};*/ 
//[0]:mu rate
//[1]:n_dc
//[2]:n_li9he8
//[3]:t_li9
//[4]:t_he9
//[5]:ratio_lihe
//[6]:n_b12
//[7]:t_b12


	double r_mu_tag = par[0];
	double r_mu_atag = par[1];
	double n_dc = fabs(par[2]);
	double n_lihe = fabs(par[3]);
	double eps_lihe = par[4];
	double ratio_lihe = par[5];
	double tau_li = par[6];
	double tau_he = par[7];
	double n_bo = fabs(par[8]);
	double eps_bo = par[9];
	double tau_bo = par[10];

	double _par[3][npars_single];

	for(int i = 0; i < 3; ++i){
		_par[i][3] = tau_li;
		_par[i][4] = tau_he;
		_par[i][5] = ratio_lihe;
		_par[i][7] = tau_bo;
	}

	_par[0][0] = r_mu_tag + r_mu_atag;
	_par[0][1] = n_dc;
	_par[0][2] = n_lihe;
	_par[0][6] = n_bo;

	_par[1][0] = r_mu_tag;
	_par[1][1] = n_dc + (1-eps_lihe) * n_lihe + (1-eps_bo) * n_bo;
	_par[1][2] = eps_lihe * n_lihe;
	_par[1][6] = eps_bo * n_bo;

	_par[2][0] = r_mu_atag;
	_par[2][1] = n_dc + eps_lihe * n_lihe + eps_bo * n_bo;
	_par[2][2] = (1-eps_lihe) * n_lihe;
	_par[2][6] = (1-eps_bo) * n_bo;

//	double pull_term = pull(par[7], tau_li, tau_li_err) + pull(par[8], tau_b12, tau_b12_err);

	double pull_term = 0;
	if(enable_eps_pull){
		pull_term += pull(eps_lihe, eps_pull[0][0], eps_pull[0][1]);
		if(eps_pull[1][1] > 0)
			pull_term += pull(eps_bo, eps_pull[1][0], eps_pull[1][1]);
	}
//	double pull_term = pull(eps_li, 0.5, 0.5) + pull(eps_bo, 0.5, 0.5) + pull(eps_he, 0.5, 0.5);
	//for(int i=0;i<npars_single;++i)
	//	printf("%2d: %10e %10e %10e\n",i,par_0[i],par_1[i],par_2[i]);
	//cout << (*fPNLL_1)(par_0) << " " << (*fPNLL_2)(par_1) << " " << (*fPNLL_3)(par_2) << endl;
	result = 0;
	for(int i = 0; i < 3; ++i)
		result += (*fPNLL[i])(_par[i]);
	result += pull_term;
}


void fillPars(double *par, TF1 *f[3]){

	double r_mu_tag = par[0];
	double r_mu_atag = par[1];
	double n_dc = fabs(par[2]);
	double n_lihe = fabs(par[3]);
	double eps_lihe = par[4];
	double ratio_lihe = par[5];
	double tau_li = par[6];
	double tau_he = par[7];
	double n_bo = fabs(par[8]);
	double eps_bo = par[9];
	double tau_bo = par[10];

	double _par[3][npars_single];

	for(int i = 0; i < 3; ++i){
		_par[i][3] = tau_li;
		_par[i][4] = tau_he;
		_par[i][5] = ratio_lihe;
		_par[i][7] = tau_bo;
	}

	_par[0][0] = r_mu_tag + r_mu_atag;
	_par[0][1] = n_dc;
	_par[0][2] = n_lihe;
	_par[0][6] = n_bo;

	_par[1][0] = r_mu_tag;
	_par[1][1] = n_dc + (1-eps_lihe) * n_lihe + (1-eps_bo) * n_bo;
	_par[1][2] = eps_lihe * n_lihe;
	_par[1][6] = eps_bo * n_bo;

	_par[2][0] = r_mu_atag;
	_par[2][1] = n_dc + eps_lihe * n_lihe + eps_bo * n_bo;
	_par[2][2] = (1-eps_lihe) * n_lihe;
	_par[2][6] = (1-eps_bo) * n_bo;

	for(int i=0;i<3;++i)
		f[i]->SetParameters(_par[i]);

}

void plotHists(int site, int range, TH1 *h[3], TF1 *f[3]){
	char *dir = "./plots/cfits";
	char buf[255];
	
	for(int i=0;i<3;++i){
		sprintf(buf, "%s/cfit_%d_%d_%d.png", dir, site, range, i);
		h[i]->Draw("E1");
		f[i]->Draw("same");
		gPad->SetLogy();
		gPad->SetLogx();
		c1->SaveAs(buf);
	}
	
}

struct Config{
	double fitMin;
	double fitMax;

	int fix_lifetime;
	int fix_rmu;
	int fix_B12;
	int fix_He8;

	int bound_eps;

	int use_eps_pull;
};

void parseInputs(Config &cfg){

	cin >> cfg.fitMin;
	cout << "fitMin: " << cfg.fitMin << endl;
	fitMin = cfg.fitMin;

	cin >> cfg.fitMax;
	cout << "fitMax: " << cfg.fitMax << endl;
	fitMax = cfg.fitMax;	

	cin >> cfg.fix_B12;
	cout << "fix B12: " << cfg.fix_B12 << endl;

	cin >> cfg.fix_He8;
	cout << "fix He8: " << cfg.fix_He8 << endl; 

	cin >> cfg.bound_eps;
	cout << "Bound tagging efficiencies: " << cfg.bound_eps << endl;

	cin >> cfg.fix_lifetime;
	cout << "Fix isotope lifetimes: " << cfg.fix_lifetime << endl;

	cin >> cfg.fix_rmu;
	cout << "Fix muon rates: " << cfg.fix_rmu << endl;

	cin >> cfg.use_eps_pull;
	cout << "Use efficiency pulls: " << cfg.use_eps_pull << endl;
}

void fitterParInit(int site, int range, TFitter &minuit, Config &cfg){
/*const char *par_names[npars] = {
	"r_mu_tag", "r_mu_atag", "N_DC",
	"N_Li/He", "EPS_Li/He", "ratio_Li/He", "tau_Li", "tau_He",
	"N_Bo", "EPS_Bo", "tau_Bo",
}; */

	double step_ratio = 1e-2;

	if(cfg.fix_rmu == 0){
		minuit.SetParameter(0, par_names[0], r_mu_measured[site][range][1], r_mu_measured[site][range][1] * step_ratio, 0, 0.1 );
		minuit.SetParameter(1, par_names[1], r_mu_measured[site][range][2], r_mu_measured[site][range][2] * step_ratio, 0, 0.1 );
	}else{
		minuit.SetParameter(0, par_names[0], r_mu_measured[site][range][1], 0, 0, 0 );
		minuit.FixParameter(0);
		minuit.SetParameter(1, par_names[1], r_mu_measured[site][range][2], 0, 0, 0 );
		minuit.FixParameter(1);
	}

	double _entries = h[0]->GetEntries();

	minuit.SetParameter(2, par_names[2], _entries, _entries * step_ratio, 0, 0);

	minuit.SetParameter(3, par_names[3], 0.01 * _entries, 0.01 * _entries * step_ratio, 0, 0);

	for(int i=4;i<npars;++i){
		if(cfg.fix_lifetime == 1 && strstr(par_names[i],"tau") != NULL){
			minuit.SetParameter(i,par_names[i],par_init[i][0],0,0,0);
			minuit.FixParameter(i);
		}else if( (cfg.fix_B12 == 1 && strstr(par_names[i],"Bo") != NULL) ){
			minuit.SetParameter(i,par_names[i],0,0,0,0);
			minuit.FixParameter(i);
		}else if( (cfg.fix_He8 == 1 && strstr(par_names[i],"ratio") != NULL) ){
			minuit.SetParameter(i,par_names[i], 1, 0, 0, 0);
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

void getMuRates(){
	ifstream ifile("murates");
	for(int s=0;s<3;++s){
		for(int r=0;r<n_range;++r){
			for(int t=0;t<3;++t){
				ifile >> r_mu_measured[s][r][t];
			}			
		}
	}	

}

void combinedFit(){

	Config cfg;
	getMuRates();
	parseInputs(cfg);

	gStyle->SetOptFit(1);
	char _f_sum[512];
	for(int i=0;i<n_pdfs;++i)
		if(i==0)
			sprintf(_f_sum,"%s",_pdfs[i]);
		else
			sprintf(_f_sum,"%s+%s",_f_sum,_pdfs[i]);
//	sprintf(_f_sum,"%s+%s+%s+%s",_f_dc,_f_li,_f_bo,_f_he);
//
	cout << "INITIALIZING FUNCTION WRAPPERS" << endl;
	char buf[255];	
	ROOT::Math::WrappedMultiTF1 *wf[3];
	for(int i=0;i<3;++i){
		sprintf(buf,"func%d",i);
		func[i] = new TF1(buf, _f_sum, 0, 5000);
		wf[i] = new ROOT::Math::WrappedMultiTF1(*func[i], 1);
	}
	
	double results[3][n_range][npars][2];
	
	ROOT::Fit::DataOptions opt;
	opt.fUseEmpty = true;
	ROOT::Fit::DataRange rangeD;
	rangeD.SetRange(cfg.fitMin, cfg.fitMax);

	cout << "INITIALIZING MINUIT" << endl;
	TFitter minuit(npars);
	double p0 = 1;
	minuit.ExecuteCommand("SET PRINTOUT",&p0,1);
	p0 = 2;
	minuit.ExecuteCommand("SET STRATEGY",&p0,1);
	p0 = 0.5;
	minuit.ExecuteCommand("SET ERRORDEF",&p0,1);
	//p0 = 1e-16;
	//minuit.ExecuteCommand("SET EPSMACHINE",&p0,1);
	double eps_measured[3][2][2];
	double eps_minmax[2][2] = {
		{1, 0},
		{1, 0},
	};
	enable_eps_pull = false;
	for(int site=0;site<3;++site){
		do_fit(site, n_range - 1, minuit, cfg, wf, opt, rangeD, results);
		eps_measured[site][0][0] = minuit.GetParameter(4);
		eps_measured[site][0][1] = minuit.GetParError(4);
		eps_measured[site][1][0] = minuit.GetParameter(9);
		eps_measured[site][1][1] = minuit.GetParError(9);
		for(int iso=0;iso<2;++iso){
			eps_minmax[iso][0] = min(eps_minmax[iso][0], eps_measured[site][iso][0]);
			eps_minmax[iso][1] = max(eps_minmax[iso][1], eps_measured[site][iso][0]);
			eps_pull[iso][0] += eps_measured[site][iso][0] / 3;
			eps_pull[iso][1] += eps_measured[site][iso][1] * eps_measured[site][iso][1] / 9;
		}
	}
	for(int iso=0;iso<2;++iso){
		double __tmp = eps_minmax[iso][1] - eps_minmax[iso][0];
		eps_pull[iso][1] += __tmp * __tmp;
		eps_pull[iso][1] = sqrt(eps_pull[iso][1]);
	}

	enable_eps_pull = (cfg.use_eps_pull == 1);

	for(int site=0;site<3;++site){
		for(int range=0;range<n_range-1;++range){
			do_fit(site, range, minuit, cfg, wf, opt, rangeD, results);
		}
	}
	printLatexTable(results);
	for(int iso=0;iso<2;++iso)
		printf("iso%d pull: %10.3f %10.3f\n", iso, eps_pull[iso][0], eps_pull[iso][1]);
	printResults(results);
}

void do_fit(int site, int range, TFitter &minuit, Config &cfg, ROOT::Math::WrappedMultiTF1 *wf[3], ROOT::Fit::DataOptions &opt, ROOT::Fit::DataRange &rangeD, double results[3][n_range][npars][2]){

	char buf[255];
	sprintf(buf,"./hists/EH%d_dtlSH_%d.root",site+1,range);
	TFile *f1 = new TFile(buf,"READ");
	h[0] = (TH1*)f1->Get("h");
	ROOT::Fit::BinData data1(opt,rangeD);
	ROOT::Fit::FillData(data1, h[0]);
	ROOT::Fit::PoissonLLFunction PNLL1(data1,*wf[0]);
	fPNLL[0] = &PNLL1;

	sprintf(buf,"./hists/EH%d_dtlSH_%d_tag.root",site+1,range);
	TFile *f2 = new TFile(buf,"READ");
	h[1] = (TH1*)f2->Get("h");
	ROOT::Fit::BinData data2(opt,rangeD);
	ROOT::Fit::FillData(data2, h[1]);
	ROOT::Fit::PoissonLLFunction PNLL2(data2,*wf[1]);
	fPNLL[1] = &PNLL2;

	sprintf(buf,"./hists/EH%d_dtlSH_%d_atag.root",site+1,range);
	TFile *f3 = new TFile(buf,"READ");
	h[2] = (TH1*)f3->Get("h");
	ROOT::Fit::BinData data3(opt,rangeD);
	ROOT::Fit::FillData(data3, h[2]);
	ROOT::Fit::PoissonLLFunction PNLL3(data3,*wf[2]);
	fPNLL[2] = &PNLL3;

	for(size_t _bin = 0;_bin<h[0]->GetNbinsX();++_bin){
		for(int _i = 0; _i < 3; ++_i)
			h_cache[_i][_bin] = h[_i]->GetBinContent(_bin);
	}
	minuit.SetFCN(wrap);
	fitterParInit(site, range, minuit, cfg);

	minuit.ExecuteCommand("SIMPLEX", simplex_args, 2);
	minuit.ExecuteCommand("MINIMIZE", minimizer_args, 2);
	minuit.ExecuteCommand("MINOS", minos_args, 1);
	//fix_and_release(minuit, cfg);
	//
	double _pars[npars];
	for(int _i=0;_i<npars;++_i){
		_pars[_i] = minuit.GetParameter(_i);
	}
	fillPars(_pars, func);

	plotHists(site, range, h, func);
	
	for(int _i=0;_i<npars;++_i){
		results[site][range][_i][0] = minuit.GetParameter(_i);
		results[site][range][_i][1] = minuit.GetParError(_i);
	}
}

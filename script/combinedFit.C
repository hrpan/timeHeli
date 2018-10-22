#include "util.C"
#include "range.h"

const char *hists_prefix = "./hists";

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

const int npars = 6 * n_range + 2;
/*
 * r_mu_tag, r_mu_atag : 2 * n_range
 * eps_bo, eps_lihe : 2 * n_range
 * n_lihe, n_bo : 2 * n_range
 * ratio_lihe : 1 
 * n_dc : 1
 */
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

double eps_pull[2][2] = {0};

bool enable_eps_pull = false;

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

ROOT::Math::IMultiGenFunction *fPNLL[n_range][3];

TF1 *func[n_range][3];

TH1 *h[n_range][3];

void parTrans(double *par, double _par[n_range][3][npars_single]){
	double r_mu_tag[n_range];
	double r_mu_atag[n_range];
	double n_lihe[n_range];
	double eps_lihe[n_range];
	double n_bo[n_range];
	double eps_bo[n_range];
	double n_dc;
	double ratio_lihe;

	int p_idx = 0;
	for(int i = 0; i < n_range; ++i){
		r_mu_tag[i] = par[p_idx++];
		r_mu_atag[i] = par[p_idx++];
		n_lihe[i] = fabs(par[p_idx++]);
		eps_lihe[i] = par[p_idx++];
		n_bo[i] = fabs(par[p_idx++]);
		eps_bo[i] = par[p_idx++];
	}
	n_dc = fabs(par[p_idx++]);
	r_lihe = par[p_idx++];

	for(int i=0;i<n_range;++i){
		double n_dc_fake = 0;
		for(int j=0;j<n_range;++j){
			if(i==j) continue;
			n_dc_fake += n_lihe[j] + n_bo[j];
		}
		for(int j=0;j<3;++j){
			if(j==0){
				_par[i][j][0] = r_mu_tag + r_mu_atag;
				_par[i][j][1] = n_dc + n_dc_fake;
				_par[i][j][2] = n_lihe[i];
				_par[i][j][6] = n_bo[i];
			}else if(j==1){
				_par[i][j][0] = r_mu_tag;
				_par[i][j][1] = n_dc + n_dc_fake + ( 1 - eps_lihe[i] ) * n_lihe[i] + ( 1 - eps_bo[i] ) * n_bo[i];
				_par[i][j][2] = eps_lihe[i] * n_lihe[i];
				_par[i][j][6] = eps_bo[i] * n_bo[i];
			}else{
				_par[i][j][0] = r_mu_atag;
				_par[i][j][1] = n_dc + n_dc_fake + eps_lihe[i] * n_lihe[i] + eps_bo[i] * n_bo[i];
				_par[i][j][2] = ( 1 - eps_lihe[i] ) * n_lihe[i];
				_par[i][j][6] = ( 1 - eps_bo[i] ) * n_bo[i];
				
			}
			_par[i][j][3] = tau_li9;
			_par[i][j][4] = tau_he8;
			_par[i][j][5] = r_lihe;
			_par[i][j][7] = tau_b12;
		}
	}
}

void wrap(int &npar, double *g, double &result, double *par, int flag){

	double _par[n_range][3][npars_single];

	parTrans(par, _par);

	double pull_term = 0;
	if(enable_eps_pull){
		//cout << par[4] << " " << eps_pull[0][0] << " " << eps_pull[0][1] << " " << pull(par[4], eps_pull[0][0], eps_pull[0][1]) << endl; 
		pull_term += pull(par[4], eps_pull[0][0], eps_pull[0][1]);
		if(eps_pull[1][1] > 0)
			pull_term += pull(par[9], eps_pull[1][0], eps_pull[1][1]);
	}
	//for(int i=0;i<npars_single;++i)
	//	printf("%2d: %10e %10e %10e\n",i,par_0[i],par_1[i],par_2[i]);
	//cout << (*fPNLL_1)(par_0) << " " << (*fPNLL_2)(par_1) << " " << (*fPNLL_3)(par_2) << endl;
	result = 0;
	for(int i = 0; i < 3; ++i)
		result += (*fPNLL[i])(_par[i]);
	result += pull_term;
}


void fillPars(double *par, TF1 *f[n_range][3]){

	double _par[n_range][3][npars_single];

	parTrans(par, _par);

	for(int i=0;i<n_range;++i)
		for(int j=0;j<3;++j)
			f[i][j]->SetParameters(_par[i][j]);

}

void plotHists(int site, TH1 *h[n_range][3], TF1 *f[n_range][3]){
	char *dir = "./plots/cfits";
	char buf[255];
	
	for(int i=0;i<3;++i){
		sprintf(buf, "%s/cfit_%d_%d_%d.png", dir, site, range, i);
		h[i]->Draw("E1");
		f[i]->Draw("same");
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
		double rmu_init = 1e-4;
		minuit.SetParameter(0, par_names[0], rmu_init, rmu_init * step_ratio, 0, 0.1 );
		minuit.SetParameter(1, par_names[1], rmu_init, rmu_init * step_ratio, 0, 0.1 );
	}else{
		minuit.SetParameter(0, par_names[0], r_mu_measured[site][range][1], 0, 0, 0 );
		minuit.FixParameter(0);
		minuit.SetParameter(1, par_names[1], r_mu_measured[site][range][2], 0, 0, 0 );
		minuit.FixParameter(1);
	}

	double _entries = h[0]->GetEntries();
	
	
	minuit.SetParameter(2, par_names[2], _entries, _entries * step_ratio, 0, 0);
	double bkg_ratio = 1e-2;
	minuit.SetParameter(3, par_names[3], bkg_ratio * _entries, bkg_ratio * _entries * step_ratio, 0, 0);

	for(int i=4;i<npars;++i){
		if(cfg.fix_lifetime == 1 && strstr(par_names[i], "tau") != NULL){
			minuit.SetParameter(i, par_names[i], par_init[i][0], 0, 0, 0);
			minuit.FixParameter(i);
		}else if( (cfg.fix_B12 == 1 && strstr(par_names[i], "Bo") != NULL) ){
			minuit.SetParameter(i, par_names[i], 0, 0, 0, 0);
			minuit.FixParameter(i);
		}else if( (cfg.fix_He8 == 1 && strstr(par_names[i], "ratio") != NULL) ){
			minuit.SetParameter(i, par_names[i], 1, 0, 0, 0);
			minuit.FixParameter(i);
		}else if( (cfg.bound_eps == 0 && strstr(par_names[i],"EPS") != NULL) ){
			minuit.SetParameter(i, par_names[i], par_init[i][0], par_init[i][0] * step_ratio, 0, 0);
		}else if( (strstr(par_names[i], "ratio") != NULL) ){
			minuit.SetParameter(i, par_names[i], 0.95, 0, 0, 0);
		}else{
			minuit.SetParameter(i, par_names[i], par_init[i][0], par_init[i][0] * step_ratio, par_init[i][1], par_init[i][2]);
		}

	}
}

void combinedFit(){

	Config cfg;
	parseInputs(cfg);

	gStyle->SetOptFit(1);
	char _f_sum[512];
	for(int i=0;i<n_pdfs;++i)
		if(i==0)
			sprintf(_f_sum,"%s",_pdfs[i]);
		else
			sprintf(_f_sum,"%s+%s",_f_sum,_pdfs[i]);

	cout << "INITIALIZING FUNCTION WRAPPERS" << endl;
	char buf[255];	
	ROOT::Math::WrappedMultiTF1 *wf[n_range][3];
	for(int r=0;r<n_range;++r){
		for(int i=0;i<3;++i){
			sprintf(buf,"func%d_%d",r,i);
			func[r][i] = new TF1(buf, _f_sum, 0, 5000);
			wf[r][i] = new ROOT::Math::WrappedMultiTF1(*func[r][i], 1);
		}
	}
	double results[3][npars][2];
	
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
			//eps_pull[iso][0] += eps_measured[site][iso][0] / 3;
			//eps_pull[iso][1] += eps_measured[site][iso][1] * eps_measured[site][iso][1] / 9;
		}
	}
	/*
	for(int iso=0;iso<2;++iso){
		double __tmp = eps_minmax[iso][1] - eps_minmax[iso][0];
		eps_pull[iso][1] += __tmp * __tmp;
		eps_pull[iso][1] = sqrt(eps_pull[iso][1]);
	}
	*/
	enable_eps_pull = (cfg.use_eps_pull == 1);

	for(int site=0;site<3;++site){
		if(enable_eps_pull){
			for(int iso=0;iso<2;++iso){
				eps_pull[iso][0] = eps_measured[site][iso][0];
				double __tmp = eps_minmax[iso][1] - eps_minmax[iso][0];
				eps_pull[iso][1] = sqrt(eps_measured[site][iso][1] * eps_measured[site][iso][1] + __tmp * __tmp);
			}
		}
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

	char *_suffix[3] = { "", "_tag", "_atag" };

	ROOT::Fit::BinData *data[n_range][3];
	for(int r=0;r<n_range;++r){
		for(int s=0;s<3;++s){
			sprintf(buf,"%s/EH%d_dtlSH_%d%s.root", hists_prefix, site+1, r, _suffix[s]);
			TFile *f = new TFile(buf, "READ");
			h[r][s] = (TH1*) f->Get("h");
			data[r][s] = new ROOT::Fit::BinData(opt, rangeD);
			ROOT::Fit::FillData(*data[r][s], h[r][s]);
			fPNLL[r][s] = new ROOT::Fit::PoissonLikelihoodFCN<ROOT::Math::IBaseFunctionMultiDim>(*data[r][s], *wf[r][s]);	
		}
	}

	minuit.SetFCN(wrap);

	fitterParInit(site, range, minuit, cfg);

	fit_procedure(minuit);

	double _pars[npars];
	for(int _i=0;_i<npars;++_i){
		_pars[_i] = minuit.GetParameter(_i);
		results[site][range][_i][0] = minuit.GetParameter(_i);
		results[site][range][_i][1] = minuit.GetParError(_i);
	}

	fillPars(_pars, func);

	plotHists(site, range, h, func);
	
	
}

void fit_procedure(TFitter &minuit){
	/*const char *par_names[npars] = {
		"r_mu_tag", "r_mu_atag", "N_DC",
		"N_Li/He", "EPS_Li/He", "ratio_Li/He", "tau_Li", "tau_He",
		"N_Bo", "EPS_Bo", "tau_Bo",
	}; */

	minuit.FixParameter(4);
	minuit.FixParameter(5);
	minuit.ExecuteCommand("MINIMIZE", minimizer_args, 2);

	minuit.ReleaseParameter(4);
	minuit.ExecuteCommand("MINIMIZE", minimizer_args, 2);

	//minuit.ReleaseParameter(5);
	//minuit.ExecuteCommand("MINIMIZE", minimizer_args, 2);

	minuit.ExecuteCommand("MINOS", minos_args, 1);
}

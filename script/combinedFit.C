#include "util.C"
#include "range.h"
#include "slice.h"

const char *hists_prefix = "./hists";

const double tau_li9 = 256.366;
const double tau_li9_err = 0.866;

const double tau_he8 = 171.17;
const double tau_he8_err = 2.31;

const double tau_b12 = 29.142;
const double tau_b12_err = 0.0288;

const double tau_n12 = 15.9;

double fitMin = 1.5;
double fitMax = 5000;

const int npars_max = 500;
const bool use_tagging = false;
const bool fix_tau_short = true;

const int slice_max = 20;

char par_names[npars_max][255]; 

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
	"[6] * ([0] + 1 / [7]) * exp(-([0] + 1 / [7]) * x )",
	"[2] * ([5] * ([0] + 1 / [3]) * exp(-([0] + 1 / [3]) * x ) + (1 - [5]) * ([0] + 1 / [4]) * exp(-([0] + 1 / [4]) * x) )",
	"[1] * [0] * exp(-[0] * x )",
};

const double livetime[3][4] = {
	{1544.17, 1737.38, 0, 0},
	{1741.12, 1554.04, 0, 0},
	{1739.26, 1739.26, 1739.26, 1551.78}
};

const double eff_muon[3][4] = {
	{0.8227, 0.8194, 0, 0},
	{0.8518, 0.8506, 0, 0},
	{0.9832, 0.9829, 0.9828, 0.9833}
};

const double eff_mult[3][4] = {
	{0.9495, 0.9502, 0, 0},
	{0.9520, 0.9518, 0, 0},
	{0.9525, 0.9525, 0.9523, 0.9524}
};

double lt;

void printResults(double results[3][npars_max][2]){
	char buf[512];
	printf("%15s %21s %21s %21s\n", "site", "EH1", "EH2", "EH3");
	for(int i=0;i<npars_max;++i){
		printf("%15s", par_names[i]);
		for(int j=0;j<3;++j)
			printf(" %10.3e/%10.3e", results[j][i][0], results[j][i][1]);
		cout << endl;
		if(par_names[i+1][0] == 0) break;
	}

}

ROOT::Math::IMultiGenFunction *fPNLL[n_range][3][slice_types][slice_max];
ROOT::Math::WrappedMultiTF1 *wf[n_range][3][slice_types][slice_max];

TF1 *func[n_range][3][slice_types][slice_max];
TH1 *h[n_range][3][slice_types][slice_max];

double _par[n_range][3][slice_types][slice_max][npars_single];

void parTrans(const double *par){

	double scale = h[0][0][0][0]->GetBinWidth(1);

	double r_mu[n_range];
	double r_mu_tag[n_range];
	double r_mu_atag[n_range];
	double n_lihe[n_range];
	double n_lihe_tag[n_range];
	double n_lihe_atag[n_range];
	double n_bo[n_range];
	double n_bo_tag[n_range];
	double n_bo_atag[n_range];
	double eps_s_dc[slice_types][slice_max];
	double eps_s_lihe[slice_types][slice_max];
	double eps_s_bo[slice_types][slice_max];
	//double eps_bo;
	double n_dc;
	double r_lihe;
	double tau_short;

	int p_idx = 0;
	for(int r = 0; r < n_range; ++r){
		if(use_tagging){
			r_mu_tag[r] = par[p_idx++];
			r_mu_atag[r] = par[p_idx++];
			n_lihe_tag[r] = fabs(par[p_idx++]) * scale * lt;
			n_lihe_atag[r] = fabs(par[p_idx++]) * scale * lt;
			n_bo_tag[r] = fabs(par[p_idx++]) * scale * lt;
			n_bo_atag[r] = fabs(par[p_idx++]) * scale * lt;
		}else{
			r_mu[r] = par[p_idx++];
			n_lihe[r] = fabs(par[p_idx++]) * scale * lt;
			n_bo[r] = fabs(par[p_idx++]) * scale * lt;
		}
	}

	for(int _type=0;_type < slice_types; ++_type){
		eps_s_dc[_type][0] = 1;
		eps_s_lihe[_type][0] = 1;
		eps_s_bo[_type][0] = 1;

		double _tmp_sum[3] = {1, 1, 1};
		int _slices = slices[_type];
		for(int _s=0;_s<_slices-1;++_s){
			eps_s_dc[_type][_s+1] = fabs(par[p_idx++]);
			_tmp_sum[0] += eps_s_dc[_type][_s+1];
			eps_s_lihe[_type][_s+1] = fabs(par[p_idx++]);
			_tmp_sum[1] += eps_s_lihe[_type][_s+1];
			eps_s_bo[_type][_s+1] = fabs(par[p_idx++]);
			_tmp_sum[2] += eps_s_bo[_type][_s+1];
		}
		for(int _s=0;_s<_slices;++_s){
			eps_s_dc[_type][_s] /= _tmp_sum[0];
			eps_s_lihe[_type][_s] /= _tmp_sum[1];
			eps_s_bo[_type][_s] /= _tmp_sum[2];
		}
	}
	r_lihe = par[p_idx++];
	n_dc = par[p_idx++] * scale * lt;
	tau_short = par[p_idx++];

	for(int r=0; r<n_range; ++r){
		double n_dc_fake_lihe = 0;
		double n_dc_fake_bo = 0;
		for(int j=0; j<n_range; ++j){
			if(r==j) continue;
			if(use_tagging){
				n_dc_fake_lihe += (n_lihe_tag[j] + n_lihe_atag[j]);
				n_dc_fake_bo += (n_bo_tag[j] + n_bo_atag[j]);
			}else{
				n_dc_fake_lihe += n_lihe[j];
				n_dc_fake_bo += n_bo[j];
			}
		}
		for(int t=0; t<3; ++t){
			for(int _type=0; _type<slice_types; ++_type){
				double *_eps_s_dc = eps_s_dc[_type];
				double *_eps_s_lihe = eps_s_lihe[_type];
				double *_eps_s_bo = eps_s_bo[_type];
				int _slices = slices[_type];
				for(int s=0; s<_slices; ++s){
					double _n_dc = n_dc * _eps_s_dc[s];
					double _n_dc_fake = n_dc_fake_lihe * _eps_s_lihe[s] + n_dc_fake_bo * _eps_s_bo[s];
					double _n_lihe = n_lihe[r] * _eps_s_lihe[s];
					double _n_bo = n_bo[r] * _eps_s_bo[s];
					if(use_tagging){
						_n_lihe = (n_lihe_tag[r] + n_lihe_atag[r]) * _eps_s_lihe[s];
						_n_bo = (n_bo_tag[r] + n_bo_atag[r]) * _eps_s_bo[s];
					}
					if(t==0){
						if(use_tagging)
							_par[r][t][_type][s][0] = r_mu_tag[r] + r_mu_atag[r];
						else
							_par[r][t][_type][s][0] = r_mu[r];
						_par[r][t][_type][s][1] = _n_dc + _n_dc_fake;
						_par[r][t][_type][s][2] = _n_lihe;			
						_par[r][t][_type][s][6] = _n_bo;
					}else if(t==1){
						_par[r][t][_type][s][0] = r_mu_tag[r];
						_par[r][t][_type][s][1] = _n_dc + _n_dc_fake 
							+ _eps_s_lihe[s] * n_lihe_atag[r] 
							+ _eps_s_bo[s] * n_bo_atag[r]; 
						_par[r][t][_type][s][2] = _eps_s_lihe[s] * n_lihe_tag[r];
						_par[r][t][_type][s][6] = _eps_s_bo[s] * n_bo_tag[r];
					}else{
						_par[r][t][_type][s][0] = r_mu_atag[r];
						_par[r][t][_type][s][1] = _n_dc + _n_dc_fake
							+ _eps_s_lihe[s] * n_lihe_tag[r] 
							+ _eps_s_bo[s] * n_bo_tag[r];
						_par[r][t][_type][s][2] = _eps_s_lihe[s] * n_lihe_atag[r];
						_par[r][t][_type][s][6] = _eps_s_bo[s] * n_bo_atag[r];
					}
					_par[r][t][_type][s][3] = tau_li9;
					_par[r][t][_type][s][4] = tau_he8;
					_par[r][t][_type][s][5] = r_lihe;
					_par[r][t][_type][s][7] = tau_short;
									
	/*	
					cout << i << " " << j << " " << k << " ";
					for(int l=0;l<npars_single;++l)
						printf("%6e ", _par[i][j][k][l]);
					cout << endl;
	*/

				} //s
			} //_type
			if(!use_tagging) break;
		} //t
		
	}

}

double likelihood(const double *par){

	parTrans(par);
	
	double result = 0;
	int idx = 0;
	for(int r = 0; r < n_range; ++r){
		for(int t = 0; t < 3; ++t){
			if(use_tagging && t==0) continue;
			for(int type = 0; type < slice_types; ++type){
				for(int s = 0; s < slices[type]; ++s){
					double _tmp = (*fPNLL[r][t][type][s])(_par[r][t][type][s]);
					result += _tmp;
				}//slice
			}//slice type
			if(!use_tagging) break;
		}//tag
	}//range
	return result;
}

void fillPars(double *par, TF1 *f[n_range][3][slice_types][slice_max]){

	parTrans(par);

	for(int r=0;r<n_range;++r)
		for(int t=0;t<3;++t){
			for(int type=0;type<slice_types;++type)
				for(int s=0;s<slices[type];++s)
					f[r][t][type][s]->SetParameters(_par[r][t][type][s]);
			if(!use_tagging) break;
		}//tag

}

void plotHists(int site, double *par){
	fillPars(par, func);
	char *dir = "./plots/cfits";
	char buf[255];
	parTrans(par);
	for(int r=0;r<n_range;++r)
		for(int t=0;t<3;++t){
			for(int type=0;type<slice_types;++type){
				for(int s=0;s<slices[type];++s){
					TF1 ibd("ibd","[1] * [0] * exp(-[0] * x)",0, 5000);
					TF1 heli("heli","[1] * [0] * exp(-[0] * x) + [2] * ([5] * ([0] + 1 / [3]) * exp(-([0] + 1 / [3]) * x ) + (1 - [5]) * ([0] + 1 / [4]) * exp(-([0] + 1 / [4]) * x) )",0, 5000);
					ibd.SetParameters(_par[r][t][type][s]);
					ibd.SetLineColor(kGreen);
					heli.SetParameters(_par[r][t][type][s]);
					heli.SetLineColor(kBlue);
					sprintf(buf, "%s/cfit_%d_%d_%d_%d_%d.png", dir, site, r, t, type, s);
					h[r][t][type][s]->Draw("E1");
					func[r][t][type][s]->Draw("same");
					ibd.Draw("same");
					heli.Draw("same");
					gPad->SetLogx();
					c1->SaveAs(buf);
				}
			}
			if(!use_tagging) break;
		}
	
}

void initialize_minimizer(int site, ROOT::Math::Minimizer *minim, bool verbose){
	cout << "INITIALIZING MINIMIZER" << endl;
	minim->SetErrorDef(0.5);
	minim->SetTolerance(0.01);
	minim->SetPrintLevel(1);
	//minim->SetPrecision(1e-16);
	minim->SetStrategy(2);
	minim->SetMaxFunctionCalls(1000000);
	minim->SetMaxIterations(1000000);

	double rmu_init;
	double n_lihe_init = 1;
	double n_bo_init = 1;
	double eps_init = 0.5;
	double r_init = 0.5;

	double step_ratio = 1e-1;
	char buf[255];
	int p_idx = 0;
	char *init_str = "Initializing parameter %2d %20s %10.2e %10.2e\n";
	for(int r=0;r<n_range;++r){
		if(use_tagging){		
			h[r][1][0][0]->Fit("expo","LQN0","", 300, 5000);
			rmu_init = fabs(expo->GetParameter(1));
			sprintf(buf, "rmu_tag_%d",r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, rmu_init, rmu_init * step_ratio);
			if(verbose)
				printf(init_str, p_idx-1, par_names[p_idx-1], rmu_init, rmu_init * step_ratio);

			h[r][2][0][0]->Fit("expo","LQN0","", 300, 5000);
			rmu_init = fabs(expo->GetParameter(1));
			sprintf(buf, "rmu_atag_%d", r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, rmu_init, rmu_init * step_ratio);
			if(verbose)
				printf(init_str, p_idx-1, par_names[p_idx-1], rmu_init, rmu_init * step_ratio);

			sprintf(buf, "n_lihe_tag_%d", r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, n_lihe_init, n_lihe_init * step_ratio);
			if(verbose)
				printf(init_str, p_idx-1, par_names[p_idx-1], n_lihe_init, n_lihe_init * step_ratio);

			sprintf(buf, "n_lihe_atag_%d", r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, n_lihe_init, n_lihe_init * step_ratio);
			if(verbose)
				printf(init_str, p_idx-1, par_names[p_idx-1], n_lihe_init, n_lihe_init * step_ratio);

			sprintf(buf, "n_bo_tag_%d", r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, n_bo_init, n_bo_init * step_ratio);
			if(verbose)
				printf(init_str, p_idx-1, par_names[p_idx-1], n_bo_init, n_bo_init * step_ratio);

			sprintf(buf, "n_bo_atag_%d", r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, n_bo_init, n_bo_init * step_ratio);
			if(verbose)
				printf(init_str, p_idx-1, par_names[p_idx-1], n_bo_init, n_bo_init * step_ratio);

		}else{
			h[r][0][0][0]->Fit("expo","LQN0","", 300, 5000);
			rmu_init = fabs(expo->GetParameter(1));
			sprintf(buf, "rmu_%d",r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, rmu_init, rmu_init * step_ratio);
			if(verbose)
				printf(init_str, p_idx-1, par_names[p_idx-1], rmu_init, rmu_init * step_ratio);

			sprintf(buf, "n_lihe_%d",r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, n_lihe_init, n_lihe_init * step_ratio);
			if(verbose)
				printf(init_str, p_idx-1, par_names[p_idx-1], n_lihe_init, n_lihe_init * step_ratio);

			sprintf(buf, "n_bo_%d",r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, n_bo_init, n_bo_init * step_ratio);
			if(verbose)
				printf(init_str, p_idx-1, par_names[p_idx-1], n_bo_init, n_bo_init * step_ratio);

		}
	}

	for(int type=0;type<slice_types;++type){
		for(int s=0;s<slices[type]-1;++s){
			eps_init = h[0][0][type][s+1]->GetEntries() / h[0][0][type][0]->GetEntries();
			sprintf(buf, "eps_s_dc_%d_%d", type, s);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, eps_init, eps_init * step_ratio);
			if(verbose)
				printf(init_str, p_idx-1, par_names[p_idx-1], eps_init, eps_init * step_ratio);

			sprintf(buf, "eps_s_lihe_%d_%d", type, s);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, eps_init, eps_init * step_ratio);
			if(verbose)
				printf(init_str, p_idx-1, par_names[p_idx-1], eps_init, eps_init * step_ratio);

			eps_init = 1;
			sprintf(buf, "eps_s_bo_%d_%d", type, s);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, eps_init, eps_init * step_ratio);
			if(verbose)
				printf(init_str, p_idx-1, par_names[p_idx-1], eps_init, eps_init * step_ratio);
		}
	}
	
	sprintf(par_names[p_idx], "r_lihe");
	minim->SetLimitedVariable(p_idx++, "r_lihe", r_init, r_init * step_ratio, 0, 1);
	//minim->SetFixedVariable(p_idx++, "r_lihe", 0.95);
	if(verbose)
		printf(init_str, p_idx-1, par_names[p_idx-1], r_init, r_init * step_ratio);

	double n_dc_init = 0;
	for(int s=0;s<slices[0];++s)
		n_dc_init += h[0][0][0][s]->GetEntries();
	n_dc_init /= lt;	

	sprintf(par_names[p_idx], "n_dc");
	minim->SetVariable(p_idx++, "n_dc", n_dc_init, n_dc_init * step_ratio);
	if(verbose)
		printf(init_str, p_idx-1, par_names[p_idx-1], n_dc_init, n_dc_init * step_ratio);

	double tau_init = 10;
	sprintf(par_names[p_idx], "tau_short");
	if(fix_tau_short)
		minim->SetFixedVariable(p_idx++, "tau_short", tau_b12 / 2.0);
	else
		minim->SetLimitedVariable(p_idx++, "tau_short", tau_init, tau_init * step_ratio, 1, 50);
	if(verbose)
		printf(init_str, p_idx-1, par_names[p_idx-1], tau_init, tau_init * step_ratio);

	par_names[p_idx][0] = 0;
	cout << endl;
}

void combinedFit(){

	gStyle->SetOptFit(1);

	char _f_sum[512];
	for(int i=0;i<n_pdfs;++i)
		if(i==0)
			sprintf(_f_sum,"%s",_pdfs[i]);
		else
			sprintf(_f_sum,"%s+%s",_f_sum,_pdfs[i]);

	cout << "INITIALIZING FUNCTION WRAPPERS" << endl;
	char buf[255];	
	for(int r=0;r<n_range;++r){
		for(int t=0;t<3;++t){
			for(int type=0;type<slice_types;++type)
				for(int s=0;s<slices[type];++s){
					sprintf(buf,"func%d_%d_%d_%d",r,t,type,s);
					func[r][t][type][s] = new TF1(buf, _f_sum, 0, 5000);
					wf[r][t][type][s] = new ROOT::Math::WrappedMultiTF1(*func[r][t][type][s], 1);
				}
			if(!use_tagging) break;
		}
	}
	double results[3][npars_max][2];
	double daily_rates[3][2];
	ROOT::Fit::DataOptions opt;
	opt.fUseEmpty = true;
	for(int site=0;site<3;++site){
		lt = 0;
		for(int ad=0;ad<4;++ad)
			lt += livetime[site][ad] * eff_muon[site][ad] * eff_mult[site][ad];

		do_fit(site, opt, results, daily_rates);
	}
	printResults(results);
	for(int s=0;s<3;++s){
		printf("EH%d: %10.3e %10.3e\n", s+1, daily_rates[s][0], daily_rates[s][1]);
	}	
	cout << endl;
}

void do_fit(int site, ROOT::Fit::DataOptions &opt, double results[3][npars_max][2], double daily_rates[3][2]){

	char buf[255];

	char *_suffix[3] = { "", "_tag", "_atag" };

	ROOT::Fit::BinData *data[n_range][3][slice_types][slice_max];
	ROOT::Fit::DataRange rangeD(fitMin, fitMax);
	for(int r=0;r<n_range;++r){
		for(int t=0;t<3;++t){
			for(int type=0;type<slice_types;++type){
				for(int s=0;s<slices[type];++s){
					sprintf(buf,"%s/EH%d_dtlSH_%d%s_%d_%d.root", hists_prefix, site+1, r, _suffix[t], type, s);
					cout << "READING FILE:" << buf << endl;
					TFile *f = new TFile(buf, "READ");
					h[r][t][type][s] = (TH1*) f->Get("h");
					data[r][t][type][s] = new ROOT::Fit::BinData(opt, rangeD);
					ROOT::Fit::FillData(*data[r][t][type][s], h[r][t][type][s]);
					fPNLL[r][t][type][s] = new ROOT::Fit::PoissonLikelihoodFCN<ROOT::Math::IBaseFunctionMultiDim>(*data[r][t][type][s], *wf[r][t][type][s]);	
				}
			}
			if(!use_tagging) break;
		}
	}

	int npars;
	int slice_pars = 0;
	for(int type=0;type<slice_types;++type)
		slice_pars += slices[type] - 1;
	slice_pars *= 3;
	if(use_tagging)
		npars = 6 * n_range + slice_pars + 3;
	else
		npars = 3 * n_range + slice_pars + 3;

	ROOT::Math::Functor f_tmp(likelihood, npars);
	ROOT::Math::Minimizer *minim = ROOT::Math::Factory::CreateMinimizer("Minuit2", "migrad");
	ROOT::Math::Minimizer *minim_simplex = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
	minim->SetFunction(f_tmp);
	minim_simplex->SetFunction(f_tmp);
	initialize_minimizer(site, minim, true);
	initialize_minimizer(site, minim_simplex, false);
	
	//minim_simplex->Minimize();
	//minim->SetVariableValues(minim_simplex->X());
	for(int i=0;i<5;++i){
		minim->Hesse();
		minim->Minimize();
	}
	minim->Hesse();
	cout << "COV STATUS: " << minim->CovMatrixStatus() << endl;
	minim->PrintResults();

	double *x = minim->X();
	double *x_err = minim->Errors();
	for(int i=0;i<npars;++i){
		results[site][i][0]	= x[i];
		results[site][i][1] = x_err[i];
	}
	
	int offset = use_tagging? 2:1;
	int step = use_tagging? 6:3;
	daily_rates[site][0] = 0;
	daily_rates[site][1] = 0;
	for(int i=0;i<n_range;++i){
		double _n = fabs(x[offset + i * step]);
		if(use_tagging)
			_n += fabs(x[offset + i * step + 1]);
		if(i==n_range-1)
			daily_rates[site][0] += 0.211 * _n;
		else
			daily_rates[site][0] += _n;
		if(use_tagging){
			for(int _i=0;_i<2;++_i)
				for(int j=0;j<n_range;++j){
					for(int _j=0;_j<2;++_j){
						double _cov = minim->CovMatrix(offset + step * i + _i, offset + step * j + _j);
						if(i==n_range-1)
							_cov *= 0.211;
						if(j==n_range-1)
							_cov *= 0.211;
						daily_rates[site][1] += _cov;
					}
				}
		}else{
			for(int j=0;j<n_range;++j){
				double _cov = minim->CovMatrix(offset + step * i, offset + step * j);
				if(i==n_range-1)
					_cov *= 0.211;
				if(j==n_range-1)
					_cov *= 0.211;
				daily_rates[site][1] += _cov;
			}
		}
	}
	daily_rates[site][1] = sqrt(daily_rates[site][1]);
	plotHists(site, x);
	plotSlice(site, minim);

}

void plotSlice(int site, ROOT::Math::Minimizer *minim){
	double *par = minim->X();
	double *par_err = minim->Errors();

	const int ncomps = 3;
	
	int offset=0;
	for(;offset<npars_max;++offset)
		if(strstr(par_names[offset],"eps_s")!=NULL)	break;

	for(int type=0;type<slice_types;++type){
		double x[ncomps][slice_max], x_err[ncomps][slice_max];
		double y[ncomps][slice_max], y_err[ncomps][slice_max];
		double p[ncomps][slice_max], p_err[ncomps][slice_max];
		double y_sum[ncomps], y_sum_err[ncomps];

		for(int i=0;i<ncomps;++i){
			y_sum[i] = 0;
			y_sum_err[i] = 0;
		}
		double _start = slice_range[type][0];
		double _end = slice_range[type][1];
		double _step = (_end-_start) / slices[type];

		for(int s=0;s<slices[type];++s){
			for(int j=0;j<ncomps;++j){
				x[j][s] = _start + _step * (s+0.5) + j * 0.1 * _step;
				if(s==0){
					y[j][s] = 1;
					y_err[j][s] = 0;
				}else{
					int idx = offset + ncomps * (s-1) + j;
					y[j][s] = fabs(par[idx]);
					//y_err[j][s] = par_err[idx];
					y_err[j][s] = sqrt(minim->CovMatrix(idx, idx));
				}
				y_sum[j] += y[j][s];
				if(s>0)
					for(int s1=0;s1<slices[type]-1;++s1)
						y_sum_err[j] += minim->CovMatrix(offset + ncomps * (s-1) + j, offset + ncomps * s1 + j);
			}
		}
		
		for(int i=0;i<ncomps;++i){
			y_sum_err[i] = sqrt(y_sum_err[i]);
			for(int j=0;j<slices[type];++j){
				p[i][j] = y[i][j] / y_sum[i];
				double _tmp_cov=0;
				if(j!=0)
					for(int k=1;k<slices[type];++k){
						_tmp_cov += minim->CovMatrix(offset + ncomps * (k-1) + i, offset + ncomps * (j-1) + i);		
						//cout << i << " " << j << " " << k << " " << minim->VariableName(offset + ncomps * (k-1) + i) << " " << minim->VariableName(offset + ncomps * (j-1) + i) << " " << minim->CovMatrix(offset + ncomps * (k-1) + i, offset + ncomps * (j-1) + i) << endl;
						//if(j==k)
						//	cout << par_err[offset + ncomps * (j-1) + i] << endl;
					}
					
				p_err[i][j] = p[i][j] * sqrt( pow(y_err[i][j] / y[i][j], 2) + pow(y_sum_err[i] / y_sum[i], 2) - 2 * _tmp_cov / (y[i][j] * y_sum[i]) );
				//if(type>1)
				//	cout << j << " " << p[i][j] << " " << p_err[i][j] << " " << y[i][j] << " " << y_sum[i] << " " << _tmp_cov << " " << y_err[i][j] / y[i][j] << " " << y_sum_err[i] / y_sum[i] << " " << _tmp_cov / (y[i][j] * y_sum[i]) << " " << pow(y_err[i][j] / y[i][j], 2) + pow(y_sum_err[i] / y_sum[i], 2) - 2 * _tmp_cov / (y[i][j] * y_sum[i]) << endl;	
			}
		}

		TMultiGraph *mg = new TMultiGraph();
		TGraph *gr[ncomps];
		for(int i=0;i<ncomps;++i){
			gr[i] = new TGraphErrors(slices[type], x[i], p[i], 0, p_err[i]);
			gr[i]->SetMarkerStyle(20);
			gr[i]->SetMarkerColor(i+2);
			gr[i]->SetFillStyle(0);
			gr[i]->SetFillColor(0);
			if(i==0)
				gr[i]->SetTitle("IBD");
			else if(i==1)
				gr[i]->SetTitle("Li9");
			else
				gr[i]->SetTitle("B12");
			mg->Add(gr[i]);
		}
		char buf[255];
		sprintf(buf, "EH%d %s distribution", site+1, slice_vars[type]);	
		mg->SetTitle(buf);
		mg->SetMinimum(0);
		//mg->SetMaximum(0.55);
		mg->Draw("APC");
		c1->BuildLegend(0.7, 0.7, 0.9, 0.9);
		if(strstr(slice_vars[type],"ep") != NULL){
			TFile *file_spec = new TFile("./data/toyli9spec_BCWmodel_v1.root","READ");
			TH1 *h_spec = file_spec->Get("h_eVisAllSmeared");
			vector<double> bins;	
			bins.push_back(slice_range[type][0]);
			for(int i=0;i<slices[type];++i){
				double edge = (slice_range[type][1] - slice_range[type][0]) / slices[type] * (i+1) + slice_range[type][0];
				bins.push_back(edge);
			}
			TH1 *h_spec_rebin = h_spec->Rebin(bins.size()-1, "h_spec_rebin", &bins[0]);
			h_spec_rebin->SetLineColor(3);
			h_spec_rebin->SetLineStyle(9);
			h_spec_rebin->DrawNormalized("same");
		}
		char buf[255];
		char *dir = "./plots/cfits_slice";
		sprintf(buf,"%s/cfit_%s_%d.png",dir,slice_vars[type],site);
		gPad->SetLogx(0);
		gPad->SetLogy(0);
		c1->SaveAs(buf);
		offset += ncomps * (slices[type] - 1);
	}//type
}

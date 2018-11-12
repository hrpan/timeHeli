#include "util.C"
#include "range.h"
#include "slice.h"

const char *hists_prefix = "./hists";

const double tau_li9 = 256.366;
const double tau_li9_err = 0.866;

const double tau_he8 = 171.68;

const double tau_b12 = 29.142;
const double tau_b12_err = 0.0288;

const double tau_n12 = 15.9;

double fitMin = 1.5;
double fitMax = 5000;

const int npars_max = 500;
const bool use_tagging = false;

char par_names[npars_max][255]; 

const int npars_single = 10;
//[0]:mu rate
//[1]:n_dc
//[2]:n_li9he8
//[3]:t_li9
//[4]:t_he9
//[5]:ratio_lihe
//[6]:n_b12
//[7]:t_b12
//[8]:t_n12
//[9]:r_b12n12
const int n_pdfs = 3;
const char *_pdfs[n_pdfs] = {
	"[6] * ([9] * ([0] + 1 / [7]) * exp(-([0] + 1 / [7]) * x ) + (1 - [9]) * ([0] + 1 / [8]) * exp(-([0] + 1 / [8]) * x) )",
	"[2] * ([5] * ([0] + 1 / [3]) * exp(-([0] + 1 / [3]) * x ) + (1 - [5]) * ([0] + 1 / [4]) * exp(-([0] + 1 / [4]) * x) )",
	"[1] * [0] * exp(-[0] * x )",
};

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
	double n_bo[n_range];
	double eps_t_lihe[n_range];
	double eps_t_bo[n_range];
	double eps_s_dc[slice_types][slice_max];
	double eps_s_lihe[slice_types][slice_max];
	double eps_s_bo[slice_types][slice_max];
	//double eps_bo;
	double n_dc;
	double r_lihe;
	double r_boni;

	int p_idx = 0;
	for(int r = 0; r < n_range; ++r){
		if(use_tagging){
			r_mu_tag[r] = par[p_idx++];
			r_mu_atag[r] = par[p_idx++];
			n_lihe[r] = par[p_idx++] * scale;
			eps_t_lihe[r] = par[p_idx++];
			n_bo[r] = par[p_idx++] * scale;
			eps_t_bo[r] = par[p_idx++];
		}else{
			r_mu[r] = par[p_idx++];
			n_lihe[r] = par[p_idx++] * scale;
			n_bo[r] = par[p_idx++] * scale;
		}
	}

	for(int _type=0;_type < slice_types; ++_type){
		eps_s_dc[_type][0] = 1;
		eps_s_lihe[_type][0] = 1;
		eps_s_bo[_type][0] = 1;

		double _tmp_sum[3] = {1, 1, 1};
		int _slices = slices[_type];
		for(int _s=0;_s<_slices-1;++_s){
			eps_s_dc[_type][_s+1] = par[p_idx++];
			_tmp_sum[0] += eps_s_dc[_type][_s+1];
			eps_s_lihe[_type][_s+1] = par[p_idx++];
			_tmp_sum[1] += eps_s_lihe[_type][_s+1];
			eps_s_bo[_type][_s+1] = par[p_idx++];
			_tmp_sum[2] += eps_s_bo[_type][_s+1];
		}
		for(int _s=0;_s<_slices;++_s){
			eps_s_dc[_type][_s] /= _tmp_sum[0];
			eps_s_lihe[_type][_s] /= _tmp_sum[1];
			eps_s_bo[_type][_s] /= _tmp_sum[2];
		}
	}
	r_lihe = par[p_idx++];
	r_boni = par[p_idx++];
	n_dc = par[p_idx++] * scale;

	for(int r=0; r<n_range; ++r){
		double n_dc_fake_lihe = 0;
		double n_dc_fake_bo = 0;
		for(int j=0;j<n_range;++j){
			if(r==j) continue;
			n_dc_fake_lihe += n_lihe[j];
			n_dc_fake_bo += n_bo[j];
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
							+ (( 1 - eps_t_lihe[r] ) * _n_lihe + ( 1 - eps_t_bo[r] ) * _n_bo); 
						_par[r][t][_type][s][2] = eps_t_lihe[r] * _n_lihe;
						_par[r][t][_type][s][6] = eps_t_bo[r] * _n_bo;
					}else{
						_par[r][t][_type][s][0] = r_mu_atag[r];
						_par[r][t][_type][s][1] = _n_dc + _n_dc_fake
							+ ((eps_t_lihe[r] * _n_lihe) + (eps_t_bo[r] * _n_bo));
						_par[r][t][_type][s][2] = ( 1 - eps_t_lihe[r] ) * _n_lihe;
						_par[r][t][_type][s][6] = ( 1 - eps_t_bo[r] ) * _n_bo;
					}
					_par[r][t][_type][s][3] = tau_li9;
					_par[r][t][_type][s][4] = tau_he8;
					_par[r][t][_type][s][5] = r_lihe;
					//_par[r][t][_type][s][7] = tau_b12;
					_par[r][t][_type][s][7] = tau_b12;
					_par[r][t][_type][s][8] = tau_n12;
					_par[r][t][_type][s][9] = r_boni;
					
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
					ibd.SetParameters(_par[r][t][type][s]);
					ibd.SetLineColor(kGreen);
					sprintf(buf, "%s/cfit_%d_%d_%d_%d_%d.png", dir, site, r, t, type, s);
					h[r][t][type][s]->Draw("E1");
					func[r][t][type][s]->Draw("same");
					ibd.Draw("same");
					gPad->SetLogx();
					c1->SaveAs(buf);
				}
			}
			if(!use_tagging) break;
		}
	
}

void initialize_minimizer(int site, ROOT::Math::Minimizer *minim){
	cout << "INITIALIZING MINIMIZER" << endl;
	minim->SetErrorDef(0.5);
	minim->SetTolerance(1e-1);
	minim->SetPrintLevel(1);
	//minim->SetPrecision(1e-16);
	minim->SetStrategy(2);
	minim->SetMaxFunctionCalls(1000000);
	minim->SetMaxIterations(1000000);

	double rmu_init;
	double n_lihe_init = 1e3;
	double n_bo_init = 1e3;
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
			printf(init_str, p_idx-1, par_names[p_idx-1], rmu_init, rmu_init * step_ratio);

			h[r][2][0][0]->Fit("expo","LQN0","", 300, 5000);
			rmu_init = fabs(expo->GetParameter(1));
			sprintf(buf, "rmu_atag_%d", r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, rmu_init, rmu_init * step_ratio);
			printf(init_str, p_idx-1, par_names[p_idx-1], rmu_init, rmu_init * step_ratio);

			sprintf(buf, "n_lihe_%d", r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetLowerLimitedVariable(p_idx++, buf, n_lihe_init, n_lihe_init * step_ratio, 0);
			//minim->SetVariable(p_idx++, buf, n_lihe_init, n_lihe_init * step_ratio);
			printf(init_str, p_idx-1, par_names[p_idx-1], n_lihe_init, n_lihe_init * step_ratio);

			sprintf(buf, "eps_tag_lihe_%d", r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0, 1);
			printf(init_str, p_idx-1, par_names[p_idx-1], eps_init, eps_init * step_ratio);

			sprintf(buf, "n_bo_%d", r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetLowerLimitedVariable(p_idx++, buf, n_bo_init, n_bo_init * step_ratio, 0);
			printf(init_str, p_idx-1, par_names[p_idx-1], n_bo_init, n_bo_init * step_ratio);
			
			sprintf(buf, "eps_tag_bo_%d", r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0, 1);
			printf(init_str, p_idx-1, par_names[p_idx-1], eps_init, eps_init * step_ratio);

		}else{
			h[r][0][0][0]->Fit("expo","LQN0","", 300, 5000);
			rmu_init = fabs(expo->GetParameter(1));
			sprintf(buf, "rmu_%d",r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, rmu_init, rmu_init * step_ratio);
			printf(init_str, p_idx-1, par_names[p_idx-1], rmu_init, rmu_init * step_ratio);

			sprintf(buf, "n_lihe_%d",r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetLowerLimitedVariable(p_idx++, buf, n_lihe_init, n_lihe_init * step_ratio, 0);
			printf(init_str, p_idx-1, par_names[p_idx-1], n_lihe_init, n_lihe_init * step_ratio);

			sprintf(buf, "n_bo_%d",r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetLowerLimitedVariable(p_idx++, buf, n_bo_init, n_bo_init * step_ratio, 0);
			printf(init_str, p_idx-1, par_names[p_idx-1], n_bo_init, n_bo_init * step_ratio);

		}
	}

	for(int type=0;type<slice_types;++type){
		for(int s=0;s<slices[type]-1;++s){
			eps_init = h[0][0][type][s+1]->GetEntries() / h[0][0][type][0]->GetEntries();
			sprintf(buf, "eps_s_dc_%d_%d", type, s);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetLowerLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0);
			printf(init_str, p_idx-1, par_names[p_idx-1], eps_init, eps_init * step_ratio);

			sprintf(buf, "eps_s_lihe_%d_%d", type, s);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetLowerLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0);
			//minim->SetLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0, 1);
			printf(init_str, p_idx-1, par_names[p_idx-1], eps_init, eps_init * step_ratio);

			sprintf(buf, "eps_s_bo_%d_%d", type, s);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetLowerLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0);
			//minim->SetLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0, 1);
			printf(init_str, p_idx-1, par_names[p_idx-1], eps_init, eps_init * step_ratio);
		}
	}
	
	sprintf(par_names[p_idx], "r_lihe");
	minim->SetLimitedVariable(p_idx++, "r_lihe", r_init, r_init * step_ratio, 0, 1);
	printf(init_str, p_idx-1, par_names[p_idx-1], r_init, r_init * step_ratio);

	sprintf(par_names[p_idx], "r_boni");
	minim->SetLimitedVariable(p_idx++, "r_boni", r_init, r_init * step_ratio, 0, 1);
	printf(init_str, p_idx-1, par_names[p_idx-1], r_init, r_init * step_ratio);
	
	double n_dc_init = 0;
	for(int s=0;s<slices[0];++s){
		n_dc_init += h[0][0][0][s]->GetEntries();
	}
	sprintf(par_names[p_idx], "n_dc");
	minim->SetVariable(p_idx++, "n_dc", n_dc_init, n_dc_init * step_ratio);
	printf(init_str, p_idx-1, par_names[p_idx-1], n_dc_init, n_dc_init * step_ratio);


	par_names[p_idx][0] = 0;
}

void combinedFit(){


	gStyle->SetOptFit(1);
	char _f_sum[512];
	for(int i=0;i<n_pdfs;++i)
		if(i==0)
			sprintf(_f_sum,"%s",_pdfs[i]);
		else
			sprintf(_f_sum,"%s+%s",_f_sum,_pdfs[i]);
	cout << _f_sum << endl;
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
	
	ROOT::Fit::DataOptions opt;
	opt.fUseEmpty = true;
	for(int site=0;site<3;++site){
		do_fit(site, opt, results);
	}
	printResults(results);
}

void do_fit(int site, ROOT::Fit::DataOptions &opt, double results[3][npars_max][2]){

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
		npars = 5 * n_range + slice_pars + 2;
	else
		npars = 3 * n_range + slice_pars + 3;

	ROOT::Math::Functor f_tmp(likelihood, npars);
	ROOT::Minuit2::Minuit2Minimizer *minim = ROOT::Math::Factory::CreateMinimizer("Minuit2");
	minim->SetFunction(f_tmp);
	initialize_minimizer(site, minim);

	minim->Hesse();
	minim->Minimize();
	minim->Hesse();
	
	minim->PrintResults();
	/*
	minim->SetPrintLevel(1);
	for(int i=0;i<npars;++i){
		double low, up;
		if(strstr(par_names[i], "n_lihe") != NULL){
			minim->GetMinosError(i, low, up);
			cout << par_names[i] << " " << minim->X()[i] << " " << low << " " << up << endl;
		}
	}
	*/
	double *x = minim->X();
	double *x_err = minim->Errors();
	for(int i=0;i<npars;++i){
		results[site][i][0]	= x[i];
		results[site][i][1] = x_err[i];
	}

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
					y[j][s] = par[idx];
					y_err[j][s] = par_err[idx];
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
				p_err[i][j] = sqrt( pow(y_err[i][j] / y_sum[i], 2) + pow(y[i][j] * y_sum_err[i] / pow(y_sum[i], 2), 2));
			}
		}

		TMultiGraph *mg = new TMultiGraph();
		TGraph *gr[ncomps];
		for(int i=0;i<ncomps;++i){
			gr[i] = new TGraphErrors(slices[type], x[i], p[i], 0, p_err[i]);
			gr[i]->SetMarkerStyle(20);
			gr[i]->SetMarkerColor(i+2);
			mg->Add(gr[i]);
		}
		mg->SetMinimum(0);
		mg->Draw("APC");
		char buf[255];
		char *dir = "./plots/cfits_slice";
		sprintf(buf,"%s/cfit_%s_%d.png",dir,slice_vars[type],site);
		gPad->SetLogx(0);
		gPad->SetLogy(0);
		c1->SaveAs(buf);
		offset += ncomps * (slices[type] - 1);
	}//type
}

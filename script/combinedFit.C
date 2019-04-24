#include "util.C"
#include "range.h"
#include "../hists/slice.h"
#include "TH1.h"
#include "TCanvas.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TROOT.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <iostream>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TFile.h"
#include "TStyle.h"
#include "Fit/LogLikelihoodFCN.h"
#include "Fit/PoissonLikelihoodFCN.h"
#include "omp.h"
//#include <Eigen/Core>
#include "LBFGS.h"
#include "cppoptlib/meta.h"
#include "cppoptlib/boundedproblem.h"
#include "cppoptlib/solver/lbfgsbsolver.h"
#include "cppoptlib/solver/bfgssolver.h"
#include "cppoptlib/solver/conjugatedgradientdescentsolver.h"

using namespace std;
using namespace cppoptlib;
using namespace LBFGSpp;
//using Eigen::VectorXd;

const char *hists_prefix = "./hists";

const double tau_li9 = 256.366;
const double tau_li9_err = 0.866;

const double tau_he8 = 171.17;
const double tau_he8_err = 2.31;

const double tau_b12 = 29.142;
const double tau_b12_err = 0.0288;

const double tau_n12 = 15.9;

const double _tau_short = tau_n12;

const double rmu_scale = 1e6;

double fitMin = 1.5;
double fitMax = 5000;

const int npars_max = 200;
int npars;
const bool use_tagging = false;
const bool fix_tau_short = true;
bool use_predicted_li9 = true;

const bool use_isotope[3] = { 1, 1, 1 };
const bool use_slice[4] = { 1, 1, 1, 1 };

char par_names[npars_max][255]; 

const int npars_single = 11;
//[0]:mu rate
//[1]:n_dc
//[2]:n_li9he8
//[3]:t_li9
//[4]:n_b12
//[5]:t_b12
//[6]:n_n12
//[7]:t_n12
//[8]:scale
const int n_pdfs = 4;
const char *_pdfs[4] = {
	"[1] * [0] * exp(-[0] * x )",
	"[2] * ([0] + 1 / [3]) * exp(-([0] + 1 / [3]) * x )",
	"[4] * ([0] + 1 / [5]) * exp(-([0] + 1 / [5]) * x )",
	"[6] * ([0] + 1 / [7]) * exp(-([0] + 1 / [7]) * x )",
};

const int colors[6] = {
	2, 3, 4, 6, 7, 8
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

void do_fit(int site, ROOT::Fit::DataOptions &opt, double results[3][npars_max][2], double daily_rates[3][2]);
void plotSlice(int site, ROOT::Math::Minimizer *minim);

void printResults(double results[3][npars_max][2]){
	printf("%15s %21s %21s %21s\n", "site", "EH1", "EH2", "EH3");
	for(int i=0;i<npars_max;++i){
		printf("%15s", par_names[i]);
		for(int j=0;j<3;++j)
			printf(" %10.3e/%10.3e", results[j][i][0], results[j][i][1]);
		cout << endl;
		if(par_names[i+1][0] == 0) break;
	}

}

ROOT::Fit::PoissonLikelihoodFCN<ROOT::Math::IBaseFunctionMultiDim> *fPNLL[n_range][3][slice_types][slice_max];
ROOT::Math::IMultiGenFunction *fChi2[n_range][3][slice_types][slice_max];
ROOT::Math::WrappedMultiTF1 *wf[n_range][3][slice_types][slice_max];

TF1 *func[n_range][3][slice_types][slice_max];
TH1 *h[n_range][3][slice_types][slice_max];

//double _par[n_range][3][slice_types][slice_max][npars_single];

void parTrans(const double *par, double _par[n_range][3][slice_types][slice_max][npars_single]){

	double scale = h[0][0][0][0]->GetBinWidth(1);

	double r_mu[n_range];
	double r_mu_tag[n_range];
	double r_mu_atag[n_range];
	double n_lihe_tag[n_range];
	double n_lihe_atag[n_range];
	double n_bo_tag[n_range];
	double n_bo_atag[n_range];
	double n_iso[n_range][n_pdfs-1];
	double eps_s[slice_types][n_pdfs][slice_max];
	double n_dc;
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
			for(int iso=0; iso<n_pdfs-1;++iso)
				//n_iso[r][iso] = fabs(par[p_idx++]) * scale * lt;
                //if(iso == 0 && r == n_range - 1)
                //    n_iso[r][iso] = (fabs(par[p_idx++]) - n_iso[0][0] - n_iso[1][0]) / 0.211;
                //else
    			//	n_iso[r][iso] = fabs(par[p_idx++]);
    			//n_iso[r][iso] = fabs(par[p_idx++]);
    			n_iso[r][iso] = par[p_idx++];
		}
	}
	for(int _type=0;_type < slice_types; ++_type){
		if(!use_slice[_type]) continue;
		for(int i=0;i<n_pdfs;++i){
            double _tmp_sum = 0;
            for(int _s=0; _s<slices[_type]; ++_s){
                eps_s[_type][i][_s] = par[p_idx++];
                _tmp_sum += eps_s[_type][i][_s];
            }
			for(int _s=0;_s<slices[_type];++_s){
				eps_s[_type][i][_s] /= _tmp_sum;
			}
		}
	}

	n_dc = par[p_idx++];
	tau_short = par[p_idx++];

	for(int r=0; r<n_range; ++r){
		double n_dc_fake[n_pdfs-1];
		for(int iso=0;iso<n_pdfs-1;++iso)
			n_dc_fake[iso] = 0;
		for(int j=0; j<n_range; ++j){
			if(r==j) continue;
			if(use_tagging){
		//		double n_dc_fake_lihe += (n_lihe_tag[j] + n_lihe_atag[j]);
		//		double n_dc_fake_bo += (n_bo_tag[j] + n_bo_atag[j]);
			}else{
				for(int iso=0;iso<n_pdfs-1;++iso)
					n_dc_fake[iso] += n_iso[j][iso];
			}
		}

		for(int t=0; t<3; ++t){
			for(int _type=0; _type<slice_types; ++_type){
				if(!use_slice[_type]) continue;
				int _slices = slices[_type];
				for(int s=0; s<_slices; ++s){
					double _n_dc = n_dc * eps_s[_type][0][s];
					double _n_dc_fake = 0;
					double _n_iso[n_pdfs-1];
					for(int iso=0;iso<n_pdfs-1;++iso){
						_n_dc_fake += n_dc_fake[iso] * eps_s[_type][iso+1][s];
						_n_iso[iso] = n_iso[r][iso] * eps_s[_type][iso+1][s];
					}
					/*
					if(use_tagging){
						_n_lihe = (n_lihe_tag[r] + n_lihe_atag[r]) * _eps_s_lihe[s];
						_n_bo = (n_bo_tag[r] + n_bo_atag[r]) * _eps_s_bo[s];
					}
					*/
					if(t==0){
						if(use_tagging)
							_par[r][t][_type][s][0] = r_mu_tag[r] + r_mu_atag[r];
						else
							_par[r][t][_type][s][0] = r_mu[r];
						_par[r][t][_type][s][1] = _n_dc + _n_dc_fake;
						//_par[r][t][_type][s][2] = _n_iso[0];
						//_par[r][t][_type][s][6] = _n_iso[1];
						//_par[r][t][_type][s][8] = _n_iso[2];
                        for(int iso=0;iso<n_pdfs-1;++iso)
    						_par[r][t][_type][s][2*iso+2] = _n_iso[iso];
					}else if(t==1){
                    /*
						_par[r][t][_type][s][0] = r_mu_tag[r];
						_par[r][t][_type][s][1] = _n_dc + _n_dc_fake 
							+ eps_s_lihe[s] * n_lihe_atag[r] 
							+ eps_s_bo[s] * n_bo_atag[r]; 
						_par[r][t][_type][s][2] = eps_s_lihe[s] * n_lihe_tag[r];
						_par[r][t][_type][s][6] = eps_s_bo[s] * n_bo_tag[r];
                    
					}else{
						_par[r][t][_type][s][0] = r_mu_atag[r];
						_par[r][t][_type][s][1] = _n_dc + _n_dc_fake
							+ eps_s_lihe[s] * n_lihe_tag[r] 
							+ eps_s_bo[s] * n_bo_tag[r];
						_par[r][t][_type][s][2] = eps_s_lihe[s] * n_lihe_atag[r];
						_par[r][t][_type][s][6] = eps_s_bo[s] * n_bo_atag[r];
                    */
					}
					_par[r][t][_type][s][3] = tau_li9;
					_par[r][t][_type][s][5] = tau_b12;
					_par[r][t][_type][s][7] = tau_n12;
		            _par[r][t][_type][s][8] = scale * lt;				
									
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

    static double _par[n_range][3][slice_types][slice_max][npars_single];

	parTrans(par, _par);
	double result = 0;
	for(int r = 0; r < n_range; ++r){
		for(int t = 0; t < 3; ++t){
			if(use_tagging && t==0) continue;
			for(int type = 0; type < slice_types; ++type){
				if(!use_slice[type]) continue;
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

double likelihood_derivative(const double *par, int coord){

    static double _par[n_range][3][slice_types][slice_max][npars_single];

	parTrans(par, _par);
    static double par_cache[npars_max];
    static double g_cache[n_range][3][slice_types][slice_max][npars_single];
    static double eps_cache[n_pdfs][slice_types][slice_max];
    static double eps_sum_cache[n_pdfs][slice_types];
    static double p_cache[n_pdfs][slice_types][slice_max];
    bool recalculate = false;
    for(int i=0;i<npars;++i){
        if(par[i] != par_cache[i]){
            par_cache[i] = par[i];
            recalculate = true;
        }
    }

    if(recalculate){
        for(int r = 0; r < n_range; ++r){
            for(int type = 0; type < slice_types; ++type){
                if(!use_slice[type]) continue;
                for(int s = 0; s < slices[type]; ++s){
                    fPNLL[r][0][type][s]->Gradient(_par[r][0][type][s], g_cache[r][0][type][s]);
                }
            }
        }
        int _offset = n_pdfs * n_range;
        for(int type = 0; type < slice_types; ++type){
            if(!use_slice[type]) continue;
            for(int iso = 0; iso < n_pdfs; ++iso){
                double _tmp_sum = 0;
                for(int s = 0; s < slices[type]; ++s){
                    eps_cache[iso][type][s] = par[_offset + s];
                    _tmp_sum += eps_cache[iso][type][s];
                }
                for(int s = 0; s < slices[type]; ++s){
                    p_cache[iso][type][s] = eps_cache[iso][type][s] / _tmp_sum;
                }
                eps_sum_cache[iso][type] = _tmp_sum;
                _offset += slices[type];
            }
        }
    }

    double _g = 0;
    double *g;
    if(coord < n_range * n_pdfs){
        int r = coord / n_pdfs;
        int d = coord % n_pdfs;
        if(d != 0){
            for(int _r = 0; _r < n_range; ++_r){
                for(int type = 0; type < slice_types; ++type){
                    if(!use_slice[type]) continue;
                    for(int s = 0; s < slices[type]; ++s){
                        g = g_cache[_r][0][type][s];
                        if(_r == r)
                            _g += g[2 * d] * p_cache[d][type][s]; 
                        else
                            _g += g[1] * p_cache[d][type][s]; 
                    }
                }

            }
            
        }else{
            for(int type = 0; type < slice_types; ++type){
                if(!use_slice[type]) continue;
                for(int s = 0; s < slices[type]; ++s){
                    g = g_cache[r][0][type][s];
                    _g += g[0]; 
                }
            }
            _g;
        }
        
    }else if(strncmp(par_names[coord], "eps_s", 5) == 0){
        int _c;
        if(strstr(par_names[coord], "dc") != NULL)
            _c = 0;
        else if(strstr(par_names[coord], "li") != NULL)
            _c = 1;
        else if(strstr(par_names[coord], "bo") != NULL)
            _c = 2;
        else if(strstr(par_names[coord], "ni") != NULL)
            _c = 3;
        
        int type = par_names[coord][9] - '0';
        int _s = par_names[coord][11] - '0'; 

        for(int r = 0; r < n_range; ++r){
            double _n = 0;
            if(_c > 0)
                for(int r_1 = 0; r_1 < n_range; ++r_1)
                    if(r_1 != r) 
                        _n += fabs(par[n_pdfs * r_1 + _c]);
                    

            for(int s = 0; s < slices[type]; ++s){
                g = g_cache[r][0][type][s];
                int _idx;
                double __n = 0;
                if(_c == 0){
                    _idx = 1;
                    for(int _i=0;_i<npars_max;++_i)
                        if(strncmp(par_names[_i],"n_dc",4) == 0){
                            __n = par[_i];
                            break;
                        }
                }else{
                    _idx = _c * 2;
                    __n = par[n_pdfs * r + _c];
                }
                double p_der = (s == _s)? ( 1 - p_cache[_c][type][s] ) / eps_sum_cache[_c][type] : - p_cache[_c][type][s] / eps_sum_cache[_c][type]; 
                //double p_der = (s == _s)? p_s_0 * ( 1 - p_s_0 ): -p_s_0 * p_s_1;
                _g += g[_idx] * __n * p_der;
                if( _c > 0 )
                    _g += g[1] * _n * p_der;
            }
        }   
    }else if(strncmp(par_names[coord], "n_dc", 4) == 0){
        for(int r = 0; r < n_range; ++r){
            int _offset = n_range * n_pdfs;
            for(int t = 0;t < slice_types; ++t){
                for(int s = 0; s < slices[t]; ++s){
                    g = g_cache[r][0][t][s];
                    //_g += g[1] * _eps_tmp[s]; 
                    //_g += g[1] * eps_cache[0][t][s] / eps_sum_cache[0][t]; 
                    _g += g[1] * p_cache[0][t][s]; 
                }
                _offset += slices[t] * n_pdfs;
            }
        }
    }
//    cout << coord << " " << par[coord] << " " << _g << endl;
    return _g;
    
}

/*
class MyFunc{
    public:
    double operator()(const VectorXd &x, VectorXd &grad){

        double _x[1000], _grad[1000];
        for(int i=0;i<npars;++i)
            _x[i] = x[i];
        for(int i=0;i<npars;++i)
            grad[i] = likelihood_derivative(_x, i); 
        double f = likelihood(_x);
        printf("%.10f\n", f);
        fflush(stdout);
        return f;
    } 
};
*/

namespace cppoptlib{

template<typename T>
class MyFunc2 : public BoundedProblem<T>{
    double _x[1000];
    public:
        using Superclass = BoundedProblem<T>;
        using typename Superclass::TVector;

    public:
        MyFunc2(int npars) :
            Superclass(npars){}

        T value(const TVector &x) {
            
            for(int i=0;i<npars;++i)
                _x[i] = x[i];
            double f = likelihood(_x);
            //printf("%.10f\n", f);
            //fflush(stdout);
            return f;
            
        }

        void gradient(const TVector &x, TVector &grad) {
            
            for(int i=0;i<npars;++i)
                _x[i] = x[i];
            for(int i=0;i<npars;++i)
                grad[i] = likelihood_derivative(_x, i);
            
        }

};

}


double chi2(const double *par){

    static double _par[n_range][3][slice_types][slice_max][npars_single];

	parTrans(par, _par);
	
	double result = 0;
	for(int r = 0; r < n_range; ++r){
		for(int t = 0; t < 3; ++t){
			if(use_tagging && t==0) continue;
			for(int type = 0; type < slice_types; ++type){
				if(!use_slice[type]) continue;
				for(int s = 0; s < slices[type]; ++s){
					double _tmp = (*fChi2[r][t][type][s])(_par[r][t][type][s]);
					result += _tmp;
				}//slice
			}//slice type
			if(!use_tagging) break;
		}//tag
	}//range
	return result;

}

void fillPars(const double *par, TF1 *f[n_range][3][slice_types][slice_max]){

    static double _par[n_range][3][slice_types][slice_max][npars_single];
	parTrans(par, _par);

	for(int r=0;r<n_range;++r)
		for(int t=0;t<3;++t){
			for(int type=0;type<slice_types;++type){
				if(!use_slice[type]) continue;
				for(int s=0;s<slices[type];++s)
					f[r][t][type][s]->SetParameters(_par[r][t][type][s]);
			}
			if(!use_tagging) break;
		}//tag

}

void plotHists(int site, const double *par){

	fillPars(par, func);
	const char *dir = "./plots/cfits";
	char buf[255];
    static double _par[n_range][3][slice_types][slice_max][npars_single];
	parTrans(par, _par);
	
	char _fs[1000];
	TF1 *_tf1s[n_pdfs];
	for(int i=0;i<n_pdfs;++i){
        _fs[0] = 0;
		for(int j=0;j<i+1;++j){
            if(j!=0)
                strcat(_fs, " + ");
            strcat(_fs, _pdfs[j]);
        }
        char _final[1000];
        sprintf(_final, "[8] * (%s) ", _fs);
		sprintf(buf,"comp%d",i);
		_tf1s[i] = new TF1(buf, _final, 0, 5000);
		_tf1s[i]->SetLineColor(colors[i]);
	}

	for(int r=0;r<n_range;++r)
		for(int t=0;t<3;++t){
			for(int type=0;type<slice_types;++type){
				if(!use_slice[type]) continue;
				for(int s=0;s<slices[type];++s){
                    TCanvas *c1 = new TCanvas("c1");
					sprintf(buf, "%s/cfit_%d_%d_%d_%d_%d.png", dir, site, r, t, type, s);
					h[r][t][type][s]->Draw("E1");
					func[r][t][type][s]->Draw("same");
					for(int f=0;f<n_pdfs;++f){
						_tf1s[f]->SetParameters(_par[r][t][type][s]);
						_tf1s[f]->Draw("same");
					}
					gPad->SetLogx();
					c1->SaveAs(buf);
                    delete c1;
				}
			}
			if(!use_tagging) break;
		}
	
}

void initialize_minimizer(int site, ROOT::Math::Minimizer *minim, bool verbose){
	cout << "INITIALIZING MINIMIZER" << endl;
	minim->SetErrorDef(0.5);
	//minim->SetTolerance(0.001);
	minim->SetTolerance(1);
    if(verbose)
    	minim->SetPrintLevel(2);
    else
    	minim->SetPrintLevel(1);
	minim->SetStrategy(2);
	minim->SetMaxFunctionCalls(1000000);
	minim->SetMaxIterations(1000000);

	double eps_init = 0.5;

	const char *_pre[5] = {
		"dc",
		"li",
        "bo",
		"ni",
		"ni"
	};
    
    const double rmu_init[3][3] = {
        { 1.178e-2, 8.613e-3, 1.441e-4 },
        { 8.118e-3, 7.128e-3, 1.412e-4 },
        { 5.825e-4, 4.615e-4, 1.373e-5 },
    };
    
    /*
      rmu_0  1.178e-02/ 6.739e-06  8.118e-03/ 4.500e-06  5.825e-04/ 8.595e-07
      rmu_1  8.613e-03/ 5.307e-06  7.128e-03/ 4.544e-06  4.615e-04/ 8.604e-07
      rmu_2  1.441e-04/ 1.743e-07  1.412e-04/ 1.760e-07  1.373e-05/ 6.071e-08
    */
    const double n_iso_init[3][3][3] = {
        {
            { 0.01, 1.58, 0.84},
            { 1.65, 1.66, 0.11},
            { 4.60, 1.36, 0.92}
        },
        {
            { 0.01, 1.14, 0.05},
            { 1.14, 1.30, 0.05},
            { 3.71, 1.29, 0.47}
        },
        {
            { 0.01, 0.09, 0.01},
            { 0.15, 0.13, 0.01},
            { 0.45, 0.17, 0.09}
        }
    };
    /*
    n_li_0  5.440e-04/ 2.485e+00  2.913e-08/ 2.994e-01  1.061e-13/ 9.430e-03
    n_bo_0  1.580e+00/ 1.404e-01  1.136e+00/ 1.004e-01  8.671e-02/ 9.803e-03
    n_ni_0  8.458e-02/ 1.302e-01  4.518e-02/ 7.161e-02  2.723e-03/ 5.413e-03
    n_li_1  1.641e+00/ 8.394e-01  1.138e+00/ 6.248e-01  1.468e-01/ 3.046e-02
    n_bo_1  1.660e+00/ 1.318e-01  1.295e+00/ 1.044e-01  1.254e-01/ 8.885e-03
    n_ni_1  1.145e-01/ 1.202e-01  5.101e-02/ 7.329e-02  5.169e-09/ 1.523e-02
    n_li_2  4.597e+00/ 1.367e-01  3.712e+00/ 1.246e-01  4.509e-01/ 1.230e-02
    n_bo_2  1.359e+00/ 7.931e-02  1.294e+00/ 7.512e-02  1.726e-01/ 1.546e-02
    n_ni_2  9.173e-01/ 7.183e-02  4.700e-01/ 6.310e-02  9.028e-02/ 1.308e-02
    */
	double step_ratio = 1e-2;
	char buf[255];
	int p_idx = 0;
	const char *init_str = "Initializing parameter %2d %20s %10.2e %10.2e\n";
	for(int r=0;r<n_range;++r){
		if(use_tagging){	
            /*
            TF1 *f_tmp = new TF1("f_tmp","expo", 300, fitMax);	
			h[r][1][0][0]->Fit(f_tmp,"LQN0","", 300, fitMax);
			rmu_init = fabs(f_tmp->GetParameter(1));
			sprintf(buf, "rmu_tag_%d",r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, rmu_init, rmu_init * step_ratio);

			h[r][2][0][0]->Fit(f_tmp,"LQN0","", 300, fitMax);
			rmu_init = fabs(f_tmp->GetParameter(1));
			sprintf(buf, "rmu_atag_%d", r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, rmu_init, rmu_init * step_ratio);

			for(int i=1;i<n_pdfs;++i){
				for(int t=0;t<2;++t){
					sprintf(buf, "n_%s_%d_%d", _pre[i], r, t);
					sprintf(par_names[p_idx], "%s", buf);
					minim->SetVariable(p_idx++, buf, n_iso_init, n_iso_init * step_ratio);
				}
			}
	        */
		}else{
			sprintf(buf, "rmu_%d",r);
			sprintf(par_names[p_idx], "%s", buf);
            double _tmp = rmu_init[site][r];
			minim->SetVariable(p_idx++, buf, rmu_init[site][r], rmu_init[site][r] * step_ratio);
			//minim->SetLimitedVariable(p_idx++, buf, _tmp, step_ratio * _tmp, 0.9 * _tmp, 1.1 * _tmp);

			for(int i=1;i<n_pdfs;++i){
				sprintf(buf, "n_%s_%d",_pre[i], r);
				sprintf(par_names[p_idx], "%s", buf);
    			//minim->SetVariable(p_idx++, buf, n_iso_init, n_iso_init * step_ratio);
                double n_tmp = n_iso_init[site][r][i-1];    
                if(!use_isotope[i-1])
                    minim->SetFixedVariable(p_idx++, buf, 0);
                else
    			    minim->SetLowerLimitedVariable(p_idx++, buf, n_tmp, n_tmp * step_ratio, 0);
			}
			
		}
	}

	for(int type=0;type<slice_types;++type){
		if(!use_slice[type]) continue;
		for(int i=0;i<n_pdfs;++i){

            if(i == 1 && type==0){
        		TFile *file_spec = new TFile("./data/toyli9spec_BCWmodel_v1.root","READ");
                TH1 *h_spec = (TH1*)file_spec->Get("h_eVisAllSmeared");
                vector<double> bins;	
                bins.push_back(slice_range[type][0]);
                for(int s=0;s<slices[type];++s){
                    double edge = (slice_range[type][1] - slice_range[type][0]) / slices[type] * (s+1) + slice_range[type][0];
                    bins.push_back(edge);
                }
                TH1 *h_spec_rebin = h_spec->Rebin(bins.size()-1, "h_spec_rebin", &bins[0]);
                h_spec_rebin->Scale(1.0/h_spec_rebin->GetBinContent(1));
                for(int s=0;s<slices[type];++s){
                    sprintf(buf, "eps_s_%s_%d_%d", _pre[i], type, s);
				    sprintf(par_names[p_idx], "%s", buf);
                    double _tmp = h_spec_rebin->GetBinContent(s+1);
                    if(use_predicted_li9 || s==0)
        				minim->SetFixedVariable(p_idx++, buf, _tmp);
                    else
        				minim->SetLowerLimitedVariable(p_idx++, buf, _tmp, _tmp * step_ratio, 0);
    				//minim->SetFixedVariable(p_idx++, buf, log(h_spec_rebin->GetBinContent(s+1)/h_spec_rebin->GetBinContent(1)));
                }
                file_spec->Close();
            }else if( i>0 && !use_isotope[i-1] ){ 
                for(int s=0;s<slices[type];++s){
                    sprintf(buf, "eps_s_%s_%d_%d", _pre[i], type, s);
                    sprintf(par_names[p_idx], "%s", buf);
                    minim->SetFixedVariable(p_idx++, buf, 1);
                }
            }else{
                for(int s=0;s<slices[type];++s){
                    if(i == 0 || i == 1)
                        eps_init = h[0][0][type][s]->GetEntries() / h[0][0][type][0]->GetEntries();
                    else 
                        eps_init = 1;
                        
                    sprintf(buf, "eps_s_%s_%d_%d", _pre[i], type, s);
                    sprintf(par_names[p_idx], "%s", buf);
                    //minim->SetVariable(p_idx++, buf, eps_init, eps_init * step_ratio);
                    if(s==0)
                        minim->SetFixedVariable(p_idx++, buf, 1);
                        //minim->SetFixedVariable(p_idx++, buf, 0);
                    else
                        minim->SetLowerLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0);
                        //minim->SetLimitedVariable(p_idx++, buf, log(eps_init), 0.1, -10, 10);
                }

            }
		}
	}

	double n_dc_init = 0;
	for(int s=0;s<slices[0];++s)
		n_dc_init += h[0][0][0][s]->GetEntries();
	n_dc_init /= lt;	

	sprintf(par_names[p_idx], "n_dc");
	minim->SetVariable(p_idx++, "n_dc", n_dc_init, n_dc_init * step_ratio);
	
	double tau_init = 10;
	sprintf(par_names[p_idx], "tau_short");
	if(fix_tau_short)
		minim->SetFixedVariable(p_idx++, "tau_short", _tau_short);
	else
		minim->SetLimitedVariable(p_idx++, "tau_short", tau_init, tau_init * step_ratio, 1, 50);

	par_names[p_idx][0] = 0;

	if(verbose)
		for(int i=0;i<p_idx;++i)
			printf(init_str, i, par_names[i], minim->X()[i], minim->Errors()[i]);
	cout << endl;
}

int main(int argc, char **argv){

    if(argc < 2 || strcmp(argv[1],"fix") == 0)
        use_predicted_li9 = true;
    else if(strcmp(argv[1],"free") == 0)
        use_predicted_li9 = false;
        
	gStyle->SetOptFit(1);
    
	char _f_sum[512] = "";
    char _f_tmp[512] = "";
	for(int i=0;i<n_pdfs;++i){
        if(i!=0) 
            strcat(_f_tmp, " + ");
        strcat(_f_tmp, _pdfs[i]);
    }
    sprintf(_f_sum,"[8] * ( %s )", _f_tmp);
	cout << "INITIALIZING FUNCTION WRAPPERS" << endl;
	char buf[255];	
	for(int r=0;r<n_range;++r){
		for(int t=0;t<3;++t){
			for(int type=0;type<slice_types;++type){
				if(!use_slice[type]) continue;
				for(int s=0;s<slices[type];++s){
					sprintf(buf,"func%d_%d_%d_%d",r,t,type,s);
					func[r][t][type][s] = new TF1(buf, _f_sum, 0, 5000);
					wf[r][t][type][s] = new ROOT::Math::WrappedMultiTF1(*func[r][t][type][s], 1);
                    wf[r][t][type][s]->SetDerivPrecision(1e-7);
				}
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
    return 0;
}

void gradient_descent(const double x[1000]){

    double alpha = 1e-15;
    double beta = 0.9;

    double _x[1000];
    for(int i=0;i<npars;++i)
        _x[i] = x[i];

    double _g[1000], _m[1000];
    for(int i=0;i<npars;++i)
        _m[i] = 0;
    for(int iter=0;iter<100000;++iter){
        for(int i=0;i<npars;++i)
            _g[i] = likelihood_derivative(_x, i);
        double g2 = 0;
        for(int i=0;i<npars;++i)
            g2 += _g[i] * _g[i];
        for(int i=0;i<npars;++i)
            _m[i] = beta * _m[i] + alpha * _g[i];
        for(int i=0;i<npars;++i){
            double _tmp = _x[i] - _m[i];
            if(_tmp < 0){
                _m[i] = _x[i];
                _x[i] = 0;
            }else{
                _x[i] = _tmp;
            }
        }
        printf("Iter %d %.10f g2: %.10f x[5]: %.10f\n", iter, likelihood(_x), g2, _x[5]);
        fflush(stdout);
    }
    for(int i=0;i<npars;++i)
        cout << par_names[i] << " " << _x[i] << endl;

}

void do_fit(int site, ROOT::Fit::DataOptions &opt, double results[3][npars_max][2], double daily_rates[3][2]){

	char buf[255];

	const char *_suffix[3] = { "", "_tag", "_atag" };

	ROOT::Fit::BinData *data[n_range][3][slice_types][slice_max];
	ROOT::Fit::DataRange rangeD(fitMin, fitMax);
	ROOT::Fit::BinData *data_chi2[n_range][3][slice_types][slice_max];
	ROOT::Fit::DataOptions opt_chi2;
	int ndf = 0;
	for(int r=0;r<n_range;++r){
		for(int t=0;t<3;++t){
			for(int type=0;type<slice_types;++type){
				if(!use_slice[type]) continue;
				for(int s=0;s<slices[type];++s){
					sprintf(buf,"%s/EH%d_dtlSH_%d%s_%d_%d.root", hists_prefix, site+1, r, _suffix[t], type, s);
					cout << "READING FILE:" << buf << endl;
					TFile *f = new TFile(buf, "READ");
					h[r][t][type][s] = (TH1*) f->Get("h");
					data[r][t][type][s] = new ROOT::Fit::BinData(opt, rangeD);
					ROOT::Fit::FillData(*data[r][t][type][s], h[r][t][type][s]);
					data_chi2[r][t][type][s] = new ROOT::Fit::BinData(opt_chi2, rangeD);
					ROOT::Fit::FillData(*data_chi2[r][t][type][s], h[r][t][type][s]);
					ndf += data_chi2[r][t][type][s]->NPoints();
					fPNLL[r][t][type][s] = new ROOT::Fit::PoissonLikelihoodFCN<ROOT::Math::IBaseFunctionMultiDim>(*data[r][t][type][s], *wf[r][t][type][s]);	
					fChi2[r][t][type][s] = new ROOT::Fit::Chi2FCN<ROOT::Math::IBaseFunctionMultiDim>(*data_chi2[r][t][type][s], *wf[r][t][type][s]);	
				}
			}
			if(!use_tagging) break;
		}
	}

	int slice_pars = 0;
	for(int type=0;type<slice_types;++type){
		if(!use_slice[type]) continue;
		slice_pars += slices[type];
	}
	slice_pars *= n_pdfs;
	if(use_tagging)
		npars = 2 * n_pdfs * n_range + slice_pars + 3;
	else
		npars = n_pdfs * n_range + slice_pars + 2;
	//ROOT::Math::Functor f_tmp(&likelihood, npars);
	ROOT::Math::GradFunctor f_tmp(&likelihood, &likelihood_derivative, npars);
	ROOT::Math::Minimizer *minim = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
	ROOT::Math::Minimizer *minim_simplex = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
	minim->SetFunction(f_tmp);
	minim_simplex->SetFunction(f_tmp);
	initialize_minimizer(site, minim, true);
	initialize_minimizer(site, minim_simplex, false);
    minim_simplex->SetTolerance(1);
	minim_simplex->Minimize();
	minim->SetVariableValues(minim_simplex->X());
	minim->Hesse();
	minim->Minimize();

    /*
    LBFGSParam<double> param;
    param.epsilon=1e-6;
    param.max_iterations = 100000;
    LBFGSSolver<double> solver(param);
    MyFunc mf;
    VectorXd v_x = VectorXd::Zero(npars);
    for(int i=0;i<npars;++i)
        v_x[i] = minim->X()[i];
    double fx;
    solver.minimize(mf, v_x, fx);
    cout << fx << v_x << endl;
    */
     
    typedef double T;
    typedef typename MyFunc2<T>::TVector TVector;
    MyFunc2<double> mf(npars);
    TVector v_x(npars);
    for(int i=0;i<npars;++i)
        v_x[i] = minim->X()[i];
    mf.setLowerBound(TVector::Zero(npars));
    LbfgsbSolver<MyFunc2<double>> solver;
    Criteria<double> crit = Criteria<double>::defaults();
    crit.iterations = 10000000;
    crit.gradNorm = 1e-8;
    solver.setHistorySize(3);
    solver.setStopCriteria(crit); 
    solver.setDebug(DebugLevel::High);
    solver.minimize(mf, v_x);
    for(int i=0;i<npars;++i)
        cout << par_names[i] << " " << v_x[i] << endl;
    cout << "f: " << mf(v_x) << endl;
    
    double _x[1000];
    for(int i=0;i<npars;++i){
        _x[i] = minim->X()[i];
    }
    for(int i=0;i<npars;++i){
        double _eps = _x[i] * 0.0001;
        double l0 = likelihood(_x);
        _x[i] += _eps;
        double l1 = likelihood(_x);
        _x[i] -= _eps;
        cout << par_names[i] << " " << (l1-l0) / _eps << " " << likelihood_derivative(_x, i) << endl;
    }
    
	minim_simplex->Minimize();
	minim->SetVariableValues(minim_simplex->X());
	minim->Hesse();
	minim->Minimize();
    getchar();
    //gradient_descent(minim->X());
    /* 
	for(int i=0;i<3;++i){
		minim->Hesse();
	    minim->Minimize();
	}
    */
    for(int i=0;i<npars;++i)
        cout << par_names[i] << " " << likelihood_derivative(minim->X(), i) << endl;
	minim->Hesse();
	cout << "COV STATUS: " << minim->CovMatrixStatus() << endl;
	minim->PrintResults();

    ndf -= minim->NFree();
	cout << "CHI2/NDF: " << chi2(minim->X()) << "/" << ndf << endl;

	const double *x = minim->X();
	const double *x_err = minim->Errors();
	for(int i=0;i<npars;++i){
		results[site][i][0]	= x[i];
		results[site][i][1] = x_err[i];
	}
	daily_rates[site][0] = 0;
	daily_rates[site][1] = 0;

	int indices[n_range];
	for(int idx=0;idx<npars;++idx){
		if(strstr(par_names[idx], "n_li") != NULL){
			int r_idx = int(par_names[idx][5]) - 48;
			indices[r_idx] = idx;
			if(r_idx==n_range-1)
				daily_rates[site][0] += 0.211 * fabs(x[idx]);
			else
				daily_rates[site][0] += fabs(x[idx]);
		}
	}
	for(int i=0;i<n_range;++i){
		for(int j=0;j<n_range;++j){
			double _cov = minim->CovMatrix(indices[i], indices[j]);
			if(i == n_range-1)
				_cov *= 0.211;
			if(j == n_range-1)
				_cov *= 0.211;
			daily_rates[site][1] += _cov;
		}
	}
	daily_rates[site][1] = sqrt(daily_rates[site][1]);

	TH2D *h_cov = new TH2D("h_cov","h_cov",npars,0,npars,npars,0,npars);
	for(int i=0;i<npars;++i)
		for(int j=0;j<npars;++j)
			h_cov->SetBinContent(i,j,minim->CovMatrix(i,j));
	sprintf(buf, "h_cov_%d.png", site);	
	gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1");
	h_cov->Draw("colz");
	c1->SaveAs(buf);
    delete c1;
	gStyle->SetOptStat(1);
	plotHists(site, x);
	plotSlice(site, minim);
}

void plotSlice(int site, ROOT::Math::Minimizer *minim){
	const double *par = minim->X();
	
	int offset=0;
	for(;offset<npars_max;++offset)
		if(strstr(par_names[offset],"eps_s")!=NULL)	break;

	for(int type=0;type<slice_types;++type){
		if(!use_slice[type]) continue;
		double x[n_pdfs][slice_max], x_err[n_pdfs][slice_max];
		double y[n_pdfs][slice_max], y_err[n_pdfs][slice_max];
		double p[n_pdfs][slice_max], p_err[n_pdfs][slice_max];
		double y_sum[n_pdfs], y_sum_err[n_pdfs];

		for(int i=0;i<n_pdfs;++i){
			y_sum[i] = 0;
			y_sum_err[i] = 0;
		}
		double _start = slice_range[type][0];
		double _end = slice_range[type][1];
		double _step = (_end-_start) / slices[type];
		for(int j=0;j<n_pdfs;++j){
			//int _off = offset + (slices[type]-1) * j;
			int _off = offset + slices[type] * j;
			for(int s=0;s<slices[type];++s){
				x[j][s] = _start + _step * (s+0.5) + j * 0.05 * _step;
				y[j][s] = par[_off + s];
				y_err[j][s] = sqrt(minim->CovMatrix(_off + s, _off + s));

				y_sum[j] += y[j][s];

				for(int s1=0;s1<slices[type];++s1){
					y_sum_err[j] += minim->CovMatrix(_off + s, _off + s1);
                }
			}
		}
		
		for(int i=0;i<n_pdfs;++i){
			//int _off = offset + (slices[type]-1) * i;
			int _off = offset + slices[type] * i;
			y_sum_err[i] = sqrt(y_sum_err[i]);
			for(int j=0;j<slices[type];++j){
				p[i][j] = y[i][j] / y_sum[i];
				double _tmp_cov=0;
				for(int k=0;k<slices[type];++k)
					_tmp_cov += minim->CovMatrix(_off + k, _off + j);		
				p_err[i][j] = p[i][j] * sqrt( pow(y_err[i][j] / y[i][j], 2) + pow(y_sum_err[i] / y_sum[i], 2) - 2 * _tmp_cov / (y[i][j] * y_sum[i]) );
			}
		}

		TMultiGraph *mg = new TMultiGraph();
		TGraph *gr[n_pdfs];
        const char *titles[4] = {
            "Uncorr.",
            "Li9/He8",
            "B12",
            "N12"
        };
		for(int i=0;i<n_pdfs;++i){
            if(i>0 && !use_isotope[i-1]) continue;
			gr[i] = new TGraphErrors(slices[type], x[i], p[i], 0, p_err[i]);
			gr[i]->SetMarkerStyle(20);
			gr[i]->SetMarkerColor(colors[i]);
			gr[i]->SetFillStyle(0);
			gr[i]->SetFillColor(0);
            gr[i]->SetTitle(titles[i]);
			mg->Add(gr[i]);
		}
		char buf[255];
		sprintf(buf, "EH%d %s distribution", site+1, slice_vars[type]);	
		mg->SetTitle(buf);
		mg->SetMinimum(0);
		//mg->SetMaximum(0.55);
        TCanvas *c1 = new TCanvas("c1");
		mg->Draw("APC");
		c1->BuildLegend(0.7, 0.7, 0.9, 0.9);
		if(strstr(slice_vars[type],"ep") != NULL){
			TFile *file_spec = new TFile("./data/toyli9spec_BCWmodel_v1.root","READ");
			TH1 *h_spec = (TH1*)file_spec->Get("h_eVisAllSmeared");
			vector<double> bins;	
			bins.push_back(slice_range[type][0]);
			for(int i=0;i<slices[type];++i){
				double edge = (slice_range[type][1] - slice_range[type][0]) / slices[type] * (i+1) + slice_range[type][0];
				bins.push_back(edge);
			}
			TH1 *h_spec_rebin = h_spec->Rebin(bins.size()-1, "h_spec_rebin", &bins[0]);
			h_spec_rebin->SetLineColor(colors[1]);
			h_spec_rebin->SetLineStyle(9);
			h_spec_rebin->DrawNormalized("same");
		}
		const char *dir = "./plots/cfits_slice";
		sprintf(buf,"%s/cfit_%s_%d.png",dir,slice_vars[type],site);
		gPad->SetLogx(0);
		gPad->SetLogy(0);
		c1->SaveAs(buf);
        delete c1;
		//offset += n_pdfs * (slices[type] - 1);
		offset += n_pdfs * (slices[type]);
	}//type
}

#include<sstream>
char *dirname = "results/release_noeff/";

void parseFile(char *filename, double rates[3][2]){
    ifstream ifile;
    ifile.open(filename);
    if(ifile.is_open()){
        vector<string> lines;
        string line;
        
        while( getline(ifile, line) )
            lines.push_back(line);
        ifile.close();
        for(int site=0;site<3;++site){
            istringstream ss(lines[lines.size()-4+site]);
            string tmp;
            double mean, err;
            ss >> tmp >> mean >> err;
            rates[site][0] = mean;
            rates[site][1] = err;
        }
    }
}

double stddev(double *x, int n){
    double _x = 0;
    double _x2 = 0;
    for(int i=0;i<n;++i){
        _x += x[i];
        _x2 += x[i] * x[i];
    }
    return sqrt((_x2 - (_x * _x) / n) / (n - 1));
}
const int start = 4;
const int bins = 7;

void plotResults(){

    char buf[500];

    char *analysis[3] = {
        "nGd",
        "nH",
        "unified"
    };
    double x[bins], _x[bins];
    for(int i=0;i<bins;++i){
        x[i] = start+i;
        _x[i] = start+i+0.1;
    }
    for(int ana=2;ana<3;++ana){
        double rate[2][3][bins];
        double rate_err[2][3][bins];
        double rate_avg[2][3];
        double rate_err_avg[2][3];
        for(int i=0;i<2;++i){
            for(int j=0;j<3;++j){
                    rate_avg[i][j] = 0;
                    rate_err_avg[i][j] = 0;
                for(int k=0;k<bins;++k){
                    rate[i][j][k] = 0;
                    rate_err[i][j][k] = 0;
                }
            }
        }
        for(int i=0;i<bins;++i){
            sprintf(buf,"%s/%s/hists_4var_%dslice_fix/fit_result",dirname,analysis[ana],i+start);
            //cout << buf << endl;
            double rates[3][2];
            parseFile(buf, rates);
            for(int site=0;site<3;++site){
                rate[0][site][i] = rates[site][0];
                rate_err[0][site][i] = rates[site][1];
                rate_avg[0][site] += rates[site][0] / bins;
                rate_err_avg[0][site] += rates[site][1] / bins;
            }
            sprintf(buf,"%s/%s/hists_4var_%dslice_free/fit_result",dirname,analysis[ana],i+start);
            //cout << buf << endl;
            double rates[3][2];
            parseFile(buf, rates);
            for(int site=0;site<3;++site){
                rate[1][site][i] = rates[site][0];
                rate_err[1][site][i] = rates[site][1];
                rate_avg[1][site] += rates[site][0] / bins;
                rate_err_avg[1][site] += rates[site][1] / bins;
            }
            
        }
        for(int site=0;site<3;++site){
            TF1 *f_fix = new TF1("f_fix","pol0");
            f_fix->SetLineColor(kBlue);
            TGraphErrors *g_fix = new TGraphErrors(bins, x, rate[0][site], 0, rate_err[0][site]);
            g_fix->SetMarkerColor(kBlue);
            g_fix->SetMarkerStyle(20);
            g_fix->SetFillStyle(0);
            g_fix->SetTitle("Fixed");
            //g_fix->Fit(f_fix);
            TF1 *f_free = new TF1("f_free","pol0");
            f_free->SetLineColor(kRed);
            TGraphErrors *g_free = new TGraphErrors(bins, _x, rate[1][site], 0, rate_err[1][site]);
            g_free->SetMarkerColor(kRed);
            g_free->SetMarkerStyle(20);
            g_free->SetFillStyle(0);
            g_free->SetTitle("Free");
            //g_free->Fit(f_free);
            TMultiGraph *g = new TMultiGraph();
            g->Add(g_fix);
            g->Add(g_free);
            sprintf(buf,"EH%d fit results;IBD slices;Li9/He8 per day",site+1);
            g->SetTitle(buf);            
            g->Draw("AP");
        
            c1->Update();
            TLine *lfix = new TLine(c1->GetUxmin(), rate_avg[0][site], c1->GetUxmax(), rate_avg[0][site]);
            lfix->SetLineColor(kBlue);
            lfix->SetLineWidth(2);
            lfix->Draw("same");
            TLine *lfree = new TLine(c1->GetUxmin(), rate_avg[1][site], c1->GetUxmax(), rate_avg[1][site]);
            lfree->SetLineColor(kRed);
            lfree->SetLineWidth(2);
            lfree->Draw("same");

            TLegend *leg = new TLegend(0.7,0.8,0.9,0.9);
            leg->AddEntry(g_fix);
            leg->AddEntry(g_free);
            leg->Draw("same");
            c1->Update();
            sprintf(buf,"./plots/cfits_result/EH%d_%s_results.png", site+1, analysis[ana]);
            c1->SaveAs(buf);
            sprintf(buf,"./plots/cfits_result/EH%d_%s_results.pdf", site+1, analysis[ana]);
            c1->SaveAs(buf);

            printf("EH%d fixed: %.3f %.3f(%.1f) %.3f(%.1f) %.3f(%.1f) %.3f\n", 
                site+1, 
                rate_avg[0][site], 
                rate_err_avg[0][site], 
                rate_err_avg[0][site] * 100 / rate_avg[0][site], 
                stddev(rate[0][site], bins), 
                stddev(rate[0][site], bins) * 100 / rate_avg[0][site],
                fabs(rate_avg[0][site]-rate_avg[1][site]),
                fabs(rate_avg[0][site]-rate_avg[1][site]) * 100 / rate_avg[0][site],
                sqrt(pow(stddev(rate[0][site], bins), 2.0) + pow(fabs(rate_avg[0][site] - rate_avg[1][site]), 2.0)));
     //           sqrt(pow(rate_err_avg[0][site], 2.0) + pow(stddev(rate[0][site], bins), 2.0) + pow(fabs(rate_avg[0][site] - rate_avg[1][site]), 2.0)));
            printf("EH%d free: %.3f %.3f(%.1f) %.3f(%.1f) %.3f(%.1f) %.3f\n", 
                site+1, 
                rate_avg[1][site], 
                rate_err_avg[1][site], 
                rate_err_avg[1][site] * 100 / rate_avg[1][site],
                stddev(rate[1][site], bins), 
                stddev(rate[1][site], bins) * 100 / rate_avg[1][site], 
                fabs(rate_avg[0][site]-rate_avg[1][site]),
                fabs(rate_avg[0][site]-rate_avg[1][site]) * 100 / rate_avg[1][site],
                sqrt(pow(stddev(rate[1][site], bins), 2.0) + pow(fabs(rate_avg[0][site] - rate_avg[1][site]), 2.0)));
   //             sqrt(pow(rate_err_avg[1][site], 2.0) + pow(stddev(rate[1][site], bins), 2.0) + pow(fabs(rate_avg[0][site] - rate_avg[1][site]), 2.0)));
            cout << endl;
        }
    }
    
}

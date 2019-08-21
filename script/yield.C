double dc_nd[2] = {5.51 ,0.52};
double dc_nd_emu[2] = {32.1, 2.0};
double dc_fd[2] = {7.90 ,0.51};
double dc_fd_emu[2] = {63.7, 5.5};

double kamland[2] = {22, 2};
double kamland_emu[2] = {260, 8};

double boxerino[2] = {29, 3};
double boxerino_emu[2] = {283, 19};

double p_mu[3] = {.6236, .6240, .6242};
double l_avg[3] = {204.1, 204.5, 204.9};
double delta_pl[3] = {5.53/100, 5.58/100, 3.80/100};
double rho = 0.86;
double li9_branch = 0.50;

double dyb_emu[3][2] = {
    {63.9, 3.8},
    {64.7, 3.9},
    {143.0, 8.6}
};

double dyb_rmu[3][2] = {
    {20.54, 8.17e-3},
    {15.39, 5.94e-3},
    {1.06, 1.23e-3}
};
double dyb_li9[3][2] = {
    {5.87, 1.04},
    {5.09, 0.85},
    {0.60, 0.05} 
};

double ed_eff[2] = {0.9271, 0.97/100};
/*
EH1 fixed: 12.136 1.724(14.2) 0.549(4.5) 0.431(3.6) 0.698
EH1 free: 12.567 1.815(14.4) 0.486(3.9) 0.431(3.4) 0.650

EH2 fixed: 8.972 1.209(13.5) 0.328(3.7) 0.656(7.3) 0.733
EH2 free: 9.627 1.317(13.7) 0.327(3.4) 0.656(6.8) 0.733

EH3 fixed: 0.945 0.054(5.7) 0.012(1.3) 0.035(3.7) 0.037
EH3 free: 0.910 0.057(6.2) 0.011(1.3) 0.035(3.9) 0.037
*/
double dyb_li9_unified[3][2] = {
    {12.567, 1.927},
    {9.627, 1.507},
    {0.910, 0.005}
};
double l_avg_unified[2] = {258, 25.8};
double gd_cap = 0.854;
double cut_eff_unified = 0.89 * 0.78; //eff dt, eff dist
void yield(){

    TMultiGraph *mg = new TMultiGraph();
    TLegend *leg = new TLegend(0.6, 0.15, 0.88, 0.6);

    TGraphErrors *ge_dc_nd = new TGraphErrors(1, dc_nd_emu, dc_nd, (dc_nd_emu+1), (dc_nd+1));
    ge_dc_nd->SetMarkerStyle(20);
    ge_dc_nd->SetFillStyle(0);
    ge_dc_nd->SetTitle("DC Near");
    TGraphErrors *ge_dc_fd = new TGraphErrors(1, dc_fd_emu, dc_fd, (dc_fd_emu+1), (dc_fd+1));
    ge_dc_fd->SetMarkerStyle(21);
    ge_dc_fd->SetFillStyle(0);
    ge_dc_fd->SetTitle("DC Far");

    TGraphErrors *ge_kamland = new TGraphErrors(1, kamland_emu, kamland, (kamland_emu+1), (kamland+1));
    ge_kamland->SetMarkerStyle(22);
    ge_kamland->SetFillStyle(0);
    ge_kamland->SetTitle("KamLAND");
    TGraphErrors *ge_boxerino = new TGraphErrors(1, boxerino_emu, boxerino, (boxerino_emu+1), (boxerino+1));
    ge_boxerino->SetMarkerStyle(23);
    ge_boxerino->SetFillStyle(0);
    ge_boxerino->SetTitle("Borexino");


    double yields[3];
    double yields_err[3];
    double x[3];
    double x_err[3];
    TGraphErrors *ge_dyb[3];
    for(int i=0;i<3;++i){

        yields[i] = (dyb_li9[i][0] / 86400 / li9_branch / ed_eff[0]) / (dyb_rmu[i][0] * p_mu[i] * l_avg[i] * rho * gd_cap) * 1e8; 
        yields_err[i] = yields[i] * sqrt( pow(dyb_li9[i][1] / dyb_li9[i][0], 2.) + pow(dyb_rmu[i][1] / dyb_rmu[i][0], 2.) + pow(delta_pl[i], 2.) + pow(ed_eff[1], 2.) ); 
        cout << yields[i] << " " << yields_err[i] << endl;
        x[i] = dyb_emu[i][0];
        x_err[i] = dyb_emu[i][1];
        ge_dyb[i] = new TGraphErrors(1, x+i, yields+i, x_err+i, yields_err+i);
        ge_dyb[i]->SetMarkerStyle(24+i);
        ge_dyb[i]->SetFillStyle(0);
        ge_dyb[i]->SetTitle(TString::Format("DYB EH%d",i+1));
        mg->Add(ge_dyb[i]);
        leg->AddEntry(ge_dyb[i], "", "p");
    }
    /*
    for(int i=0;i<3;++i){

        yields[i] = (dyb_li9_unified[i][0] / 86400 / li9_branch / cut_eff_unified) / (dyb_rmu[i][0] * l_avg_unified[0] * rho) * 1e8; 
        yields_err[i] = yields[i] * sqrt( pow(dyb_li9_unified[i][1] / dyb_li9_unified[i][0], 2.) + pow(dyb_rmu[i][1] / dyb_rmu[i][0], 2.) + pow(0.1, 2.) ); 
        cout << yields[i] << " " << yields_err[i] << endl;
        ge_dyb[i] = new TGraphErrors(1, x+i, yields+i, x_err+i, yields_err+i);
        ge_dyb[i]->SetMarkerStyle(35+i);
        ge_dyb[i]->SetFillStyle(0);
        ge_dyb[i]->SetTitle(TString::Format("DYB (Uni.) EH%d",i+1));
        mg->Add(ge_dyb[i]);
        leg->AddEntry(ge_dyb[i], "", "p");
    }
    */


    TF1 *power = new TF1("power", "[0] * x ^ [1]", 0, 1000);
    power->SetParameter(0, 1);    
    power->SetParameter(1, 1);    

    mg->Add(ge_dc_nd);
    mg->Add(ge_dc_fd);
    mg->Add(ge_kamland);
    mg->Add(ge_boxerino);
    mg->Fit(power);
    mg->Draw("AP");
    mg->GetHistogram()->SetTitle("Lithium-9 Yield");
    mg->GetXaxis()->SetTitle("Mean Muon Energy [GeV]");
//    mg->GetYaxis()->SetTitle("\\text{Yield} [10^{-8}\\mu^{-1}\\text{g}^{-1}\\text{cm}^2]");
    mg->GetYaxis()->SetTitle("Yield [10^{-8}#mu^{-1}g^{-1}cm^{2}]");
    mg->GetYaxis()->SetMoreLogLabels();
    gPad->SetLogy();
    power->Draw("same");
    leg->AddEntry(ge_dc_nd, "", "p");
    leg->AddEntry(ge_dc_fd, "", "p");
    leg->AddEntry(ge_kamland, "", "p");
    leg->AddEntry(ge_boxerino, "", "p");
    leg->AddEntry(power, "Best Fit", "l");
    leg->SetBorderSize(0);
    leg->Draw("same");
    c1->SaveAs("yield.pdf");



//B12
    mg = new TMultiGraph();
    leg = new TLegend(0.6, 0.15, 0.88, 0.5);
    leg->SetBorderSize(0);

    kamland[0] = 42.9;
    kamland[1] = 3.3;
    ge_kamland = new TGraphErrors(1, kamland_emu, kamland, (kamland_emu+1), (kamland+1));
    ge_kamland->SetMarkerStyle(22);
    ge_kamland->SetFillStyle(0);
    ge_kamland->SetTitle("KamLAND");

    boxerino[0] = 56;
    boxerino[1] = 3;
    ge_boxerino = new TGraphErrors(1, boxerino_emu, boxerino, (boxerino_emu+1), (boxerino+1));
    ge_boxerino->SetMarkerStyle(23);
    ge_boxerino->SetFillStyle(0);
    ge_boxerino->SetTitle("Borexino");



    double rates_b12[3][2] = {
        { 4.24, 0.20 },
        { 3.54, 0.24 },
        { 0.39, 0.03 }
    };

    double rp[3] = {57, 54, 53};
    double eps_mult = 0.97;
    
    for(int i=0;i<3;++i){
        double _poisson_p = exp(-rp[i] * 200e-6) * rp[i] * 200e-6;
        yields[i] = (rates_b12[i][0] / 86400) / (eps_mult * dyb_rmu[i][0] * 0.54 * rho * l_avg_unified[0] * _poisson_p ) * 1e7; 
        yields_err[i] = yields[i] * sqrt( pow(rates_b12[i][1] / rates_b12[i][0], 2.) + pow(dyb_rmu[i][1] / dyb_rmu[i][0], 2.) + 0.01); 
        cout << yields[i] << " " << yields_err[i] << endl;
        ge_dyb[i] = new TGraphErrors(1, x+i, yields+i, x_err+i, yields_err+i);
        ge_dyb[i]->SetMarkerStyle(24+i);
        ge_dyb[i]->SetFillStyle(0);
        ge_dyb[i]->SetTitle(TString::Format("DYB EH%d",i+1));
        mg->Add(ge_dyb[i]);
        leg->AddEntry(ge_dyb[i], "", "p");
    }

    power->SetParameter(0, 1);    
    power->SetParameter(1, 1);    
    mg->Add(ge_kamland);
    mg->Add(ge_boxerino);
    leg->AddEntry(ge_kamland, "", "p");
    leg->AddEntry(ge_boxerino, "", "p");
    mg->Fit(power);
    mg->Draw("AP");
    power->Draw("same");
    leg->Draw("same");
    mg->GetHistogram()->SetTitle("Boron-12 Yield");
    mg->GetXaxis()->SetTitle("Mean Muon Energy [GeV]");
    mg->GetYaxis()->SetTitle("Yield [10^{-7}#mu^{-1}g^{-1}cm^{2}]");
    mg->GetYaxis()->SetMoreLogLabels();
    c1->SaveAs("yield_b12.pdf");
}

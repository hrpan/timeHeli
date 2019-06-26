/******************************************************************************
*Tree    :Heli      :                                                        *
*Entries :    96041 : Total =       494629221 bytes  File  Size =  157132299 *
*        :          : Tree compression factor =   3.15                       *
******************************************************************************
*Br    0 :site      : site/I                                                 *
*Entries :    96041 : Total  Size=     385546 bytes  File Size  =       3147 *
*Baskets :       12 : Basket Size=      52224 bytes  Compression= 122.34     *
*............................................................................*
*Br    1 :detector  : detector/I                                             *
*Entries :    96041 : Total  Size=     385610 bytes  File Size  =      40150 *
*Baskets :       12 : Basket Size=      52224 bytes  Compression=   9.59     *
*............................................................................*
*Br    2 :date      : date/I                                                 *
*Entries :    96041 : Total  Size=     385546 bytes  File Size  =       3169 *
*Baskets :       12 : Basket Size=      52224 bytes  Compression= 121.49     *
*............................................................................*
*Br    3 :ep        : ep/F                                                   *
*Entries :    96041 : Total  Size=     385514 bytes  File Size  =     339498 *
*Baskets :       12 : Basket Size=      52224 bytes  Compression=   1.13     *
*............................................................................*
*Br    4 :ed        : ed/F                                                   *
*Entries :    96041 : Total  Size=     385514 bytes  File Size  =     338746 *
*Baskets :       12 : Basket Size=      52224 bytes  Compression=   1.14     *
*............................................................................*
*Br    5 :dt        : dt/F                                                   *
*Entries :    96041 : Total  Size=     385514 bytes  File Size  =     304999 *
*Baskets :       12 : Basket Size=      52224 bytes  Compression=   1.26     *
*............................................................................*
*Br    6 :xp        : xp/F                                                   *
*Entries :    96041 : Total  Size=     385514 bytes  File Size  =     352064 *
*Baskets :       12 : Basket Size=      52224 bytes  Compression=   1.09     *
*............................................................................*
*Br    7 :yp        : yp/F                                                   *
*Entries :    96041 : Total  Size=     385514 bytes  File Size  =     351844 *
*Baskets :       12 : Basket Size=      52224 bytes  Compression=   1.09     *
*............................................................................*
*Br    8 :zp        : zp/F                                                   *
*Entries :    96041 : Total  Size=     385514 bytes  File Size  =     353744 *
*Baskets :       12 : Basket Size=      52224 bytes  Compression=   1.09     *
*............................................................................*
*Br    9 :xd        : xd/F                                                   *
*Entries :    96041 : Total  Size=     385514 bytes  File Size  =     351522 *
*Baskets :       12 : Basket Size=      52224 bytes  Compression=   1.10     *
*............................................................................*
*Br   10 :yd        : yd/F                                                   *
*Entries :    96041 : Total  Size=     385514 bytes  File Size  =     351325 *
*Baskets :       12 : Basket Size=      52224 bytes  Compression=   1.10     *
*............................................................................*
*Br   11 :zd        : zd/F                                                   *
*Entries :    96041 : Total  Size=     385514 bytes  File Size  =     353342 *
*Baskets :       12 : Basket Size=      52224 bytes  Compression=   1.09     *
*............................................................................*
*Br   12 :dist      : dist/F                                                 *
*Entries :    96041 : Total  Size=     385546 bytes  File Size  =     342164 *
*Baskets :       12 : Basket Size=      52224 bytes  Compression=   1.13     *
*............................................................................*
*Br   13 :dtlSH     : vector<double>                                         *
*Entries :    96041 : Total  Size=   84154612 bytes  File Size  =   65705246 *
*Baskets :      495 : Basket Size=   25600000 bytes  Compression=   1.28     *
*............................................................................*
*Br   14 :nPESum    : vector<double>                                         *
*Entries :    96041 : Total  Size=   84155111 bytes  File Size  =   23741172 *
*Baskets :      495 : Basket Size=   25600000 bytes  Compression=   3.54     */

#include "TH2D.h"
#include "TChain.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include<iostream>
char *gcut = "dt<200&&dt>1&&ed>6";
char *hcut = "dt<400&&dist<500&&ed>1.9&&ed<2.7&&ep>1.5";
char *ucut = "dt<400&&dist<500&&ed>1.5";

int bins1d = 1000;
int bins2d = 100;
void plotIBD(){
    gStyle->SetOptStat(0);
    TChain *chain = new TChain("Heli");
    chain->Add("../p17b/data_heli/*.root");
    //chain->Add("../p17b/data_heli/21221.root");

    float _ep, _ed, _dt, _dist, _xp, _xd, _yp, _yd, _zp, _zd;
    int _site, _det;

    chain->SetBranchAddress("ep",&_ep);
    chain->SetBranchAddress("ed",&_ed);
    chain->SetBranchAddress("dt",&_dt);
    chain->SetBranchAddress("dist",&_dist);
    chain->SetBranchAddress("xp",&_xp);
    chain->SetBranchAddress("xd",&_xd);
    chain->SetBranchAddress("yp",&_yp);
    chain->SetBranchAddress("yd",&_yd);
    chain->SetBranchAddress("zp",&_zp);
    chain->SetBranchAddress("zd",&_zd);
    chain->SetBranchAddress("site",&_site);
    chain->SetBranchAddress("detector",&_det);

    TFile *hists = new TFile("data/hists.root","RECREATE");
    hists->cd();
    //ep, ed, dt, dist
    TH1 *h_var1d[4][4][3];

    char buf[500];

    char *var1d_unit[4] = {
        "[MeV]",
        "[MeV]",
        "[us]",
        "[mm]"
    };

    char *var1d_name[4] = {
        "prompt energy",
        "delayed energy",
        "capture time",
        "p-d distance",
    };

    double var1d_bounds[4][2] = {
        {0.7, 12},
        {1.5, 12},
        {1, 400},
        {0, 5000}
    };

    for(int var=0; var<4; ++var){
        for(int ana=0; ana<4; ++ana){
            for(int site=0; site<3; ++site){
                sprintf(buf, "h_var1d_%d_%d_%d", var, ana, site);                
                char title[500];
                sprintf(title, "EH%d %s;%s", site+1, var1d_name[var], var1d_unit[var]);
                h_var1d[var][ana][site] = new TH1D(buf, title, bins1d, var1d_bounds[var][0], var1d_bounds[var][1]);
            }
        }
    }

    //r2z, ep-d
    TH2 *h_var2d[2][4][3];
   
    char *var2d_unit[2][2] = {
        {"[mm^2]", "[mm]"},
        {"Ed[MeV]", "Ep[MeV]"}
    };

    char *var2d_name[4] = {
        "R2-Z",
        "prompt energy vs. delayed energy"
    };

    double var2d_bounds[2][2][2] = {
        {
            {0, 5000000},
            {-3000, 3000}
        },
        {
            {1.5, 12},
            {0.7, 12},
        }
    };

    for(int var=0; var<2; ++var){
        for(int ana=0; ana<4; ++ana){
            for(int site=0; site<3; ++site){
                sprintf(buf, "h_var2d_%d_%d_%d", var, ana, site);                
                char title[500];
                sprintf(title, "EH%d %s;%s;%s", site+1, var2d_name[var], var2d_unit[var][0], var2d_unit[var][1]);
                h_var2d[var][ana][site] = new TH2D(buf, title, bins2d, var2d_bounds[var][0][0], var2d_bounds[var][0][1], bins2d, var2d_bounds[var][1][0], var2d_bounds[var][1][1]);
            }
        }
    }
    
    size_t events = chain->GetEntries();
    for(size_t i=0; i<events; ++i){

        if(i % (events/10000) == 0){
            printf("%.2f\%\n", 100.0 * i / events);
            fflush(stdout);
        }

        chain->GetEntry(i);

        int site = _site==4? 2:_site-1;
        bool _gd = (_ep>0.7 && _ep<12) && (_ed>6 && _ed<12) && (_dt>1 && _dt<200);
        bool _h = (_ep>1.5 && _ep<12) && (_ed>1.9 && _ed<2.7) && (_dist < 500) && (_dt>1 && _dt<400);
        bool _uni = (_dist < 500);

        bool ana_check[4] = {true, _gd, _h, _uni};

        for(int ana=0; ana<4; ++ana){
            if(!ana_check[ana]) continue;
            h_var1d[0][ana][site]->Fill(_ep);
            h_var1d[1][ana][site]->Fill(_ed);
            h_var1d[2][ana][site]->Fill(_dt);
            h_var1d[3][ana][site]->Fill(_dist);
            h_var2d[0][ana][site]->Fill((_xp*_xp+_yp*_yp), _zp);
            h_var2d[1][ana][site]->Fill(_ed, _ep);
        } 
    }
    hists->Write();
    hists->Close();
}

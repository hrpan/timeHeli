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
*Baskets :      495 : Basket Size=   25600000 bytes  Compression=   3.54     *
*............................................................................*
*Br   15 :nTag_n    : vector<int>                                            *
*Entries :    96041 : Total  Size=   42750953 bytes  File Size  =    6121397 *
*Baskets :      257 : Basket Size=    5808128 bytes  Compression=   6.98     *
*............................................................................*
*Br   16 :nTag_e    : vector<float>                                          *
*Entries :    96041 : Total  Size=   34916352 bytes  File Size  =   14294807 *
*Baskets :      212 : Basket Size=    4742144 bytes  Compression=   2.44     *
*............................................................................*
*Br   17 :nTag_dt   : vector<float>                                          *
*Entries :    96041 : Total  Size=   34916568 bytes  File Size  =    8948382 *
*Baskets :      212 : Basket Size=    4742144 bytes  Compression=   3.90     *
*............................................................................*
*Br   18 :wpTag_n   : vector<int>                                            *
*Entries :    96041 : Total  Size=   42751214 bytes  File Size  =    2540909 *
*Baskets :      257 : Basket Size=    5808640 bytes  Compression=  16.82     *
*............................................................................*
*Br   19 :wpTag_nHit : vector<int>                                           *
*Entries :    96041 : Total  Size=   82986344 bytes  File Size  =   19591522 *
*Baskets :      488 : Basket Size=   25600000 bytes  Compression=   4.24     *
*............................................................................*
*Br   20 :wpTag_dt  : vector<float>                                          *
*Entries :    96041 : Total  Size=   82985360 bytes  File Size  =   12674815 *
*Baskets :      488 : Basket Size=   25600000 bytes  Compression=   6.55     *
*............................................................................*/

void splitIBD(){
    TChain *chain = new TChain("Heli");
    //chain->Add("../p17b/data_heli/*.root");
    chain->Add("../p17b/data_heli/2*.root");

    float ep, ed, dt, dist, xp, xd, yp, yd, zp, zd;
    int det, site;
    chain->SetBranchAddress("ep", &ep);
    chain->SetBranchAddress("ed", &ed);
    chain->SetBranchAddress("dt", &dt);
    chain->SetBranchAddress("dist", &dist);
    chain->SetBranchAddress("xp", &xp);
    chain->SetBranchAddress("xd", &xd);
    chain->SetBranchAddress("yp", &yp);
    chain->SetBranchAddress("yd", &yd);
    chain->SetBranchAddress("zp", &zp);
    chain->SetBranchAddress("zd", &zd);
    chain->SetBranchAddress("detector", &det);
    chain->SetBranchAddress("site", &site);
    
    char buf[500];

    TFile *f_out[3];
    TTree *tr[3];
    for(int s=0;s<3;++s){
        sprintf(buf, "./data/ibd%d.root", s);
        f_out[s] = new TFile(buf, "RECREATE");
        f_out[s]->cd();
        tr[s] = new TTree("tr","tr");
        tr[s]->Branch("ep", &ep);
        tr[s]->Branch("ed", &ed);
        tr[s]->Branch("dt", &dt);
        tr[s]->Branch("dist", &dist);
        tr[s]->Branch("xp", &xp);
        tr[s]->Branch("xd", &xd);
        tr[s]->Branch("yp", &yp);
        tr[s]->Branch("yd", &yd);
        tr[s]->Branch("zp", &zp);
        tr[s]->Branch("zd", &zd);
        tr[s]->Branch("det", &det);
    }

    size_t entries = chain->GetEntries();
    for(size_t i=0;i<entries;++i){
        if(i% (entries/10000) == 0){
            printf("%.3f\%\n", i * 100.0 / entries);
            fflush(stdout);
        }
        chain->GetEntry(i);

        int s = site==4? 2 : site-1; 

        tr[s]->Fill(); 

    }
}

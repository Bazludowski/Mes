#pragma once
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
struct node {
    double x = 0, y = 0, t0 = 0, t1 =0, BC = 0;
    node()
    {
        x = 0.0;
        y = 0.0;
        t0 = 0.0;
        t1 = 0.0;
        BC = 0.0;
    }
};

struct element {
    int id[4]; //kazdy element ma 4 id dla 4 wierzcholkow
    double HL[4][4] = {0.0};
    double CL[4][4] = {0.0};
    double HBCL[4][4] = {0.0};
    double PL[4] = {0.0};
};

struct GlobalData {
    double w = 0.1;
    double h = 0.1;
    int n_w = 31;
    int n_h = 31;
    int n_n = n_w*n_h;
    int n_e = (n_w - 1)*(n_h - 1);
    int ro = 7800;
    int c = 700;
    double st = 100;
    double dt = 1;
    double Talfa = 1200;
    double **HG;
    double **CG;
    double *PG;
    double **HGZ;
    double *PGZ;
    double k = 25.0;
    double alfa = 300;
    double *t0;
    GlobalData (){
        HG = new double *[n_n];
        CG = new double *[n_n];
        PG = new double [n_n];
        HGZ = new double *[n_n];
        PGZ = new double [n_n];
        t0 = new double [n_n];
        for (int i = 0; i < n_n; i++){
            HG[i] = new double [n_n];
            CG[i] = new double [n_n];
            HGZ[i] = new double [n_n];
            for(int j = 0; j < n_n; j++){
                HG[i][j] = 0.0;
                CG[i][j] = 0.0;
                HGZ[i][j] = 0.0;
            }
            PG[i] = 0.0;
            PGZ[i] = 0.0;
            t0[i] = 100.0;
        }
    }
};
struct elem4 {
    double ksi2[4]={0.0};
    double eta2[4]={0.0};
    double ksi3[9]={0.0};
    double eta3[9]={0.0};
    double ksi4[16]={0.0};
    double eta4[16]={0.0};
    double pc1 = 1.0 / sqrt(3);//1 punkt calkowania
    double pc2 = sqrt(3)/sqrt(5);
    double pc41 = 0.339981;
    double pc42 = 0.861136;
    double w41 = 0.652145; // z plusem
    double w42 = 0.347854; // z minusem
    double waga[4] = {0.0};
    elem4()
    {
        ksi2[0] = -pc1;
        ksi2[1] = pc1;
        ksi2[2] = pc1;
        ksi2[3] = -pc1;

        eta2[0] = -pc1;
        eta2[1] = -pc1;
        eta2[2] = pc1;
        eta2[3] = pc1;

        ksi3[0] = -pc2;
        ksi3[1] = 0.0;
        ksi3[2] = pc2;
        ksi3[3] = -pc2;
        ksi3[4] = 0.0;
        ksi3[5] = pc2;
        ksi3[6] = -pc2;
        ksi3[7] = 0;
        ksi3[8] = pc2;

        eta3[0] = -pc2;
        eta3[1] = -pc2;
        eta3[2] = -pc2;
        eta3[3] = 0.0;
        eta3[4] = 0.0;
        eta3[5] = 0.0;
        eta3[6] = pc2;
        eta3[7] = pc2;
        eta3[8] = pc2;
//        waga[0] = w42;
//        waga[1] = w41;
//        waga[2] = w41;
//        waga[3] = w42;

        ksi4[0] = -pc42;
        ksi4[1] = -pc41;
        ksi4[2] = pc41;
        ksi4[3] = pc42;
        ksi4[4] = -pc42;
        ksi4[5] = -pc41;
        ksi4[6] = pc41;
        ksi4[7] = pc42;
        ksi4[8] = -pc42;
        ksi4[9] = -pc41;
        ksi4[10] = pc41;
        ksi4[11] = pc42;
        ksi4[12] = -pc42;
        ksi4[13] = -pc41;
        ksi4[14] = pc41;
        ksi4[15] = pc42;

        eta4[0] = -pc42;
        eta4[1] = -pc42;
        eta4[2] = -pc42;
        eta4[3] = -pc42;
        eta4[4] = -pc41;
        eta4[5] = -pc41;
        eta4[6] = -pc41;
        eta4[7] = -pc41;
        eta4[8] = pc41;
        eta4[9] = pc41;
        eta4[10] = pc41;
        eta4[11] = pc41;
        eta4[12] = pc42;
        eta4[13] = pc42;
        eta4[14] = pc42;
        eta4[15] = pc42;
    }

};

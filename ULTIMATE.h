#pragma once
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include "struct.h"
#include <iomanip>
using namespace std;

void HBC(element *elementy, node *wierzcholki, int bok, int flaga){
    double pitagoras = pow(wierzcholki[elementy->id[bok]-1].x-wierzcholki[elementy->id[(bok+1)%4]-1].x,2.0)+pow(wierzcholki[elementy->id[bok]-1].y-wierzcholki[elementy->id[(bok+1)%4]-1].y,2.0);
//zmienna bok mowi dla ktorego boku kwadratu bedziemy liczyc, pitagoras
    double detJ = sqrt(pitagoras)/2.0; //dzielimy przez dwa bo w ukladzie lokalnym (kwadracie) bok ma 2 //DetJ to stosunek tego w rzeczywistosci do tego w ukladzie lokalnym
    double n[4][4] = {0.0};
    elem4 el;
    GlobalData data;
    int flagahbc = flaga*flaga;
    double ksix, etax = 0.0;
    switch(flaga){
        case 2:{
            el.waga[0] = 1.0;
            el.waga[1] = 1.0;
            for(int i = 0; i <2; i++){
                if(bok==0){
                    etax = -1.0;
                    ksix = el.ksi2[i];
                }
                else if (bok == 1){
                    etax = el.eta2[i+bok];
                    ksix = 1.0;
                }
                else if (bok ==2){
                    etax = 1.0;
                    ksix = el.ksi2[i+bok];
                }
                else if (bok == 3){
                    etax = el.eta2[(i+bok)%4];
                    ksix = -1.0;
                }
                n[i][0] = (1.0-ksix)*(1.0 - etax)/4.0;
                n[i][1] = (1.0+ksix)*(1.0 - etax)/4.0;
                n[i][2] = (1.0+ksix)*(1.0 + etax)/4.0;
                n[i][3] = (1.0-ksix)*(1.0 + etax)/4.0;

            }
//        for(int ksi2 = 0; ksi2 < 2; ksi2++){
//            for (int eta2 = 0; eta2<4; eta2++){
//                cout<<n[ksi2][eta2]<< " ";
//            }
//            cout<<endl;
//        }
//        cout <<endl;
            break;
        }
        case 3:{
            el.waga[0] = 5.0 / 9.0;
            el.waga[1] = 8.0 / 9.0;
            el.waga[2] = 5.0 / 9.0;
            for(int i = 0; i <3; i++){
                if(bok==0){
                    etax = -1.0;
                    ksix = el.ksi3[i];
                }
                else if (bok == 1){
                    etax = el.ksi3[i];
                    ksix = 1.0;
                }
                else if (bok ==2){
                    etax = 1.0;
                    ksix = el.ksi3[flaga-1-i]; //jak masz kwadrat i masz dolna podstawe, dajesz ksi[i]
                }
                else if (bok == 3){
                    etax = el.ksi3[flaga-1-i];
                    ksix = -1.0;
                }
                n[i][0] = (1.0-ksix)*(1.0 - etax)/4.0;
                n[i][1] = (1.0+ksix)*(1.0 - etax)/4.0;
                n[i][2] = (1.0+ksix)*(1.0 + etax)/4.0;
                n[i][3] = (1.0-ksix)*(1.0 + etax)/4.0;
            }
            break;
        }
        case 4:{
            el.waga[0] = 0.347854;
            el.waga[1] = 0.652145;
            el.waga[2] = 0.652145;
            el.waga[3] = 0.347854;
            for(int i = 0; i <4; i++){

                if(bok==0){
                    etax = -1.0;
                    ksix = el.ksi4[i];
                }
                else if (bok == 1){
                    etax = el.ksi4[i];
                    ksix = 1.0;
                }
                else if (bok ==2){
                    etax = 1.0;
                    ksix = el.ksi4[flaga-1-i];
                }
                else if (bok == 3){
                    etax = el.ksi4[flaga-1-i];
                    ksix = -1.0;
                }
                n[i][0] = (1.0-ksix)*(1.0 - etax)/4.0;
                n[i][1] = (1.0+ksix)*(1.0 - etax)/4.0;
                n[i][2] = (1.0+ksix)*(1.0 + etax)/4.0;
                n[i][3] = (1.0-ksix)*(1.0 + etax)/4.0;
            }
            break;
        }

    }
//cout << detJ<<endl;
    for(int hbc = 0; hbc < flaga; hbc++){
        for (int hbc2 = 0; hbc2 < 4; hbc2++) {
            for (int hbc3 = 0; hbc3 < 4; hbc3++) {
                elementy->HBCL[hbc2][hbc3] += detJ*data.alfa*n[hbc][hbc2]*n[hbc][hbc3]*el.waga[hbc];
//            cout<<detJ*data.alfa*n[hbc][hbc2]*n[hbc][hbc3]<<" ";
            }
            elementy->PL[hbc2]+=-data.alfa*detJ*n[hbc][hbc2]*el.waga[hbc]*data.Talfa;
        }
//    cout<<endl;
    }
}

double* gauss (double **HZast, double *PZast, int size) {

    int i,j,k,pivot;

    double *x=new double [size];
    double **a= new double *[size];
    for(i=0; i<size; i++) {
        a[i] = new double[size+1];
        x[i] = 0.0;
        for(j=0; j<size; j++) {
            a[i][j] = HZast[i][j];
        }
        a[i][size] = PZast[i];
    }
    for(i=0; i<size; i++) {
        for(pivot=i+1; pivot<size; pivot++) {
            if(a[i][i]<a[pivot][i]) {
                for(k=0; k<=size; k++) {
                    double temp = a[pivot][k];
                    a[pivot][k] = a[i][k];
                    a[i][k] = temp;
                }
            }
        }
    }
    for (i = 0; i<size-1; i++){
        for(pivot = i+1; pivot < size; pivot++){
            double temp = a[pivot][i]/a[i][i];
            for(k = 0; k<=size; k++){
                a[pivot][k] -= temp*a[i][k];
            }
        }
    }
    for(i=size-1; i>=0; i--) {
        x[i]=a[i][size];
        for(k=0; k < size; k++){
            if(k!=i){
                x[i]-=a[i][k]*x[k];
            }
        }
        x[i] = x[i]/a[i][i];
    }


    return x;
}

void Macierz(element *elementy, node *wierzcholki, GlobalData*dane) {
    int flaga = 3;
    int flag = flaga * flaga;
    GlobalData data = *dane;
    element ele;
    for (int elementid = 0; elementid < data.n_e; elementid++) {
        elem4 el;
        double x[4];
        double y[4];
        for (int j = 0; j < 4; j++) {
            int w = elementy[elementid].id[j] - 1;
            x[j] = wierzcholki[w].x;
            y[j] = wierzcholki[w].y;
        }
        double N_ksi[16][4] = {0.0};
        double N_eta[16][4] = {0.0};
        double dx_dksi[16] = {0.0};
        double dy_dksi[16] = {0.0};
        double dx_deta[16] = {0.0};
        double dy_deta[16] = {0.0};
        double Ni[16][4] = {0.0};
        double NTransponowane[16][4][4] = {0.0};
        double J[16][4][4] = {0.0};
        double det_J[16] = {0.0};
        for(int bok = 0; bok < 4; bok++)
        {
            if(wierzcholki[elementy[elementid].id[bok]-1].BC==1)  {
//                cout<<"Pierwszy if: " << elementy[elementid].id[bok]-1 << " " << wierzcholki[elementy[elementid].id[bok]-1].x<< " " << wierzcholki[elementy[elementid].id[bok]-1].BC<<endl;
                if(wierzcholki[elementy[elementid].id[((bok+1)%4)]-1].BC==1){
//                    cout <<"Drugi if: "<< elementy[elementid].id[((bok+1)%4)]-1 << " " << wierzcholki[elementy[elementid].id[((bok+1)%4)]-1].x << " " << wierzcholki[elementy[elementid].id[((bok+1)%4)]-1].BC<<endl;
                    HBC(&elementy[elementid], wierzcholki, bok, flaga);
                }
            }
        }

        switch (flaga) {
            case 2: {
                el.waga[0] = 1.0;
                el.waga[1] = 1.0;
                for (int o = 0; o < 4; o++) {
                    N_ksi[o][0] = -1.0 / 4.0 * (1 - el.eta2[o]);
                    N_ksi[o][1] = 1.0 / 4.0 * (1 - el.eta2[o]);
                    N_ksi[o][2] = 1.0 / 4.0 * (1 + el.eta2[o]);
                    N_ksi[o][3] = -1.0 / 4.0 * (1 + el.eta2[o]);

                    N_eta[o][0] = -1.0 / 4.0 * (1 - el.ksi2[o]);
                    N_eta[o][1] = -1.0 / 4.0 * (1 + el.ksi2[o]);
                    N_eta[o][2] = 1.0 / 4.0 * (1 + el.ksi2[o]);
                    N_eta[o][3] = 1.0 / 4.0 * (1 - el.ksi2[o]);

                    Ni[o][0] = 0.25 * (1 - el.ksi2[o]) * (1 - el.eta2[o]);
                    Ni[o][1] = 0.25 * (1 + el.ksi2[o]) * (1 - el.eta2[o]);
                    Ni[o][2] = 0.25 * (1 + el.ksi2[o]) * (1 + el.eta2[o]);
                    Ni[o][3] = 0.25 * (1 - el.ksi2[o]) * (1 + el.eta2[o]);
                }
//                for(int ksi2 = 0; ksi2 < 4; ksi2++){
//                    for (int eta2 = 0; eta2<4; eta2++){
//                        cout << "ETA: "<<N_eta[ksi2][eta2]<<endl;
//                    }
//                    cout<<endl;
//                }
//                for(int ksi2 = 0; ksi2 < 4; ksi2++){
//                    for (int eta2 = 0; eta2<4; eta2++){
//                        cout << "KSI: "<<N_ksi[ksi2][eta2]<<endl;
//                    }
//                }
//                cout<<endl;
                break;
            }

            case 3: {
                el.waga[0] = 5.0 / 9.0;
                el.waga[1] = 8.0 / 9.0;
                el.waga[2] = 5.0 / 9.0;
                for (int fi = 0; fi < 9; fi++) {
                    N_ksi[fi][0] = -1.0 / 4.0 * (1 - el.eta3[fi]);
                    N_ksi[fi][1] = 1.0 / 4.0 * (1 - el.eta3[fi]);
                    N_ksi[fi][2] = 1.0 / 4.0 * (1 + el.eta3[fi]);
                    N_ksi[fi][3] = -1.0 / 4.0 * (1 + el.eta3[fi]);

                    N_eta[fi][0] = -1.0 / 4.0 * (1 - el.ksi3[fi]);
                    N_eta[fi][1] = -1.0 / 4.0 * (1 + el.ksi3[fi]);
                    N_eta[fi][2] = 1.0 / 4.0 * (1 + el.ksi3[fi]);
                    N_eta[fi][3] = 1.0 / 4.0 * (1 - el.ksi3[fi]);

                    Ni[fi][0] = 0.25 * (1 - el.ksi3[fi]) * (1 - el.eta3[fi]);
                    Ni[fi][1] = 0.25 * (1 + el.ksi3[fi]) * (1 - el.eta3[fi]);
                    Ni[fi][2] = 0.25 * (1 + el.ksi3[fi]) * (1 + el.eta3[fi]);
                    Ni[fi][3] = 0.25 * (1 - el.ksi3[fi]) * (1 + el.eta3[fi]);
                }
//                for(int ksi2 = 0; ksi2 < 9; ksi2++){
//                    for (int eta2 = 0; eta2<4; eta2++){
//                        cout << "ETA: "<<N_eta[ksi2][eta2];
//                    }
//                    cout<<endl;
//                }
//                for(int ksi2 = 0; ksi2 < 9; ksi2++){
//                    for (int eta2 = 0; eta2<4; eta2++){
//                        cout << "KSI: "<<N_ksi[ksi2][eta2];
//                    }
//                    cout<<endl;
//                }
                break;
            }

            case 4: {
                el.waga[0] = 0.347854;
                el.waga[1] = 0.652145;
                el.waga[2] = 0.652145;
                el.waga[3] = 0.347854;
                for (int fi = 0; fi < 16; fi++) {
                    N_ksi[fi][0] = -1.0 / 4.0 * (1 - el.eta4[fi]);
                    N_ksi[fi][1] = 1.0 / 4.0 * (1 - el.eta4[fi]);
                    N_ksi[fi][2] = 1.0 / 4.0 * (1 + el.eta4[fi]);
                    N_ksi[fi][3] = -1.0 / 4.0 * (1 + el.eta4[fi]);

                    N_eta[fi][0] = -1.0 / 4.0 * (1 - el.ksi4[fi]);
                    N_eta[fi][1] = -1.0 / 4.0 * (1 + el.ksi4[fi]);
                    N_eta[fi][2] = 1.0 / 4.0 * (1 + el.ksi4[fi]);
                    N_eta[fi][3] = 1.0 / 4.0 * (1 - el.ksi4[fi]);

                    Ni[fi][0] = 0.25 * (1 - el.ksi4[fi]) * (1 - el.eta4[fi]);
                    Ni[fi][1] = 0.25 * (1 + el.ksi4[fi]) * (1 - el.eta4[fi]);
                    Ni[fi][2] = 0.25 * (1 + el.ksi4[fi]) * (1 + el.eta4[fi]);
                    Ni[fi][3] = 0.25 * (1 - el.ksi4[fi]) * (1 + el.eta4[fi]);
                }
                break;
            }
        }

        for (int ip = 0; ip < flag; ip++) {
            for (int hi = 0; hi < 4; hi++) {
                dx_dksi[ip] += N_ksi[ip][hi] * x[hi]; //przypisanie wartosci do macierzy
                dy_dksi[ip] += N_ksi[ip][hi] * y[hi];
                dx_deta[ip] += N_eta[ip][hi] * x[hi];
                dy_deta[ip] += N_eta[ip][hi] * y[hi];
            }
        }

        for (int ip = 0; ip < flag; ip++) {
            J[ip][0][0] = dy_deta[ip];
            J[ip][1][0] = -dy_dksi[ip];
            J[ip][0][1] = -dx_deta[ip]; //okreslanie Jakobianu na odwrotnosci macierzy
            J[ip][1][1] = dx_dksi[ip];

            det_J[ip] = J[ip][0][0] * J[ip][1][1] - J[ip][0][1] * J[ip][1][0]; //liczenie wyznacznika
        }
//        for (int i = 0; i < 4; i++) {
//            cout << det_J[i] << endl;
//        }

        double dN_dx[16][4];
        double dN_dy[16][4];

        for (int ip = 0; ip < flag; ip++) {
            for (int j = 0; j < 4; j++) {
                dN_dx[ip][j] = 1.0 / det_J[ip] *
                               (N_ksi[ip][j] * dy_deta[ip] + N_eta[ip][j] * (-dy_dksi[ip])); // liczenie ze wzoru dN/dX
                dN_dy[ip][j] = 1.0 / det_J[ip] * (N_ksi[ip][j] * (-dx_deta[ip]) + N_eta[ip][j] * dx_dksi[ip]);
            }
        }
//        for (int i=0; i<4;i++){
//            for(int j = 0; j<4; j++) {
//                cout <<"dN_dx: "<< dN_dx[i][j]<< " ";
//            }
//        cout<<endl;
//        }
//        for (int i=0; i<4;i++){
//            for(int j = 0; j<4; j++) {
//                cout <<"dN_dy: "<< dN_dy[i][j]<< " ";
//            }
//            cout<<endl;
//        }
        double dN_dx_dN_dx_T[16][4][4];
        double dN_dy_dN_dy_T[16][4][4];
        double funkcja[16][4][4];

        for (int ip = 0; ip < flag; ip++) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    dN_dx_dN_dx_T[ip][i][j] = dN_dx[ip][i] * dN_dx[ip][j];
                    dN_dy_dN_dy_T[ip][i][j] = dN_dy[ip][i] * dN_dy[ip][j];
                    NTransponowane[ip][i][j] = Ni[ip][i] * Ni[ip][j];
                    funkcja[ip][i][j] = data.k * (dN_dx_dN_dx_T[ip][i][j] + dN_dy_dN_dy_T[ip][i][j]) * det_J[ip]* el.waga[ip / flaga] *
                                        el.waga[ip % flaga];
                    elementy[elementid].HL[i][j] += funkcja[ip][i][j]; //przypisanie wynikow do macierzy HL
                    elementy[elementid].CL[i][j] +=
                            NTransponowane[ip][i][j] * data.c * data.ro * det_J[ip] * el.waga[ip / flaga] *
                            el.waga[ip % flaga];
                }
            }
        }

//        for (int  id = 0; id < 4; id++) {
//            for (int j = 0; j < 4; j++) {
//                cout << elementy[elementid].CL[id][j] << setw(9);
//
//            }
//            cout << endl;
//        }
//        cout<<"Macierz HBC:"<<endl;
//        for(int wyp = 0; wyp<4; wyp++){
//            for (int wyp2 = 0; wyp2 < 4; wyp2++){
//
//                cout<<setw(9)<<elementy[elementid].HBCL[wyp][wyp2];
//            }
//            cout<<endl;
//        }

    }
    for (int hgp = 0; hgp < data.n_e; hgp++) {
        for (int hgi = 0; hgi < 4; hgi++) {
            for (int hgj = 0; hgj < 4; hgj++) {
                data.HG[elementy[hgp].id[hgi] - 1][elementy[hgp].id[hgj] - 1] += elementy[hgp].HL[hgi][hgj]+elementy[hgp].HBCL[hgi][hgj];
                data.CG[elementy[hgp].id[hgi] - 1][elementy[hgp].id[hgj] - 1] += elementy[hgp].CL[hgi][hgj];
            }
            data.PG[elementy[hgp].id[hgi]-1]+=elementy[hgp].PL[hgi];
        }
    }
    for (int zas1 = 0; zas1 < data.n_n; zas1++){
        data.PGZ[zas1] = -1 * data.PG[zas1];
        for (int zas2 = 0; zas2<data.n_n; zas2++){
            data.HGZ[zas1][zas2] = data.HG[zas1][zas2] + (data.CG[zas1][zas2]/data.dt);
            data.PGZ[zas1] += (data.CG[zas1][zas2]/data.dt) * data.t0[zas2];
        }
    }
//    for (int cgy = 0; cgy < 16; cgy++) {
//        for (int cgl = 0; cgl < 16; cgl++) {
//            cout<<setw(7) << data.CG[cgy][cgl] << " ";
//        }
//        cout << endl;
//
//    }
//    for (int cgy = 0; cgy < 16; cgy++) {
//        cout << data.PGZ[cgy];
//        cout << endl;
//    }
//for(int px = 0; px < 16; px++){
//    cout<<data.PG[px]<<" "<<endl;
//}

    dane->t0 = gauss(data.HGZ, data.PGZ, data.n_n);


    for (int i = 0; i < data.n_n; i++){
        data.PG[i] = 0.0;
        data.PGZ[i] = 0.0;
        for (int j = 0; j< data.n_n; j++){
            data.HG[i][j] = 0.0;
            data.CG[i][j] = 0.0;
            data.HGZ[i][j] = 0.0;
        }
    }
    for (int i = 0; i<data.n_e; i++){
        for (int j = 0; j<4; j++){
            elementy[i].PL[j] = 0.0;
            for (int k = 0; k<4; k++){
                elementy[i].HL[j][k] = 0.0;
                elementy[i].HBCL[j][k] = 0.0;
                elementy[i].CL[j][k] = 0.0;
            }
        }
    }
}

void Siatka() {
    GlobalData data;
    double dx = data.w / ((double)data.n_w - 1);
    double dy = data.h / ((double)data.n_h - 1);
    node* nodes = new node[data.n_n];
    int n = 0;
    for (int i = 0; i < data.n_w; i++) {
        for (int j = 0; j < data.n_h; j++) {
            nodes[i * data.n_h + j].x = i * dx;
            nodes[i * data.n_h + j].y = j * dy;
            n++;
//sprawdzamy, czy node[(i*data.n_h +j)] 1)czy sa najbardziej po lewej 2)czy sa najbardziej po prawej itp
            if((i*data.n_h +j)/data.n_h == 0 ||(i*data.n_h +j)/data.n_h == data.n_w-1 || (i*data.n_h +j)%data.n_h == 0 || (i*data.n_h +j)%data.n_h == data.n_h-1){
                nodes[(i*data.n_h +j)].BC = 1;
            }
//            cout << i*data.n_h+j << " " << nodes[(i*data.n_h +j)].x << " " << nodes[(i*data.n_h +j)].y << " " << nodes[(i*data.n_h +j)].BC<< endl;
        }
    }
    auto* el = new element[data.n_e];
    int licznik = 0;
    for (int i=0; i< data.n_e; i++) {
        el[i].id[0] = licznik + 1;
        el[i].id[1] = licznik + data.n_h + 1;
        el[i].id[2] = licznik + data.n_h + 2;
        el[i].id[3] = licznik + 2;
        licznik++;
        if((i+1)%(data.n_h-1)==0){
            licznik++;
        }
//     cout << "Element " << i + 1 << ": \t" << el[i].id[0] << " " << el[i].id[1] << " " << el[i].id[2] << " "<< el[i].id[3] << endl;
    }
    for(int i = 0; i<=ceil(data.st/data.dt); i++) {
        Macierz(el, nodes,&data);
        double tmax = 0;
        double tmin = 1000000000000;
        for (int t00 = 0; t00 < data.n_n; t00++){
            if(data.t0[t00]>tmax){
                tmax = data.t0[t00];
            }
            if(data.t0[t00]<tmin){
                tmin = data.t0[t00];
            }
//        cout << data.t0[t00] << endl;
        }
    cout << "tmax = "<<tmax<<endl;
    cout << "tmin = "<<tmin<<"\n\n";
    }
}


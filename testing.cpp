#include<iostream>
#include"testing.h"
#include"PseudoHB.h"
#include"HBRandom.h"
#include"Sommer.h"

using namespace std;


//basic test
void test1(){

SU3Grid mygrid(true);
SU3Grid mygrid2(true);

std::cout<<mygrid.GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();
std::cout<<mygrid2.GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();
mygrid2=mygrid;

std::cout<<mygrid.GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();
std::cout<<mygrid2.GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();

SU3Grid mygrid3=mygrid;
std::cout<<mygrid.GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();
std::cout<<mygrid2.GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();
    std::cout<<mygrid3.GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();


Modell mymodell(true);
std::cout<<"modell:"<<std::endl;
std::cout<<mymodell.GetModellGrid()(0).GetGrid()(0,0,0,0)<<std::endl;
cout<<mygrid2.GetGrid()(0,0,0,1)<<endl;
mymodell.ModifyLink(0,0,0,0,0,mygrid2.GetGrid()(0,0,0,1));
std::cout<<mymodell.GetModellGrid()(0).GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();

for(int i=0;i<10;i++){
cout<<"random Bernoulli: "<<Flip(0.5)<<endl;
}

for(int i=0;i<10;i++){
cout<<"random real: "<<GetRealRandom(-1,5)<<endl;
}

for(int i=0;i<10;i++){
    vector<double> randsph=RandOnSphere(3);
    cout<<"random on Sphere: "<<endl;
    for(vector<double>::iterator it=randsph.begin();it!=randsph.end();it++){
        cout<<*it<<endl;
    }
}

}//test 1

void WilsonAVGtest(){
    arma::cx_mat id3d(3,3,arma::fill::eye);
    Modell mymodell("p400config");

    ScaleSetV scaler(mymodell,3,3,0,0,0,0,1);

    std::ofstream resultfile;
    resultfile.precision(6);
    resultfile.open("wilsonloopAVGsomm.dat",std::ios::out);

    std::ofstream resultidfile;
    resultidfile.precision(6);
    resultidfile.open("isitid.dat",std::ios::out);

    for(int i=0;i<150;i++){
        scaler.WilsonAVG();
        arma::cx_mat isitid=(scaler.GetCorrelT().t()*scaler.GetCorrelT().i());
        resultfile<<real(scaler.GetCorrelT()(0,0))<<'\t'
                  <<real(scaler.GetCorrelT()(0,1))<<'\t'
                  <<real(scaler.GetCorrelT()(1,0))<<'\t'
                  <<real(scaler.GetCorrelT()(1,1))<<endl;
        resultidfile<<real(isitid(0,0))<<'\t'
                  <<real(isitid(0,1))<<'\t'
                  <<real(isitid(1,0))<<'\t'
                  <<real(isitid(1,1))<<endl;

        mymodell.HeatBathSweep();
    }
    resultfile.close();
    resultidfile.close();
   // std::cin.ignore();
    //std::cin.get();

}//wilsonavg test

void CorrelAVGtest(){

    arma::cx_mat id3d(3,3,arma::fill::eye);
    arma::cx_mat rescorr0(2,2,arma::fill::zeros);
    arma::cx_mat rescorrT(2,2,arma::fill::zeros);
    arma::cx_mat rescorrT1(2,2,arma::fill::zeros);

    Modell mymodell("p400config");

    ScaleSetV scaler(mymodell,3,3,0,0,0,0,1);
    scaler.CorrelMAVG(50,rescorr0,rescorrT,rescorrT1,"./");

}//correlavg test


void WilsonPot(){
    arma::cx_mat id3d(3,3,arma::fill::eye);
    Modell mymodell(id3d);
    for(int mcrun=0;mcrun<200;mcrun++){
        mymodell.HeatBathSweep();
    }

    for(int i=1;i<8;i++){
        ScaleSetV scaler(mymodell,3,i,0,0,0,0,1);
        scaler.BuildCorrelM();
        std::cout<<scaler.GetCorrelT()<<std::endl;
        std::cin.ignore();
        std::cin.get();
    }
}//wilsonpot

void TestPot(){
    arma::cx_mat id3d(3,3,arma::fill::eye);
    Modell mymodell("p400config");
    ScaleSetV scaler(mymodell,6,3,0,0,0,0,1);
    double potential=0;
    potential=scaler.CountV();
    std::ofstream resultfile;
    resultfile.precision(6);
    resultfile.open("potential.dat",std::ios::out);
    resultfile<<potential<<std::endl;
    resultfile.close();
}//testpot

void wilsontest(){
arma::cx_mat id3d(3,3,arma::fill::eye);

Modell mymodell(id3d);

string dir="./";

mymodell.writeToFileModell((dir+"initconfig").c_str());

//reach eq.
//Monte Carlo Run - 400 sweep
for(int mcrun=0;mcrun<400;mcrun++){
    mymodell.HeatBathSweep();
}

mymodell.writeToFileModell((dir+"p400config").c_str());

std::ofstream resultfile;
resultfile.precision(6);
resultfile.open("wilsonloop.dat",std::ios::out);

std::ofstream resultavgf;
resultavgf.precision(6);
resultavgf.open("wilsonloopAVGs.dat",std::ios::out);

for(int mcrun=0;mcrun<400;mcrun++){
    complex<double> wilson;
    wilson=mymodell.WilsonLoop(3,3,0,0,0,0,1);
    resultfile<<real(wilson)<<'\t'<<imag(wilson)<<endl;

    complex<double> wilsonavg(0,0);
    const int maxtdim=SU3Grid::GetTDim();
    const int maxdim=SU3Grid::GetDim();
    int counter=0;
    for(int i=0;i<maxtdim;i++){
        for(int j=0;j<maxdim;j++){
            for(int k=0;k<maxdim;k++){
                for(int l=0;l<maxdim;l++){
                    wilsonavg+=mymodell.WilsonLoop(3,3,i,j,k,l,1);
                    counter++;
                }
            }
        }
    }
    wilsonavg/=counter;
    resultavgf<<real(wilsonavg)<<'\t'<<imag(wilsonavg)<<endl;
    mymodell.HeatBathSweep();
}

}//wilsontest

void wilsontest2(){
arma::cx_mat id3d(3,3,arma::fill::eye);

Modell mymodell("p400config");

string dir="./";

mymodell.writeToFileModell((dir+"initconfig").c_str());

//reach eq.
//Monte Carlo Run - 400 sweep
//for(int mcrun=0;mcrun<400;mcrun++){
//    mymodell.HeatBathSweep();
//}

//mymodell.writeToFileModell((dir+"p400config").c_str());

std::ofstream resultfile;
resultfile.precision(6);
resultfile.open("wilsonloop.dat",std::ios::out);

std::ofstream resultavgf;
resultavgf.precision(6);
resultavgf.open("wilsonloopAVGs.dat",std::ios::out);

for(int mcrun=0;mcrun<400;mcrun++){
    complex<double> wilson;
    wilson=mymodell.WilsonLoop(3,3,0,0,0,0,1);
    resultfile<<real(wilson)<<'\t'<<imag(wilson)<<endl;

    complex<double> wilsonavg(0,0);
    int counter=0;
    for(int i=0;i<7;i++){
        for(int j=0;j<7;j++){
            for(int k=0;k<7;k++){
                for(int l=0;l<7;l++){
                    wilsonavg+=mymodell.WilsonLoop(3,3,i,j,k,l,1);
                    counter++;
                }
            }
        }
    }
    wilsonavg/=counter;
    resultavgf<<real(wilsonavg)<<'\t'<<imag(wilsonavg)<<endl;
    mymodell.HeatBathSweep();
}

}//wilsontest2

void isUnitary(){
string config;
cout<<"kerem a konfigfajlt"<<endl;
cin>>config;
Modell mymodell((config).c_str());
int maxtdim=SU3Grid::GetTDim();
int maxdim=SU3Grid::GetDim();
for(int g=0;g<4;g++){
    for(int i=0;i<maxtdim;i++){
        for(int j=0;j<maxdim;j++){
            for(int k=0;k<maxdim;k++){
                for(int l=0;l<maxdim;l++){
                    cout<<"is it Id?"<<endl;
                    cout<<mymodell.GetLink(g,i,j,k,l)*mymodell.GetLink(g,i,j,k,l).t()<<endl;
                }//for l space3
            }//for k space2
        }//for space (j)
    }//for time

}//for grid

}

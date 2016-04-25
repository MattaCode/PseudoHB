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

}


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

}

void CorrelAVGtest(){

    arma::cx_mat id3d(3,3,arma::fill::eye);
    arma::cx_mat rescorr0(2,2,arma::fill::zeros);
    arma::cx_mat rescorrT(2,2,arma::fill::zeros);
    arma::cx_mat rescorrT1(2,2,arma::fill::zeros);

    Modell mymodell("p400config");

    ScaleSetV scaler(mymodell,3,3,0,0,0,0,1);
    scaler.CorrelMAVG(50,rescorr0,rescorrT,rescorrT1,"./");

}


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
}

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
}

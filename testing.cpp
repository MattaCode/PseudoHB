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

void testsymmcorrelM(){
    arma::cx_mat id3d(3,3,arma::fill::eye);
    Modell mymodell(id3d);
    //reach eq.
    //Monte Carlo Run - 400 sweep
    for(int mcrun=0;mcrun<200;mcrun++){
        mymodell.HeatBathSweep();
    }
    ScaleSetV scaler(mymodell,3,3,0,0,0,0,1);
    scaler.isitsymm();
    scaler.CorrelMAVG(20);

}

#include<iostream>
#include<vector>
#include"PseudoHB.h"
#include"HBRandom.h"
#include"testing.h"
#include"algorithm.h"

using namespace std;

//AutoCorrMain
void AutoCorrMain(Modell & mymodell,const int & maxstep,string dir){

const int gridmax=SU3Grid::GetDim();
const int tgridmax=SU3Grid::GetTDim();

int tmin=0;
cout<<"tmin? "<<endl;
cin>>tmin;

int tincr=0;
cout<<"tincr?"<<endl;
cin>>tincr;

//int avgnum=0;
//cout<<"avgnum? "<<endl;
//cin>>avgnum;

int tnum=0;
cout<<"how many t-s? "<<endl;
cin>>tnum;

std::ofstream infofile;
infofile.open((dir+"info.dat").c_str(),std::ios::out);
infofile<<"AutoCorr of Polyakov loop averaged for space"<<endl;
infofile<<"dimensions: timelike: "<<tgridmax<<", spacelike: "<<gridmax<<endl;
infofile<<"beta: "<<mymodell.GetBeta()<<endl;
infofile<<"tmin: "<<tmin<<endl;
infofile<<"maxtime: "<<maxstep<<endl;
if(tmin==0){
    tmin=1;
    infofile<<"tmin changed to "<<tmin<<endl;
}
infofile<<"tincr: "<<tincr<<endl;
//infofile<<"avgnum: "<<avgnum<<endl;
infofile<<"number of t-s: "<<tnum<<endl;
infofile.close();
AutoCorrel autcorr(tnum,tmin,tincr);

complex<double> result=0;

//Monte Carlo Run
for(int mctime=0;mctime<maxstep;mctime++){
    result=mymodell.PolyakovLoopAVG();
    autcorr.InsertData(real(result));
    mymodell.HeatBathSweep();
}


std::ofstream resultfile;
resultfile.precision(6);
resultfile.open((dir+"NewAutCorr.dat").c_str(),std::ios::out);

for(int row=0;row<tnum;row++){
    resultfile<<(tmin+row*tincr)<<'\t'<<autcorr.GetResult(row)<<'\t'<<autcorr.GetCorr(row)<<'\t'<<autcorr.GetLoop0AVG()<<
    '\t'<<autcorr.GetAVG(row)<<endl;
}

resultfile.close();

}


int main(){

const int mcmaxtime=10;
arma::cx_mat id3d(3,3,arma::fill::eye);
Modell mymodell(id3d);

string dir="./";

//reach eq.
//Monte Carlo Run - 400 sweep
//for(int mcrun=0;mcrun<400;mcrun++){
//    mymodell.HeatBathSweep();
//}

//after 400 sweep we measure mcmaxtime step

//AutoCorrMain(mymodell,mcmaxtime,dir);

//debug
//Modell::GetPauli();

std::ofstream resultfile;
resultfile.precision(6);
resultfile.open("REALofPolyaAVG.dat",std::ios::out);

std::ofstream energyout;
energyout.precision(6);
energyout.open("MeanOfPlaqEn.dat",std::ios::out);

std::ofstream configuration;
configuration.precision(6);
configuration.open("Configuration.dat",std::ios::out);

complex<double> result=0;
for(int t=0;t<mcmaxtime;t++){

    result=mymodell.PolyakovLoopAVG();
    cout<<"*********REAL of Polya.AVG RESULT: "<<real(result)<<endl;
    resultfile<<t<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
    energyout<<t<<'\t'<<mymodell.CountMeanEnergyDens()<<endl;
    arma::cx_mat config(3,3,arma::fill::eye);
for(int i=0;i<4;i++){
    configuration<<(mymodell.GetModellGrid())(0).GetGrid()(i,0,0,0)<<endl;
    config*=(mymodell.GetModellGrid())(0).GetGrid()(i,0,0,0);
}

configuration<<trace(config)<<endl;
        mymodell.HeatBathSweep();

}
resultfile.close();
energyout.close();





configuration.close();

return 0;
}

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
//    result=mymodell.PolyakovLoopAVG();
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

const int mcmaxtime=2;
arma::cx_mat id3d(3,3,arma::fill::eye);
try{
Modell mymodell(id3d);

string dir="./";

mymodell.writeToFileModell((dir+"initconfig").c_str());

//reach eq.
//Monte Carlo Run - 400 sweep
for(int mcrun=0;mcrun<2;mcrun++){
    mymodell.HeatBathSweep();
}

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

    ofstream polyaloop;
    polyaloop.precision(6);

    ofstream energydens;
    energydens.precision(6);


complex<double> result=0;
for(int t=0;t<mcmaxtime;t++){
    ostringstream convert;
    convert<<t;
    polyaloop.open("PolyaLOOP"+convert.str()+".dat",ios::out);
    result=mymodell.PolyakovLoopAVG(polyaloop);
    polyaloop.close();
    cout<<"*********REAL of Polya.AVG RESULT: "<<real(result)<<endl;
    resultfile<<t<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
    energydens.open("EnergyDens"+convert.str()+".dat",ios::out);
    energyout<<t<<'\t'<<mymodell.CountMeanEnergyDens(energydens)<<endl;
    energydens.close();
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

cin.ignore();
cin.get();

mymodell.writeToFileModell((dir+"midconfig").c_str());

cin.ignore();
cin.get();


Modell mymod2((dir+"midconfig").c_str());
mymod2.writeToFileModell((dir+"afterreadconfig").c_str());

cin.ignore();
cin.get();


//reach eq.
//Monte Carlo Run - 400 sweep
for(int mcrun=0;mcrun<2;mcrun++){
    mymod2.HeatBathSweep();
    cin.ignore();
cin.get();
    mymodell.HeatBathSweep();
}

mymod2.writeToFileModell((dir+"finalconfig").c_str());
mymodell.writeToFileModell((dir+"originfinalconfig").c_str());

configuration.close();

}catch(const char * a){
    std::cerr<<"error detected: "<<a<<endl;
    cin.ignore();
    cin.get();

}

return 0;
}

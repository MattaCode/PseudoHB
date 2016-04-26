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

void MeasureHistosPolya(){
    const int mcmaxtime=200;
arma::cx_mat id3d(3,3,arma::fill::eye);

Modell mymodell(id3d);

string dir="./";

mymodell.writeToFileModell((dir+"initconfig").c_str());

//reach eq.
//Monte Carlo Run - 400 sweep
//for(int mcrun=0;mcrun<600;mcrun++){
//    mymodell.HeatBathSweep();
//}

//mymodell.writeToFileModell((dir+"p600config").c_str());

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


    ofstream polyaloop;
    polyaloop.precision(6);

    ofstream energydens;
    energydens.precision(6);


complex<double> result=0;
for(int t=0;t<200;t++){
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

        mymodell.HeatBathSweep();

}


mymodell.writeToFileModell((dir+"200config").c_str());

//reach eq.
//Monte Carlo Run - 400 sweep
for(int mcrun=0;mcrun<400;mcrun++){
    mymodell.HeatBathSweep();
}

mymodell.writeToFileModell((dir+"p600config").c_str());

result=0;
for(int t=600;t<800;t++){
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

        mymodell.HeatBathSweep();

}
resultfile.close();
energyout.close();

mymodell.writeToFileModell((dir+"finalonfig").c_str());
}

void Wilson(){
string dir;
string initconfig;
int r=0;
int t=0;
cout<<"kerem a konyvtarat: "<<endl;
cin>>dir;
cout<<dir<<endl;
cout<<"kerem a kezdeti konfiguraciot!"<<endl;
cin>>initconfig;
cout<<initconfig<<endl;
cout<<"wilson loop terszeru hossza"<<endl;
cin>>r;
cout<<r<<endl;
cout<<"wilson loop idoszeru hossza"<<endl;
cin>>t;
cout<<t<<endl;

std::cout<<(dir+"info.dat").c_str()<<std::endl;
    std::ofstream infofile;
    infofile.open((dir+"info.dat").c_str(),std::ios::out);
    if(!infofile.is_open()) throw "can not open";
    infofile<<"tdim: "<<SU3Grid::GetTDim()<<endl;
    infofile<<"space dim: "<<SU3Grid::GetDim()<<endl;
    infofile<<"beta: "<<Modell::GetBeta()<<endl;
    infofile<<"initconfig: "<<initconfig<<endl;
    infofile<<"wilson loop R: "<<r<<endl;
    infofile<<"wilson loop T: "<<t<<endl;
    infofile.close();

    Modell mymodell((initconfig).c_str());
    std::complex<double> wilsonavg;
    std::ofstream resultfile;
    resultfile.precision(6);
    resultfile.open((dir+"wilsonavg.dat").c_str(),std::ios::out);
    for(int i=0;i<20;i++){
    wilsonavg=mymodell.WilsonAvg(r,t,1);
    resultfile<<real(wilsonavg)<<'\t'<<imag(wilsonavg)<<std::endl;
    for(int j=0;j<10;j++){
        mymodell.HeatBathSweep();
    }
    }
    resultfile.close();

}

int main(){
try{

wilsontest2();

}catch(const char * a){
    std::cerr<<"error detected: "<<a<<endl;
    cin.ignore();
    cin.get();

}

return 0;
}

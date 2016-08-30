#include<iostream>
#include<vector>
#include"PseudoHB.h"
#include"HBRandom.h"
#include"testing.h"
#include"algorithm.h"
#include"Sommer.h"
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

void MeasureHistoPolya(){
//    const int mcmaxtime=200;
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
    polyaloop.open(("PolyaLOOP"+convert.str()+".dat").c_str(),ios::out);
    result=mymodell.PolyakovLoopAVG(polyaloop);
    polyaloop.close();
    cout<<"*********REAL of Polya.AVG RESULT: "<<real(result)<<endl;
    resultfile<<t<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
    energydens.open(("EnergyDens"+convert.str()+".dat").c_str(),ios::out);
    energyout<<t<<'\t'<<mymodell.CountMeanEnergyDens(energydens,false)<<endl;
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
    polyaloop.open(("PolyaLOOP"+convert.str()+".dat").c_str(),ios::out);
    result=mymodell.PolyakovLoopAVG(polyaloop);
    polyaloop.close();
    cout<<"*********REAL of Polya.AVG RESULT: "<<real(result)<<endl;
    resultfile<<t<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
    energydens.open(("EnergyDens"+convert.str()+".dat").c_str(),ios::out);
    energyout<<t<<'\t'<<mymodell.CountMeanEnergyDens(energydens,false)<<endl;
    energydens.close();

        mymodell.HeatBathSweep();

}
resultfile.close();
energyout.close();

mymodell.writeToFileModell((dir+"finalonfig").c_str());
}

//SommerScale potential
void SommerPot(string dir,string initconfig,const int r,const int t,const int maxsmear){


    Modell mymodell((initconfig).c_str());
    ScaleSetV scaler(mymodell,r,t,0,0,0,0,1,maxsmear);
    double potential=0;
    potential=scaler.CountV(dir);
    std::ofstream resultfile;
    resultfile.precision(6);
    resultfile.open((dir+"potential.dat").c_str(),std::ios::out);
    resultfile<<potential<<std::endl;
    resultfile.close();
mymodell.writeToFileModell((dir+"Finalconfig").c_str());
}
void SommerPotMain(){
string dir;
string initconfig;
int r=0;
int t=0;
int maxsmear=0;
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
cout<<"maxsmear"<<endl;
cin>>maxsmear;
std::cout<<"maxsmear: "<<maxsmear<<std::endl;
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
infofile<<"maxsmear: "<<maxsmear<<endl;
    infofile.close();
    SommerPot(dir,initconfig,r,t,maxsmear);

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
mymodell.writeToFileModell((dir+"finalconfig").c_str());
}
//Just HeatBath
void HBRun(){
string dir;
string initconfig="none";
 string outconfig;
 int maxtime=400;
 bool isIDconfig=false;
cout<<"kerem a konyvtarat: "<<endl;
cin>>dir;
cout<<dir<<endl;
 cout<<"identity init config? <1,0>"<<endl;
 cin>>isIDconfig;
 if(isIDconfig){
   arma::cx_mat id3d(3,3,arma::fill::eye);
   Modell initmodell(id3d);
   initmodell.writeToFileModell((dir+"initIDconfig").c_str());
   initconfig=(dir+"initIDconfig");
}
 else{
cout<<"kerem a kezdeti konfiguraciot!"<<endl;
cin>>initconfig;
 }
cout<<initconfig<<endl;
 cout<<"kerem a max idot"<<endl;
 cin>>maxtime;
 cout<<"veg konfig neve: "<<endl;
 cin>>outconfig;

std::cout<<(dir+"info.dat").c_str()<<std::endl;
    std::ofstream infofile;
    infofile.open((dir+"info.dat").c_str(),std::ios::out);
    if(!infofile.is_open()) throw "can not open";
    infofile<<"tdim: "<<SU3Grid::GetTDim()<<endl;
    infofile<<"space dim: "<<SU3Grid::GetDim()<<endl;
    infofile<<"beta: "<<Modell::GetBeta()<<endl;
    infofile<<"initconfig: "<<initconfig<<endl;
    infofile<<"futasido: "<<maxtime<<endl;
    infofile<<"outconfig: "<<outconfig<<endl;
 infofile.close();

    Modell mymodell((initconfig).c_str());

    for(int i=0;i<maxtime;i++){
        mymodell.HeatBathSweep();
    }
   mymodell.writeToFileModell((dir+outconfig).c_str());
}

void isUnitary(){
string config;
cout<<"kerem a konfigfajlt"<<endl;
cin>>config;
 cout<<"valasztott konfigfajl: "<<endl;
 cout<<config<<endl;
Modell mymodell((config).c_str());
int maxtdim=SU3Grid::GetTDim();
int maxdim=SU3Grid::GetDim();
for(int g=0;g<4;g++){
    for(int i=0;i<maxtdim;i++){
        for(int j=0;j<maxdim;j++){
            for(int k=0;k<maxdim;k++){
                for(int l=0;l<maxdim;l++){
                    cout<<"is it Id?"<<endl;
                    cout<<det(mymodell.GetLink(g,i,j,k,l)*mymodell.GetLink(g,i,j,k,l).t())<<endl;
                }//for l space3
            }//for k space2
        }//for space (j)
    }//for time

}//for grid

}

void CoarseEnergyHist(){
    string dir;
    string initconfig="none";
    string outconfig;
    int maxtime=400;
    bool isIDconfig=false;
    bool istimelike=false;
    cout<<"kerem a konyvtarat: "<<endl;
    cin>>dir;
    cout<<dir<<endl;
    cout<<"identity init config? <1,0>"<<endl;
    cin>>isIDconfig;
    if(isIDconfig){
    arma::cx_mat id3d(3,3,arma::fill::eye);
    Modell initmodell(id3d);
    initmodell.writeToFileModell((dir+"initIDconfig").c_str());
    initconfig=(dir+"initIDconfig");
    }
    else{
    cout<<"kerem a kezdeti konfiguraciot!"<<endl;
    cin>>initconfig;
    }
    cout<<initconfig<<endl;
    cout<<"kerem a max idot"<<endl;
    cin>>maxtime;
    cout<<"veg konfig neve: "<<endl;
    cin>>outconfig;
    cout<<"timelike only? <y=1, n=0> "<<endl;
    cin>>istimelike;

    std::cout<<(dir+"info.dat").c_str()<<std::endl;
    std::ofstream infofile;
    infofile.open((dir+"info.dat").c_str(),std::ios::out);
    if(!infofile.is_open()) throw "can not open";
    infofile<<"BoxEnHist meres!"<<endl;
    infofile<<"tdim: "<<SU3Grid::GetTDim()<<endl;
    infofile<<"space dim: "<<SU3Grid::GetDim()<<endl;
    infofile<<"beta: "<<Modell::GetBeta()<<endl;
    infofile<<"initconfig: "<<initconfig<<endl;
    infofile<<"futasido: "<<maxtime<<endl;
    infofile<<"outconfig: "<<outconfig<<endl;
    infofile<<"istimelike: "<<istimelike<<endl;
    infofile.close();

    Modell mymodell((initconfig).c_str());

    const int timedim=SU3Grid::GetTDim();
    const int spacedim=SU3Grid::GetDim();
    const int mindim=min(timedim,spacedim);
    //const int maxlimit=mindim/2+1;
    const int maxlimit=2;
    std::ofstream energyout;
    energyout.precision(6);

std::ofstream energydens;
energydens.precision(6);

for(int t=0;t<maxtime;t++){
    ostringstream convert;
    convert<<t;
    for (int boxs=2;boxs<maxlimit;boxs++){
    ostringstream conv2;
    conv2<<boxs;
    energyout.open((dir+"t"+convert.str()+"box"+conv2.str()+"BoxEnHisto.dat").c_str(),std::ios::out);
    mymodell.BoxEnHisto(boxs,energyout,istimelike);
    energyout.close();
    }
    energydens.open((dir+"EnergyDens"+convert.str()+".dat").c_str(),ios::out);
    energyout<<t<<'\t'<<mymodell.CountMeanEnergyDens(energydens,istimelike)<<endl;
    energydens.close();
    mymodell.HeatBathSweep();
}

    mymodell.writeToFileModell((dir+outconfig).c_str());

}


int main(){
try{
int switcher=0;
//menu
do{
	cout<<"1: HeatBathRun"<<endl;
	cout<<"3: SommerPot run from initconfig"<<endl;
	cout<<"4: WilsonPot run from initconfig"<<endl;
	cout<<"5: AutoCorr"<<endl;
	cout<<"6: MeasureHistoPolya"<<endl;
	cout<<"8: Is it Unitary?"<<endl;
	cout<<"9: Box Histogram"<<endl;
	cout<<"7: exit"<<endl;
	cout<<endl;
	cout<<"Choose your destiny: "<<endl;
	cin>>switcher;
	switch(switcher){
	case 1:
	  HBRun();

	break;
	case 2:
	//TO DO
	break;
	case 3:
		SommerPotMain();
	break;
	case 4:
		Wilson();

	break;
	case 5:
		//AutoCorrMain();
		//cin.get();
	break;
	case 6:
		MeasureHistoPolya();
	break;
	case 7:
	  break;
    case 8:
        isUnitary();
    break;
    case 9:
        CoarseEnergyHist();
    break;
	}//switch

}while(switcher!=7);


}catch(const char * a){
    std::cerr<<"error detected: "<<a<<endl;
    cin.ignore();
    cin.get();

}

return 0;
}

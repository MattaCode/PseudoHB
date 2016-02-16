#include<iostream>
#include<vector>
#include"PseudoHB.h"
#include"HBRandom.h"
#include"testing.h"

using namespace std;




int main(){

const int mcmaxtime=100;
arma::cx_mat id3d(3,3,arma::fill::eye);
Modell mymodell(id3d);

//debug
//Modell::GetPauli();

std::ofstream resultfile;
resultfile.precision(6);
resultfile.open("REALofPolyaAVG.dat",std::ios::out);
complex<double> result=0;
for(int t=0;t<mcmaxtime;t++){
    mymodell.HeatBathSweep();
    result=mymodell.PolyakovLoopAVG();
    cout<<"*********REAL of Polya.AVG RESULT: "<<real(result)<<endl;
    resultfile<<t<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;

}
resultfile.close();
return 0;
}

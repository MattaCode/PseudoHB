#include<iostream>
#include"testing.h"
#include"PseudoHB.h"
#include"HBRandom.h"

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

}

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

}

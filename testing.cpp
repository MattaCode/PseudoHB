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

//void testsymmcorrelM(){
//    arma::cx_mat id3d(3,3,arma::fill::eye);
//    Modell mymodell(id3d);
//    //reach eq.
//    //Monte Carlo Run - 400 sweep
//    for(int mcrun=0;mcrun<1;mcrun++){
//        mymodell.HeatBathSweep();
//    }
//    ScaleSetV scaler(mymodell,1,1,0,0,0,0,1);
//    scaler.isitsymm();
//    std::cin.ignore();
//    std::cin.get();
//
//    for(int mcrun=0;mcrun<50;mcrun++){
//        mymodell.HeatBathSweep();
//    }
//
//    scaler.CorrelMAVG(30);
//
//}

//void testsmear(){
//    arma::cx_mat id3d(3,3,arma::fill::eye);
//    Modell mymodell(id3d);
//mymodell.HeatBathSweep();
//    ScaleSetV scaler(mymodell,3,3,0,0,0,0,1);
//    for(int smearlevel=0;smearlevel<10;smearlevel++){
//        std::cout<<"smear level: "<<smearlevel<<std::endl;
//        std::cout<<scaler.GetSpace0Grid()(0)<<std::endl;
//        std::cout<<scaler.GetSpaceTGrid()(0)<<std::endl;
//        std::cin.ignore();
//        std::cin.get();
//        scaler.Smearing0();
//       // scaler.SmearingT();
//        std::cin.ignore();
//        std::cin.get();
//        std::cout<<"scaler.getspacegrid:"<<std::endl;
//        std::cout<<scaler.GetSpace0Grid()(0)<<std::endl;
//        std::cout<<scaler.GetSpaceTGrid()(0)<<std::endl;
//        std::cin.ignore();
//        std::cin.get();
//        std::cout<<"mymodell element"<<std::endl;
//        std::cout<<mymodell.GetModellGrid()(1).GetGrid()(0,0,0,0)<<std::endl;
//
//    }
//
//}

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
//
//void correltest(){
//    arma::cx_mat id3d(3,3,arma::fill::eye);
//    Modell mymodell(id3d);
//
//
//    ScaleSetV scaler(mymodell,3,3,0,0,0,0,1);
//    arma::cx_mat correl(2,2,arma::fill::zeros);
//    mymodell.HeatBathSweep();
//    //for(int i=0;i<200;i++){
//
//    arma::cx_mat spline0(3,3,arma::fill::eye);
//    arma::cx_mat splineT(3,3,arma::fill::eye);
//    arma::cx_mat tlineup(3,3,arma::fill::eye);
//    arma::cx_mat tlinedown(3,3,arma::fill::eye);
//    scaler.CountTimeLineDown(tlinedown,0,0,0,0+3-1,3);
//    scaler.CountTimeLineUp(tlineup,0+3,0,0,0,3);
//    scaler.InitSpaceLikeTInv(3,scaler.ModSpaceTGrid());
//for(int i=0;i<2;i++){
//    scaler.InitSpaceLikeT(0,scaler.ModSpace0Grid());
//    for(int j=0;j<2;j++){
//        scaler.CountSpaceLine0(spline0);
//        scaler.CountSpaceLineT(splineT);
//        //IT should BE Changed!
//        correl(i,j)+=real(trace(spline0*tlineup*splineT*tlinedown));
//        scaler.Smearing0();
//
//    }//for j smear
//    scaler.SmearingT(3);
//}//for i smear
////mymodell.HeatBathSweep();
////}
//
//std::cout<<correl<<std::endl;
//
//}

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
    ScaleSetV scaler(mymodell,3,3,0,0,0,0,1);
    double potential=0;
    potential=scaler.CountV();
    std::ofstream resultfile;
    resultfile.precision(6);
    resultfile.open("potential.dat",std::ios::out);
    resultfile<<potential<<std::endl;
    resultfile.close();
}

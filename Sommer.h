#ifndef SOMMER_H
#define SOMMER_H
#include<complex>
#include<armadillo>
#include"fftw++-1.13/Array.h"
#include"PseudoHB.h"

class ScaleSetV{
Modell & mymodell;

const int R; //spacelike length of Wilson loop
const int T; //timelike length of Wilson loop
const int initt; //init idx timelike
const int initx;
const int inity;
const int initz;
const int spacegrididx; //spacelike direction of Wilson loop
const double alpha; //smear parameter
const int maxsmearlevel;

Array::array1<arma::cx_mat> spacelike_0; //smeared spacelike edges at t=0
Array::array1<arma::cx_mat> spacelike_T; //smeared spacelike edges at t=T-1

arma::cx_mat correlT; //correl matrix at T (and fix R)
arma::cx_mat correlT1;//correl matrix at T+1 (and fix R)



public:

//ctr
ScaleSetV(Modell &,const int,const int,const int,const int,const int,const int,const int);

//to do
//copy

//to do
//assignment

//helper for ctr
void InitSpaceLikeT(const int,Array::array1<arma::cx_mat> );

//for one link
void SmearPt1(arma::cx_mat &,int,int,int,const int);
void SmearPt2(arma::cx_mat &);
void Smearing0();
void SmearingT();

void CountTimeLineUp(arma::cx_mat &,const int,const int,const int,const int,const int);
void CountTimeLineDown(arma::cx_mat &,const int,const int,const int,const int,const int);
void CountSpaceLine0(arma::cx_mat &);
void CountSpaceLineT(arma::cx_mat &);
void BuildCorrelM();

void CorrelMAVG();

double CountV();




~ScaleSetV();


};

#endif

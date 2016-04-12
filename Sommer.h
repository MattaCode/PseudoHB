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
int initt; //init idx timelike
int initx;
int inity;
int initz;
const int oinitt; //init idx timelike
const int oinitx;
const int oinity;
const int oinitz;
const int spacegrididx; //spacelike direction of Wilson loop
const double alpha; //smear parameter
const int maxsmearlevel;

Array::array1<arma::cx_mat> spacelike_0; //smeared spacelike edges at t=0
/****/
/*inverse-aka-.t() links*/
/****/
Array::array1<arma::cx_mat> spacelike_T; //smeared spacelike edges at t=T-1

arma::cx_mat correlT; //correl matrix at T (and fix R)
arma::cx_mat correlT1;//correl matrix at T+1 (and fix R)
arma::cx_mat correl0;//correl matrix at t0+1


public:

//ctr
ScaleSetV(Modell &,const int,const int,const int,const int,const int,const int,const int);

//to do
//copy

//to do
//assignment

//helper for ctr
void InitSpaceLikeT(const int,Array::array1<arma::cx_mat> &);
void InitSpaceLikeTInv(const int,Array::array1<arma::cx_mat> & );

void TriplForw(arma::cx_mat &,int,int,int,const int,const int);
void TriplBack(arma::cx_mat &,int,int,int,const int,const int);
//for one link
void SmearPt1Forw(arma::cx_mat &,int,int,int,const int);
void SmearPt1Back(arma::cx_mat &,int,int,int,const int);
void SmearPt2(arma::cx_mat &);
void Smearing0();
void SmearingT(const int);

void CountTimeLineUp(arma::cx_mat &,const int,const int,const int,const int,const int);
void CountTimeLineDown(arma::cx_mat &,const int,const int,const int,const int,const int);
void CountSpaceLine0(arma::cx_mat &);
void CountSpaceLineT(arma::cx_mat &);
void BuildCorrelM();

void CorrelMAVG(const int);
void WilsonAVG();

double CountV();

void isitsymm();

const Array::array1<arma::cx_mat> & GetSpace0Grid()const;
const Array::array1<arma::cx_mat> & GetSpaceTGrid()const;
Array::array1<arma::cx_mat> & ModSpace0Grid();
Array::array1<arma::cx_mat> & ModSpaceTGrid();


~ScaleSetV();


};

#endif

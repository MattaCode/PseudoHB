#ifndef SOMMER_H
#define SOMMER_H
#include<complex>
#include<armadillo>
#include<string>
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

arma::cx_mat correlT; //correl matrix at T (and fix R)
arma::cx_mat correlT1;//correl matrix at T+1 (and fix R)
arma::cx_mat correl0;//correl matrix at t0+1

Array::array2<SU3Grid> smeared;
Array::array2<SU3Grid> smearedinv;

public:

//ctr
ScaleSetV(Modell &,const int,const int,const int,const int,const int,const int,const int);
ScaleSetV(Modell &,const int,const int,const int,const int,const int,const int,const int,const int);

//to do
//copy

//to do
//assignment


void InitSm();



void TriplForw(arma::cx_mat &,int,int,int,const int,const int,const int,const int);
void TriplBack(arma::cx_mat &,int,int,int,const int,const int,const int,const int);
//for one link
void SmearPt1Forw(arma::cx_mat &,int,int,int,const int,const int,const int);
void SmearPt1Back(arma::cx_mat &,int,int,int,const int,const int,const int);
void SmearPt2(arma::cx_mat &);

void Smear();

void CountTimeLineUp(arma::cx_mat &,const int,const int,const int,const int,const int);
void CountTimeLineDown(arma::cx_mat &,const int,const int,const int,const int,const int);
void CountSpaceLine0(arma::cx_mat &,const int);
void CountSpaceLineT(arma::cx_mat &,const int,const int);
void BuildCorrelM();

void CorrelMAVG(const int,arma::cx_mat &,arma::cx_mat &,arma::cx_mat &,std::string);
void WilsonAVG(std::ofstream &,std::ofstream &,std::ofstream &);
void WilsonAVG();

double CountV();
double CountV(std::string);

void isitsymm();

void Symmetrize(arma::mat &);

const Array::array1<arma::cx_mat> & GetSpace0Grid()const;
const Array::array1<arma::cx_mat> & GetSpaceTGrid()const;
Array::array1<arma::cx_mat> & ModSpace0Grid();
Array::array1<arma::cx_mat> & ModSpaceTGrid();

const arma::cx_mat & GetCorrelT()const;

~ScaleSetV();


};

#endif

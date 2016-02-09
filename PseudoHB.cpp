#include<iostream>
#include<cmath>
#include<armadillo>
#include"PseudoHB.h"

using namespace arma;
using namespace std;

//random su3 matrix
cx_mat RandSU3(cx_mat su3){
    //generate random matrix - gauss(0,1)
    arma_rng::set_seed_random();  // set the seed to a random value
    mat A(3,3,fill::randn);
    mat B(3,3,fill::randn);
    cx_mat temp(A,B);
    su3=temp;
    //qr decomp
    cx_mat Q,R;
    bool qr_succ=true;
    qr_succ=qr(Q,R,su3);
    if(!qr_succ) throw "constructor qr fails";
    //qr postproc to get unique result
    cx_mat Lambda(3,3,fill::zeros);
    for(int i=0;i<3;i++){
        Lambda(i,i)=R(i,i)/abs(R(i,i));
    }
    su3=Q*Lambda;
    //let it be SU(3) instead of U(3)
    complex<double> determinant=det(su3);
    su3=su3/determinant;
    return su3;
}


SU3Grid::SU3Grid():ei(tdim,dim,dim,dim){}

//construct from random
SU3Grid::SU3Grid(bool flag):ei(tdim,dim,dim,dim){
    for(int i=0;i<tdim;i++){
        for(int j=0;j<dim;j++){
            for(int k=0;k<dim;k++){
                for(int l=0;l<dim;l++){
                    cx_mat su3(3,3);
                    ei(i,j,k,l)=RandSU3(su3);
                }//for l spacelike
            }//for k spacelike
        }//for j spacelike
    }//for i timelike
}

SU3Grid::~SU3Grid(){}

//get value of dimension
const int SU3Grid::GetDim(){
    return dim;
}

//get value of dimension
const int SU3Grid::GetTDim(){
    return tdim;
}

//GetGrid for reading
const Array::array4<arma::cx_mat>& SU3Grid::GetGrid()const{
    return ei;
}

//GetGrid for modifing a link
Array::array4<cx_mat>& SU3Grid::ModifyGrid(){
    return ei;
}

//default constr
Modell::Modell():grid(4){}

//RandomInit
void Modell::RandomInit(){
    int tgridmax=SU3Grid::GetTDim();
    int gridmax=SU3Grid::GetDim();
    //sweep cycle
        for(int ei=0;ei<4;ei++){
            for(int i=0;i<tgridmax;i++){
                for(int j=0;j<gridmax;j++){
                    for(int k=0;k<gridmax;k++){
                        for(int l=0;l<gridmax;l++){
                            cx_mat su3;
                            grid(ei).ModifyGrid()(i,j,k,l)=RandSU3(su3);
                        }//for l
                    }
                }//for j
            }//for timelike
        }//for grid
    }

//construct from random
Modell::Modell(bool flag):grid(4){
    RandomInit();
}

//access for reading
const Array::array1<SU3Grid>& Modell::GetModellGrid()const{
    return grid;
}

Modell::~Modell(){}

const double Modell::beta=6;
const int SU3Grid::dim=4;
const int SU3Grid::tdim=8;

#include<iostream>
#include<complex>
#include<cmath>
#include<vector>
#include<armadillo>
#include"PseudoHB.h"
#include "HBRandom.h"

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

/*******************/
/******SU3GRID******/
/*******************/

SU3Grid::SU3Grid():ei(tdim,dim,dim,dim){
//debug
cout<<"su3grid def. ctr. call"<<endl;
}

//construct from random
SU3Grid::SU3Grid(bool flag):ei(tdim,dim,dim,dim){
//debug
cout<<"su3grid ctr from rand call"<<endl;
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



//assignment operator
SU3Grid& SU3Grid::operator=(const SU3Grid& fromsu3grid){
//debug
cout<<"su3grid assignment call"<<endl;
    if(this!=&fromsu3grid){
        ei=fromsu3grid.ei;
    }
    return *this;
}

//copy
SU3Grid::SU3Grid(const SU3Grid& su3grid):ei(tdim,dim,dim,dim){
//debug
cout<<"su3grid copy call"<<endl;
    *this=su3grid;
}

//TO DO
    //write SU3Grid to os
    //void writeSU3Grid(std::ostream &){}

    //read SU3Grid from is
    //void readSU3Grid(std::istream&){}

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

SU3Grid::~SU3Grid(){}
/*******************************************************/


/*********************/
/***** MODELL ********/
/*********************/

//default constr
Modell::Modell():grid(4),su3staple(3,3,fill::zeros),su2staple(2,2,fill::zeros),su2strootdet(0){}

//construct from random
Modell::Modell(bool flag):grid(4),su3staple(3,3,fill::zeros),su2staple(2,2,fill::zeros),su2strootdet(0){
    RandomInit();
}

//TO DO
    //construct from file
    //Modell(const char* );
    //SU3Grid(std::istream&);

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

//ASSIGNMENT
//Modell& operator=(const Modell&);

    //COPY
    //Modell(const Modell&);
    //get beta
   // static const double GetBeta();

    //write Modell to os
    //void writeModell( std::ostream& )const;

    //write Modell to file
    //void writeToFileModell( const char* )const;

    //read Modell from is
    //void readModell( std::istream& );

    //read Modell from file
    //void readFromFileModell( const char* );

    //count plaquett at idx, for two selected direction i,j
    //up initialized as Identity
    //assumes calls with existing indices idx-y-z-k and i-j
    //void CountUp(int ,int ,int ,int ,const unsigned int , const unsigned int ,arma::cx_mat & );

//count forward staple (plaquett without the selected link)
//result initialized as Identity
void Modell::TriplUForw(const unsigned int grididx, const unsigned int grididx2,
int idx, int idy, int idz, int idk, arma::cx_mat & result){
    //maxdim of SU3Grid
    const int maxdim=SU3Grid::GetDim();
    const int maxtdim=SU3Grid::GetTDim();
    if(grididx==grididx2) throw "equal direction error";
        //step along first direction: index->index+ei*a
        switch(grididx){
        case 0:
            idx++;
            //modulo for periodic BC
            idx=(idx+maxtdim)%maxtdim;
            break;
        case 1:
            idy++;
            idy=(idy+maxdim)%maxdim;
            break;
        case 2:
            idz++;
            idz=(idz+maxdim)%maxdim;
            break;
        case 3:
            idk++;
            idk=(idk+maxdim)%maxdim;
            break;
        }//switch
        //multipl by second edge (j) at index+ei*a
        result=(grid(grididx2).GetGrid())(idx,idy,idz,idk)*result;
        //step back along first direction
        switch(grididx){
        case 0:
            idx--;
            idx=(idx+maxtdim)%maxtdim;
            break;
        case 1:
            idy--;
            idy=(idy+maxdim)%maxdim;
            break;
        case 2:
            idz--;
            idz=(idz+maxdim)%maxdim;
            break;
        case 3:
            idk--;
            idk=(idk+maxdim)%maxdim;
            break;
        }//switch
        //step along second direction
        switch(grididx2){
        case 0:
            idx++;
            idx=(idx+maxtdim)%maxtdim;
            break;
        case 1:
            idy++;
            idy=(idy+maxdim)%maxdim;
            break;
        case 2:
            idz++;
            idz=(idz+maxdim)%maxdim;
            break;
        case 3:
            idk++;
            idk=(idk+maxdim)%maxdim;
            break;
        }//switch
        //multipl by third edge at index+ej*a
        result=(grid(grididx).GetGrid())(idx,idy,idz,idk).t()*result;
        //step back along second direction
        switch(grididx2){
        case 0:
            idx--;
            idx=(idx+maxtdim)%maxtdim;
            break;
        case 1:
            idy--;
            idy=(idy+maxdim)%maxdim;
            break;
        case 2:
            idz--;
            idz=(idz+maxdim)%maxdim;
            break;
        case 3:
            idk--;
            idk=(idk+maxdim)%maxdim;
            break;
        }//switch
        //multipl by fourth edge at index - closing the cycle
        result=(grid(grididx2).GetGrid())(idx,idy,idz,idk).t()*result;
}

//count backward staple (plaquett without the selected link)
//result initialized as Identity
void Modell::TriplURev(const unsigned int grididx, const unsigned int grididx2,
int idx, int idy, int idz, int idk, arma::cx_mat & result){
    //maxdim of SU3Grid
    const int maxdim=SU3Grid::GetDim();
    const int maxtdim=SU3Grid::GetTDim();
    if(grididx==grididx2) throw "equal direction error";
    //step back along second direction
    switch(grididx2){
        case 0:
            idx--;
            idx=(idx+maxtdim)%maxtdim;
            break;
        case 1:
            idy--;
            idy=(idy+maxdim)%maxdim;
            break;
        case 2:
            idz--;
            idz=(idz+maxdim)%maxdim;
            break;
        case 3:
            idk--;
            idk=(idk+maxdim)%maxdim;
            break;
        }//switch
        //step along first direction
        switch(grididx){
        case 0:
            idx++;
            //modulo for periodic BC
            idx=(idx+maxtdim)%maxtdim;
            break;
        case 1:
            idy++;
            idy=(idy+maxdim)%maxdim;
            break;
        case 2:
            idz++;
            idz=(idz+maxdim)%maxdim;
            break;
        case 3:
            idk++;
            idk=(idk+maxdim)%maxdim;
            break;
        }//switch
        result=(grid(grididx2).GetGrid())(idx,idy,idz,idk).t()*result;
        //step back along first direction
        switch(grididx){
        case 0:
            idx--;
            idx=(idx+maxtdim)%maxtdim;
            break;
        case 1:
            idy--;
            idy=(idy+maxdim)%maxdim;
            break;
        case 2:
            idz--;
            idz=(idz+maxdim)%maxdim;
            break;
        case 3:
            idk--;
            idk=(idk+maxdim)%maxdim;
            break;
        }//switch
        result=(grid(grididx).GetGrid())(idx,idy,idz,idk).t()*result;

        result=(grid(grididx2).GetGrid())(idx,idy,idz,idk)*result;

}

//count all six staple for a selected link
//sixstaple init.ed as Zero matrix
void Modell::Count6Staple(const unsigned int grididx,int idx,int idy,int idz,int idk){
    int grididx2=grididx;
    su3staple.zeros();
    cx_mat result(3,3,fill::eye);
    for(int i=1;i<4;i++){
        grididx2=grididx+i;
        grididx2=(grididx2+4)%4;
        result.eye();
        TriplUForw(grididx,grididx2,idx,idy,idz,idk,result);
        su3staple+=result;
        result.eye();
        TriplURev(grididx,grididx2,idx,idy,idz,idk,result);
        su3staple+=result;
    }
}

//find Pauli coefficients for 2x2 submatrix of su3staple
//strowcol: first row and column idx of submatrix
//negative sign because different way for building SU2 staple matrix as opposed to other su2 matrix
vector<double> Modell::CountCoeffs(const int strowcol){

vector<double> coeffs={-real((1./2)*iunit*(su3staple(strowcol,strowcol+1)+su3staple(strowcol+1,strowcol))),
                       -real((1./2)*(su3staple(strowcol+1,strowcol)-su3staple(strowcol,strowcol+1))),
                       -real((1./2)*iunit*(su3staple(strowcol,strowcol)-su3staple(strowcol+1,strowcol+1)))};
    return coeffs;
}

//generating new su2 matrix coeffs
double Modell::GenerateCoeff0(){
    double a0=GetRealRandom(exp(-2.*real(su2strootdet)),1);
//debug
cout<<su2strootdet<<endl;
cout<<"is it real? Equals this? "<<real(su2strootdet)<<endl;
    //generate a_0
    //with accept-reject
    bool accept=Flip(sqrt(1-a0*a0));
//debug
int counter=1;
    while(!accept){
        a0=GetRealRandom(exp(-2.*real(su2strootdet)),1);
        accept=Flip(sqrt(1-a0*a0));
        //debug
        counter++;
        //debug
        cout<<"num of trials: "<<counter<<endl;
    }
return a0;
}

//generate new su2 matrix coeffs 3d sphere
std::vector<double> Modell::GenerateCoeffs(){
    vector<double> coeff3d=RandOnSphere(3);
    return coeff3d;
}

// build SU2  matrix
//su2: result matrix, initialized as zeros
void Modell::BuildSU2(const double coeff0, const vector<double> & coeffs,cx_mat & su2){

//reinitialize su2staple
//su2staple.zeros();
su2.zeros(); //just in case
if (coeffs.empty()) throw "BuildSU2staple fails! empty coeffs";
if(coeffs.size()!=3) throw "Build SU2staple fails! wrong coeffs size";
//build su2 staple like matrix
su2+=coeff0*identity2;
su2+=iunit*coeffs.front()*pauli1;
su2+=iunit*coeffs[1]*pauli2;
su2+=iunit*coeffs.back()*pauli3;
}

//su2staple
//and determinant
void Modell::BuildSU2staple(const int strowcol){
    BuildSU2(real(1./2*(su3staple(strowcol,strowcol)+su3staple(strowcol+1,strowcol+1))),
    CountCoeffs(strowcol),su2staple);
    su2strootdet=sqrt(det(su2staple));
}


//Modify the selected link
void Modell::ModifyLink(int grididx, int idx, int idy, int idz, int idk, const arma::cx_mat& newlink){
//debug
cout<<"modell modifylink call"<<endl;
    grid(grididx).ModifyGrid()(idx,idy,idz,idk)=newlink;
}

//access for reading
const Array::array1<SU3Grid>& Modell::GetModellGrid()const{
    return grid;
}

Modell::~Modell(){}

/****************************************************************/

const double Modell::beta=6;
const int SU3Grid::dim=4;
const int SU3Grid::tdim=8;
const std::complex<double> Modell::iunit(0,1);
const arma::cx_mat Modell::pauli1={{{0,0},{1,0}},{{1,0},{0,0}}};
const arma::cx_mat Modell::pauli2={{{0,0},{0,1}},{{0,-1},{0,0}}};
const arma::cx_mat Modell::pauli3={{{1,0},{0,0}},{{0,0},{-1,0}}};
const arma::cx_mat Modell::identity2(2,2,fill::eye);



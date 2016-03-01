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
//debug
cout<<"original det: "<<determinant<<endl;
    su3=su3/determinant;
//debug
cout<<"det: "<<det(su3)<<endl;
    return su3;
}

/*******************/
/******SU3GRID******/
/*******************/

SU3Grid::SU3Grid():ei(tdim,dim,dim,dim){
//debug
//cout<<"su3grid def. ctr. call"<<endl;
}

SU3Grid::SU3Grid(const cx_mat & initmatrix):ei(tdim,dim,dim,dim){
    for(int i=0;i<tdim;i++){
        for(int j=0;j<dim;j++){
            for(int k=0;k<dim;k++){
                for(int l=0;l<dim;l++){
                    cx_mat su3(3,3);
                    ei(i,j,k,l)=initmatrix;
                }//for l spacelike
            }//for k spacelike
        }//for j spacelike
    }//for i timelike
}

//construct from random
SU3Grid::SU3Grid(bool flag):ei(tdim,dim,dim,dim){
//debug
//cout<<"su3grid ctr from rand call"<<endl;
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
//cout<<"su3grid assignment call"<<endl;
    if(this!=&fromsu3grid){
        ei=fromsu3grid.ei;
    }
    return *this;
}

//copy
SU3Grid::SU3Grid(const SU3Grid& su3grid):ei(tdim,dim,dim,dim){
//debug
//cout<<"su3grid copy call"<<endl;
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
//debug
//cout<<"su3 grid modifygrid call"<<endl;
    return ei;
}

SU3Grid::~SU3Grid(){}
/*******************************************************/


/*********************/
/***** MODELL ********/
/*********************/

//default constr
Modell::Modell():grid(4),su3staple(3,3,fill::zeros),su2staple(2,2,fill::zeros),su2strootdet(0){}

//same matrix for all element
Modell::Modell(cx_mat & initmatrix):grid(4),su3staple(3,3,fill::zeros),su2staple(2,2,fill::zeros),su2strootdet(0){
//debug
cout<<"Modell ctr with same matrix"<<endl;
    MatrixInit(initmatrix);
}

//construct from random
Modell::Modell(bool flag):grid(4),su3staple(3,3,fill::zeros),su2staple(2,2,fill::zeros),su2strootdet(0){
//debug
cout<<"Modell ctr from rand"<<endl;
    RandomInit();
}

//TO DO
    //construct from file
    //Modell(const char* );
    //SU3Grid(std::istream&);

//Matrix Init
void Modell::MatrixInit(arma::cx_mat & initmatrix){
 //debug
cout<<"Modell matrixinit call"<<endl;
    int tgridmax=SU3Grid::GetTDim();
    int gridmax=SU3Grid::GetDim();
    //sweep cycle
        for(int ei=0;ei<4;ei++){
            for(int i=0;i<tgridmax;i++){
                for(int j=0;j<gridmax;j++){
                    for(int k=0;k<gridmax;k++){
                        for(int l=0;l<gridmax;l++){
                            grid(ei).ModifyGrid()(i,j,k,l)=initmatrix;
                        }//for l
                    }
                }//for j
            }//for timelike
        }//for grid
}

//RandomInit
void Modell::RandomInit(){
//debug
cout<<"Modell randominit call"<<endl;
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
Modell& Modell::operator=(const Modell& frommodell){
    if(this!=&frommodell){
        //debug
        cout<<"modell assignment call"<<endl;
        grid=frommodell.grid;
        su3staple=frommodell.su3staple;
        su2staple=frommodell.su2staple;
        su2strootdet=frommodell.su2strootdet;
    }
return *this;

}

//COPY
Modell::Modell(const Modell& mymodell):grid(4){
    *this=mymodell;
}

//get beta
const double Modell::GetBeta(){
    return beta;
}

    //write Modell to os
    //void writeModell( std::ostream& )const;

    //write Modell to file
    //void writeToFileModell( const char* )const;

    //read Modell from is
    //void readModell( std::istream& );

    //read Modell from file
    //void readFromFileModell( const char* );



//count forward staple (plaquett without the selected link)
//result initialized as Identity
//OR: HELPER FOR PLAQUETT COUNTER - INITED AS THE SELECTED LINK
void Modell::TriplUForw(const unsigned int grididx, const unsigned int grididx2,
int idx, int idy, int idz, int idk, arma::cx_mat & result){
//debug
//cout<<"Modell tripluforw call"<<endl;
//cout<<"grids: "<<grididx<<", "<<grididx2<<endl;
//cout<<"idxs: "<<idx<<", "<<idy<<", "<<idz<<", "<<idk<<endl;
//cout<<"grid link matrix: "<<(grid(grididx2).GetGrid())(idx,idy,idz,idk)<<endl;
//cout<<"result: "<<result<<endl;
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
//debug
//cout<<"grids: "<<grididx<<", "<<grididx2<<endl;
//cout<<"idxs: "<<idx<<", "<<idy<<", "<<idz<<", "<<idk<<endl;
//cout<<"grid link matrix: "<<(grid(grididx2).GetGrid())(idx,idy,idz,idk)<<endl;
        //multipl by second edge (j) at index+ei*a
        result=(grid(grididx2).GetGrid())(idx,idy,idz,idk)*result;
//cout<<"result: "<<result<<endl;
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
//debug
//cout<<"grids: "<<grididx<<", "<<grididx2<<endl;
//cout<<"idxs: "<<idx<<", "<<idy<<", "<<idz<<", "<<idk<<endl;
//cout<<"grid link matrix: "<<(grid(grididx).GetGrid())(idx,idy,idz,idk)<<endl;
//cout<<"grid link matrix herm.adj: "<<(grid(grididx).GetGrid())(idx,idy,idz,idk).t()<<endl;
        //multipl by third edge at index+ej*a
        result=(grid(grididx).GetGrid())(idx,idy,idz,idk).t()*result;
//cout<<"result: "<<result<<endl;
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
//debug
//cout<<"grids: "<<grididx<<", "<<grididx2<<endl;
//cout<<"idxs: "<<idx<<", "<<idy<<", "<<idz<<", "<<idk<<endl;
//cout<<"grid link matrix: "<<(grid(grididx2).GetGrid())(idx,idy,idz,idk)<<endl;
//cout<<"grid link matrix herm.adj: "<<(grid(grididx2).GetGrid())(idx,idy,idz,idk).t()<<endl;
        //multipl by fourth edge at index - closing the cycle
        result=(grid(grididx2).GetGrid())(idx,idy,idz,idk).t()*result;
//cout<<"result: "<<result<<endl;
}

//count plaquett at idx, for two selected direction i,j
//result initialized as Identity
//assumes calls with existing indices idx-y-z-k and i-j
void Modell::CountUp(int idx,int idy,int idz,int idk,const unsigned int grididx1, const unsigned int grididx2,arma::cx_mat & result){
    //result should be initialized as Identity
    result.eye();
    //result is the selected link
    result=(grid(grididx1).GetGrid())(idx,idy,idz,idk)*result;
    //multiplication with the other 3 links
    TriplUForw(grididx1,grididx2,idx,idy,idz,idk,result);

}

//count action for a plaquett up
 double Modell::CountPlaqEnergy(const cx_mat & plaquett){
        double spup;
        spup=(1-1/3.0*real(trace(plaquett)));
        return spup;
}

//count mean for plaquett energy on lattice
double Modell::CountMeanEnergyDens(ofstream & file){
    const int timedim=SU3Grid::GetTDim();
    const int spacedim=SU3Grid::GetDim();
    double meanEDens=0;
    cx_mat plaquett(3,3,fill::eye);
    int counter=0;
    //evergy plaquett counted 4 times!!!
    //cycle on selected links, select all links
    for(int grid=0;grid<4;grid++){
     for(int grid2=grid+1;grid2<4;grid2++){
        for(int i=0;i<timedim;i++){
            for(int j=0;j<spacedim;j++){
                for(int k=0;k<spacedim;k++){
                    for(int l=0;l<spacedim;l++){
                        CountUp(i,j,k,l,grid,grid2,plaquett);
                        file<<grid<<'\t'<<grid2<<'\t'<<i<<'\t'<<j<<'\t'<<k<<'\t'<<l<<'\t'<<CountPlaqEnergy(plaquett)<<endl;
                        meanEDens+=CountPlaqEnergy(plaquett);
                        counter++;
                    }//for l
                }//for k
            }//for j
        }//for i
     }//grid2
    }//for grid

    meanEDens/=counter;
    return meanEDens;
}

//count backward staple (plaquett without the selected link)
//result initialized as Identity
//except when count plaquett - in this case it is the actual link variable
void Modell::TriplURev(const unsigned int grididx, const unsigned int grididx2,
int idx, int idy, int idz, int idk, arma::cx_mat & result){
//debug
//cout<<"Modell triplurev call"<<endl;
//just in case - reinit.
result.eye();
//cout<<"grids: "<<grididx<<", "<<grididx2<<endl;
//cout<<"idxs: "<<idx<<", "<<idy<<", "<<idz<<", "<<idk<<endl;
//cout<<"grid link matrix: "<<(grid(grididx2).GetGrid())(idx,idy,idz,idk)<<endl;
//cout<<"result: "<<result<<endl;
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
//cout<<"grids: "<<grididx<<", "<<grididx2<<endl;
//cout<<"idxs: "<<idx<<", "<<idy<<", "<<idz<<", "<<idk<<endl;
//cout<<"grid link matrix: "<<(grid(grididx2).GetGrid())(idx,idy,idz,idk)<<endl;
        result=(grid(grididx2).GetGrid())(idx,idy,idz,idk).t()*result;
//cout<<"result: "<<result<<endl;
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
//cout<<"grids: "<<grididx<<", "<<grididx2<<endl;
//cout<<"idxs: "<<idx<<", "<<idy<<", "<<idz<<", "<<idk<<endl;
//cout<<"grid link matrix: "<<(grid(grididx2).GetGrid())(idx,idy,idz,idk)<<endl;
        result=(grid(grididx).GetGrid())(idx,idy,idz,idk).t()*result;
//cout<<"result: "<<result<<endl;

        result=(grid(grididx2).GetGrid())(idx,idy,idz,idk)*result;
//cout<<"result: "<<result<<endl;
}

//count all six staple for a selected link
//sixstaple init.ed as Zero matrix
void Modell::Count6Staple(const unsigned int grididx,int idx,int idy,int idz,int idk){
//debug
//cout<<"Modell count6staple call"<<endl;
    int grididx2=grididx;
    su3staple.zeros();
    cx_mat result(3,3,fill::eye);
//debug
//cout<<"su3staple: "<<su3staple<<endl;
//cout<<"result "<<result<<endl;
    for(int i=1;i<4;i++){
        grididx2=grididx+i;
        grididx2=(grididx2+4)%4;
        result.eye();
//cout<<"idxs: "<<idx<<", "<<idy<<", "<<idz<<", "<<idk<<endl;
        TriplUForw(grididx,grididx2,idx,idy,idz,idk,result);
        su3staple+=result;
//cout<<"result "<<result<<endl;
//cout<<"su3staple: "<<su3staple<<endl;
        result.eye();
        //debug
        //cout<<result<<endl;
//cout<<"idxs: "<<idx<<", "<<idy<<", "<<idz<<", "<<idk<<endl;
        TriplURev(grididx,grididx2,idx,idy,idz,idk,result);
        su3staple+=result;
//cout<<"result "<<result<<endl;
//cout<<"su3staple: "<<su3staple<<endl;
        //debug
        //cout<<"su3staple"<<su3staple<<endl;
    }
}

//find Pauli coefficients for 2x2 submatrix of su3staple
//strowcol: first row and column idx of submatrix
//negative sign because different way for building SU2 staple matrix as opposed to other su2 matrix
vector<double> Modell::CountCoeffs(const int strowcol,const cx_mat & ustaple){
//debug
//cout<<"Modell countcoeffs call"<<endl;
//cout<<"iunit: "<<iunit<<endl;
//cout<<"coeffs: "<<-real((1./2)*iunit*(su3staple(strowcol,strowcol+1)+su3staple(strowcol+1,strowcol)))<<", "
 //               <<-real((1./2)*(su3staple(strowcol+1,strowcol)-su3staple(strowcol,strowcol+1)))<<", "
 //               <<-real((1./2)*iunit*(su3staple(strowcol,strowcol)-su3staple(strowcol+1,strowcol+1)))<<endl;
vector<double> coeffs={-real((1./2)*iunit*(ustaple(strowcol,strowcol+1)+ustaple(strowcol+1,strowcol))),
                       -real((1./2)*(ustaple(strowcol+1,strowcol)-ustaple(strowcol,strowcol+1))),
                       -real((1./2)*iunit*(ustaple(strowcol,strowcol)-ustaple(strowcol+1,strowcol+1)))};
//cout<<"coeffs vector: "<<coeffs.front()<<", "<<coeffs[1]<<", "<<coeffs.back()<<endl;
//debug
//cout<<"coeff1"<<coeffs.front()<<endl;
//cout<<"coeff2"<<coeffs[1]<<endl;
//cout<<"coeff3"<<coeffs.back()<<endl;
    return coeffs;
}

//generating new su2 matrix coeffs
double Modell::GenerateCoeff0(){
//debug
//cout<<"Modell gen.coeff0 call"<<endl;
//cout<<"su2strootdet: "<<real(su2strootdet)<<endl;

double a0=GetRealRandom(exp(-2.*Modell::beta*real(su2strootdet)*2./3),1);
    a0=1+1./(Modell::beta*2./3*real(su2strootdet))*log(a0);
//cout<<"a0: "<<a0<<endl;

//debug
//cout<<su2strootdet<<endl;
//cout<<"is it real? Equals this? "<<real(su2strootdet)<<endl;
    //generate a_0
    //with accept-reject
    bool accept=Flip(sqrt(1-a0*a0));
//cout<<"flip? "<<accept<<endl;
//debug
int counter=1;
//debug
//cout<<"a0 "<<a0<<"success: "<<sqrt(1-a0*a0)<<endl;
//debug
//cout<<"su2strootdet"<<real(su2strootdet)<<"lower lim: "<<exp(-2.*real(su2strootdet))<<endl;
    while(!accept){
        a0=GetRealRandom(exp(-2.*Modell::beta*real(su2strootdet)*2./3),1);
//cout<<"su2strootdet: "<<real(su2strootdet)<<endl;
        a0=1+1./(Modell::beta*2./3*real(su2strootdet))*log(a0);
//cout<<"a0: "<<a0<<endl;
        accept=Flip(sqrt(1-a0*a0));
//cout<<"flip? "<<accept<<endl;
//debug
counter++;
//debug
//cout<<"num of trials: "<<counter<<endl;
    }
//debug
cout<<"num of trials: "<<counter<<" a0: "<<a0<<endl;
return a0;
}

//generate new su2 matrix coeffs 3d sphere
std::vector<double> Modell::GenerateCoeffs(double a0){
//debug modell
//cout<<"Modell gen.coeffs call"<<endl;
    vector<double> coeff3d=RandOnSphere(3);
    double factor=sqrt(1-a0*a0);
//cout<<"generated vector: "<<coeff3d.front()<<", "<<coeff3d[1]<<", "<<coeff3d.back();
//cout<<"factor: "<<factor<<endl;
    coeff3d[0]*=factor;
    coeff3d[1]*=factor;
    coeff3d[2]*=factor;
//cout<<"generated vector: "<<coeff3d.front()<<", "<<coeff3d[1]<<", "<<coeff3d.back();
    return coeff3d;
}

// build SU2  matrix
//su2: result matrix, initialized as zeros
void Modell::BuildSU2(const double coeff0, const vector<double> & coeffs,cx_mat & su2){
//debug
//cout<<"Modell build su2 call"<<endl;
//reinitialize su2staple
//su2staple.zeros();
su2.zeros(); //just in case
//cout<<"su2 to build: "<<su2<<endl;
if (coeffs.empty()) throw "BuildSU2staple fails! empty coeffs";
if(coeffs.size()!=3) throw "Build SU2staple fails! wrong coeffs size";

//debug
//cout<<"coeff0: "<<coeff0<<endl;
//debug
//cout<<"coeffs 0(1)"<<coeffs.front()<<endl;
//debug
//cout<<"coeffs 1(2)"<<coeffs[1]<<endl;
//debug
//cout<<"coeffs 2(3)"<<coeffs.back()<<endl;
//debug
//cout<<"equals? "<<coeffs[2]<<endl;

//build su2 staple like matrix
su2+=coeff0*identity2;
//debug
//cout<<su2<<endl;
su2+=iunit*coeffs.front()*pauli1;
//cout<<su2<<endl;
su2+=iunit*coeffs[1]*pauli2;
//cout<<su2<<endl;
su2+=iunit*coeffs.back()*pauli3;
//cout<<su2<<endl;
//debug
//cout<<"determinant of build: "<<det(su2)<<endl;
}

//su2staple
//and determinant
void Modell::BuildSU2staple(const int strowcol,const unsigned int grididx,const int idx,const int idy,const int idz,const int idk){
//debug
//cout<<"Modell buildsu2staple call"<<endl;
//cout<<"strowcol: "<<strowcol<<endl;
//cout<<"su2staple: "<<su2staple<<endl;

//in the SU3 case instead of S_l one need U_l*S_l
cx_mat ustaple=(grid(grididx).GetGrid())(idx,idy,idz,idk)*su3staple;

    BuildSU2(real(1./2*(ustaple(strowcol,strowcol)+ustaple(strowcol+1,strowcol+1))),
    CountCoeffs(strowcol,ustaple),su2staple);
    su2strootdet=sqrt(det(su2staple));
//debug
//cout<<"su2staple"<<su2staple<<endl;
//cout<<"rootdet: "<<su2strootdet<<endl;
//cout<<"rootdet*invstaple"<<su2strootdet*inv(su2staple)<<endl;
//cout<<"determinant of this MUST equal 1: "<<det(su2strootdet*inv(su2staple))<<endl;

}

//build generated and transformed su3 matrix
//refresh link with it
void Modell::RefreshLinkpart(const int grididx, const int idx, const int idy, const int idz, const int idk, const int strowcol){
//debug
//cout<<"Modell refreshlinkpart call"<<endl;
    //count su2 from generated coeffs
    cx_mat alphamat(2,2,fill::zeros);
    double a0=GenerateCoeff0();
    BuildSU2(a0,GenerateCoeffs(a0),alphamat);
//debug
//cout<<"generated su2 det "<<det(alphamat)<<endl;
//cout<<"generated su2: "<<alphamat<<endl;
//cout<<"su2staple: "<<su2staple<<endl;
//cout<<"su2staple inv: "<<su2staple.i()<<endl;
//cout<<"is this Id? "<<su2staple*su2staple.i()<<endl;
    //transfom alpha
    alphamat=alphamat*su2strootdet*su2staple.i();


//debug
//cout<<"generated transformed su2 det "<<det(alphamat)<<endl;
//cout<<"gen.d transf.d. su2: "<<alphamat<<endl;
    //build refresher matrix
    cx_mat refresher(3,3,fill::eye);
    refresher(strowcol,strowcol)=alphamat(0,0);
    refresher(strowcol+1,strowcol)=alphamat(1,0);
    refresher(strowcol,strowcol+1)=alphamat(0,1);
    refresher(strowcol+1,strowcol+1)=alphamat(1,1);
//debug
//cout<<"alphamat: "<<alphamat<<endl;
//debug
//cout<<"refresher"<<refresher<<endl;
//cout<<"refresher det"<<det(refresher)<<endl;

//debug - original det
//cout<<"orig. det: "<<det(grid(grididx).GetGrid()(idx,idy,idz,idk))<<endl;
//cout<<"orig matrix "<<grid(grididx).GetGrid()(idx,idy,idz,idk)<<endl;
//    //modify link
    grid(grididx).ModifyGrid()(idx,idy,idz,idk)=refresher*grid(grididx).GetGrid()(idx,idy,idz,idk);
//debug - new det
cout<<"new det: "<<det(grid(grididx).GetGrid()(idx,idy,idz,idk))<<endl;
//cout<<"new matrix "<<grid(grididx).GetGrid()(idx,idy,idz,idk)<<endl;

}

//Modify the selected link
void Modell::ModifyLink(int grididx, int idx, int idy, int idz, int idk, const arma::cx_mat& newlink){
//debug
//cout<<"modell modifylink call"<<endl;
    grid(grididx).ModifyGrid()(idx,idy,idz,idk)=newlink;
}

//access for reading
const Array::array1<SU3Grid>& Modell::GetModellGrid()const{
//debug
//cout<<"Modell getgrid call"<<endl;
    return grid;
}

//heat bath step
void Modell::HeatBathStep(const unsigned int grididx,int idx,int idy, int idz, int idk){
//debug
//cout<<"Modell heatbathstep call"<<endl;
//cout<<"su3staple: "<<su3staple<<endl;
    //count staple
    Count6Staple(grididx,idx,idy,idz,idk);
    //for: fist and second sub SU2
    for(int strawcol=0;strawcol<2;strawcol++){
//cout<<"strowcol: "<<strawcol<<endl;
//cout<<"su2staple (old): "<<su2staple<<endl;
//cout<<"su3staple: "<<su3staple<<endl;

        //build su2 analog matrix
        BuildSU2staple(strawcol,grididx,idx,idy,idz,idk);
//cout<<"su2staple new: "<<su2staple<<endl;
        //generate coeffs and refresh
        RefreshLinkpart(grididx,idx,idy,idz,idk,strawcol);
    }
}

//heat bath sweep
void Modell::HeatBathSweep(){
//debug
//cout<<"Modell heatbathsweep call"<<endl;
    const int tdim=SU3Grid::GetTDim();
    const int dim=SU3Grid::GetDim();
    for(int grididx=0;grididx<4;grididx++){
        for(int i=0;i<tdim;i++){
            for(int j=0;j<dim;j++){
                for(int k=0;k<dim;k++){
                    for(int l=0;l<dim;l++){
                        HeatBathStep(grididx,i,j,k,l);
cout<<"*********************"<<endl;
cout<<"****IS IT UNITARY?***"<<endl;
cout<<"*********************"<<endl;
cout<<grid(grididx).GetGrid()(i,j,k,l)*grid(grididx).GetGrid()(i,j,k,l).t()<<endl;
//cout<<grid(grididx).GetGrid()(i,j,k,l).t()<<endl;
//cout<<grid(grididx).GetGrid()(i,j,k,l).i()<<endl;
                    }//for space 3
                }//for space 2
            }//for space1
        }//for time
    }//for grid
}

//Polyakov loop
void Modell::PolyakovMatrix(const int x,const int y,const int z,arma::cx_mat & result){
    result.eye(); //reinitialise just in case and actually it is needed
    const int maxtdim=SU3Grid::GetTDim();
    for(int i=0;i<maxtdim;i++){
            result=(grid(0).GetGrid())(i,x,y,z)*result;
        }
}

//Polyakov space avg
complex<double> Modell::PolyakovLoopAVG(std::ofstream & polyaloop){
    const int maxdim=SU3Grid::GetDim();



//debug
complex<double> cmplavg(0,0);
    cx_mat polyamatrix(3,3,fill::eye);
    //double polyaavg=0;
    for(int i=0;i<maxdim;i++){
        for(int j=0;j<maxdim;j++){
            for(int k=0;k<maxdim;k++){
                PolyakovMatrix(i,j,k,polyamatrix);
                //polyakov loop trace is not real
                //polyaavg+=real(trace(polyamatrix));
                polyaloop<<i<<'\t'<<j<<'\t'<<k<<'\t'<<trace(polyamatrix)<<'\t'<<real(trace(polyamatrix))<<'\t'<<imag(trace(polyamatrix))<<endl;
                cmplavg+=trace(polyamatrix);

            }//for k
        }//for j
    }//for i
    //polyaavg/=(maxdim*maxdim*maxdim);

cmplavg/=(maxdim*maxdim*maxdim);
//debug
cout<<cmplavg<<endl;
return cmplavg;
}

    //debug
    void Modell::GetPauli(){
        cout<<"Identity"<<endl;
        cout<<Modell::identity2<<endl;
        cout<<"Pauli1"<<endl;
        cout<<Modell::pauli1<<endl;
        cout<<"Pauli2"<<endl;
        cout<<Modell::pauli2<<endl;
        cout<<"Pauli3"<<endl;
        cout<<Modell::pauli3<<endl;
    }

Modell::~Modell(){}

/****************************************************************/

const double Modell::beta=6;
const int SU3Grid::dim=10;
const int SU3Grid::tdim=10;
const std::complex<double> Modell::iunit(0,1);
const arma::cx_mat Modell::pauli1={{{0,0},{1,0}},{{1,0},{0,0}}};
const arma::cx_mat Modell::pauli2={{{0,0},{0,-1}},{{0,1},{0,0}}};
const arma::cx_mat Modell::pauli3={{{1,0},{0,0}},{{0,0},{-1,0}}};
const arma::cx_mat Modell::identity2(2,2,fill::eye);



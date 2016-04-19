#include<cmath>
#include<armadillo>
#include<complex>
#include<string>
#include"Sommer.h"

void ScaleSetV::InitSpaceLikeT(const int time,Array::array1<arma::cx_mat> & matarray){
const int maxdim=SU3Grid::GetDim();
const int maxtdim=SU3Grid::GetTDim();
for(int i=0;i<R;i++){


    switch(spacegrididx){
        case 1:
            {int idx=(initx+i+maxdim)%maxdim;
            matarray(i)=this->mymodell.GetModellGrid()(spacegrididx).GetGrid()((initt+time+maxtdim)%maxtdim,idx,inity,initz);}
            break;
        case 2:
            {int idy=(inity+i+maxdim)%maxdim;
            matarray=this->mymodell.GetModellGrid()(spacegrididx).GetGrid()((initt+time+maxtdim)%maxtdim,initx,idy,initz);}
            break;
        case 3:
            {int idz=(initz+i+maxdim)%maxdim;
            matarray=this->mymodell.GetModellGrid()(spacegrididx).GetGrid()((initt+time+maxtdim)%maxtdim,initx,inity,idz);}
            break;
    }//switch


}//for

}

void ScaleSetV::InitSpaceLikeTInv(const int time,Array::array1<arma::cx_mat> & matarray){
const int maxdim=SU3Grid::GetDim();
const int maxtdim=SU3Grid::GetTDim();
for(int i=0;i<R;i++){


    switch(spacegrididx){
        case 1:
            {int idx=(initx+i+maxdim)%maxdim;
            matarray(i)=this->mymodell.GetModellGrid()(spacegrididx).GetGrid()((initt+time+maxtdim)%maxtdim,idx,inity,initz).t();}
            break;
        case 2:
            {int idy=(inity+i+maxdim)%maxdim;
            matarray=this->mymodell.GetModellGrid()(spacegrididx).GetGrid()((initt+time+maxtdim)%maxtdim,initx,idy,initz).t();}
            break;
        case 3:
            {int idz=(initz+i+maxdim)%maxdim;
            matarray=this->mymodell.GetModellGrid()(spacegrididx).GetGrid()((initt+time+maxtdim)%maxtdim,initx,inity,idz).t();}
            break;
    }//switch


}//for

}

//CTR
ScaleSetV::ScaleSetV(Modell & mymodell,
const int r,const int t,const int tidx,
const int xidx,const int yidx,const int zidx,const int grididx):
mymodell(mymodell),R(r),T(t),oinitt(tidx),oinitx(xidx),oinity(yidx),
oinitz(zidx),spacegrididx(grididx),alpha(0.5),maxsmearlevel(2),spacelike_0(R),spacelike_T(R),
correlT(maxsmearlevel,maxsmearlevel,arma::fill::zeros),correlT1(maxsmearlevel,maxsmearlevel,arma::fill::zeros),
correl0(maxsmearlevel,maxsmearlevel,arma::fill::zeros){
initt=tidx;
initx=xidx;
inity=yidx;
initz=zidx;
InitSpaceLikeT(0,spacelike_0);
InitSpaceLikeTInv(T,spacelike_T);

}

//smear pt1 helper
void ScaleSetV::TriplForw(arma::cx_mat & produc,int actidx,int actidy,int actidz,const int timeidx,const int spacidx){
    produc.eye();
        const int maxdim=SU3Grid::GetDim();
    if(spacidx==spacegrididx) throw "spacidx equals wilson space direction";
        produc=mymodell.GetModellGrid()(spacidx).GetGrid()(timeidx,actidx,actidy,actidz)*produc;
    //step along second direction
        switch(spacidx){
        case 1:
            actidx++;
            actidx=(actidx+maxdim)%maxdim;
            break;
        case 2:
            actidy++;
            actidy=(actidy+maxdim)%maxdim;
            break;
        case 3:
            actidz++;
            actidz=(actidz+maxdim)%maxdim;
            break;
        }//switch
        produc=mymodell.GetModellGrid()(spacegrididx).GetGrid()(timeidx,actidx,actidy,actidz)*produc;
        //step back along second direction
        switch(spacidx){
        case 1:
            actidx--;
            actidx=(actidx+maxdim)%maxdim;
            break;
        case 2:
            actidy--;
            actidy=(actidy+maxdim)%maxdim;
            break;
        case 3:
            actidz--;
            actidz=(actidz+maxdim)%maxdim;
            break;
        }//switch
        //step along first direction
        switch(spacegrididx){
        case 1:
            actidx++;
            actidx=(actidx+maxdim)%maxdim;
            break;
        case 2:
            actidy++;
            actidy=(actidy+maxdim)%maxdim;
            break;
        case 3:
            actidz++;
            actidz=(actidz+maxdim)%maxdim;
            break;
        }//switch
        produc=mymodell.GetModellGrid()(spacidx).GetGrid()(timeidx,actidx,actidy,actidz).t()*produc;
//restore indices to original
        //step back along first direction
        switch(spacegrididx){
        case 1:
            actidx--;
            actidx=(actidx+maxdim)%maxdim;
            break;
        case 2:
            actidy--;
            actidy=(actidy+maxdim)%maxdim;
            break;
        case 3:
            actidz--;
            actidz=(actidz+maxdim)%maxdim;
            break;
        }//switch
}

void ScaleSetV::TriplBack(arma::cx_mat & produc,int actidx,int actidy,int actidz,const int timeidx,const int spacidx){
produc.eye();
    const int maxdim=SU3Grid::GetDim();
if(spacidx==spacegrididx) throw "spacidx equals wilson space direction";

        //step back along second direction
        switch(spacidx){
        case 1:
            actidx--;
            actidx=(actidx+maxdim)%maxdim;
            break;
        case 2:
            actidy--;
            actidy=(actidy+maxdim)%maxdim;
            break;
        case 3:
            actidz--;
            actidz=(actidz+maxdim)%maxdim;
            break;
        }//switch
            produc=mymodell.GetModellGrid()(spacidx).GetGrid()(timeidx,actidx,actidy,actidz).t()*produc;
            produc=mymodell.GetModellGrid()(spacegrididx).GetGrid()(timeidx,actidx,actidy,actidz)*produc;
        //step along first direction
        switch(spacegrididx){
        case 1:
            actidx++;
            actidx=(actidx+maxdim)%maxdim;
            break;
        case 2:
            actidy++;
            actidy=(actidy+maxdim)%maxdim;
            break;
        case 3:
            actidz++;
            actidz=(actidz+maxdim)%maxdim;
            break;
        }//switch
        produc=mymodell.GetModellGrid()(spacidx).GetGrid()(timeidx,actidx,actidy,actidz)*produc;
        //set back original indices
                //step along second direction
        switch(spacidx){
        case 1:
            actidx++;
            actidx=(actidx+maxdim)%maxdim;
            break;
        case 2:
            actidy++;
            actidy=(actidy+maxdim)%maxdim;
            break;
        case 3:
            actidz++;
            actidz=(actidz+maxdim)%maxdim;
            break;
        }//switch
                //step back along first direction
        switch(spacegrididx){
        case 1:
            actidx--;
            actidx=(actidx+maxdim)%maxdim;
            break;
        case 2:
            actidy--;
            actidy=(actidy+maxdim)%maxdim;
            break;
        case 3:
            actidz--;
            actidz=(actidz+maxdim)%maxdim;
            break;
        }//switch
}

//for one link
void ScaleSetV::SmearPt1Forw(arma::cx_mat & link,int actidx,int actidy,int actidz,const int timeidx){
    arma::cx_mat produc(3,3,arma::fill::eye);
    //for space dim.s
    for(int i=1;i<4;i++){
        produc.eye();
        if(i!=spacegrididx){
           TriplForw(produc,actidx,actidy,actidz,timeidx,i);
            link+=alpha*produc;
            TriplBack(produc,actidx,actidy,actidz,timeidx,i);
            link+=alpha*produc;
        }//if
        else{}
    }//for spacedim.s

}

//for one link
void ScaleSetV::SmearPt1Back(arma::cx_mat & link,int actidx,int actidy,int actidz,const int timeidx){
    arma::cx_mat produc(3,3,arma::fill::eye);
    //for space dim.s
    for(int i=1;i<4;i++){
        produc.eye();
        if(i!=spacegrididx){
           TriplForw(produc,actidx,actidy,actidz,timeidx,i);
            link+=(alpha*produc).i();
            TriplBack(produc,actidx,actidy,actidz,timeidx,i);
            link+=(alpha*produc).i();
        }//if
        else{}
    }//for spacedim.s

}

//project back to su3
void ScaleSetV::SmearPt2(arma::cx_mat & link){
    arma::cx_mat diagonal(3,3,arma::fill::eye);
    arma::cx_mat eigvecs;
    arma::vec eigvals;
    arma::eig_sym(eigvals,eigvecs,link.t()*link,"std");
    diagonal(0,0)=eigvals(0)*diagonal(0,0);
    diagonal(1,1)=eigvals(1)*diagonal(1,1);
    diagonal(2,2)=eigvals(2)*diagonal(2,2);
//    std::cout<<"eigvals"<<std::endl;
//    std::cout<<eigvals<<std::endl;
//    std::cout<<"eigvecs"<<std::endl;
//    std::cout<<eigvecs<<std::endl;
//    std::cout<<"diagonal"<<std::endl;
//    std::cout<<diagonal<<std::endl;
    arma::cx_mat sqrtmatrix=eigvecs*sqrt(diagonal)*eigvecs.i();
//    std::cout<<"original link"<<std::endl;
//    std::cout<<link<<std::endl;
    link=link*sqrtmatrix.i()*pow(arma::det(link.i()*link.t()),1./6);
//    std::cout<<"new link"<<std::endl;
//    std::cout<<link<<std::endl;
//    std::cout<<"is it identity?"<<std::endl;
//    std::cout<<link.t()*link<<std::endl;
//    std::cout<<arma::det(link)<<std::endl;
}

void ScaleSetV::Smearing0(){
int maxdim=SU3Grid::GetDim();
    for(int i=0;i<R;i++){
        //smearing a link
        switch(spacegrididx){
        case 1:
            SmearPt1Forw(spacelike_0(i),(initx+i+maxdim)%maxdim,inity,initz,initt);
            SmearPt2(spacelike_0(i));
            break;
        case 2:
            SmearPt1Forw(spacelike_0(i),initx,(inity+i+maxdim)%maxdim,initz,initt);
            SmearPt2(spacelike_0(i));
            break;
        case 3:
            SmearPt1Forw(spacelike_0(i),initx,inity,(initz+i+maxdim)%maxdim,initt);
            SmearPt2(spacelike_0(i));
            break;
        }//switch
    }
}

void ScaleSetV::SmearingT(const int tshift){
int maxdim=SU3Grid::GetDim();
int tidx=(initt+tshift+SU3Grid::GetTDim())%SU3Grid::GetTDim();
    for(int i=0;i<R;i++){
        //smearing a link
        switch(spacegrididx){
        case 1:
            SmearPt1Back(spacelike_T(i),(initx+i+maxdim)%maxdim,inity,initz,tidx);
            SmearPt2(spacelike_T(i));
            break;
        case 2:
            SmearPt1Back(spacelike_T(i),initx,(inity+i+maxdim)%maxdim,initz,tidx);
            SmearPt2(spacelike_T(i));
            break;
        case 3:
            SmearPt1Back(spacelike_T(i),initx,inity,(initz+i+maxdim)%maxdim,tidx);
            SmearPt2(spacelike_T(i));
            break;
        }//switch
    }
}

void ScaleSetV::CountTimeLineUp(arma::cx_mat & timeline,
const int actidx,const int actidy,const int actidz,const int actidt,const int maxtime){
timeline.eye();
int sumaxtime=SU3Grid::GetTDim();
for(int i=0;i<maxtime;i++){
    timeline=mymodell.GetModellGrid()(0).GetGrid()((actidt+i+sumaxtime)%sumaxtime,actidx,actidy,actidz)*timeline;
}

}

void ScaleSetV::CountTimeLineDown(arma::cx_mat & timeline,
const int actidx,const int actidy,const int actidz,const int actidt,const int maxtime){
timeline.eye();
int sumaxtime=SU3Grid::GetTDim();
for(int i=0;i<maxtime;i++){
    timeline=mymodell.GetModellGrid()(0).GetGrid()((actidt-i+sumaxtime)%sumaxtime,actidx,actidy,actidz).t()*timeline;
}
}
void ScaleSetV::CountSpaceLine0(arma::cx_mat & spaceline){
    spaceline.eye();
    for(int i=0;i<R;i++){
        spaceline=spacelike_0(i)*spaceline;
    }
}
void ScaleSetV::CountSpaceLineT(arma::cx_mat & spaceline){
    spaceline.eye();
    for(int i=0;i<R;i++){
        spaceline=spacelike_T(R-1-i)*spaceline;
    }
}
void ScaleSetV::BuildCorrelM(){
int maxtimedim=SU3Grid::GetTDim();
int maxspacedim=SU3Grid::GetDim();
int TplusOne=(T+1+maxtimedim)%maxtimedim;
InitSpaceLikeTInv(T,spacelike_T);
//debug
//for(int i=0;i<spacelike_T.Nx();i++){
//std::cout<<"spacelikeT"<<std::endl;
//std::cout<<spacelike_T(i)<<std::endl;
//std::cout<<"mymodell"<<std::endl;
////specificated for spacegrdidx==1
//std::cout<<mymodell.GetModellGrid()(spacegrididx).GetGrid()(initt+T,initx+i,inity,initz).t()<<std::endl;
//std::cin.ignore();
//std::cin.get();
//}

InitSpaceLikeT(0,spacelike_0);
//debug
//for(int i=0;i<spacelike_T.Nx();i++){
//std::cout<<"spacelike0"<<std::endl;
//std::cout<<spacelike_0(i)<<std::endl;
//std::cout<<"mymodell"<<std::endl;
////specificated for spacegrdidx==1
//std::cout<<mymodell.GetModellGrid()(spacegrididx).GetGrid()(initt,initx+i,inity,initz)<<std::endl;
//std::cin.ignore();
//std::cin.get();
//}
arma::cx_mat spline0(3,3,arma::fill::eye);
arma::cx_mat splineT(3,3,arma::fill::eye);
arma::cx_mat tlineup(3,3,arma::fill::eye);
arma::cx_mat tlinedown(3,3,arma::fill::eye);
arma::cx_mat tlineupplus(3,3,arma::fill::eye);
arma::cx_mat tlinedownplus(3,3,arma::fill::eye);
arma::cx_mat tlineup1(3,3,arma::fill::eye);
arma::cx_mat tlinedown1(3,3,arma::fill::eye);
CountTimeLineDown(tlinedown,initx,inity,initz,(initt+T-1+maxtimedim)%maxtimedim,T);
//debug
//std::cout<<"tlinedown"<<std::endl;
//std::cout<<tlinedown<<std::endl;
//std::cout<<"mymodell"<<std::endl;
////specificated for spacegrdidx==1
//std::cout<<mymodell.GetModellGrid()(0).GetGrid()(initt,initx,inity,initz).t()<<std::endl;
//std::cin.ignore();
//std::cin.get();

CountTimeLineDown(tlinedownplus,initx,inity,initz,(initt+T+maxtimedim)%maxtimedim,TplusOne);
tlinedown1=mymodell.GetModellGrid()(0).GetGrid()(initt,initx,inity,initz).t();
switch(spacegrididx){
case 1:
    CountTimeLineUp(tlineup,(initx+R+maxspacedim)%maxspacedim,inity,initz,initt,T);
//    //debug
//std::cout<<"tlineup"<<std::endl;
//std::cout<<tlineup<<std::endl;
//std::cout<<"mymodell"<<std::endl;
////specificated for spacegrdidx==1
//std::cout<<mymodell.GetModellGrid()(0).GetGrid()(initt,initx+R,inity,initz)<<std::endl;
//std::cin.ignore();
//std::cin.get();

    CountTimeLineUp(tlineupplus,(initx+R+maxspacedim)%maxspacedim,inity,initz,initt,TplusOne);
    tlineup1=mymodell.GetModellGrid()(0).GetGrid()(initt,(initx+R+maxspacedim)%maxspacedim,inity,initz);
    break;
case 2:
    CountTimeLineUp(tlineup,initx,(inity+R+maxspacedim)%maxspacedim,initz,initt,TplusOne);
    CountTimeLineUp(tlineupplus,initx,(inity+R+maxspacedim)%maxspacedim,initz,initt,TplusOne);
    tlineup1=mymodell.GetModellGrid()(0).GetGrid()(initt,initx,(inity+R+maxspacedim)%maxspacedim,initz);
    break;
case 3:
    CountTimeLineUp(tlineup,initx,inity,(initz+R+maxspacedim)%maxspacedim,initt,T);
    CountTimeLineUp(tlineupplus,initx,inity,(initz+R+maxspacedim)%maxspacedim,initt,TplusOne);
    tlineup1=mymodell.GetModellGrid()(0).GetGrid()(initt,initx,inity,(initz+R+maxspacedim)%maxspacedim);
    break;

}//switch

arma::cx_mat temp(3,3,arma::fill::eye);
InitSpaceLikeTInv(T,spacelike_T);
for(int i=0;i<maxsmearlevel;i++){
    InitSpaceLikeT(0,spacelike_0);
    for(int j=0;j<maxsmearlevel;j++){
        CountSpaceLine0(spline0);
        CountSpaceLineT(splineT);
        temp.eye();
        //debug
//        std::cout<<"spline0"<<std::endl;
//        std::cout<<spline0<<std::endl;
//        std::cout<<"splineT"<<std::endl;
//        std::cout<<splineT<<std::endl;
        temp=spline0*temp;
        temp=tlineup*temp;
        temp=splineT*temp;
        temp=tlinedown*temp;
//        std::cout<<"correl before"<<std::endl;
//        std::cout<<correlT<<std::endl;
        correlT(i,j)+=real(trace(temp));
//        std::cout<<"correl after"<<std::endl;
//        std::cout<<correlT<<std::endl;
        //correlT(i,j)+=real(trace(spline0*tlineup*splineT*tlinedown));
        Smearing0();
        //debug
//for(int k=0;k<spacelike_T.Nx();k++){
//std::cout<<"spacelike0"<<std::endl;
//std::cout<<spacelike_0(k)<<std::endl;
//std::cin.ignore();
//std::cin.get();
//}
//debug
//for(int k=0;k<spacelike_T.Nx();k++){
//std::cout<<"spacelikeT"<<std::endl;
//std::cout<<spacelike_T(k)<<std::endl;
//std::cin.ignore();
//std::cin.get();
//}
    }//for j smear
    SmearingT(T);
}//for i smear

InitSpaceLikeTInv(TplusOne,spacelike_T);
for(int i=0;i<maxsmearlevel;i++){
    InitSpaceLikeT(0,spacelike_0);
    for(int j=0;j<maxsmearlevel;j++){
        CountSpaceLine0(spline0);
        CountSpaceLineT(splineT);
        temp.eye();
        temp=spline0*temp;
        temp=tlineup*temp;
        temp=splineT*temp;
        temp=tlinedown*temp;
        correlT1(i,j)+=real(trace(temp));
        //correlT1(i,j)+=real(trace(spline0*tlineupplus*splineT*tlinedownplus));
        Smearing0();

    }//for j smear
    SmearingT(TplusOne);
}//for i smear

InitSpaceLikeTInv(1,spacelike_T);
for(int i=0;i<maxsmearlevel;i++){
    InitSpaceLikeT(0,spacelike_0);
    for(int j=0;j<maxsmearlevel;j++){
        CountSpaceLine0(spline0);
        CountSpaceLineT(splineT);
        temp.eye();
        temp=spline0*temp;
        temp=tlineup*temp;
        temp=splineT*temp;
        temp=tlinedown*temp;
        correl0(i,j)+=real(trace(temp));
        //correl0(i,j)+=real(trace(spline0*tlineup1*splineT*tlinedown1));
        Smearing0();

    }//for j smear
    SmearingT(1);
}//for i smear

}


void ScaleSetV::CorrelMAVG(const int num,arma::cx_mat & rescorr0,arma::cx_mat & rescorrT,arma::cx_mat & rescorrT1,std::string dir){
    correl0.zeros();
    correlT1.zeros();
    correlT.zeros();

     rescorr0.zeros();
     rescorrT.zeros();
     rescorrT1.zeros();
std::ofstream outcorr0;
outcorr0.precision(6);
outcorr0.open((dir+"outcorr0.dat").c_str(),std::ios::out);

std::ofstream outcorrT;
outcorrT.precision(6);
outcorrT.open((dir+"outcorrT.dat").c_str(),std::ios::out);

std::ofstream outcorrT1;
outcorrT1.precision(6);
outcorrT1.open((dir+"outcorrT1.dat").c_str(),std::ios::out);

    int counter=0;
    for(int i=0;i<num;i++){
        this->WilsonAVG(outcorr0,outcorrT,outcorrT1);
        rescorr0+=correl0;
        rescorrT+=correlT;
        rescorrT1+=correlT1;
        counter++;
        for(int j=0;j<10;j++){
            mymodell.HeatBathSweep();
        }

    }//for i

    rescorr0/=counter;
    rescorrT/=counter;
    rescorrT1/=counter;

outcorr0.close();
outcorrT.close();
outcorrT1.close();

    std::cout<<"after avg:"<<std::endl;
    std::cout<<rescorrT<<std::endl;
    std::cout<<"is it identity?"<<std::endl;
    std::cout<<rescorrT.t()*rescorrT.i()<<std::endl;
    std::cin.ignore();
    std::cin.get();
}

void ScaleSetV::WilsonAVG(std::ofstream & outcorrel0,std::ofstream & outcorrelT,std::ofstream & outcorrelT1){
if(!outcorrel0.is_open()) throw "outcorrel0 file is not open";
if(!outcorrelT.is_open()) throw "outcorrelT file is not open";
if(!outcorrelT1.is_open()) throw "outcorrelT1 file is not open";

int maxspacedim=SU3Grid::GetDim();
int maxtimedim=SU3Grid::GetTDim();

    correl0.zeros();
    correlT1.zeros();
    correlT.zeros();

int counter=0;

for(int i=0;i<maxtimedim;i++){
initt=(oinitt+i+maxtimedim)%maxtimedim;
    for(int j=0;j<maxspacedim;j++){
    initx=(oinitx+j+maxspacedim)%maxspacedim;
        for(int k=0;k<maxspacedim;k++){
        inity=(oinity+k+maxspacedim)%maxspacedim;
            for(int l=0;l<maxspacedim;l++){
            initz=(oinitz+l+maxspacedim)%maxspacedim;
                BuildCorrelM();
                counter++;
                //std::cout<<real(correlT(0,0))/counter<<std::endl;
            }//for z
        }//for y
    }//for x


}//for i time

initt=oinitt;
initx=oinitx;
inity=oinity;
initz=oinitz;

//avg-ing
correl0/=counter;
correlT/=counter;
correlT1/=counter;

outcorrel0<<real(correl0(0,0))<<'\t'
          <<real(correl0(0,1))<<'\t'
          <<real(correl0(1,0))<<'\t'
          <<real(correl0(1,1))<<'\t'
          <<imag(correl0(0,0))<<'\t'
          <<imag(correl0(0,1))<<'\t'
          <<imag(correl0(1,0))<<'\t'
          <<imag(correl0(1,1))<<std::endl;
outcorrelT<<real(correlT(0,0))<<'\t'
          <<real(correlT(0,1))<<'\t'
          <<real(correlT(1,0))<<'\t'
          <<real(correlT(1,1))<<'\t'
          <<imag(correlT(0,0))<<'\t'
          <<imag(correlT(0,1))<<'\t'
          <<imag(correlT(1,0))<<'\t'
          <<imag(correlT(1,1))<<std::endl;
outcorrelT1<<real(correlT1(0,0))<<'\t'
          <<real(correlT1(0,1))<<'\t'
          <<real(correlT1(1,0))<<'\t'
          <<real(correlT1(1,1))<<'\t'
          <<imag(correlT1(0,0))<<'\t'
          <<imag(correlT1(0,1))<<'\t'
          <<imag(correlT1(1,0))<<'\t'
          <<imag(correlT1(1,1))<<std::endl;

//debug
//std::cout<<"avg-ed for wilson loops"<<std::endl;
//std::cout<<correlT<<std::endl;
//std::cout<<"is it id?"<<std::endl;
//std::cout<<correlT.t()*correlT.i()<<std::endl;

}

void ScaleSetV::WilsonAVG(){

int maxspacedim=SU3Grid::GetDim();
int maxtimedim=SU3Grid::GetTDim();

    correl0.zeros();
    correlT1.zeros();
    correlT.zeros();

int counter=0;

for(int i=0;i<maxtimedim;i++){
initt=(oinitt+i+maxtimedim)%maxtimedim;
    for(int j=0;j<maxspacedim;j++){
    initx=(oinitx+j+maxspacedim)%maxspacedim;
        for(int k=0;k<maxspacedim;k++){
        inity=(oinity+k+maxspacedim)%maxspacedim;
            for(int l=0;l<maxspacedim;l++){
            initz=(oinitz+l+maxspacedim)%maxspacedim;
                BuildCorrelM();
                counter++;
                //std::cout<<real(correlT(0,0))/counter<<std::endl;
            }//for z
        }//for y
    }//for x


}//for i time

initt=oinitt;
initx=oinitx;
inity=oinity;
initz=oinitz;

//avg-ing
correl0/=counter;
correlT/=counter;
correlT1/=counter;

//debug
//std::cout<<"avg-ed for wilson loops"<<std::endl;
//std::cout<<correlT<<std::endl;
//std::cout<<"is it id?"<<std::endl;
//std::cout<<correlT.t()*correlT.i()<<std::endl;

}

void ScaleSetV::isitsymm(){
    correl0.zeros();
    correlT1.zeros();
    correlT.zeros();
    BuildCorrelM();
//    std::cout<<"is it symmetric? "<<std::endl;
//    std::cout<<correl0<<std::endl;
//    std::cout<<correlT<<std::endl;
//    std::cout<<correlT1<<std::endl;
//    std::cin.ignore();
//    std::cin.get();
//    std::cout<<"is it identity?"<<std::endl;
//    std::cout<<correlT.t()*correlT.i()<<std::endl;

}

//double ScaleSetV::CountV(){
//    double value=0;
//    CorrelMAVG(15);
//not works in debian testing
//    arma::cx_mat t0t=arma::sqrtmat(correl0).i()*correlT*arma::sqrtmat(correl0).i();
//    arma::cx_mat t0t1=arma::sqrtmat(correl0).i()*correlT1*arma::sqrtmat(correl0).i();
//
//    arma::cx_vec eigvals;
//    arma::eig_gen(eigvals,t0t);
//    std::complex<double> maxeigt0t=eigvals.max();
//    double eig1=real(maxeigt0t);
//    arma::eig_gen(eigvals,t0t1);
//    std::complex<double> maxeigt0t1=eigvals.max();
//    double eig2=real(maxeigt0t1);
//    value=log(eig1/eig2);
//    return value;
//}

const Array::array1<arma::cx_mat> & ScaleSetV::GetSpace0Grid()const{
    return spacelike_0;
}

const Array::array1<arma::cx_mat> & ScaleSetV::GetSpaceTGrid()const{
    return spacelike_T;
}

Array::array1<arma::cx_mat> & ScaleSetV::ModSpace0Grid(){
    return spacelike_0;
}
Array::array1<arma::cx_mat> & ScaleSetV::ModSpaceTGrid(){
    return spacelike_T;
}

const arma::cx_mat & ScaleSetV::GetCorrelT()const{
return correlT;
}


ScaleSetV::~ScaleSetV(){}

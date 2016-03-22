#include<cmath>
#include<armadillo>
#include<complex>
#include"Sommer.h"

void ScaleSetV::InitSpaceLikeT(const int time,Array::array1<arma::cx_mat> & matarray){
const int maxdim=SU3Grid::GetDim();
for(int i=0;i<R;i++){


    switch(spacegrididx){
        case 1:
            {int idx=(initx+i+maxdim)%maxdim;
            matarray(i)=this->mymodell.GetModellGrid()(spacegrididx).GetGrid()(initt+time,idx,inity,initz);}
            break;
        case 2:
            {int idy=(inity+i+maxdim)%maxdim;
            matarray=this->mymodell.GetModellGrid()(spacegrididx).GetGrid()(initt+time,initx,idy,initz);}
            break;
        case 3:
            {int idz=(initz+i+maxdim)%maxdim;
            matarray=this->mymodell.GetModellGrid()(spacegrididx).GetGrid()(initt+time,initx,inity,idz);}
            break;
    }//switch


}//for

}

//CTR
ScaleSetV::ScaleSetV(Modell & mymodell,
const int r,const int t,const int tidx,
const int xidx,const int yidx,const int zidx,const int grididx):
mymodell(mymodell),R(r),T(t),initt(tidx),initx(xidx),inity(yidx),
initz(zidx),spacegrididx(grididx),alpha(0.5),maxsmearlevel(5),spacelike_0(R),spacelike_T(R),
correlT(maxsmearlevel,maxsmearlevel,arma::fill::zeros),correlT1(maxsmearlevel,maxsmearlevel,arma::fill::zeros),
correl0(maxsmearlevel,maxsmearlevel,arma::fill::zeros){

InitSpaceLikeT(0,spacelike_0);
InitSpaceLikeT(T,spacelike_T);

}

//for one link
void ScaleSetV::SmearPt1(arma::cx_mat & link,int actidx,int actidy,int actidz,const int timeidx){
    arma::cx_mat produc(3,3,arma::fill::eye);
    const int maxdim=SU3Grid::GetDim();
    //for space dim.s
    for(int i=1;i<4;i++){
        produc.eye();
        if(i!=spacegrididx){
            produc=mymodell.GetModellGrid()(i).GetGrid()(timeidx,actidx,actidy,actidz)*produc;
        //step along second direction
        switch(i){
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
        switch(i){
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
        produc=mymodell.GetModellGrid()(i).GetGrid()(timeidx,actidx,actidy,actidz).t()*produc;

        link+=alpha*produc;

        produc.eye();
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
                //step back along second direction
        switch(i){
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
            produc=mymodell.GetModellGrid()(i).GetGrid()(timeidx,actidx,actidy,actidz).t()*produc;
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
            produc=mymodell.GetModellGrid()(i).GetGrid()(timeidx,actidx,actidy,actidz)*produc;
        link+=alpha*produc;

        }//if
        else{}
    }//for spacedim.s

}
//project back to su3
//void ScaleSetV::SmearPt2(arma::cx_mat & link){
//    double value=real(trace(link*link.t()));
//    double newvalue=0;
//    arma::cx_mat tmp(3,3,arma::fill::zeros);
//    arma::cx_mat subsu3(2,2,arma::fill::zeros);
//    std::cout<<"original link: "<<std::endl;
//    std::cout<<link<<std::endl;
//    //initialize tmp as link
//    for(int i=0;i<3;i++){
//        for(int j=0;j<3;j++){
//            tmp(i,j)=link(i,j);
//        }
//    }
//    int counter=0;
//    do{
//        counter++;
//        //firts subgroup
//        subsu3(0,0)=link(0,0);
//        subsu3(1,0)=link(1,0);
//        subsu3(0,1)=link(0,1);
//        subsu3(1,1)=link(1,1);
//        subsu3=subsu3/det(subsu3);
//        tmp(0,0)=subsu3(0,0);
//        tmp(1,0)=subsu3(1,0);
//        tmp(0,1)=subsu3(0,1);
//        tmp(1,1)=subsu3(1,1);
//        newvalue=real(trace(tmp*link.t()));
//        link(0,0)=tmp(0,0);
//        link(1,0)=tmp(1,0);
//        link(0,1)=tmp(0,1);
//        link(1,1)=tmp(1,1);
//        if(newvalue<value){
//            std::cout<<"first sub."<<std::endl;
//            std::cout<<"val: "<<value<<" newval: "<<newvalue<<std::endl;
//            std::cout<<"no more task"<<std::endl;
//        }
//        else{
//            value=newvalue;
//            //second subgroup
//        subsu3(0,0)=link(1,1);
//        subsu3(1,0)=link(2,1);
//        subsu3(0,1)=link(1,2);
//        subsu3(1,1)=link(2,2);
//        subsu3=subsu3/det(subsu3);
//        tmp(1,1)=subsu3(0,0);
//        tmp(2,1)=subsu3(1,0);
//        tmp(1,2)=subsu3(0,1);
//        tmp(2,2)=subsu3(1,1);
//        newvalue=real(trace(tmp*link.t()));
//        link(1,1)=tmp(1,1);
//        link(2,1)=tmp(2,1);
//        link(1,2)=tmp(1,2);
//        link(2,2)=tmp(2,2);
//            if(newvalue<value){
//                std::cout<<"second sub."<<std::endl;
//                std::cout<<"val: "<<value<<" newval: "<<newvalue<<std::endl;
//                std::cout<<"no more task"<<std::endl;
//            }
//            else{
//            value=newvalue;
//            //second subgroup
//            subsu3(0,0)=link(0,0);
//            subsu3(1,0)=link(2,0);
//            subsu3(0,1)=link(0,2);
//            subsu3(1,1)=link(2,2);
//            subsu3=subsu3/det(subsu3);
//            tmp(0,0)=subsu3(0,0);
//            tmp(2,0)=subsu3(1,0);
//            tmp(0,2)=subsu3(0,1);
//            tmp(2,2)=subsu3(1,1);
//            newvalue=real(trace(tmp*link.t()));
//            link(0,0)=tmp(0,0);
//            link(2,0)=tmp(2,0);
//            link(0,2)=tmp(0,2);
//            link(2,2)=tmp(2,2);
//            std::cout<<"third sub."<<std::endl;
//            std::cout<<"val: "<<value<<" newval: "<<newvalue<<std::endl;
//            }//
//        }//
//                std::cout<<"act number of iterations: "<<counter<<std::endl;
//
//        }while(newvalue>value);
//        //debug
//        std::cout<<"final number of iterations: "<<counter<<std::endl;
//    std::cout<<"final link: "<<std::endl;
//    std::cout<<link<<std::endl;
//    std::cout<<"is it identity?"<<std::endl;
//    std::cout<<link*link.t()<<std::endl;
//    std::cout<<"det: "<<arma::det(link)<<std::endl;
//}

//project back to su3
void ScaleSetV::SmearPt2(arma::cx_mat & link){
    arma::cx_mat diagonal(3,3,arma::fill::eye);
    arma::cx_mat eigvecs;
    arma::vec eigvals;
    arma::eig_sym(eigvals,eigvecs,link.t()*link,"std");
    diagonal(0,0)=eigvals(0)*diagonal(0,0);
    diagonal(1,1)=eigvals(1)*diagonal(1,1);
    diagonal(2,2)=eigvals(2)*diagonal(2,2);
    std::cout<<"eigvals"<<std::endl;
    std::cout<<eigvals<<std::endl;
    std::cout<<"eigvecs"<<std::endl;
    std::cout<<eigvecs<<std::endl;
    std::cout<<"diagonal"<<std::endl;
    std::cout<<diagonal<<std::endl;
    arma::cx_mat sqrtmatrix=eigvecs*sqrt(diagonal)*eigvecs.i();
    std::cout<<"original link"<<std::endl;
    std::cout<<link<<std::endl;
    link=link*sqrtmatrix.i()*pow(arma::det(link.i()*link.t()),1./6);
    std::cout<<"new link"<<std::endl;
    std::cout<<link<<std::endl;
    std::cout<<"is it identity?"<<std::endl;
    std::cout<<link.t()*link<<std::endl;
    std::cout<<arma::det(link)<<std::endl;
}

void ScaleSetV::Smearing0(){

    for(int i=0;i<R;i++){
        //smearing a link
        switch(spacegrididx){
        case 1:
            SmearPt1(spacelike_0(i),initx+i,inity,initz,initt);
            SmearPt2(spacelike_0(i));
            break;
        case 2:
            SmearPt1(spacelike_0(i),initx,inity+i,initz,initt);
            SmearPt2(spacelike_0(i));
            break;
        case 3:
            SmearPt1(spacelike_0(i),initx,inity,initz+i,initt);
            SmearPt2(spacelike_0(i));
            break;
        }//switch
    }
}

void ScaleSetV::SmearingT(){

    for(int i=0;i<R;i++){
        //smearing a link
        switch(spacegrididx){
        case 1:
            SmearPt1(spacelike_T(i),initx+i,inity,initz,initt);
            SmearPt2(spacelike_T(i));
            break;
        case 2:
            SmearPt1(spacelike_T(i),initx,inity+i,initz,initt);
            SmearPt2(spacelike_T(i));
            break;
        case 3:
            SmearPt1(spacelike_T(i),initx,inity,initz+i,initt);
            SmearPt2(spacelike_T(i));
            break;
        }//switch
    }
}

void ScaleSetV::CountTimeLineUp(arma::cx_mat & timeline,
const int actidx,const int actidy,const int actidz,const int actidt,const int maxtime){
timeline.eye();
for(int i=0;i<maxtime;i++){
    timeline=mymodell.GetModellGrid()(0).GetGrid()(actidt+i,actidx,actidy,actidz)*timeline;
}

}

void ScaleSetV::CountTimeLineDown(arma::cx_mat & timeline,
const int actidx,const int actidy,const int actidz,const int actidt,const int maxtime){
timeline.eye();
for(int i=0;i<maxtime;i++){
    timeline=mymodell.GetModellGrid()(0).GetGrid()(actidt-i,actidx,actidy,actidz).t()*timeline;
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
        spaceline=spacelike_T(R-1-i).t()*spaceline;
    }
}
void ScaleSetV::BuildCorrelM(){
arma::cx_mat spline0(3,3,arma::fill::eye);
arma::cx_mat splineT(3,3,arma::fill::eye);
arma::cx_mat tlineup(3,3,arma::fill::eye);
arma::cx_mat tlinedown(3,3,arma::fill::eye);
arma::cx_mat tlineupplus(3,3,arma::fill::eye);
arma::cx_mat tlinedownplus(3,3,arma::fill::eye);
arma::cx_mat tlineup1(3,3,arma::fill::eye);
arma::cx_mat tlinedown1(3,3,arma::fill::eye);
CountTimeLineDown(tlinedown,initx,inity,initz,initt+T-1,T);
CountTimeLineDown(tlinedownplus,initx,inity,initz,initt+T,T+1);
tlinedown1=mymodell.GetModellGrid()(0).GetGrid()(initt,initx,inity,initz).t();
switch(spacegrididx){
case 1:
    CountTimeLineUp(tlineup,initx+R,inity,initz,initt,T);
    CountTimeLineUp(tlineupplus,initx+R,inity,initz,initt,T+1);
    tlineup1=mymodell.GetModellGrid()(0).GetGrid()(initt,initx+R,inity,initz);
    break;
case 2:
    CountTimeLineUp(tlineup,initx,inity+R,initz,initt,T);
    CountTimeLineUp(tlineupplus,initx,inity+R,initz,initt,T+1);
    tlineup1=mymodell.GetModellGrid()(0).GetGrid()(initt,initx,inity+R,initz);
    break;
case 3:
    CountTimeLineUp(tlineup,initx,inity,initz+R,initt,T);
    CountTimeLineUp(tlineupplus,initx,inity,initz+R,initt,T+1);
    tlineup1=mymodell.GetModellGrid()(0).GetGrid()(initt,initx,inity,initz+R);
    break;

}//switch

InitSpaceLikeT(T,spacelike_T);
for(int i=0;i<maxsmearlevel;i++){
    InitSpaceLikeT(0,spacelike_0);
    for(int j=0;j<maxsmearlevel;j++){
        CountSpaceLine0(spline0);
        CountSpaceLineT(splineT);
        correlT(i,j)+=trace(spline0*tlineup*splineT*tlinedown);
        Smearing0();

    }//for j smear
    SmearingT();
}//for i smear

InitSpaceLikeT(T+1,spacelike_T);
for(int i=0;i<maxsmearlevel;i++){
    InitSpaceLikeT(0,spacelike_0);
    for(int j=0;j<maxsmearlevel;j++){
        CountSpaceLine0(spline0);
        CountSpaceLineT(splineT);
        correlT1(i,j)+=trace(spline0*tlineupplus*splineT*tlinedownplus);
        Smearing0();

    }//for j smear
    SmearingT();
}//for i smear

InitSpaceLikeT(1,spacelike_T);
for(int i=0;i<maxsmearlevel;i++){
    InitSpaceLikeT(0,spacelike_0);
    for(int j=0;j<maxsmearlevel;j++){
        CountSpaceLine0(spline0);
        CountSpaceLineT(splineT);
        correl0(i,j)+=trace(spline0*tlineup1*splineT*tlinedown1);
        Smearing0();

    }//for j smear
    SmearingT();
}//for i smear

}

void ScaleSetV::CorrelMAVG(const int num){
    correl0.zeros();
    correlT1.zeros();
    correlT.zeros();
    for(int i=0;i<num;i++){
        BuildCorrelM();
        for(int j=0;j<10;j++){
            mymodell.HeatBathSweep();
        }

    }
    correlT1/=num;
    correlT/=num;
    std::cout<<"after avg:"<<std::endl;
    std::cout<<correlT<<std::endl;
}

void ScaleSetV::isitsymm(){
    correl0.zeros();
    correlT1.zeros();
    correlT.zeros();
    BuildCorrelM();
    std::cout<<"is it symmetric? "<<std::endl;
    std::cout<<correl0<<std::endl;
    std::cout<<correlT<<std::endl;
    std::cout<<correlT1<<std::endl;

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
ScaleSetV::~ScaleSetV(){}

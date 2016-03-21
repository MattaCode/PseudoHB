#include<cmath>
#include"Sommer.h"

void ScaleSetV::InitSpaceLikeT(const int time,Array::array1<arma::cx_mat> matarray){
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
correlT(maxsmearlevel,maxsmearlevel,arma::fill::zeros),correlT1(maxsmearlevel,maxsmearlevel,arma::fill::zeros){

InitSpaceLikeT(0,spacelike_0);
InitSpaceLikeT(T-1,spacelike_T);

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
void ScaleSetV::SmearPt2(arma::cx_mat & link){
    double value=real(trace(link*link.t()));
    double newvalue=0;
    arma::cx_mat tmp(3,3,arma::fill::zeros);
    arma::cx_mat subsu3(2,2,arma::fill::zeros);

    //initialize tmp as link
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            tmp(i,j)=link(i,j);
        }
    }
    int counter=0;
    do{
        counter++;
        //firts subgroup
        subsu3(0,0)=link(0,0);
        subsu3(1,0)=link(1,0);
        subsu3(0,1)=link(0,1);
        subsu3(1,1)=link(1,1);
        subsu3=subsu3/det(subsu3);
        tmp(0,0)=subsu3(0,0);
        tmp(1,0)=subsu3(1,0);
        tmp(0,1)=subsu3(0,1);
        tmp(1,1)=subsu3(1,1);
        newvalue=real(trace(tmp*link.t()));
        link(0,0)=tmp(0,0);
        link(1,0)=tmp(1,0);
        link(0,1)=tmp(0,1);
        link(1,1)=tmp(1,1);
        if(newvalue<value){}
        else{
            value=newvalue;
            //second subgroup
        subsu3(0,0)=link(1,1);
        subsu3(1,0)=link(2,1);
        subsu3(0,1)=link(1,2);
        subsu3(1,1)=link(2,2);
        subsu3=subsu3/det(subsu3);
        tmp(1,1)=subsu3(0,0);
        tmp(2,1)=subsu3(1,0);
        tmp(1,2)=subsu3(0,1);
        tmp(2,2)=subsu3(1,1);
        newvalue=real(trace(tmp*link.t()));
        link(1,1)=tmp(1,1);
        link(2,1)=tmp(2,1);
        link(1,2)=tmp(1,2);
        link(2,2)=tmp(2,2);
            if(newvalue<value){}
            else{
            value=newvalue;
            //second subgroup
            subsu3(0,0)=link(0,0);
            subsu3(1,0)=link(3,0);
            subsu3(0,1)=link(0,3);
            subsu3(1,1)=link(3,3);
            subsu3=subsu3/det(subsu3);
            tmp(0,0)=subsu3(0,0);
            tmp(3,0)=subsu3(1,0);
            tmp(0,3)=subsu3(0,1);
            tmp(3,3)=subsu3(1,1);
            newvalue=real(trace(tmp*link.t()));
            link(0,0)=tmp(0,0);
            link(3,0)=tmp(3,0);
            link(0,3)=tmp(0,3);
            link(3,3)=tmp(3,3);
            }
        }
                std::cout<<"number of iterations: "<<counter<<std::endl;

        }while(newvalue>value);
        //debug
        std::cout<<"final number of iterations: "<<counter<<std::endl;

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
CountTimeLineDown(tlinedown,initx,inity,initz,initt+T-1,T);
CountTimeLineDown(tlinedownplus,initx,inity,initz,initt+T,T+1);

switch(spacegrididx){
case 1:
    CountTimeLineUp(tlineup,initx+R,inity,initz,initt,T);
        CountTimeLineUp(tlineupplus,initx+R,inity,initz,initt,T+1);
    break;
case 2:
    CountTimeLineUp(tlineup,initx,inity+R,initz,initt,T);
        CountTimeLineUp(tlineupplus,initx,inity+R,initz,initt,T+1);

    break;
case 3:
    CountTimeLineUp(tlineup,initx,inity,initz+R,initt,T);
        CountTimeLineUp(tlineupplus,initx,inity,initz+R,initt,T+1);

    break;

}//switch


for(int i=0;i<maxsmearlevel;i++){
    InitSpaceLikeT(0,spacelike_0);
    for(int j=0;j<maxsmearlevel;j++){
        CountSpaceLine0(spline0);
        CountSpaceLineT(splineT);
        correlT(i,j)=trace(spline0*tlineup*splineT*tlinedown);
        correlT1(i,j)=trace(spline0*tlineupplus*splineT*tlinedownplus);
        Smearing0();

    }//for j smear
    SmearingT();
}//for i smear

}

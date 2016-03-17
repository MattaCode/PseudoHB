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
initz(zidx),spacegrididx(grididx),spacelike_0(R),spacelike_T(R),
correlT(3,3,arma::fill::zeros),correlT1(3,3,arma::fill::zeros){

InitSpaceLikeT(0,spacelike_0);
InitSpaceLikeT(T-1,spacelike_T);

}


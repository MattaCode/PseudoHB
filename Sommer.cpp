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
initz(zidx),spacegrididx(grididx),alpha(0.5),spacelike_0(R),spacelike_T(R),
correlT(3,3,arma::fill::zeros),correlT1(3,3,arma::fill::zeros){

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

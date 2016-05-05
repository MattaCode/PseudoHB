#include<cmath>
#include<armadillo>
#include<complex>
#include<string>
#include"Sommer.h"

void ScaleSetV::InitSm(){
int maxtdim=SU3Grid::GetTDim();
int maxdim=SU3Grid::GetDim();
for(int i=0;i<maxtdim;i++){
    for(int j=0;j<maxdim;j++){
        for(int k=0;k<maxdim;k++){
            for(int l=0;l<maxdim;l++){
                smeared(0,0).ModifyGrid()(i,j,k,l)=mymodell.GetModellGrid()(0).GetGrid()(i,j,k,l);
                smeared(0,1).ModifyGrid()(i,j,k,l)=mymodell.GetModellGrid()(1).GetGrid()(i,j,k,l);
                smeared(0,2).ModifyGrid()(i,j,k,l)=mymodell.GetModellGrid()(2).GetGrid()(i,j,k,l);
                smeared(0,3).ModifyGrid()(i,j,k,l)=mymodell.GetModellGrid()(3).GetGrid()(i,j,k,l);
                smearedinv(0,0).ModifyGrid()(i,j,k,l)=mymodell.GetModellGrid()(0).GetGrid()(i,j,k,l).t();
                smearedinv(0,1).ModifyGrid()(i,j,k,l)=mymodell.GetModellGrid()(1).GetGrid()(i,j,k,l).t();
                smearedinv(0,2).ModifyGrid()(i,j,k,l)=mymodell.GetModellGrid()(2).GetGrid()(i,j,k,l).t();
                smearedinv(0,3).ModifyGrid()(i,j,k,l)=mymodell.GetModellGrid()(3).GetGrid()(i,j,k,l).t();
            }//for l
        }//for k
    }//for j
}//for i

}

//CTR
ScaleSetV::ScaleSetV(Modell & mymodell,
const int r,const int t,const int tidx,
const int xidx,const int yidx,const int zidx,const int grididx):
mymodell(mymodell),R(r),T(t),oinitt(tidx),oinitx(xidx),oinity(yidx),
oinitz(zidx),spacegrididx(grididx),alpha(0.5),maxsmearlevel(3),
correlT(maxsmearlevel,maxsmearlevel,arma::fill::zeros),correlT1(maxsmearlevel,maxsmearlevel,arma::fill::zeros),
correl0(maxsmearlevel,maxsmearlevel,arma::fill::zeros),smeared(maxsmearlevel,4),smearedinv(maxsmearlevel,4){
initt=tidx;
initx=xidx;
inity=yidx;
initz=zidx;
InitSm();
Smear();

}

//CTR - overload - set maxsmearlevel
ScaleSetV::ScaleSetV(Modell & mymodell,
const int r,const int t,const int tidx,
const int xidx,const int yidx,const int zidx,const int grididx,const int maxsmear):
mymodell(mymodell),R(r),T(t),oinitt(tidx),oinitx(xidx),oinity(yidx),
oinitz(zidx),spacegrididx(grididx),alpha(0.5),maxsmearlevel(maxsmear),
correlT(maxsmearlevel,maxsmearlevel,arma::fill::zeros),correlT1(maxsmearlevel,maxsmearlevel,arma::fill::zeros),
correl0(maxsmearlevel,maxsmearlevel,arma::fill::zeros),smeared(maxsmearlevel,4),smearedinv(maxsmearlevel,4){
initt=tidx;
initx=xidx;
inity=yidx;
initz=zidx;
InitSm();
Smear();

}

//smear pt1 helper
void ScaleSetV::TriplForw(arma::cx_mat & produc,int actidx,int actidy,int actidz,const int timeidx,
const int spacidx,const int wilspaceidx,const int actsmearlevel){
    produc.eye();
        const int maxdim=SU3Grid::GetDim();
    if(spacidx==wilspaceidx) throw "spacidx equals wilson space direction";
        produc=smeared(actsmearlevel,spacidx).GetGrid()(timeidx,actidx,actidy,actidz)*produc;
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
        produc=smeared(actsmearlevel,wilspaceidx).GetGrid()(timeidx,actidx,actidy,actidz)*produc;
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
        produc=smeared(actsmearlevel,spacidx).GetGrid()(timeidx,actidx,actidy,actidz).t()*produc;
//restore indices to original
        //step back along first direction
        switch(wilspaceidx){
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

void ScaleSetV::TriplBack(arma::cx_mat & produc,int actidx,int actidy,int actidz,const int timeidx,
const int spacidx,const int wilspaceidx,const int actsmearlevel){
produc.eye();
    const int maxdim=SU3Grid::GetDim();
if(spacidx==wilspaceidx) throw "spacidx equals wilson space direction";

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
            produc=smeared(actsmearlevel,spacidx).GetGrid()(timeidx,actidx,actidy,actidz).t()*produc;
            produc=smeared(actsmearlevel,wilspaceidx).GetGrid()(timeidx,actidx,actidy,actidz)*produc;
        //step along first direction
        switch(wilspaceidx){
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
        produc=smeared(actsmearlevel,spacidx).GetGrid()(timeidx,actidx,actidy,actidz)*produc;
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
        switch(wilspaceidx){
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
void ScaleSetV::SmearPt1Forw(arma::cx_mat & link,int actidx,int actidy,int actidz,const int timeidx,
const int wilspaceidx,const int actsmearlevel){
    arma::cx_mat produc(3,3,arma::fill::eye);
    //for space dim.s
    for(int i=1;i<4;i++){
        produc.eye();
        if(i!=wilspaceidx){
           TriplForw(produc,actidx,actidy,actidz,timeidx,i,wilspaceidx,actsmearlevel);
            link+=alpha*produc;
            TriplBack(produc,actidx,actidy,actidz,timeidx,i,wilspaceidx,actsmearlevel);
            link+=alpha*produc;
        }//if
        else{}
    }//for spacedim.s

}

//for one link
void ScaleSetV::SmearPt1Back(arma::cx_mat & link,int actidx,int actidy,int actidz,const int timeidx,
const int wilspaceidx,const int actsmearlevel){
    arma::cx_mat produc(3,3,arma::fill::eye);
    //for space dim.s
    for(int i=1;i<4;i++){
        produc.eye();
        if(i!=wilspaceidx){
           TriplForw(produc,actidx,actidy,actidz,timeidx,i,wilspaceidx,actsmearlevel);
            link+=(alpha*produc).i();
            TriplBack(produc,actidx,actidy,actidz,timeidx,i,wilspaceidx,actsmearlevel);
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

    arma::cx_mat sqrtmatrix=eigvecs*sqrt(diagonal)*eigvecs.i();

    link=link*sqrtmatrix.i()*pow(arma::det(link.i()*link.t()),1./6);

}

void ScaleSetV::Smear(){
int maxtdim=SU3Grid::GetTDim();
int maxdim=SU3Grid::GetDim();
for(int actsmearlevel=0;actsmearlevel<(maxsmearlevel-1);actsmearlevel++){
for(int grididx=1;grididx<4;grididx++){
for(int i=0;i<maxtdim;i++){
    for(int j=0;j<maxdim;j++){
        for(int k=0;k<maxdim;k++){
            for(int l=0;l<maxdim;l++){
                //init actsmeared+1
                smeared(actsmearlevel+1,grididx).ModifyGrid()(i,j,k,l)
                    =smeared(actsmearlevel,grididx).GetGrid()(i,j,k,l);
                //smear pt1
                SmearPt1Forw(smeared(actsmearlevel+1,grididx).ModifyGrid()(i,j,k,l),j,k,l,i,grididx,actsmearlevel);
                //smear pt2
                SmearPt2(smeared(actsmearlevel+1,grididx).ModifyGrid()(i,j,k,l));
                //init actsmeared+1
                smearedinv(actsmearlevel+1,grididx).ModifyGrid()(i,j,k,l)
                    =smearedinv(actsmearlevel,grididx).GetGrid()(i,j,k,l);
                //smear pt1
                SmearPt1Back(smearedinv(actsmearlevel+1,grididx).ModifyGrid()(i,j,k,l),j,k,l,i,grididx,actsmearlevel);
                //smear pt2
                SmearPt2(smearedinv(actsmearlevel+1,grididx).ModifyGrid()(i,j,k,l));
            }//for l
        }//for k
    }//for j
}//for i
}//for grid
}//for actsmearlevel
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
void ScaleSetV::CountSpaceLine0(arma::cx_mat & spaceline,const int actsmearlevel){
    spaceline.eye();
    int maxdim=SU3Grid::GetDim();
    switch(spacegrididx){
        case 1:
            for(int i=0;i<R;i++){
                spaceline=smeared(actsmearlevel,spacegrididx).GetGrid()(initt,(initx+i+maxdim)%maxdim,inity,initz)*spaceline;
            }
        break;
        case 2:
            for(int i=0;i<R;i++){
                spaceline=smeared(actsmearlevel,spacegrididx).GetGrid()(initt,initx,(inity+i+maxdim)%maxdim,initz)*spaceline;
            }
        break;
        case 3:
            for(int i=0;i<R;i++){
                spaceline=smeared(actsmearlevel,spacegrididx).GetGrid()(initt,initx,inity,(initz+i+maxdim)%maxdim)*spaceline;
            }
        break;

    }


}
void ScaleSetV::CountSpaceLineT(arma::cx_mat & spaceline,const int actsmearlevel,const int tshift){
    spaceline.eye();
    int maxdim=SU3Grid::GetDim();
    int maxtdim=SU3Grid::GetTDim();
    switch(spacegrididx){
        case 1:
            for(int i=0;i<R;i++){
                spaceline=smearedinv(actsmearlevel,spacegrididx).GetGrid()((initt+tshift+maxtdim)%maxtdim,(initx+R-1-i+maxdim)%maxdim,inity,initz)*spaceline;
            }
        break;
        case 2:
            for(int i=0;i<R;i++){
                spaceline=smearedinv(actsmearlevel,spacegrididx).GetGrid()((initt+tshift+maxtdim)%maxtdim,initx,(inity+R-1-i+maxdim)%maxdim,initz)*spaceline;
            }
        break;
        case 3:
            for(int i=0;i<R;i++){
                spaceline=smearedinv(actsmearlevel,spacegrididx).GetGrid()((initt+tshift+maxtdim)%maxtdim,initx,inity,(initz+R-1-i+maxdim)%maxdim)*spaceline;
            }
        break;

    }

}

void ScaleSetV::BuildCorrelM(){
int maxtimedim=SU3Grid::GetTDim();
int maxspacedim=SU3Grid::GetDim();
int TplusOne=(T+1+maxtimedim)%maxtimedim;
arma::cx_mat spline0(3,3,arma::fill::eye);
arma::cx_mat splineT(3,3,arma::fill::eye);
arma::cx_mat tlineup(3,3,arma::fill::eye);
arma::cx_mat tlinedown(3,3,arma::fill::eye);
arma::cx_mat tlineupplus(3,3,arma::fill::eye);
arma::cx_mat tlinedownplus(3,3,arma::fill::eye);
arma::cx_mat tlineup1(3,3,arma::fill::eye);
arma::cx_mat tlinedown1(3,3,arma::fill::eye);
CountTimeLineDown(tlinedown,initx,inity,initz,(initt+T-1+maxtimedim)%maxtimedim,T);

CountTimeLineDown(tlinedownplus,initx,inity,initz,(initt+T+maxtimedim)%maxtimedim,TplusOne);
tlinedown1=mymodell.GetModellGrid()(0).GetGrid()(initt,initx,inity,initz).t();
switch(spacegrididx){
case 1:
    CountTimeLineUp(tlineup,(initx+R+maxspacedim)%maxspacedim,inity,initz,initt,T);
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
for(int i=0;i<maxsmearlevel;i++){
    for(int j=0;j<maxsmearlevel;j++){
        CountSpaceLine0(spline0,i);
        CountSpaceLineT(splineT,j,T);
        temp.eye();
        temp=spline0*temp;
        temp=tlineup*temp;
        temp=splineT*temp;
        temp=tlinedown*temp;
        correlT(i,j)+=trace(temp);

    }//for j smear
}//for i smear

for(int i=0;i<maxsmearlevel;i++){
    for(int j=0;j<maxsmearlevel;j++){
        CountSpaceLine0(spline0,i);
        CountSpaceLineT(splineT,j,TplusOne);
        temp.eye();
        temp=spline0*temp;
        temp=tlineupplus*temp;
        temp=splineT*temp;
        temp=tlinedownplus*temp;
        correlT1(i,j)+=trace(temp);

    }//for j smear
}//for i smear

for(int i=0;i<maxsmearlevel;i++){
    for(int j=0;j<maxsmearlevel;j++){
        CountSpaceLine0(spline0,i);
        CountSpaceLineT(splineT,j,1);
        temp.eye();
        temp=spline0*temp;
        temp=tlineup1*temp;
        temp=splineT*temp;
        temp=tlinedown1*temp;
        correl0(i,j)+=trace(temp);

    }//for j smear
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
        this->WilsonAVG(outcorr0,outcorrT,outcorrT1);
        rescorr0+=correl0;
        rescorrT+=correlT;
        rescorrT1+=correlT1;
        counter++;

    for(int i=1;i<num;i++){
        for(int j=0;j<10;j++){
            mymodell.HeatBathSweep();
        }
        InitSm();
        Smear();
        this->WilsonAVG(outcorr0,outcorrT,outcorrT1);
        rescorr0+=correl0;
        rescorrT+=correlT;
        rescorrT1+=correlT1;
        counter++;


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

for(int mi=0;mi<maxsmearlevel;mi++){
	for(int mj=0;mj<maxsmearlevel;mj++){
		if((mi==(maxsmearlevel-1))&&(mj==(maxsmearlevel-1))){
			outcorrel0<<real(correl0(mi,mj));
			outcorrelT<<real(correlT(mi,mj));
			outcorrelT1<<real(correlT1(mi,mj));
		}//if
		else{
		outcorrel0<<real(correl0(mi,mj))<<'\t';
		outcorrelT<<real(correlT(mi,mj))<<'\t';
		outcorrelT1<<real(correlT1(mi,mj))<<'\t';
		}//else
	}//for mj
}//for mi
outcorrel0<<std::endl;
outcorrelT<<std::endl;
outcorrelT1<<std::endl;

	
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

}

void ScaleSetV::Symmetrize(arma::mat & matrix){
    matrix=(matrix+matrix.t())/2;
}

void ScaleSetV::isitsymm(){
    correl0.zeros();
    correlT1.zeros();
    correlT.zeros();
    BuildCorrelM();

}

double ScaleSetV::CountV(){
    double value=0;
    arma::cx_mat rescorr0(maxsmearlevel,maxsmearlevel,arma::fill::zeros);
    arma::cx_mat rescorrT(maxsmearlevel,maxsmearlevel,arma::fill::zeros);
    arma::cx_mat rescorrT1(maxsmearlevel,maxsmearlevel,arma::fill::zeros);

    CorrelMAVG(50,rescorr0,rescorrT,rescorrT1,"./");
//not works in debian testing
//    arma::cx_mat t0t=arma::sqrtmat(correl0).i()*correlT*arma::sqrtmat(correl0).i();
//    arma::cx_mat t0t1=arma::sqrtmat(correl0).i()*correlT1*arma::sqrtmat(correl0).i();

    arma::mat realrescorr0=real(rescorr0);
    arma::mat realrescorrT=real(rescorrT);
    arma::mat realrescorrT1=real(rescorrT1);

    Symmetrize(realrescorr0);
    Symmetrize(realrescorrT);
    Symmetrize(realrescorrT1);
    //Diagonalize c(t0) aka realrescorr0
    arma::vec eigval0;
    arma::mat eigvec0;
    eig_sym(eigval0, eigvec0, realrescorr0);
    arma::cx_mat diag0mat(maxsmearlevel,maxsmearlevel,arma::fill::eye);
    std::cout<<"eigval: "<<eigval0<<std::endl;
    std::cout<<"eigvec: "<<eigvec0<<std::endl;
    //build sqrtmat
    arma::cx_vec cpleigval0(maxsmearlevel);
    for(int i=0;i<maxsmearlevel;i++){
        cpleigval0(i)=eigval0(i);
        diag0mat(i,i)*=sqrt(cpleigval0(i));
    }
    arma::cx_mat sqrt0mat=eigvec0*diag0mat*eigvec0.i();
    std::cout<<"sqrtmat: "<<sqrt0mat<<std::endl;

    //find eigvals of sqrtinv*c*sqrtinv
    arma::cx_vec eigvals;
    arma::eig_gen(eigvals,sqrt0mat.i()*realrescorrT*sqrt0mat.i());
    std::cout<<"eigvals! "<<eigvals<<std::endl;
    double maxeigt0t=real(eigvals.max());
    std::cout<<eigvals.max()<<std::endl;
    arma::cx_vec eigvals2;
    arma::eig_gen(eigvals2,sqrt0mat.i()*realrescorrT1*sqrt0mat.i());
    double maxeigt0t1=real(eigvals2.max());
    std::cout<<eigvals2.max()<<std::endl;
    std::cin.ignore();
    std::cin.get();
    value=log(maxeigt0t/maxeigt0t1);
    return value;
}
double ScaleSetV::CountV(std::string dir){
    double value=0;
    arma::cx_mat rescorr0(maxsmearlevel,maxsmearlevel,arma::fill::zeros);
    arma::cx_mat rescorrT(maxsmearlevel,maxsmearlevel,arma::fill::zeros);
    arma::cx_mat rescorrT1(maxsmearlevel,maxsmearlevel,arma::fill::zeros);

    CorrelMAVG(50,rescorr0,rescorrT,rescorrT1,dir);
//not works in debian testing
//    arma::cx_mat t0t=arma::sqrtmat(correl0).i()*correlT*arma::sqrtmat(correl0).i();
//    arma::cx_mat t0t1=arma::sqrtmat(correl0).i()*correlT1*arma::sqrtmat(correl0).i();

    arma::mat realrescorr0=real(rescorr0);
    arma::mat realrescorrT=real(rescorrT);
    arma::mat realrescorrT1=real(rescorrT1);

    Symmetrize(realrescorr0);
    Symmetrize(realrescorrT);
    Symmetrize(realrescorrT1);
    //Diagonalize c(t0) aka realrescorr0
    arma::vec eigval0;
    arma::mat eigvec0;
    eig_sym(eigval0, eigvec0, realrescorr0);
    arma::cx_mat diag0mat(maxsmearlevel,maxsmearlevel,arma::fill::eye);
    //std::cout<<"eigval: "<<eigval0<<std::endl;
    //std::cout<<"eigvec: "<<eigvec0<<std::endl;
    //build sqrtmat
    arma::cx_vec cpleigval0(maxsmearlevel);
    for(int i=0;i<maxsmearlevel;i++){
        cpleigval0(i)=eigval0(i);
        diag0mat(i,i)*=sqrt(cpleigval0(i));
    }
    arma::cx_mat sqrt0mat=eigvec0*diag0mat*eigvec0.i();
    //std::cout<<"sqrtmat: "<<sqrt0mat<<std::endl;

    //find eigvals of sqrtinv*c*sqrtinv
    arma::cx_vec eigvals;
    arma::eig_gen(eigvals,sqrt0mat.i()*realrescorrT*sqrt0mat.i());
    //std::cout<<"eigvals! "<<eigvals<<std::endl;
    double maxeigt0t=real(eigvals.max());
    //std::cout<<eigvals.max()<<std::endl;
    arma::cx_vec eigvals2;
    arma::eig_gen(eigvals2,sqrt0mat.i()*realrescorrT1*sqrt0mat.i());
    double maxeigt0t1=real(eigvals2.max());
    //std::cout<<eigvals2.max()<<std::endl;
    //std::cin.ignore();
    //std::cin.get();
    std::ofstream eigs;
    eigs.precision(6);
    eigs.open((dir+"eigvals.dat").c_str(),std::ios::out);
    eigs<<maxeigt0t<<'\t'<<maxeigt0t1<<std::endl;
    eigs.close();
    value=log(maxeigt0t/maxeigt0t1);
    return value;
}
const arma::cx_mat & ScaleSetV::GetCorrelT()const{
return correlT;
}


ScaleSetV::~ScaleSetV(){}

#include"algorithm.h"

using namespace std;

/*******************************/
/********* AutoCorrel **********/
/*******************************/

//Constr
AutoCorrel::AutoCorrel(const int& dim,
const int & tmin,const int & tincr):dim(dim),tmin(tmin),tincr(tincr),
tmax(tmin+(dim-1)*tincr), result(dim,3),curr(dim),tmp(tmax){
if(tmin==0) throw "tmin should not be zero!";
section=0;
actidx=0;
loop0avg=0;

for(int i=0;i<dim;i++){
    curr(i)=(-i*tincr);
    for(int j=0;j<3;j++){
        result(i,j)=0;
    }//for j
}//for i

    for(int i=0;i<tmax;i++){
        tmp(i)=0;
    }

}

void AutoCorrel::InsertData(const double & newdata){
        //debud
    cout<<"curr i"<<endl;
    for(int i=0;i<curr.Nx();i++){
        cout<<curr(i)<<endl;
    }

    //only saving
    if(actidx<tmin){
        tmp(actidx)=newdata;
        loop0avg+=newdata;
        actidx++;
    }
    //saving and counting
    else if ((actidx>=tmin)&&(actidx<tmax)){
        loop0avg+=newdata;
        //cycle on result array
            for(int i=0;i<dim;i++){
                if(curr(i)>=0){
                    result(i,0)+=tmp(curr(i))*newdata; //<A_0*A_tau>
                    result(i,1)+=newdata;//<A_tau>
                    result(i,2)++;//averager counter


                }
                curr(i)+=1;
            }//for
            tmp(actidx)=newdata;
            actidx++;
    }//elsif
    else{
        loop0avg+=newdata;
        for(int i=0;i<dim;i++){
        //debug
        cout<<tmp(curr(i)%tmax)*newdata<<endl;
           result(i,0)+=tmp(curr(i)%tmax)*newdata;
           result(i,1)+=newdata;
           result(i,2)++;
           curr(i)+=1;
                               //debug
                    cout<<"result"<<endl;
                    cout<<result(i,0)<<endl;
        }//for i
        tmp(actidx%tmax)=newdata;
        actidx++;
    }//else
}

double AutoCorrel::GetLoop0AVG()const{
    return loop0avg/actidx;
}

double AutoCorrel::GetResult(int rownum)const{
    double avgedcorr1=result(rownum,0)/result(rownum,2);
    double avgedcorr2=result(rownum,1)/result(rownum,2);
    double avgedcorr3=GetLoop0AVG();
    return (avgedcorr1-avgedcorr2*avgedcorr3);
}

double AutoCorrel::GetCorr(int row)const{
    double avgedcorr1=result(row,0)/result(row,2);
    return avgedcorr1;
}

double AutoCorrel::GetAVG(int row)const{
    double avgedcorr2=result(row,1)/result(row,2);
    return avgedcorr2;
}

    //Get data
    const int AutoCorrel::GetDim()const{
        return dim;
    }
    const int AutoCorrel::GetTMin()const{
        return tmin;
    }
    const int AutoCorrel::GetTIncr()const{
        return tincr;
    }
    const int AutoCorrel::GetTMax()const{
        return tmax;
    }
/***************************************************************/

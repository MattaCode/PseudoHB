#ifndef ALGORITHM_H
#define ALGORITHM_H

#include"fftw++-1.13/Array.h"

class AutoCorrel{

    const int dim;
    const int tmin;
    const int tincr;
    const int tmax;


    Array::array2<double> result;
    Array::array1<int> curr;
    Array::array1<double> tmp;

    int section;
    int actidx;
    double loop0avg;
    public:
    AutoCorrel(const int&,const int&,const int&);
    //COPY
    AutoCorrel(const AutoCorrel&)=delete;

    //Get data
    const int GetDim()const;
    const int GetTMin()const;
    const int GetTIncr()const;
    const int GetTMax()const;
    double GetLoop0AVG()const;
    double GetResult(int)const;
    double GetCorr(int)const;
    double GetAVG(int)const;

    //insert data
    void InsertData(const double &);

};

#endif

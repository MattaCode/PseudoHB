#ifndef PSEUDOHB_H
#define PSEUDOHB_H
#include<complex>
#include<armadillo>
#include"fftw++-1.13/Array.h"

class SU3Grid{

    static const int dim;
    static const int tdim;
    Array::array4<arma::cx_mat> ei;

    public:
    //default
    SU3Grid();
    //from random
    SU3Grid(bool);

//TODO
    //construct from file
    //SU3Grid(const char * filename);

    //SU3Grid(std::istream&);

    //assignment operator
    SU3Grid& operator=(const SU3Grid& fromsu3grid);

    //copy
    SU3Grid(const SU3Grid& su3grid)=delete;

    //write SU3Grid to os
    void writeSU3Grid(std::ostream &)const;

    //read SU3Grid from is
    void readSU3Grid(std::istream&);

    //get value of dimension
    static const int GetDim();
    static const int GetTDim();
    //GetGrid for reading
    const Array::array4<arma::cx_mat>& GetGrid()const;

    //GetGrid for modifing a link
    Array::array4<arma::cx_mat>& ModifyGrid();

    //destructor
    ~SU3Grid();
};

class Modell{
    //links
    Array::array1<SU3Grid> grid;
    static const double beta;

    public:

    //default constr
    Modell();
    //construct from random
    Modell(bool);

    //construct from file
    Modell(const char* );

    void RandomInit();

    //ASSIGNMENT
    Modell& operator=(const Modell&);

    //COPY
    Modell(const Modell&);
    //get beta
    static const double GetBeta();

    //write Modell to os
    void writeModell( std::ostream& )const;

    //write Modell to file
    void writeToFileModell( const char* )const;

    //read Modell from is
    void readModell( std::istream& );

    //read Modell from file
    void readFromFileModell( const char* );

    //count plaquett at idx, for two selected direction i,j
    //up initialized as Identity
    //assumes calls with existing indices idx-y-z-k and i-j
    void CountUp(int ,int ,int ,int ,const unsigned int , const unsigned int ,arma::cx_mat & );

    //Modify the selected link
    void ModifyLink(int, int, int, int, int, const arma::cx_mat& );

    //access for reading
    const Array::array1<SU3Grid>& GetModellGrid()const;


    ~Modell();

};

#endif

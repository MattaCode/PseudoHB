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
    //with the same matrix
    SU3Grid(const arma::cx_mat & );

//TODO
    //construct from file
    //SU3Grid(const char * filename);

    //SU3Grid(std::istream&);

    //assignment operator
    SU3Grid& operator=(const SU3Grid& fromsu3grid);

    //copy
    SU3Grid(const SU3Grid& );

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
    arma::cx_mat su3staple;
    arma::cx_mat su2staple;
    std::complex<double> su2strootdet;
    static const std::complex<double> iunit;
    static const arma::cx_mat pauli1;
    static const arma::cx_mat pauli2;
    static const arma::cx_mat pauli3;
    static const arma::cx_mat identity2;



    public:

    //default constr
    Modell();
    //construct from random
    Modell(bool);
    //same matrix for all element
    Modell(arma::cx_mat &);

    //construct from file
    Modell(const char* );

    void RandomInit();

    void MatrixInit(arma::cx_mat &);

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
    //void readFromFileModell( const char* );

    //count plaquett at idx, for two selected direction i,j
    //result initialized as Identity
    //assumes calls with existing indices idx-y-z-k and i-j
    void CountUp(int ,int ,int ,int ,const unsigned int , const unsigned int ,arma::cx_mat & );

    //count action for a plaquett up
     double CountPlaqEnergy(const arma::cx_mat & );
    //count mean for plaquett energy on lattice
     double CountMeanEnergyDens(std::ofstream &,const bool,const int);
     //helper for CountMeanEnergyDens but one grididx fixed
     void CountForGrididxEnergyDens(const int,double &,int &,std::ofstream &,const bool);
    //count mean for plaquett energy on lattice without output
     double CountMeanEnergyDens(const bool,const int);
     //helper for CountMeanEnergyDens but one grididx fixed without output
     void CountForGrididxEnergyDens(const int,double &,int &,const bool);

    //count energy in a subsystem
    double CountBoxEnergy(const int,const int,const int,const int,const int,const bool,const int);
    //helper for CountBoxenergy but one grididx fixed
    void CountForGrididxBoxEnergy(const int,const int,const int,const int,const int,const int,double &,int &,const bool);
    //count boxenergy histogram
    void BoxEnHisto(const int,std::ofstream &,const bool,const int);

    //count forward staple (plaquett without the selected link)
    //result initialized as Identity
    void TriplUForw(const unsigned int, const unsigned int, int, int, int, int, arma::cx_mat &);

    //count all six staple for a selected link
    //result init.ed as Identity
    void Count6Staple(const unsigned int,int,int,int,int);

    //count backward staple (plaquett without the selected link)
    //result initialized as Identity
    void TriplURev(const unsigned int, const unsigned int, int, int, int, int, arma::cx_mat &);

    //find Pauli coefficients for 2x2 submatrix of su3staple
    std::vector<double> CountCoeffs(const int,const arma::cx_mat &);

    //generating new su2 matrix coeff0
    double GenerateCoeff0();

    //generate new su2 matrix coeffs 3d sphere
    std::vector<double> GenerateCoeffs(double);

    // build SU2 matrix
    void BuildSU2(const double, const std::vector<double> & ,arma::cx_mat &);

    //build generated and transformed su3 matrix
    //refresh link with it
    void RefreshLinkpart(const int, const int, const int, const int, const int,const int);

    //su2staple
    //and determinant
    void BuildSU2staple(const int,const unsigned int ,const int ,const int ,const int ,const int );

    //Modify the selected link
    void ModifyLink(int, int, int, int, int, const arma::cx_mat& );

    //Get selected link
    const arma::cx_mat&  GetLink(int, int, int, int, int)const;

    //access for reading
    const Array::array1<SU3Grid>& GetModellGrid()const;

    //heat bath step
    void HeatBathStep(const unsigned int ,int ,int , int , int);
    //heat bath sweep
    void HeatBathSweep();

    //Polyakov loop
    void PolyakovMatrix(const int,const int,const int,arma::cx_mat &);

    //Polyakov space avg
    std::complex<double> PolyakovLoopAVG(std::ofstream &);

    //Polyakov space avg without output
    std::complex<double> PolyakovLoopAVG();

    //Wilson loop
    std::complex<double> WilsonLoop(const int,const int,int,int,int,int,int);
    //wilson avg
    std::complex<double> WilsonAvg(const int,const int,int);
    //debug
    static void GetPauli();

    ~Modell();

};

#endif

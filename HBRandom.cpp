#include<boost/random/bernoulli_distribution.hpp>
#include<boost/random/uniform_real_distribution.hpp>
#include<boost/random/uniform_on_sphere.hpp>
#include<boost/random/mersenne_twister.hpp>
#include<cmath>
#include"HBRandom.h"

//random generator
boost::mt19937 gen(time(NULL));

//Bernoulli true with param.
//return: true - accept, false - reject
bool Flip(const double accept){
    boost::random::bernoulli_distribution<> dist(accept);
    return dist(gen);
}


//uniform real
//min. param to max param.
double GetRealRandom(const double xmin,const double xmax){
	boost::random::uniform_real_distribution<> dist(xmin, xmax);
	return dist(gen);
}

//uniform on D sphere
//param: dimension of sphere
std::vector<double> RandOnSphere(const int dim){
    boost::random::uniform_on_sphere<> dist(dim);
    return dist(gen);
}
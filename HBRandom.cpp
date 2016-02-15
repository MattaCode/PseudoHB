#include<boost/random/bernoulli_distribution.hpp>
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

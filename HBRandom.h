#ifndef HBRANDOM_H
#define HBRANDOM_H

//uniform real
//min. param to max param.
double GetRealRandom(const double ,const double);

//Bernoulli true with param.
bool Flip(const double);

//uniform on D sphere
//param: dimension of sphere
std::vector<double> RandOnSphere(const int);

#endif // HBRandom.h

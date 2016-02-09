#include<iostream>
#include"PseudoHB.h"

int main(){

SU3Grid mygrid(true);

std::cout<<mygrid.GetGrid()(0,0,0,0)<<std::endl;

Modell mymodell(true);
std::cout<<"modell:"<<std::endl;
std::cout<<mymodell.GetModellGrid()(0).GetGrid()(0,0,0,0)<<std::endl;

return 0;
}

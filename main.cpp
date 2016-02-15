#include<iostream>
#include"PseudoHB.h"
#include"HBRandom.h"

using namespace std;

int main(){

SU3Grid mygrid(true);
SU3Grid mygrid2(true);

std::cout<<mygrid.GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();
std::cout<<mygrid2.GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();
mygrid2=mygrid;

std::cout<<mygrid.GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();
std::cout<<mygrid2.GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();

SU3Grid mygrid3=mygrid;
std::cout<<mygrid.GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();
std::cout<<mygrid2.GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();
    std::cout<<mygrid3.GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();


Modell mymodell(true);
std::cout<<"modell:"<<std::endl;
std::cout<<mymodell.GetModellGrid()(0).GetGrid()(0,0,0,0)<<std::endl;
cout<<mygrid2.GetGrid()(0,0,0,1)<<endl;
mymodell.ModifyLink(0,0,0,0,0,mygrid2.GetGrid()(0,0,0,1));
std::cout<<mymodell.GetModellGrid()(0).GetGrid()(0,0,0,0)<<std::endl;
    cin.ignore();
    cin.get();

for(int i=0;i<10;i++){
cout<<"random Bernoulli: "<<Flip(0.5)<<endl;
}

return 0;
}

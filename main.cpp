#include<iostream>
#include<vector>
#include"PseudoHB.h"
#include"HBRandom.h"
#include"testing.h"

using namespace std;




int main(){

const int mcmaxtime=1;
arma::cx_mat id3d(3,3,arma::fill::eye);
Modell mymodell(id3d);

//debug
//Modell::GetPauli();


for(int t=0;t<mcmaxtime;t++){
    mymodell.HeatBathSweep();
    //mymodell.PolyakovLoopAVG();
}

return 0;
}

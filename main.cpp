#include<iostream>
#include<vector>
#include"PseudoHB.h"
#include"HBRandom.h"
#include"testing.h"

using namespace std;




int main(){

const int mcmaxtime=1;
Modell mymodell(true);

for(int t=0;t<mcmaxtime;t++){
    mymodell.HeatBathSweep();
    mymodell.PolyakovLoopAVG();
}

return 0;
}

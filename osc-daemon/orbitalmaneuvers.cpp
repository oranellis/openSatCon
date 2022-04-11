#include <iostream>
#include <vector>
#include <math.h>

#include "osctypes.hpp"
#include "planet.cpp"

namespace osc{
orbparam VNBburn(orbparam koearg, vnb VNBdeltav){
    orbparam koeret;
    //eci curpos, eci curvel=KOEtoECI(koeret);//how do?
    eci ECIdeltaV;

    newECIvel.i =   (curpos.i/sqrt(pow(curpos.i,2)+pow(curpos.j,2)+pow(curpos.k,2)))*VNBdeltav.v
                    +((curpos.j*curvel.k-curpos.k*curvel.j)
                    /sqrt(pow((curpos.j*curvel.k-curpos.k*curvel.j),2)+pow((curpos.k*curvel.i-curpos.i*curvel.k),2)+pow((curpos.i*curvel.j-curpos.j*curvel.i),2)))*VNBdeltav.n
                    +(((curvel.k*curpos.i-curvel.i*curpos.k)*curvel.k-(curvel.i*curpos.j-curvel.j*curpos.i)*curvel.j)
                    /sqrt(pow(((curvel.k*curpos.i-curvel.i*curpos.k)*curvel.k-(curvel.i*curpos.j-curvel.j*curpos.i)*curvel.j),2)
                    +pow(((curvel.i*curpos.j-curvel.j*curpos.i)*curvel.i-(curvel.j*curpos.k-curvel.k*curpos.j)*curvel.k),2)
                    +pow(((curvel.j*curpos.k-curvel.k*curpos.j)*curvel.j-(curvel.k*curpos.i-curvel.i*curpos.k)*curvel.i),2)))*VNBdeltav.b+curvel.i;

    newECIvel.j =   (curpos.j/sqrt(pow(curpos.i,2)+pow(curpos.j,2)+pow(curpos.k,2)))*VNBdeltav.v
                    +((curpos.k*curvel.i-curpos.i*curvel.k)
                    /sqrt(pow((curpos.j*curvel.k-curpos.k*curvel.j),2)+pow((curpos.k*curvel.i-curpos.i*curvel.k),2)+pow((curpos.i*curvel.j-curpos.j*curvel.i),2)))*VNBdeltav.n
                    +(((curvel.i*curpos.j-curvel.j*curpos.i)*curvel.i-(curvel.j*curpos.k-curvel.k*curpos.j)*curvel.k)
                    /sqrt(pow(((curvel.k*curpos.i-curvel.i*curpos.k)*curvel.k-(curvel.i*curpos.j-curvel.j*curpos.i)*curvel.j),2)
                    +pow(((curvel.i*curpos.j-curvel.j*curpos.i)*curvel.i-(curvel.j*curpos.k-curvel.k*curpos.j)*curvel.k),2)
                    +pow(((curvel.j*curpos.k-curvel.k*curpos.j)*curvel.j-(curvel.k*curpos.i-curvel.i*curpos.k)*curvel.i),2)))*VNBdeltav.b+curvel.j;

    newECIvel.k =   (curpos.k/sqrt(pow(curpos.i,2)+pow(curpos.j,2)+pow(curpos.k,2)))*VNBdeltav.v
                    +((curpos.i*curvel.j-curpos.j*curvel.i)
                    /sqrt(pow((curpos.j*curvel.k-curpos.k*curvel.j),2)+pow((curpos.k*curvel.i-curpos.i*curvel.k),2)+pow((curpos.i*curvel.j-curpos.j*curvel.i),2)))*VNBdeltav.n
                    +(((curvel.j*curpos.k-curvel.k*curpos.j)*curvel.j-(curvel.k*curpos.i-curvel.i*curpos.k)*curvel.i)
                    /sqrt(pow(((curvel.k*curpos.i-curvel.i*curpos.k)*curvel.k-(curvel.i*curpos.j-curvel.j*curpos.i)*curvel.j),2)
                    +pow(((curvel.i*curpos.j-curvel.j*curpos.i)*curvel.i-(curvel.j*curpos.k-curvel.k*curpos.j)*curvel.k),2)
                    +pow(((curvel.j*curpos.k-curvel.k*curpos.j)*curvel.j-(curvel.k*curpos.i-curvel.i*curpos.k)*curvel.i),2)))*VNBdeltav.b+curvel.k;

    orbparam koeret =  ECItoKOE(curpos, newECIvel);

    
};

};
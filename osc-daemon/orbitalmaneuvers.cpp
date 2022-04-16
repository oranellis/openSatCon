#include <iostream>
#include <vector>
#include <math.h>

#include "osctypes.hpp"
#include "planet.cpp"
#include "axistransforms.cpp"

namespace osc{
orbparam VNBburn(orbparam koearg, vnb VNBdV){
    orbparam koeret;
    pcs pscposvel=KOEtoPCS(koeret);
    eci posvel=PCStoECI(koeret, pscposvel);
    eci neweci;
   

    neweci.vi =   (posvel.i/sqrt(pow(posvel.i,2)+pow(posvel.j,2)+pow(posvel.k,2)))*VNBdV.v
                    +((posvel.j*posvel.vk-posvel.k*posvel.vj)
                    /sqrt(pow((posvel.j*posvel.vk-posvel.k*posvel.vj),2)+pow((posvel.k*posvel.vi-posvel.i*posvel.vk),2)+pow((posvel.i*posvel.vj-posvel.j*posvel.vi),2)))*VNBdV.n
                    +(((posvel.vk*posvel.i-posvel.vi*posvel.k)*posvel.vk-(posvel.vi*posvel.j-posvel.vj*posvel.i)*posvel.vj)
                    /sqrt(pow(((posvel.vk*posvel.i-posvel.vi*posvel.k)*posvel.vk-(posvel.vi*posvel.j-posvel.vj*posvel.i)*posvel.vj),2)
                    +pow(((posvel.vi*posvel.j-posvel.vj*posvel.i)*posvel.vi-(posvel.vj*posvel.k-posvel.vk*posvel.j)*posvel.vk),2)
                    +pow(((posvel.vj*posvel.k-posvel.vk*posvel.j)*posvel.vj-(posvel.vk*posvel.i-posvel.vi*posvel.k)*posvel.vi),2)))*VNBdV.b+posvel.vi;

    neweci.vj =   (posvel.j/sqrt(pow(posvel.i,2)+pow(posvel.j,2)+pow(posvel.k,2)))*VNBdV.v
                    +((posvel.k*posvel.vi-posvel.i*posvel.vk)
                    /sqrt(pow((posvel.j*posvel.vk-posvel.k*posvel.vj),2)+pow((posvel.k*posvel.vi-posvel.i*posvel.vk),2)+pow((posvel.i*posvel.vj-posvel.j*posvel.vi),2)))*VNBdV.n
                    +(((posvel.vi*posvel.j-posvel.vj*posvel.i)*posvel.vi-(posvel.vj*posvel.k-posvel.vk*posvel.j)*posvel.vk)
                    /sqrt(pow(((posvel.vk*posvel.i-posvel.vi*posvel.k)*posvel.vk-(posvel.vi*posvel.j-posvel.vj*posvel.i)*posvel.vj),2)
                    +pow(((posvel.vi*posvel.j-posvel.vj*posvel.i)*posvel.vi-(posvel.vj*posvel.k-posvel.vk*posvel.j)*posvel.vk),2)
                    +pow(((posvel.vj*posvel.k-posvel.vk*posvel.j)*posvel.vj-(posvel.vk*posvel.i-posvel.vi*posvel.k)*posvel.vi),2)))*VNBdV.b+posvel.vj;

    neweci.vk =   (posvel.k/sqrt(pow(posvel.i,2)+pow(posvel.j,2)+pow(posvel.k,2)))*VNBdV.v
                    +((posvel.i*posvel.vj-posvel.j*posvel.vi)
                    /sqrt(pow((posvel.j*posvel.vk-posvel.k*posvel.vj),2)+pow((posvel.k*posvel.vi-posvel.i*posvel.vk),2)+pow((posvel.i*posvel.vj-posvel.j*posvel.vi),2)))*VNBdV.n
                    +(((posvel.vj*posvel.k-posvel.vk*posvel.j)*posvel.vj-(posvel.vk*posvel.i-posvel.vi*posvel.k)*posvel.vi)
                    /sqrt(pow(((posvel.vk*posvel.i-posvel.vi*posvel.k)*posvel.vk-(posvel.vi*posvel.j-posvel.vj*posvel.i)*posvel.vj),2)
                    +pow(((posvel.vi*posvel.j-posvel.vj*posvel.i)*posvel.vi-(posvel.vj*posvel.k-posvel.vk*posvel.j)*posvel.vk),2)
                    +pow(((posvel.vj*posvel.k-posvel.vk*posvel.j)*posvel.vj-(posvel.vk*posvel.i-posvel.vi*posvel.k)*posvel.vi),2)))*VNBdV.b+posvel.vk;

    neweci.i=posvel.i; //should add the small 
    neweci.j=posvel.j; //change in time, based 
    neweci.k=posvel.k; //off calculation refresh rate

    orbparam koeret =  ECItoKOE(neweci);

    
};

};
#include <iostream>
#include <vector>
#include <math.h>

#include "osctypes.hpp"
#include "planet.cpp"

namespace osc{

    quaternion quaternionrotation(orbrot arg_i, orbrot arg_f){
        // Be aware that this does not handle the case of parallel vectors 
        // (both in the same direction or pointing in opposite directions). 
        // crossproduct will not be valid in these cases, so you first need 
        // to check dot(v1, v2) > 0.999999 and dot(v1, v2) < -0.999999, respectively, 
        // and either return an identity quat for parallel vectors, or return a 
        // 180 degree rotation (about any axis) for opposite vectors. 
        quaternion retquat;
        retquat.qw =    sqrt(pow(pow(arg_i.x,2)+pow(arg_i.y,2)+pow(arg_i.z,2),2)*(pow(pow(arg_f.x,2)+pow(arg_f.y,2)+pow(arg_f.z,2),2)));
                        +arg_i.x*arg_f.x+arg_i.y*arg_f.y+arg_i.z*arg_f.z;

        retquat.qx =    arg_i.y*arg_f.z-arg_i.z*arg_f.y;
        retquat.qy =    arg_i.z*arg_f.x-arg_i.x*arg_f.z;
        retquat.qz =    arg_i.x*arg_f.y-arg_i.y*arg_f.x;
        return retquat;
    };

    quaternion quaternionderivative(quaternion argquat, bodyframe bodyrate){
        quaternion dotquat;
        argquat.qw=sqrt(1-pow(argquat.qx,2)-pow(argquat.qy,2)-pow(argquat.qz,2));
        //dotquat.qw=0.5*(0-bodyrate.x*argquat.qx-bodyrate.y*argquat.qy-bodyrate.z*argquat.qz);
        dotquat.qx=0.5*(argquat.qw*bodyrate.x-argquat.qz*bodyrate.y+argquat.qy*bodyrate.z);
        dotquat.qy=0.5*(argquat.qz*bodyrate.x+argquat.qw*bodyrate.y-argquat.qx*bodyrate.z);
        dotquat.qz=0.5*(-argquat.qy*bodyrate.x+argquat.qx*bodyrate.y-argquat.qw*bodyrate.z);
    };
};
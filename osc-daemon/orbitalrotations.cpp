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
        retquat.qw =    sqrt(pow2(pow2(arg_i.x)+pow2(arg_i.y)+pow2(arg_i.z))*(pow2(pow2(arg_f.x)+pow2(arg_f.y)+pow2(arg_f.z))));
                        +arg_i.x*arg_f.x+arg_i.y*arg_f.y+arg_i.z*arg_f.z;

        retquat.qx =    arg_i.y*arg_f.z-arg_i.z*arg_f.y;
        retquat.qy =    arg_i.z*arg_f.x-arg_i.x*arg_f.z;
        retquat.qz =    arg_i.x*arg_f.y-arg_i.y*arg_f.x;
        return retquat;
    };

    quaternion quaternionderivative(quaternion argquat, vec3 bodyrate){
        quaternion dotquat;
        argquat.qw=sqrt(1-pow2(argquat.qx)-pow2(argquat.qy)-pow2(argquat.qz));
        //dotquat.qw=0.5*(0-bodyrate[0]*argquat.qx-bodyrate.y*argquat.qy-bodyrate.z*argquat.qz);
        dotquat.qx=0.5*(argquat.qw*bodyrate[0]-argquat.qz*bodyrate[1]+argquat.qy*bodyrate[2]);
        dotquat.qy=0.5*(argquat.qz*bodyrate[0]+argquat.qw*bodyrate[1]-argquat.qx*bodyrate[2]);
        dotquat.qz=0.5*(-argquat.qy*bodyrate[0]+argquat.qx*bodyrate[1]-argquat.qw*bodyrate[2]);
    };
};
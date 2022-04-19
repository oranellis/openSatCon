#include "../../osc-daemon/axistransforms.cpp"

double deg2rad(double deg) {
    return deg/180*M_PI;
};

double rad2deg(double rad) {
    return rad*180/M_PI;
};

int main() {
    osc::orbParam testKOE;
    testKOE.sma = 6786230;
    testKOE.ecc = 0.01;
    testKOE.inc = deg2rad(52);
    testKOE.asc = deg2rad(95);
    testKOE.truAnom = deg2rad(300);
    osc::pcs posvelPCS = osc::KOEtoPCS(testKOE);
    std::cout << posvelPCS.rPCS.data[0] << std::endl << posvelPCS.rPCS.data[1] << std::endl << posvelPCS.rPCS.data[2] << std::endl; 
    std::cout << std::endl;
    std::cout << posvelPCS.rPCS.data[0] << std::endl << posvelPCS.rPCS.data[1] << std::endl << posvelPCS.rPCS.data[2] << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    osc::eci testPosVel;
    testPosVel.rIJK.data[0]=-2981784;
    testPosVel.rIJK.data[1]=5207055;
    testPosVel.rIJK.data[2]=3161595;
    testPosVel.vIJK.data[0]=-3384;
    testPosVel.vIJK.data[1]=-4887;
    testPosVel.vIJK.data[2]=4843;
    osc::orbParam outKOE = osc::ECItoKOE(testPosVel);
    std::cout << outKOE.sma << std::endl << outKOE.ecc << std::endl << rad2deg(outKOE.inc) << std::endl; 
    std::cout << std::endl;
    std::cout << rad2deg(outKOE.asc) << std::endl << rad2deg(outKOE.aop) << std::endl << rad2deg(outKOE.truAnom) << std::endl; 
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    return 0;

}




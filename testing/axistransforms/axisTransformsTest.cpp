#include "../../osc-daemon/orbitalmechanics/axistransforms.hpp"

double deg2rad(double deg) {
    return deg/180*M_PI;
};

double rad2deg(double rad) {
    return rad*180/M_PI;
};

int main() {
    osc::orbParam testKOE;
        testKOE.sma     = 6786230;
        testKOE.ecc     = 0.01;
        testKOE.inc     = deg2rad(52);
        testKOE.asc     = deg2rad(95);
        testKOE.aop     = deg2rad(93);
        testKOE.truAnom = deg2rad(300);

    osc::ecef testECEF;
        testECEF.rXYZ.data[0] = 10766080.3;
        testECEF.rXYZ.data[1] = 14143607;
        testECEF.rXYZ.data[2] = 33992388;
        testECEF.vXYZ.data[0] = 3832;
        testECEF.vXYZ.data[1] = -4024;
        testECEF.vXYZ.data[2] = 4837;

    osc::eci testECI;
        testECI.rIJK.data[0] = -2981784;
        testECI.rIJK.data[1] = 5207055;
        testECI.rIJK.data[2] = 3161595;
        testECI.vIJK.data[0] = -3384;
        testECI.vIJK.data[1] = -4887;
        testECI.vIJK.data[2] = 4843;

    osc::lla testLLA; //glasgow uni position
        // testLLA.lat = deg2rad(55.8724);
        // testLLA.lon = deg2rad(4.2900);
        // testLLA.alt = 38;
        testLLA.lat = deg2rad(45.9132);
        testLLA.lon = deg2rad(36.7484);
        testLLA.alt = 1877753;

    osc::pcs outPCS = osc::KOEtoPCS(testKOE);
        std::cout << "testing KOE to PCS:" << std::endl;
        std::cout << outPCS.rPCS.data[0] << std::endl << outPCS.rPCS.data[1] << std::endl << outPCS.rPCS.data[2] << std::endl; 
        std::cout << std::endl;
        std::cout << outPCS.vPCS.data[0] << std::endl << outPCS.vPCS.data[1] << std::endl << outPCS.vPCS.data[2] << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    
    osc::orbParam outKOE = osc::ECItoKOE(testECI);
        std::cout << "testing ECI to KOE:" << std::endl;
        std::cout << "a=" << outKOE.sma << "m" << std::endl << "e=" << outKOE.ecc << std::endl << "inc=" << rad2deg(outKOE.inc) << "deg" << std::endl; 
        std::cout << std::endl;
        std::cout << "RAAN=" << rad2deg(outKOE.asc) << "deg" << std::endl << "argp=" << rad2deg(outKOE.aop) << "deg" << std::endl << "theta=" << rad2deg(outKOE.truAnom) << "deg" << std::endl; 
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;

    double sideRealAngle = osc::greenwichSiderealAngle();
        std::cout << "testing Sidereal Angle:" << std::endl << rad2deg(sideRealAngle) << "rad" << std::endl << std::endl << std::endl;

    osc::ecef outECEF = osc::ECItoECEF(testECI, 0.4971);
            std::cout << "testing ECI to ECEF:" << std::endl;
            std::cout << outECEF.rXYZ.data[0] << std::endl << outECEF.rXYZ.data[1] << std::endl << outECEF.rXYZ.data[2] << std::endl; 
            std::cout << std::endl;
            std::cout << outECEF.vXYZ.data[0] << std::endl << outECEF.vXYZ.data[1] << std::endl << outECEF.vXYZ.data[2] << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;

    
    osc::eci outECI = osc::ECEFtoECI(outECEF, 0.4971);
            std::cout << "testing ECEF to ECI:" << std::endl;
            std::cout << outECI.rIJK.data[0] << std::endl << outECI.rIJK.data[1] << std::endl << outECI.rIJK.data[2] << std::endl; 
            std::cout << std::endl;
            std::cout << outECI.vIJK.data[0] << std::endl << outECI.vIJK.data[1] << std::endl << outECI.vIJK.data[2] << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;

    osc::eci outECI2 = osc::PCStoECI(testKOE, outPCS);
            std::cout << "testing PCS to ECI:" << std::endl;
            std::cout << outECI2.rIJK.data[0] << std::endl << outECI2.rIJK.data[1] << std::endl << outECI2.rIJK.data[2] << std::endl; 
            std::cout << std::endl;
            std::cout << outECI2.vIJK.data[0] << std::endl << outECI2.vIJK.data[1] << std::endl << outECI2.vIJK.data[2] << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;

    osc::ecef outECEF2 = osc::LLAtoECEF(testLLA);
            std::cout << "testing LLA to ECEF:" << std::endl;
            std::cout << outECEF2.rXYZ.data[0] << std::endl << outECEF2.rXYZ.data[1] << std::endl << outECEF2.rXYZ.data[2] << std::endl; 
            std::cout << std::endl;

    osc::lla outLLA = osc::ECEFtoLLA(testECEF);
            std::cout << "testing ECEF to LLA:" << std::endl;
            std::cout << rad2deg(outLLA.lat) << std::endl << rad2deg(outLLA.lon) << std::endl << outLLA.alt << std::endl; 
            std::cout << std::endl;

    osc::ned outNED = osc::ECEFtoNED(testECEF, testLLA);
            std::cout << "testing ECEF to NED:" << std::endl;
            std::cout << outNED.rNED.data[0] << std::endl << outNED.rNED.data[1] << std::endl << outNED.rNED.data[2] << std::endl; 
            std::cout << std::endl;

    osc::ecef outECEF3 = osc::NEDtoECEF(outNED, testLLA);
            std::cout << "testing NED to ECEF:" << std::endl;
            std::cout << outECEF3.rXYZ.data[0] << std::endl << outECEF3.rXYZ.data[1] << std::endl << outECEF3.rXYZ.data[2] << std::endl; 
            std::cout << std::endl;

    osc::ear outEAR = osc::ECEFtoEAR(testECEF, testLLA);
            std::cout << "testing ECEF to EAR:" << std::endl;
            std::cout << rad2deg(outEAR.e) << std::endl << rad2deg(outEAR.a) << std::endl << outEAR.r << std::endl; 
            std::cout << std::endl;

    osc::ecef outECEF4 = osc::EARtoECEF(outEAR, testLLA);
            std::cout << "testing EAR to ECEF:" << std::endl;
            std::cout << outECEF4.rXYZ.data[0] << std::endl << outECEF4.rXYZ.data[1] << std::endl << outECEF4.rXYZ.data[2] << std::endl; 
            std::cout << std::endl;

    osc::enu outENU = osc::ECEFtoENU(testECEF, testLLA);
            std::cout << "testing ECEF to ENU:" << std::endl;
            std::cout << outENU.rENU.data[0] << std::endl << outENU.rENU.data[1] << std::endl << outENU.rENU.data[2] << std::endl; 
            std::cout << std::endl;

    osc::ecef outECEF5 = osc::EARtoECEF(outEAR, testLLA);
            std::cout << "testing ENU to ECEF:" << std::endl;
            std::cout << outECEF5.rXYZ.data[0] << std::endl << outECEF5.rXYZ.data[1] << std::endl << outECEF5.rXYZ.data[2] << std::endl; 
            std::cout << std::endl;
            
    return 0;
}




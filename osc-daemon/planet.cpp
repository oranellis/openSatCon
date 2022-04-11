#include <vector>
#include <array>
#include <math.h>

//have this editable so we can orbit Mars as well
struct celestial{
    double gravity = 9.81;
    double semimajoraxis = 6378137;
    double flattening = 1/298.257223563;
    double semiminoraxis = semimajoraxis-semimajoraxis*flattening;
    double eccentricity = sqrt(2*flattening-pow(flattening,2));
};
celestial planet;
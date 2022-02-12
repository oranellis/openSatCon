#include <iostream>
#include <vector>
#include <string>
#include <math.h>



#define PI=3.14159265
int main(){

  double x, y, result;
  x = -10.0;
  y = 10.0;
  result = atan2 (y,x) * 180 / 3.14159265;
  printf ("The arc tangent for (x=%f, y=%f) is %f degrees\n", x, y, result );
  return 0;


}
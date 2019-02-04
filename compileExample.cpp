
// compileExample.cpp

#include "include.h"

int main2(int argc, const char* argv[] )
{
//unsigned int old_cw;
//fpu_fix_start(&old_cw);

  cout.precision(60);

//qd_real readTest;
//cin >> readTest;
//cout << "readTest = " << readTest << endl;

  qd_real x = "1.0";
  x /= 3.0;

  qd_real y;
  y = pow( qd_real(2.0) , 3);

  cout << "y = " << y << endl;
  cout << "x = " << x << endl;


  qd_real a;
  qd_real b = qd_real("0.1");

  a = sqrt(b);

  cout << " sqrt(0.1) = " << a << endl;
  cout << " sqrt(0.1) * sqrt(0.1) = " << a * a << endl;

//fpu_fix_end(&old_cw);

  return 0;
}



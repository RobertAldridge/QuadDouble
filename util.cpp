
// util.cpp

#include "include.h"

void append_expn(std::string& str, int expn)
{
  int k = 0;

  str += (expn < 0 ? '-' : '+');

  expn = std::abs(expn);

  if(expn >= 100)
  {
    k = (expn / 100);

    str += '0' + k;

    expn -= 100 * k;
  }

  k = (expn / 10);

  str += '0' + k;

  expn -= 10 * k;

  str += '0' + expn;
}

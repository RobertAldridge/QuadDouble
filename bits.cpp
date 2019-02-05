
// bits.cpp

#include "include.h"

int get_double_expn(double x)
{
  if(x == 0.0)
    return INT_MIN;
  if(QD_ISINF(x) || QD_ISNAN(x) )
    return INT_MAX;

  double y = std::abs(x);
  int i = 0;
  if(y < 1.0)
  {
    while(y < 1.0)
    {
      y *= 2.0;
      i++;
    }
    return -i;
  }
  else if(y >= 2.0)
  {
    while(y >= 2.0)
    {
      y *= 0.5;
      i++;
    }
    return i;
  }
  return 0;
}

void print_double_info(std::ostream& os, double x)
{
  std::streamsize old_prec = os.precision(19);
  std::ios_base::fmtflags old_flags  = os.flags();
  os << std::scientific;

  os << setw(27) << x << ' ';
  if(QD_ISNAN(x) || QD_ISINF(x) || (x == 0.0) )
  {
    os << "                                                           ";
  }
  else
  {
    x = std::abs(x);
    int expn = get_double_expn(x);
    double d = std::ldexp(1.0, expn);
    os << setw(5) << expn << " ";
    for(int i = 0; i < 53; i++)
    {
      if(x >= d)
      {
        x -= d;
        os << '1';
      }
      else
      {
        os << '0';
      }

      d *= 0.5;
    }

    if(x != 0.0)
    {
      os << " +trailing stuff";
    }
  }

  os.precision(old_prec);
  os.flags(old_flags);
}

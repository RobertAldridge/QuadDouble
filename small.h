
// small.h

#define _QD_SPLITTER 134217729.0

#define _QD_SPLIT_THRESH 6.69692879491417e+299

namespace qd
{
static const double _d_nan = std::numeric_limits<double>::quiet_NaN();

static const double _d_inf = std::numeric_limits<double>::infinity();

double quick_two_sum(double a, double b, double& err)
{
  double s = a + b;

  err = b - (s - a);

  return s;
}

double quick_two_diff(double a, double b, double& err)
{
  double s = a - b;

  err = (a - s) - b;

  return s;
}

double two_sum(double a, double b, double& err)
{
  double s = a + b;

  double bb = s - a;

  err = (a - (s - bb) ) + (b - bb);

  return s;
}

double two_diff(double a, double b, double& err)
{
  double s = a - b;

  double bb = s - a;

  err = (a - (s - bb) ) - (b + bb);

  return s;
}

#if !defined(QD_FMS)
void split(double a, double& hi, double& lo)
{
  double temp = 0.0;

  if(a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH)
  {
    a *= 3.7252902984619140625e-09;

    temp = _QD_SPLITTER * a;

    hi = temp - (temp - a);

    lo = a - hi;

    hi *= 268435456.0;

    lo *= 268435456.0;
  }
  else
  {
    temp = _QD_SPLITTER * a;

    hi = temp - (temp - a);

    lo = a - hi;
  }
}
#endif

double two_prod(double a, double b, double& err)
{
#if defined(QD_FMS)
  double p = a * b;

  err = QD_FMS(a, b, p);

  return p;
#else
  double a_hi = 0.0;

  double a_lo = 0.0;

  double b_hi = 0.0;

  double b_lo = 0.0;

  double p = a * b;

  split(a, a_hi, a_lo);

  split(b, b_hi, b_lo);

  err = ( (a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;

  return p;
#endif
}

double two_sqr(double a, double& err)
{
#if defined(QD_FMS)
  double p = a * a;

  err = QD_FMS(a, a, p);

  return p;
#else
  double hi = 0.0;

  double lo = 0.0;

  double q = a * a;

  split(a, hi, lo);

  err = ( (hi * hi - q) + 2.0 * hi * lo) + lo * lo;

  return q;
#endif
}

double nint(double d)
{
  if(d == std::floor(d) )
  {
    return d;
  }

  return std::floor(d + 0.5);
}

double aint(double d)
{
  return (d >= 0.0) ? std::floor(d) : std::ceil(d);
}

void sincosh(double t, double& sinh_t, double& cosh_t)
{
  sinh_t = std::sinh(t);

  cosh_t = std::cosh(t);
}

double sqr(double t)
{
  return t * t;
}

double to_double(double a)
{
  return a;
}

int to_int(double a)
{
  return static_cast<int>(a);
}

}

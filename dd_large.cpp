
// dd_large.cpp

#include "include.h"

void dd_real::error(const char* msg)
{
  if(msg)
  {
    cerr << "ERROR " << msg << endl;
  }
}

dd_real sqrt(const dd_real& a)
{
  if(a.is_zero() )
  {
    return 0.0;
  }

  if(a.is_negative() )
  {
    dd_real::error("(dd_real::sqrt): Negative argument.");

    return dd_real::_nan;
  }

  double x = 1.0 / std::sqrt(a.x[0] );

  double ax = a.x[0] * x;

  return dd_real::add(ax, (a - dd_real::sqr(ax) ).x[0] * (x * 0.5) );
}

dd_real dd_real::sqrt(double d)
{
  return /*?illuminati?*/::sqrt(dd_real(d) );
}

dd_real nroot(const dd_real& a, int n)
{
  if(n <= 0)
  {
    dd_real::error("(dd_real::nroot): N must be positive.");

    return dd_real::_nan;
  }

  if( (n % 2) == 0 && a.is_negative() )
  {
    dd_real::error("(dd_real::nroot): Negative argument.");

    return dd_real::_nan;
  }

  if(n == 1)
  {
    return a;
  }

  if(n == 2)
  {
    return sqrt(a);
  }

  if(a.is_zero() )
  {
    return 0.0;
  }

  dd_real r = abs(a);

  dd_real x = std::exp(-std::log(r.x[0] ) / n);

  x += x * (1.0 - r * npwr(x, n) ) / static_cast<double>(n);

  if(a.x[0] < 0.0)
  {
    x = -x;
  }

  return 1.0 / x;
}

dd_real npwr(const dd_real& a, int n)
{
  if(n == 0)
  {
    if(a.is_zero() )
    {
      dd_real::error("(dd_real::npwr): Invalid argument.");

      return dd_real::_nan;
    }

    return 1.0;
  }

  dd_real r = a;

  dd_real s = 1.0;

  int N = std::abs(n);

  if(N > 1)
  {
    while(N > 0)
    {
      if(N % 2 == 1)
      {
        s *= r;
      }

      N /= 2;

      if(N > 0)
      {
        r = sqr(r);
      }
    }
  }
  else
  {
    s = r;
  }

  if(n < 0)
  {
    return (1.0 / s);
  }

  return s;
}

dd_real pow(const dd_real& a, int n)
{
  return npwr(a, n);
}

dd_real pow(const dd_real& a, const dd_real& b)
{
  return exp(b * log(a) );
}

static const int n_inv_fact = 15;

static const double inv_fact[n_inv_fact][2] =
{
  {
    1.66666666666666657e-01,
    9.25185853854297066e-18
  },
  {
    4.16666666666666644e-02,
    2.31296463463574266e-18
  },
  {
    8.33333333333333322e-03,
    1.15648231731787138e-19
  },
  {
    1.38888888888888894e-03,
    -5.30054395437357706e-20
  },
  {
    1.98412698412698413e-04,
    1.72095582934207053e-22
  },
  {
    2.48015873015873016e-05,
    2.15119478667758816e-23
  },
  {
    2.75573192239858925e-06,
    -1.85839327404647208e-22
  },
  {
    2.75573192239858883e-07,
    2.37677146222502973e-23
  },
  {
    2.50521083854417202e-08,
    -1.44881407093591197e-24
  },
  {
    2.08767569878681002e-09,
    -1.20734505911325997e-25
  },
  {
    1.60590438368216133e-10,
    1.25852945887520981e-26
  },
  {
    1.14707455977297245e-11,
    2.06555127528307454e-28
  },
  {
    7.64716373181981641e-13,
    7.03872877733453001e-30
  },
  {
    4.77947733238738525e-14,
    4.39920548583408126e-31
  },
  {
    2.81145725434552060e-15,
    1.65088427308614326e-31
  }
};

dd_real exp(const dd_real& a)
{
  const double k = 512.0;

  const double inv_k = 1.0 / k;

  if(a.x[0] <= -709.0)
  {
    return 0.0;
  }

  if(a.x[0] >=  709.0)
  {
    return dd_real::_inf;
  }

  if(a.is_zero() )
  {
    return 1.0;
  }

  if(a.is_one() )
  {
    return dd_real::_e;
  }

  double m = std::floor(a.x[0] / dd_real::_log2.x[0] + 0.5);

  dd_real r = mul_pwr2(a - dd_real::_log2 * m, inv_k);

  dd_real s, t, p;

  p = sqr(r);

  s = r + mul_pwr2(p, 0.5);

  p *= r;

  t = p * dd_real(inv_fact[0][0], inv_fact[0][1] );

  int i = 0;

  do
  {
    s += t;

    p *= r;

    ++i;

    t = p * dd_real(inv_fact[i][0], inv_fact[i][1] );
  }
  while(std::abs(to_double(t) ) > inv_k * dd_real::_eps && i < 5);

  s += t;

  s = mul_pwr2(s, 2.0) + sqr(s);

  s = mul_pwr2(s, 2.0) + sqr(s);

  s = mul_pwr2(s, 2.0) + sqr(s);

  s = mul_pwr2(s, 2.0) + sqr(s);

  s = mul_pwr2(s, 2.0) + sqr(s);

  s = mul_pwr2(s, 2.0) + sqr(s);

  s = mul_pwr2(s, 2.0) + sqr(s);

  s = mul_pwr2(s, 2.0) + sqr(s);

  s = mul_pwr2(s, 2.0) + sqr(s);

  s += 1.0;

  return ldexp(s, static_cast<int>(m) );
}

dd_real log(const dd_real& a)
{
  if(a.is_one() )
  {
    return 0.0;
  }

  if(a.x[0] <= 0.0)
  {
    dd_real::error("(dd_real::log): Non-positive argument.");

    return dd_real::_nan;
  }

  dd_real x = std::log(a.x[0] );

  x = x + a * exp(-x) - 1.0;

  return x;
}

dd_real log10(const dd_real& a)
{
  return log(a) / dd_real::_log10;
}

static const dd_real _pi16 = dd_real(1.963495408493620697e-01, 7.654042494670957545e-18);

static const double sin_table [4][2] =
{
  {
    1.950903220161282758e-01,
    -7.991079068461731263e-18
  },
  {
    3.826834323650897818e-01,
    -1.005077269646158761e-17
  },
  {
    5.555702330196021776e-01,
    4.709410940561676821e-17
  },
  {
    7.071067811865475727e-01,
    -4.833646656726456726e-17
  }
};

static const double cos_table [4][2] =
{
  {
    9.807852804032304306e-01,
    1.854693999782500573e-17
  },
  {
    9.238795325112867385e-01,
    1.764504708433667706e-17
  },
  {
    8.314696123025452357e-01,
    1.407385698472802389e-18
  },
  {
    7.071067811865475727e-01,
    -4.833646656726456726e-17
  }
};

static dd_real sin_taylor(const dd_real& a)
{
  const double thresh = 0.5 * std::abs(to_double(a) ) * dd_real::_eps;

  dd_real r;

  dd_real s;

  dd_real t;

  dd_real x;

  if(a.is_zero() )
  {
    return 0.0;
  }

  int i = 0;

  x = -sqr(a);

  s = a;

  r = a;

  do
  {
    r *= x;

    t = r * dd_real(inv_fact[i][0], inv_fact[i][1] );

    s += t;

    i += 2;
  }
  while(i < n_inv_fact && std::abs(to_double(t) ) > thresh);

  return s;
}

static dd_real cos_taylor(const dd_real& a)
{
  const double thresh = 0.5 * dd_real::_eps;

  dd_real r;

  dd_real s;

  dd_real t;

  dd_real x;

  if(a.is_zero() )
  {
    return 1.0;
  }

  x = -sqr(a);

  r = x;

  s = 1.0 + mul_pwr2(r, 0.5);

  int i = 1;

  do
  {
    r *= x;

    t = r * dd_real(inv_fact[i][0], inv_fact[i][1] );

    s += t;

    i += 2;
  }
  while(i < n_inv_fact && std::abs(to_double(t) ) > thresh);

  return s;
}

static void sincos_taylor(const dd_real& a, dd_real& sin_a, dd_real& cos_a)
{
  if(a.is_zero() )
  {
    sin_a = 0.0;

    cos_a = 1.0;

    return;
  }

  sin_a = sin_taylor(a);

  cos_a = sqrt(1.0 - sqr(sin_a) );
}


dd_real sin(const dd_real& a)
{
  if(a.is_zero() )
  {
    return 0.0;
  }

  dd_real z = nint(a / dd_real::_2pi);

  dd_real r = a - dd_real::_2pi * z;

  dd_real t;

  double q = std::floor(r.x[0] / dd_real::_pi2.x[0] + 0.5);

  t = r - dd_real::_pi2 * q;

  int j = static_cast<int>(q);

  q = std::floor(t.x[0] / _pi16.x[0] + 0.5);

  t -= _pi16 * q;

  int k = static_cast<int>(q);

  int abs_k = std::abs(k);

  if(j < -2 || j > 2)
  {
    dd_real::error("(dd_real::sin): Cannot reduce modulo pi/2.");

    return dd_real::_nan;
  }

  if(abs_k > 4)
  {
    dd_real::error("(dd_real::sin): Cannot reduce modulo pi/16.");

    return dd_real::_nan;
  }

  if(k == 0)
  {
    switch(j)
    {

      case 0:
      {
        return sin_taylor(t);
      }
      break;

      case 1:
      {
        return cos_taylor(t);
      }
      break;

      case -1:
      {
        return -cos_taylor(t);
      }
      break;

      default:
      {
        return -sin_taylor(t);
      }
      break;

    }
  }

  dd_real u(cos_table[abs_k-1][0], cos_table[abs_k-1][1] );

  dd_real v(sin_table[abs_k-1][0], sin_table[abs_k-1][1] );

  dd_real sin_t, cos_t;

  sincos_taylor(t, sin_t, cos_t);

  if(j == 0)
  {
    if(k > 0)
    {
      r = u * sin_t + v * cos_t;
    }
    else
    {
      r = u * sin_t - v * cos_t;
    }
  }
  else if(j == 1)
  {
    if(k > 0)
    {
      r = u * cos_t - v * sin_t;
    }
    else
    {
      r = u * cos_t + v * sin_t;
    }
  }
  else if(j == -1)
  {
    if(k > 0)
    {
      r = v * sin_t - u * cos_t;
    }
    else if(k < 0)
    {
      r = -u * cos_t - v * sin_t;
    }
  }
  else
  {
    if(k > 0)
    {
      r = -u * sin_t - v * cos_t;
    }
    else
    {
      r = v * cos_t - u * sin_t;
    }
  }

  return r;
}

dd_real cos(const dd_real& a)
{
  if(a.is_zero() )
  {
    return 1.0;
  }

  dd_real z = nint(a / dd_real::_2pi);

  dd_real r = a - z * dd_real::_2pi;

  dd_real t;

  double q = std::floor(r.x[0] / dd_real::_pi2.x[0] + 0.5);

  t = r - dd_real::_pi2 * q;

  int j = static_cast<int>(q);

  q = std::floor(t.x[0] / _pi16.x[0] + 0.5);

  t -= _pi16 * q;

  int k = static_cast<int>(q);

  int abs_k = std::abs(k);

  if(j < -2 || j > 2)
  {
    dd_real::error("(dd_real::cos): Cannot reduce modulo pi/2.");

    return dd_real::_nan;
  }

  if(abs_k > 4)
  {
    dd_real::error("(dd_real::cos): Cannot reduce modulo pi/16.");

    return dd_real::_nan;
  }

  if(k == 0)
  {
    switch(j)
    {

    case 0:
    {
      return cos_taylor(t);
    }
    break;

    case 1:
    {
      return -sin_taylor(t);
    }
    break;

    case -1:
    {
      return sin_taylor(t);
    }
    break;

    default:
    {
      return -cos_taylor(t);
    }
    break;

    }
  }

  dd_real sin_t;

  dd_real cos_t;

  sincos_taylor(t, sin_t, cos_t);

  dd_real u(cos_table[abs_k - 1][0], cos_table[abs_k - 1][1] );

  dd_real v(sin_table[abs_k - 1][0], sin_table[abs_k - 1][1] );

  if(j == 0)
  {
    if(k > 0)
    {
      r = u * cos_t - v * sin_t;
    }
    else
    {
      r = u * cos_t + v * sin_t;
    }
  }
  else if(j == 1)
  {
    if(k > 0)
    {
      r = - u * sin_t - v * cos_t;
    }
    else
    {
      r = v * cos_t - u * sin_t;
    }
  }
  else if(j == -1)
  {
    if(k > 0)
    {
      r = u * sin_t + v * cos_t;
    }
    else
    {
      r = u * sin_t - v * cos_t;
    }
  }
  else
  {
    if(k > 0)
    {
      r = v * sin_t - u * cos_t;
    }
    else
    {
      r = - u * cos_t - v * sin_t;
    }
  }

  return r;
}

void sincos(const dd_real& a, dd_real& sin_a, dd_real& cos_a)
{
  if(a.is_zero() )
  {
    sin_a = 0.0;

    cos_a = 1.0;

    return;
  }

  dd_real z = nint(a / dd_real::_2pi);

  dd_real r = a - dd_real::_2pi * z;

  dd_real t;

  double q = std::floor(r.x[0] / dd_real::_pi2.x[0] + 0.5);

  t = r - dd_real::_pi2 * q;

  int j = static_cast<int>(q);

  int abs_j = std::abs(j);

  q = std::floor(t.x[0] / _pi16.x[0] + 0.5);

  t -= _pi16 * q;

  int k = static_cast<int>(q);

  int abs_k = std::abs(k);

  if(abs_j > 2)
  {
    dd_real::error("(dd_real::sincos): Cannot reduce modulo pi/2.");

    cos_a = sin_a = dd_real::_nan;

    return;
  }

  if(abs_k > 4)
  {
    dd_real::error("(dd_real::sincos): Cannot reduce modulo pi/16.");

    cos_a = sin_a = dd_real::_nan;

    return;
  }

  dd_real sin_t;

  dd_real cos_t;

  dd_real s;

  dd_real c;

  sincos_taylor(t, sin_t, cos_t);

  if(abs_k == 0)
  {
    s = sin_t;

    c = cos_t;
  }
  else
  {
    dd_real u(cos_table[abs_k - 1][0], cos_table[abs_k - 1][1] );

    dd_real v(sin_table[abs_k - 1][0], sin_table[abs_k - 1][1] );

    if(k > 0)
    {
      s = u * sin_t + v * cos_t;

      c = u * cos_t - v * sin_t;
    }
    else
    {
      s = u * sin_t - v * cos_t;

      c = u * cos_t + v * sin_t;
    }
  }

  if(abs_j == 0)
  {
    sin_a = s;

    cos_a = c;
  }
  else if(j == 1)
  {
    sin_a = c;

    cos_a = -s;
  }
  else if(j == -1)
  {
    sin_a = -c;

    cos_a = s;
  }
  else
  {
    sin_a = -s;

    cos_a = -c;
  }
}

dd_real atan(const dd_real& a)
{
  return atan2(a, dd_real(1.0) );
}

dd_real atan2(const dd_real& y, const dd_real& x)
{
  if(x.is_zero() )
  {
    if(y.is_zero() )
    {
      dd_real::error("(dd_real::atan2): Both arguments zero.");

      return dd_real::_nan;
    }

    return y.is_positive()? dd_real::_pi2 : -dd_real::_pi2;
  }
  else if(y.is_zero() )
  {
    return x.is_positive()? dd_real(0.0) : dd_real::_pi;
  }

  if(x == y)
  {
    return y.is_positive()? dd_real::_pi4 : -dd_real::_3pi4;
  }

  if(x == -y)
  {
    return y.is_positive()? dd_real::_3pi4 : -dd_real::_pi4;
  }

  dd_real r = sqrt(sqr(x) + sqr(y) );

  dd_real xx = x / r;

  dd_real yy = y / r;

  dd_real z = std::atan2(to_double(y), to_double(x) );

  dd_real sin_z, cos_z;

  if(std::abs(xx.x[0] ) > std::abs(yy.x[0] ) )
  {
    sincos(z, sin_z, cos_z);

    z += (yy - sin_z) / cos_z;
  }
  else
  {
    sincos(z, sin_z, cos_z);

    z -= (xx - cos_z) / sin_z;
  }

  return z;
}

dd_real tan(const dd_real& a)
{
  dd_real s;

  dd_real c;

  sincos(a, s, c);

  return s / c;
}

dd_real asin(const dd_real& a)
{
  dd_real abs_a = abs(a);

  if(abs_a > 1.0)
  {
    dd_real::error("(dd_real::asin): Argument out of domain.");

    return dd_real::_nan;
  }

  if(abs_a.is_one() )
  {
    return a.is_positive()? dd_real::_pi2 : -dd_real::_pi2;
  }

  return atan2(a, sqrt(1.0 - sqr(a) ) );
}

dd_real acos(const dd_real& a)
{
  dd_real abs_a = abs(a);

  if(abs_a > 1.0)
  {
    dd_real::error("(dd_real::acos): Argument out of domain.");

    return dd_real::_nan;
  }

  if(abs_a.is_one() )
  {
    return a.is_positive()? dd_real(0.0) : dd_real::_pi;
  }

  return atan2(sqrt(1.0 - sqr(a) ), a);
}

dd_real sinh(const dd_real& a)
{
  if(a.is_zero() )
  {
    return 0.0;
  }

  if(abs(a) > 0.05)
  {
    dd_real ea = exp(a);

    return mul_pwr2(ea - inv(ea), 0.5);
  }

  dd_real s = a;

  dd_real t = a;

  dd_real r = sqr(t);

  double m = 1.0;

  double thresh = std::abs( (to_double(a) ) * dd_real::_eps);

  do
  {
    m += 2.0;

    t *= r;

    t /= (m - 1) * m;

    s += t;
  }
  while(abs(t) > thresh);

  return s;

}

dd_real cosh(const dd_real& a)
{
  if(a.is_zero() )
  {
    return 1.0;
  }

  dd_real ea = exp(a);

  return mul_pwr2(ea + inv(ea), 0.5);
}

dd_real tanh(const dd_real& a)
{
  if(a.is_zero() )
  {
    return 0.0;
  }

  if(std::abs(to_double(a) ) > 0.05)
  {
    dd_real ea = exp(a);

    dd_real inv_ea = inv(ea);

    return (ea - inv_ea) / (ea + inv_ea);
  }
  else
  {
    dd_real s, c;

    s = sinh(a);

    c = sqrt(1.0 + sqr(s) );

    return s / c;
  }
}

void sincosh(const dd_real& a, dd_real& s, dd_real& c)
{
  if(std::abs(to_double(a) ) <= 0.05)
  {
    s = sinh(a);

    c = sqrt(1.0 + sqr(s) );
  }
  else
  {
    dd_real ea = exp(a);

    dd_real inv_ea = inv(ea);

    s = mul_pwr2(ea - inv_ea, 0.5);

    c = mul_pwr2(ea + inv_ea, 0.5);
  }
}

dd_real asinh(const dd_real& a)
{
  return log(a + sqrt(sqr(a) + 1.0) );
}

dd_real acosh(const dd_real& a)
{
  if(a < 1.0)
  {
    dd_real::error("(dd_real::acosh): Argument out of domain.");

    return dd_real::_nan;
  }

  return log(a + sqrt(sqr(a) - 1.0) );
}

dd_real atanh(const dd_real& a)
{
  if(abs(a) >= 1.0)
  {
    dd_real::error("(dd_real::atanh): Argument out of domain.");

    return dd_real::_nan;
  }

  return mul_pwr2(log( (1.0 + a) / (1.0 - a) ), 0.5);
}

dd_real fmod(const dd_real& a, const dd_real& b)
{
  dd_real n = aint(a / b);

  return (a - b * n);
}

dd_real ddrand()
{
  static const double m_const = 4.6566128730773926e-10;

  double m = m_const;

  dd_real r = 0.0;

  double d;

  for(int i = 0; i < 4; i++, m *= m_const)
  {
/*d = lrand48() * m;*/

    d = std::rand() * m;

    r += d;
  }

  return r;
}

dd_real polyeval(const dd_real* c, int n, const dd_real& x)
{
  dd_real r = c[n];

  for(int i = n-1; i >= 0; i--)
  {
    r *= x;

    r += c[i];
  }

  return r;
}

dd_real polyroot(const dd_real* c, int n, const dd_real& x0, int max_iter, double thresh)
{
  dd_real x = x0;

  dd_real f;

  dd_real* d = new dd_real[n];

  bool conv = false;

  int i = 0;

  double max_c = std::abs(to_double(c[0] ) );

  double v = 0.0;

  if(thresh == 0.0)
  {
    thresh = dd_real::_eps;
  }

  for(i = 1; i <= n; i++)
  {
    v = std::abs(to_double(c[i] ) );

    if(v > max_c)
    {
      max_c = v;
    }

    d[i-1] = c[i] * static_cast<double>(i);
  }
  thresh *= max_c;

  for(i = 0; i < max_iter; i++)
  {
    f = polyeval(c, n, x);

    if(abs(f) < thresh)
    {
      conv = true;

      break;
    }

    x -= (f / polyeval(d, n-1, x) );
  }

  delete[] d;

  if( !conv)
  {
    dd_real::error("(dd_real::polyroot): Failed to converge.");

    return dd_real::_nan;
  }

  return x;
}

dd_real::dd_real(const char* s)
{
  if(dd_real::read(s, *this) )
  {
    dd_real::error("(dd_real::dd_real): INPUT ERROR.");

    *this = dd_real::_nan;
  }
}

dd_real& dd_real::operator=(const char* s)
{
  if(dd_real::read(s, *this) )
  {
    dd_real::error("(dd_real::operator=): INPUT ERROR.");

    *this = dd_real::_nan;
  }
  return *this;
}

ostream& operator<<(ostream& os, const dd_real& dd)
{
  bool showpos = ( (os.flags() & ios_base::showpos) != 0);

  bool uppercase =  ( (os.flags() & ios_base::uppercase) != 0);

  return os << dd.to_string( (int)os.precision(), (int)os.width(), os.flags(), showpos, uppercase, os.fill() );
}

istream& operator>>(istream& s, dd_real& a)
{
  char str[255] = {0};

  s >> str;

  a = dd_real(str);

  return s;
}

void dd_real::to_digits(char* s, int& expn, int precision) const
{
  int D = precision + 1;

  dd_real r = abs(*this);

  int e = 0;

  int i = 0;

  int d = 0;

  if(x[0] == 0.0)
  {
    expn = 0;

    for(i = 0; i < precision; i++)
    {
      s[i] = '0';
    }

    return;
  }

  e = to_int(std::floor(std::log10(std::abs(x[0] ) ) ) );

  if(e < -300)
  {
    r *= dd_real(10.0) ^ 300;

    r /= dd_real(10.0) ^ (e + 300);
  }
  else if(e > 300)
  {
    r = ldexp(r, -53);

    r /= dd_real(10.0) ^ e;

    r = ldexp(r, 53);
  }
  else
  {
    r /= dd_real(10.0) ^ e;
  }

  if(r >= 10.0)
  {
    r /= 10.0;

    e++;
  }
  else if(r < 1.0)
  {
    r *= 10.0;

    e--;
  }

  if(r >= 10.0 || r < 1.0)
  {
    dd_real::error("(dd_real::to_digits): can't compute exponent.");
    return;
  }

  for(i = 0; i < D; i++)
  {
    d = static_cast<int>(r.x[0] );

    r -= d;

    r *= 10.0;

    s[i] = static_cast<char>(d + '0');
  }

  for(i = D - 1; i > 0; i--)
  {
    if(s[i] < '0')
    {
      s[i - 1]--;

      s[i] += 10;
    }
    else if(s[i] > '9')
    {
      s[i - 1]++;

      s[i] -= 10;
    }
  }

  if(s[0] <= '0')
  {
    dd_real::error("(dd_real::to_digits): non-positive leading digit.");

    return;
  }

  if(s[D - 1] >= '5')
  {
    s[D - 2]++;

    i = D - 2;

    while(i > 0 && s[i] > '9')
    {
      s[i] -= 10;

      s[--i]++;
    }
  }

  if(s[0] > '9')
  {
    e++;

    for(i = precision; i >= 2; i--)
    {
      s[i] = s[i - 1];
    }

    s[0] = '1';

    s[1] = '0';
  }

  s[precision] = 0;

  expn = e;
}

void dd_real::write(char* s, int len, int precision, bool showpos, bool uppercase) const
{
  string str = to_string(precision, 0, ios_base::scientific, showpos, uppercase);

  std::strncpy(s, str.c_str(), len-1);

  s[len-1] = 0;
}


void round_string(char* s, int precision, int* offset)
{
  int i = 0;

  int D = precision;

  if(D > 0 && s[D] >= '5')
  {
    s[D - 1]++;

    i = D - 1;

    while(i > 0 && s[i] > '9')
    {
      s[i] -= 10;

      s[--i]++;
    }
  }

  if(s[0] > '9')
  {
/*e++;*/

    for(i = precision; i >= 1; i--)
    {
      s[i + 1] = s[i];
    }

    s[0] = '1';

    s[1] = '0';

    (*offset)++;

    precision++;
  }

  s[precision] = 0;
}

string dd_real::to_string(int precision, int width, ios_base::fmtflags fmt, bool showpos, bool uppercase, char fill) const
{
  string s = "";

  bool fixed = ( (fmt & ios_base::fixed) != 0);

  bool sgn = true;

  int i = 0;

  int e = 0;

  if(isnan() )
  {
    s = uppercase? "NAN" : "nan";

    sgn = false;
  }
  else
  {
    if(*this < 0.0)
    {
      s += '-';
    }
    else if(showpos)
    {
      s += '+';
    }
    else
    {
      sgn = false;
    }

    if(isinf() )
    {
      s += uppercase ? "INF" : "inf";
    }
    else if(*this == 0.0)
    {
      s += '0';

      if(precision > 0)
      {
        s += '.';

        s.append(precision, '0');
      }
    }
    else
    {
      int off = (fixed ? (1 + to_int(floor(log10(abs(*this) ) ) ) ) : 1);

      int d = precision + off;

      int d_with_extra = d;

      if(fixed)
      {
        d_with_extra = std::max(60, d);
      }

      if(fixed && (precision == 0) && (abs(*this) < 1.0) )
      {
        if(abs(*this) >= 0.5)
        {
          s += '1';
        }
        else
        {
          s += '0';
        }

        return s;
      }

      if(fixed && d <= 0)
      {
        s += '0';

        if(precision > 0)
        {
          s += '.';

          s.append(precision, '0');
        }
      }
      else
      {
        char* t = nullptr;/*= new char[d + 1];*/

        int j = 0;

        if(fixed)
        {
          t = new char[d_with_extra+1];

          to_digits(t, e, d_with_extra);
        }
        else
        {
          t = new char[d + 1];

          to_digits(t, e, d);
        }

        off = e + 1;

        if(fixed)
        {
          round_string(t, d, &off);

          if(off > 0)
          {
            for(i = 0; i < off; i++)
            {
              s += t[i];
            }

            if(precision > 0)
            {
              s += '.';

              for(j = 0; j < precision; j++, i++)
              {
                s += t[i];
              }
            }
          }
          else
          {
            s += "0.";

            if(off < 0)
            {
              s.append(-off, '0');
            }

            for(i = 0; i < d; i++)
            {
              s += t[i];
            }
          }
        }
        else
        {
          s += t[0];

          if(precision > 0)
          {
            s += '.';
          }

          for(i = 1; i <= precision; i++)
          {
            s += t[i];
          }

        }
    delete[] t;
      }
    }

    if(fixed && (precision > 0) )
    {
      double from_string = atof(s.c_str() );

      if(fabs( from_string / this->x[0] ) > 3.0 )
      {
/*int point_position = 0;*/

/*char temp = 0;*/

        for(i = 1; i < s.length(); i++)
        {
          if(s[i] == '.')
          {
            s[i] = s[i - 1];

            s[i - 1] = '.';

            break;
          }
        }

          from_string = atof(s.c_str() );

          if( fabs( from_string / this->x[0] ) > 3.0 )
          {
            dd_real::error("Re-rounding unsuccessful in large number fixed point trap.") ;
          }
      }
    }


    if( !fixed && !isinf() )
    {
      s += uppercase ? 'E' : 'e';

      append_expn(s, e);
    }
  }

  int len = (int)s.length();

  if(len < width)
  {
    int delta = width - len;

    if(fmt & ios_base::internal)
    {
      if(sgn)
      {
        s.insert(static_cast<string::size_type>(1), delta, fill);
      }
      else
      {
        s.insert(static_cast<string::size_type>(0), delta, fill);
      }
    }
    else if(fmt & ios_base::left)
    {
      s.append(delta, fill);
    }
    else
    {
      s.insert(static_cast<string::size_type>(0), delta, fill);
    }
  }

  return s;
}

int dd_real::read(const char* s, dd_real& a)
{
  const char* p = s;

  char ch = 0;

  int sign = 0;

  int point = -1;

  int nd = 0;

  int e = 0;

  bool done = false;

  dd_real r = 0.0;

  int nread = 0;

  while(*p == ' ')
  {
    p++;
  }

  while( !done && (ch = *p) != '\0')
  {
    if(ch >= '0' && ch <= '9')
    {
      int d = ch - '0';

      r *= 10.0;

      r += static_cast<double>(d);

      nd++;
    }
    else
    {

      switch(ch)
      {

      case '.':
      {
        if(point >= 0)
        {
          return -1;
        }

        point = nd;
      }
      break;

      case '-':
      {
        if(sign != 0 || nd > 0)
        {
          return -1;
        }

        sign = (ch == '-') ? -1 : 1;
      }
      break;

      case '+':
      {
        if(sign != 0 || nd > 0)
        {
          return -1;
        }

        sign = (ch == '-') ? -1 : 1;
      }
      break;

      case 'E':
      {
        nread = std::sscanf(p + 1, "%d", &e);

        done = true;

        if(nread != 1)
        {
          return -1;
        }
      }
      break;

      case 'e':
      {
        nread = std::sscanf(p + 1, "%d", &e);

        done = true;

        if(nread != 1)
        {
          return -1;
        }
      }
      break;

      default:
      {
        return -1;
      }
      break;

      }
    }

    p++;
  }

  if(point >= 0)
  {
    e -= (nd - point);
  }

  if(e != 0)
  {
    r *= (dd_real(10.0) ^ e);
  }

  a = (sign == -1) ? -r : r;

  return 0;
}

void dd_real::dump(const string& name, std::ostream& os) const
{
  std::ios_base::fmtflags old_flags = os.flags();

  std::streamsize old_prec = os.precision(19);

  os << std::scientific;

  if(name.length() > 0)
  {
    os << name << " = ";
  }

  os << "[ " << setw(27) << x[0] << ", " << setw(27) << x[1] << " ]" << endl;

  os.precision(old_prec);

  os.flags(old_flags);
}

void dd_real::dump_bits(const string& name, std::ostream& os) const
{
  string::size_type len = name.length();

  if(len > 0)
  {
    os << name << " = ";

    len +=3;
  }

  os << "[ ";

  len += 2;

  print_double_info(os, x[0] );

  os << endl;

  for(string::size_type i = 0; i < len; i++)
  {
    os << ' ';
  }

  print_double_info(os, x[1] );

  os << " ]" << endl;
}

dd_real dd_real::debug_rand()
{
  if(std::rand() % 2 == 0)
  {
    return ddrand();
  }

  int expn = 0;

  dd_real a = 0.0;

  double d;

  for(int i = 0; i < 2; i++)
  {
    d = std::ldexp(static_cast<double>(std::rand() ) / RAND_MAX, -expn);

    a += d;

    expn = expn + 54 + std::rand() % 200;
  }

  return a;
}

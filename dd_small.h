
// dd_small.h

dd_real dd_real::add(double a, double b)
{
  double s = 0.0;

  double e = 0.0;

  s = qd::two_sum(a, b, e);

  return dd_real(s, e);
}

dd_real operator+(const dd_real& a, double b)
{
  double s1 = 0.0;

  double s2 = 0.0;

  s1 = qd::two_sum(a.x[0], b, s2);

  s2 += a.x[1];

  s1 = qd::quick_two_sum(s1, s2, s2);

  return dd_real(s1, s2);
}

dd_real dd_real::ieee_add(const dd_real& a, const dd_real& b)
{
  double s1 = 0.0;

  double s2 = 0.0;

  double t1 = 0.0;

  double t2 = 0.0;

  s1 = qd::two_sum(a.x[0], b.x[0], s2);

  t1 = qd::two_sum(a.x[1], b.x[1], t2);

  s2 += t1;

  s1 = qd::quick_two_sum(s1, s2, s2);

  s2 += t2;

  s1 = qd::quick_two_sum(s1, s2, s2);

  return dd_real(s1, s2);
}

dd_real dd_real::sloppy_add(const dd_real& a, const dd_real& b)
{
  double s = 0.0;

  double e = 0.0;

  s = qd::two_sum(a.x[0], b.x[0], e);

  e += (a.x[1] + b.x[1] );

  s = qd::quick_two_sum(s, e, e);

  return dd_real(s, e);
}

dd_real operator+(const dd_real& a, const dd_real& b)
{
#if defined(QD_IEEE_ADD)
  return dd_real::ieee_add(a, b);
#else
  return dd_real::sloppy_add(a, b);
#endif
}

dd_real operator+(double a, const dd_real& b)
{
  return (b + a);
}

dd_real& dd_real::operator+=(double a)
{
  double s1 = 0.0;

  double s2 = 0.0;

  s1 = qd::two_sum(x[0], a, s2);

  s2 += x[1];

  x[0] = qd::quick_two_sum(s1, s2, x[1] );

  return *this;
}

dd_real& dd_real::operator+=(const dd_real& a)
{
#if defined(QD_IEEE_ADD)
  double s1 = 0.0;

  double s2 = 0.0;

  double t1 = 0.0;

  double t2 = 0.0;

  s1 = qd::two_sum(x[0], a.x[0], s2);

  t1 = qd::two_sum(x[1], a.x[1], t2);

  s2 += t1;

  s1 = qd::quick_two_sum(s1, s2, s2);

  s2 += t2;

  x[0] = qd::quick_two_sum(s1, s2, x[1] );

  return *this;
#else
  double s = 0.0;

  double e = 0.0;

  s = qd::two_sum(x[0], a.x[0], e);

  e += x[1];

  e += a.x[1];

  x[0] = qd::quick_two_sum(s, e, x[1] );

  return *this;
#endif
}

dd_real dd_real::sub(double a, double b)
{
  double s = 0.0;

  double e = 0.0;

  s = qd::two_diff(a, b, e);

  return dd_real(s, e);
}

dd_real operator-(const dd_real& a, double b)
{
  double s1 = 0.0;

  double s2 = 0.0;

  s1 = qd::two_diff(a.x[0], b, s2);

  s2 += a.x[1];

  s1 = qd::quick_two_sum(s1, s2, s2);

  return dd_real(s1, s2);
}

dd_real operator-(const dd_real& a, const dd_real& b)
{
#if defined(QD_IEEE_ADD)
  double s1 = 0.0;

  double s2 = 0.0;

  double t1 = 0.0;

  double t2 = 0.0;

  s1 = qd::two_diff(a.x[0], b.x[0], s2);

  t1 = qd::two_diff(a.x[1], b.x[1], t2);

  s2 += t1;

  s1 = qd::quick_two_sum(s1, s2, s2);

  s2 += t2;

  s1 = qd::quick_two_sum(s1, s2, s2);

  return dd_real(s1, s2);
#else
  double s = 0.0;

  double e = 0.0;

  s = qd::two_diff(a.x[0], b.x[0], e);

  e += a.x[1];

  e -= b.x[1];

  s = qd::quick_two_sum(s, e, e);

  return dd_real(s, e);
#endif
}

dd_real operator-(double a, const dd_real& b)
{
  double s1 = 0.0;

  double s2 = 0.0;

  s1 = qd::two_diff(a, b.x[0], s2);

  s2 -= b.x[1];

  s1 = qd::quick_two_sum(s1, s2, s2);

  return dd_real(s1, s2);
}

dd_real& dd_real::operator-=(double a)
{
  double s1 = 0.0;

  double s2 = 0.0;

  s1 = qd::two_diff(x[0], a, s2);

  s2 += x[1];

  x[0] = qd::quick_two_sum(s1, s2, x[1] );

  return *this;
}

dd_real& dd_real::operator-=(const dd_real& a)
{
#if defined(QD_IEEE_ADD)
  double s1 = 0.0;

  double s2 = 0.0;

  double t1 = 0.0;

  double t2 = 0.0;

  s1 = qd::two_diff(x[0], a.x[0], s2);

  t1 = qd::two_diff(x[1], a.x[1], t2);

  s2 += t1;

  s1 = qd::quick_two_sum(s1, s2, s2);

  s2 += t2;

  x[0] = qd::quick_two_sum(s1, s2, x[1] );

  return *this;
#else
  double s = 0.0;

  double e = 0.0;

  s = qd::two_diff(x[0], a.x[0], e);

  e += x[1];

  e -= a.x[1];

  x[0] = qd::quick_two_sum(s, e, x[1] );

  return *this;
#endif
}

dd_real dd_real::operator-() const
{
  return dd_real(-x[0], -x[1] );
}

dd_real dd_real::mul(double a, double b)
{
  double p = 0.0;

  double e = 0.0;

  p = qd::two_prod(a, b, e);

  return dd_real(p, e);
}

dd_real ldexp(const dd_real& a, int exp)
{
  return dd_real(std::ldexp(a.x[0], exp), std::ldexp(a.x[1], exp) );
}

dd_real mul_pwr2(const dd_real& a, double b)
{
  return dd_real(a.x[0] * b, a.x[1] * b);
}

dd_real operator*(const dd_real& a, double b)
{
  double p1 = 0.0;

  double p2 = 0.0;

  p1 = qd::two_prod(a.x[0], b, p2);

  p2 += (a.x[1] * b);

  p1 = qd::quick_two_sum(p1, p2, p2);

  return dd_real(p1, p2);
}

dd_real operator*(const dd_real& a, const dd_real& b)
{
  double p1 = 0.0;

  double p2 = 0.0;

  p1 = qd::two_prod(a.x[0], b.x[0], p2);

  p2 += (a.x[0] * b.x[1] + a.x[1] * b.x[0] );

  p1 = qd::quick_two_sum(p1, p2, p2);

  return dd_real(p1, p2);
}

dd_real operator*(double a, const dd_real& b)
{
  return (b * a);
}

dd_real& dd_real::operator*=(double a)
{
  double p1 = 0.0;

  double p2 = 0.0;

  p1 = qd::two_prod(x[0], a, p2);

  p2 += x[1] * a;

  x[0] = qd::quick_two_sum(p1, p2, x[1] );

  return *this;
}

dd_real& dd_real::operator*=(const dd_real& a)
{
  double p1 = 0.0;

  double p2 = 0.0;

  p1 = qd::two_prod(x[0], a.x[0], p2);

  p2 += a.x[1] * x[0];

  p2 += a.x[0] * x[1];

  x[0] = qd::quick_two_sum(p1, p2, x[1] );

  return *this;
}

dd_real dd_real::div(double a, double b)
{
  double q1 = 0.0;

  double q2 = 0.0;

  double p1 = 0.0;

  double p2 = 0.0;

  double s = 0.0;

  double e = 0.0;

  q1 = a / b;

  p1 = qd::two_prod(q1, b, p2);

  s = qd::two_diff(a, p1, e);

  e -= p2;

  q2 = (s + e) / b;

  s = qd::quick_two_sum(q1, q2, e);

  return dd_real(s, e);
}

dd_real operator/(const dd_real& a, double b)
{
  double q1 = 0.0;

  double q2 = 0.0;

  double p1 = 0.0;

  double p2 = 0.0;

  double s = 0.0;

  double e = 0.0;

  dd_real r;

  q1 = a.x[0] / b;

  p1 = qd::two_prod(q1, b, p2);

  s = qd::two_diff(a.x[0], p1, e);

  e += a.x[1];

  e -= p2;

  q2 = (s + e) / b;

  r.x[0] = qd::quick_two_sum(q1, q2, r.x[1] );

  return r;
}

dd_real dd_real::sloppy_div(const dd_real& a, const dd_real& b)
{
  double s1 = 0.0;

  double s2 = 0.0;

  double q1 = 0.0;

  double q2 = 0.0;

  dd_real r;

  q1 = a.x[0] / b.x[0];

  r = b * q1;

  s1 = qd::two_diff(a.x[0], r.x[0], s2);

  s2 -= r.x[1];

  s2 += a.x[1];

  q2 = (s1 + s2) / b.x[0];

  r.x[0] = qd::quick_two_sum(q1, q2, r.x[1] );

  return r;
}

dd_real dd_real::accurate_div(const dd_real& a, const dd_real& b)
{
  double q1 = 0.0;

  double q2 = 0.0;

  double q3 = 0.0;

  dd_real r;

  q1 = a.x[0] / b.x[0];

  r = a - q1 * b;

  q2 = r.x[0] / b.x[0];

  r -= (q2 * b);

  q3 = r.x[0] / b.x[0];

  q1 = qd::quick_two_sum(q1, q2, q2);

  r = dd_real(q1, q2) + q3;

  return r;
}

dd_real operator/(const dd_real& a, const dd_real& b)
{
#if defined(QD_SLOPPY_DIV)
  return dd_real::sloppy_div(a, b);
#else
  return dd_real::accurate_div(a, b);
#endif
}

dd_real operator/(double a, const dd_real& b)
{
  return dd_real(a) / b;
}

dd_real inv(const dd_real& a)
{
  return 1.0 / a;
}

dd_real& dd_real::operator/=(double a)
{
  *this = *this / a;

  return *this;
}

dd_real& dd_real::operator/=(const dd_real& a)
{
  *this = *this / a;

  return *this;
}

dd_real drem(const dd_real& a, const dd_real& b)
{
  dd_real n = nint(a / b);

  return (a - n * b);
}

dd_real divrem(const dd_real& a, const dd_real& b, dd_real& r)
{
  dd_real n = nint(a / b);

  r = a - n * b;

  return n;
}

dd_real sqr(const dd_real& a)
{
  double p1 = 0.0;

  double p2 = 0.0;

  double s1 = 0.0;

  double s2 = 0.0;

  p1 = qd::two_sqr(a.x[0], p2);

  p2 += 2.0 * a.x[0] * a.x[1];

  p2 += a.x[1] * a.x[1];

  s1 = qd::quick_two_sum(p1, p2, s2);

  return dd_real(s1, s2);
}

dd_real dd_real::sqr(double a)
{
  double p1 = 0.0;

  double p2 = 0.0;

  p1 = qd::two_sqr(a, p2);

  return dd_real(p1, p2);
}

dd_real dd_real::operator^(int n)
{
  return npwr(*this, n);
}

dd_real& dd_real::operator=(double a)
{
  x[0] = a;

  x[1] = 0.0;

  return *this;
}

bool operator==(const dd_real& a, double b)
{
  return (a.x[0] == b && a.x[1] == 0.0);
}

bool operator==(const dd_real& a, const dd_real& b)
{
  return (a.x[0] == b.x[0] && a.x[1] == b.x[1] );
}

bool operator==(double a, const dd_real& b)
{
  return (a == b.x[0] && b.x[1] == 0.0);
}

bool operator>(const dd_real& a, double b)
{
  return (a.x[0] > b || (a.x[0] == b && a.x[1] > 0.0) );
}

bool operator>(const dd_real& a, const dd_real& b)
{
  return (a.x[0] > b.x[0] || (a.x[0] == b.x[0] && a.x[1] > b.x[1] ) );
}

bool operator>(double a, const dd_real& b)
{
  return (a > b.x[0] || (a == b.x[0] && b.x[1] < 0.0) );
}

bool operator<(const dd_real& a, double b)
{
  return (a.x[0] < b || (a.x[0] == b && a.x[1] < 0.0) );
}

bool operator<(const dd_real& a, const dd_real& b)
{
  return (a.x[0] < b.x[0] || (a.x[0] == b.x[0] && a.x[1] < b.x[1] ) );
}

bool operator<(double a, const dd_real& b)
{
  return (a < b.x[0] || (a == b.x[0] && b.x[1] > 0.0) );
}

bool operator>=(const dd_real& a, double b)
{
  return (a.x[0] > b || (a.x[0] == b && a.x[1] >= 0.0) );
}

bool operator>=(const dd_real& a, const dd_real& b)
{
  return (a.x[0] > b.x[0] || (a.x[0] == b.x[0] && a.x[1] >= b.x[1] ) );
}

bool operator>=(double a, const dd_real& b)
{
  return (b <= a);
}

bool operator<=(const dd_real& a, double b)
{
  return (a.x[0] < b || (a.x[0] == b && a.x[1] <= 0.0) );
}

bool operator<=(const dd_real& a, const dd_real& b)
{
  return (a.x[0] < b.x[0] || (a.x[0] == b.x[0] && a.x[1] <= b.x[1] ) );
}

bool operator<=(double a, const dd_real& b)
{
  return (b >= a);
}

bool operator!=(const dd_real& a, double b)
{
  return (a.x[0] != b || a.x[1] != 0.0);
}

bool operator!=(const dd_real& a, const dd_real& b)
{
  return (a.x[0] != b.x[0] || a.x[1] != b.x[1] );
}

bool operator!=(double a, const dd_real& b)
{
  return (a != b.x[0] || b.x[1] != 0.0);
}

bool dd_real::is_zero() const
{
  return (x[0] == 0.0);
}

bool dd_real::is_one() const
{
  return (x[0] == 1.0 && x[1] == 0.0);
}

bool dd_real::is_positive() const
{
  return (x[0] > 0.0);
}

bool dd_real::is_negative() const
{
  return (x[0] < 0.0);
}

dd_real abs(const dd_real& a)
{
  return (a.x[0] < 0.0) ? -a : a;
}

dd_real fabs(const dd_real& a)
{
  return abs(a);
}

dd_real nint(const dd_real& a)
{
  double hi = qd::nint(a.x[0] );

  double lo = 0.0;

  if(hi == a.x[0] )
  {
    lo = qd::nint(a.x[1] );

    hi = qd::quick_two_sum(hi, lo, lo);
  }
  else
  {
    lo = 0.0;

    if(std::abs(hi - a.x[0] ) == 0.5 && a.x[1] < 0.0)
    {
      hi -= 1.0;
    }
  }

  return dd_real(hi, lo);
}

dd_real floor(const dd_real& a)
{
  double hi = std::floor(a.x[0] );

  double lo = 0.0;

  if(hi == a.x[0] )
  {
    lo = std::floor(a.x[1] );

    hi = qd::quick_two_sum(hi, lo, lo);
  }

  return dd_real(hi, lo);
}

dd_real ceil(const dd_real& a)
{
  double hi = std::ceil(a.x[0] );

  double lo = 0.0;

  if(hi == a.x[0] )
  {
    lo = std::ceil(a.x[1] );

    hi = qd::quick_two_sum(hi, lo, lo);
  }

  return dd_real(hi, lo);
}

dd_real aint(const dd_real& a)
{
  return (a.x[0] >= 0.0) ? floor(a) : ceil(a);
}

double to_double(const dd_real& a)
{
  return a.x[0];
}

int to_int(const dd_real& a)
{
  return static_cast<int>(a.x[0] );
}

dd_real dd_real::rand()
{
  return ddrand();
}

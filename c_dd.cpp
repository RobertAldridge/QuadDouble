
// c_dd.cpp

#include "include.h"

void c_dd_add(const double* a, const double* b, double* c)
{
  dd_real cc;

  cc = dd_real(a) + dd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_add_dd_d(const double* a, double b, double* c)
{
  dd_real cc;

  cc = dd_real(a) + b;

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_add_d_dd(double a, const double* b, double* c)
{
  dd_real cc;

  cc = a + dd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_sub(const double* a, const double* b, double* c)
{
  dd_real cc;

  cc = dd_real(a) - dd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_sub_dd_d(const double* a, double b, double* c)
{
  dd_real cc;

  cc = dd_real(a) - b;

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_sub_d_dd(double a, const double* b, double* c)
{
  dd_real cc;

  cc = a - dd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_mul(const double* a, const double* b, double* c)
{
  dd_real cc;

  cc = dd_real(a) * dd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_mul_dd_d(const double* a, double b, double* c)
{
  dd_real cc;

  cc = dd_real(a) * b;

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_mul_d_dd(double a, const double* b, double* c)
{
  dd_real cc;

  cc = a * dd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_div(const double* a, const double* b, double* c)
{
  dd_real cc;

  cc = dd_real(a) / dd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_div_dd_d(const double* a, double b, double* c)
{
  dd_real cc;

  cc = dd_real(a) / b;

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_div_d_dd(double a, const double* b, double* c)
{
  dd_real cc;

  cc = a / dd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_copy(const double* a, double* b)
{
  b[0] = a[0];

  b[1] = a[1];
}

void c_dd_copy_d(double a, double* b)
{
  b[0] = a;

  b[1] = 0.0;
}

void c_dd_sqrt(const double* a, double* b)
{
  dd_real bb;

  bb = sqrt(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_sqr(const double* a, double* b)
{
  dd_real bb;

  bb = sqr(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_abs(const double* a, double* b)
{
  dd_real bb;

  bb = abs(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_npwr(const double* a, int n, double* b)
{
  dd_real bb;

  bb = npwr(dd_real(a), n);

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_nroot(const double* a, int n, double* b)
{
  dd_real bb;

  bb = nroot(dd_real(a), n);

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_nint(const double* a, double* b)
{
  dd_real bb;

  bb = nint(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_aint(const double* a, double* b)
{
  dd_real bb;

  bb = aint(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_floor(const double* a, double* b)
{
  dd_real bb;

  bb = floor(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_ceil(const double* a, double* b)
{
  dd_real bb;

  bb = ceil(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_log(const double* a, double* b)
{
  dd_real bb;

  bb = log(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_log10(const double* a, double* b)
{
  dd_real bb;

  bb = log10(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_exp(const double* a, double* b)
{
  dd_real bb;

  bb = exp(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_sin(const double* a, double* b)
{
  dd_real bb;

  bb = sin(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_cos(const double* a, double* b)
{
  dd_real bb;

  bb = cos(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_tan(const double* a, double* b)
{
  dd_real bb;

  bb = tan(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_asin(const double* a, double* b)
{
  dd_real bb;

  bb = asin(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_acos(const double* a, double* b)
{
  dd_real bb;

  bb = acos(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_atan(const double* a, double* b)
{
  dd_real bb;

  bb = atan(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_atan2(const double* a, const double* b, double* c)
{
  dd_real cc;

  cc = atan2(dd_real(a), dd_real(b) );

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_sinh(const double* a, double* b)
{
  dd_real bb;

  bb = sinh(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_cosh(const double* a, double* b)
{
  dd_real bb;

  bb = cosh(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_tanh(const double* a, double* b)
{
  dd_real bb;

  bb = tanh(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_asinh(const double* a, double* b)
{
  dd_real bb;

  bb = asinh(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_acosh(const double* a, double* b)
{
  dd_real bb;

  bb = acosh(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_atanh(const double* a, double* b)
{
  dd_real bb;

  bb = atanh(dd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];
}

void c_dd_sincos(const double* a, double* s, double* c)
{
  dd_real ss;

  dd_real cc;

  sincos(dd_real(a), ss, cc);

  s[0] = ss.x[0];

  s[1] = ss.x[1];

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_sincosh(const double* a, double* s, double* c)
{
  dd_real ss;

  dd_real cc;

  sincosh(dd_real(a), ss, cc);

  s[0] = ss.x[0];

  s[1] = ss.x[1];

  c[0] = cc.x[0];

  c[1] = cc.x[1];
}

void c_dd_read(const char* s, double* a)
{
  dd_real aa(s);

  a[0] = aa.x[0];

  a[1] = aa.x[1];
}

void c_dd_swrite(const double* a, int precision, char* s, int len)
{
  dd_real(a).write(s, len, precision);
}

void c_dd_write(const double* a)
{
  std::cout << dd_real(a).to_string(dd_real::_ndigits) << std::endl;
}

void c_dd_neg(const double* a, double* b)
{
  b[0] = -a[0];

  b[1] = -a[1];
}

void c_dd_rand(double* a)
{
  dd_real aa;

  aa = ddrand();

  a[0] = aa.x[0];

  a[1] = aa.x[1];
}

void c_dd_comp(const double* a, const double* b, int* result)
{
  dd_real aa(a);

  dd_real bb(b);

  if(aa < bb)
  {
    *result = -1;
  }
  else if(aa > bb)
  {
    *result = 1;
  }
  else
  {
    *result = 0;
  }
}

void c_dd_comp_dd_d(const double* a, double b, int* result)
{
  dd_real aa(a);

  dd_real bb(b);

  if(aa < bb)
  {
    *result = -1;
  }
  else if(aa > bb)
  {
    *result = 1;
  }
  else
  {
    *result = 0;
  }
}

void c_dd_comp_d_dd(double a, const double* b, int* result)
{
  dd_real aa(a);

  dd_real bb(b);

  if(aa < bb)
  {
    *result = -1;
  }
  else if(aa > bb)
  {
    *result = 1;
  }
  else
  {
    *result = 0;
  }
}

void c_dd_pi(double* a)
{
  a[0] = dd_real::_pi.x[0];

  a[1] = dd_real::_pi.x[1];
}


// c_qd.cpp

#include "include.h"

void c_qd_add(const double* a, const double* b, double* c)
{
  qd_real cc;

  cc = qd_real(a) + qd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_add_qd_dd(const double* a, const double* b, double* c)
{
  qd_real cc;

  cc = qd_real(a) + dd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_add_dd_qd(const double* a, const double* b, double* c)
{
  qd_real cc;

  cc = dd_real(a) + qd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_add_qd_d(const double* a, double b, double* c)
{
  qd_real cc;

  cc = qd_real(a) + b;

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_add_d_qd(double a, const double* b, double* c)
{
  qd_real cc;

  cc = a + qd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_sub(const double* a, const double* b, double* c)
{
  qd_real cc;

  cc = qd_real(a) - qd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_sub_qd_dd(const double* a, const double* b, double* c)
{
  qd_real cc;

  cc = qd_real(a) - dd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_sub_dd_qd(const double* a, const double* b, double* c)
{
  qd_real cc;

  cc = dd_real(a) - qd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_sub_qd_d(const double* a, double b, double* c)
{
  qd_real cc;

  cc = qd_real(a) - b;

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_sub_d_qd(double a, const double* b, double* c)
{
  qd_real cc;

  cc = a - qd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_mul(const double* a, const double* b, double* c)
{
  qd_real cc;

  cc = qd_real(a) * qd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_mul_qd_dd(const double* a, const double* b, double* c)
{
  qd_real cc;

  cc = qd_real(a) * dd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_mul_dd_qd(const double* a, const double* b, double* c)
{
  qd_real cc;

  cc = dd_real(a) * qd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_mul_qd_d(const double* a, double b, double* c)
{
  qd_real cc;

  cc = qd_real(a) * b;

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_mul_d_qd(double a, const double* b, double* c)
{
  qd_real cc;

  cc = a * qd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_div(const double* a, const double* b, double* c)
{
  qd_real cc;

  cc = qd_real(a) / qd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_div_qd_dd(const double* a, const double* b, double* c)
{
  qd_real cc;

  cc = qd_real(a) / dd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_div_dd_qd(const double* a, const double* b, double* c)
{
  qd_real cc;

  cc = dd_real(a) / qd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_div_qd_d(const double* a, double b, double* c)
{
  qd_real cc;

  cc = qd_real(a) / b;

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_div_d_qd(double a, const double* b, double* c)
{
  qd_real cc;

  cc = a / qd_real(b);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_selfadd(const double* a, double* b)
{
  qd_real bb(b);

  bb += qd_real(a);

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_selfadd_dd(const double* a, double* b)
{
  qd_real bb(b);

  bb += dd_real(a);

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_selfadd_d(double a, double* b)
{
  qd_real bb(b);

  bb += a;

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_selfsub(const double* a, double* b)
{
  qd_real bb(b);

  bb -= qd_real(a);

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_selfsub_dd(const double* a, double* b)
{
  qd_real bb(b);

  bb -= dd_real(a);

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_selfsub_d(double a, double* b)
{
  qd_real bb(b);

  bb -= a;

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_selfmul(const double* a, double* b)
{
  qd_real bb(b);

  bb *= qd_real(a);

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_selfmul_dd(const double* a, double* b)
{
  qd_real bb(b);

  bb *= dd_real(a);

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_selfmul_d(double a, double* b)
{
  qd_real bb(b);

  bb *= a;

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_selfdiv(const double* a, double* b)
{
  qd_real bb(b);

  bb /= qd_real(a);

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_selfdiv_dd(const double* a, double* b)
{
  qd_real bb(b);

  bb /= dd_real(a);

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_selfdiv_d(double a, double* b)
{
  qd_real bb(b);

  bb /= a;

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_copy(const double* a, double* b)
{
  b[0] = a[0];

  b[1] = a[1];

  b[2] = a[2];

  b[3] = a[3];
}

void c_qd_copy_dd(const double* a, double* b)
{
  b[0] = a[0];

  b[1] = a[1];

  b[2] = 0.0;

  b[3] = 0.0;
}

void c_qd_copy_d(double a, double* b)
{
  b[0] = a;

  b[1] = 0.0;

  b[2] = 0.0;

  b[3] = 0.0;
}

void c_qd_sqrt(const double* a, double* b)
{
  qd_real bb;

  bb = sqrt(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_sqr(const double* a, double* b)
{
  qd_real bb;

  bb = sqr(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_abs(const double* a, double* b)
{
  qd_real bb;

  bb = abs(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_npwr(const double* a, int n, double* b)
{
  qd_real bb;

  bb = npwr(qd_real(a), n);

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_nroot(const double* a, int n, double* b)
{
  qd_real bb;

  bb = nroot(qd_real(a), n);

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_nint(const double* a, double* b)
{
  qd_real bb;

  bb = nint(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_aint(const double* a, double* b)
{
  qd_real bb;

  bb = aint(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_floor(const double* a, double* b)
{
  qd_real bb;

  bb = floor(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_ceil(const double* a, double* b)
{
  qd_real bb;

  bb = ceil(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_log(const double* a, double* b)
{
  qd_real bb;

  bb = log(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_log10(const double* a, double* b)
{
  qd_real bb;

  bb = log10(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_exp(const double* a, double* b)
{
  qd_real bb;

  bb = exp(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_sin(const double* a, double* b)
{
  qd_real bb;

  bb = sin(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_cos(const double* a, double* b)
{
  qd_real bb;

  bb = cos(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_tan(const double* a, double* b)
{
  qd_real bb;

  bb = tan(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_asin(const double* a, double* b)
{
  qd_real bb;

  bb = asin(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_acos(const double* a, double* b)
{
  qd_real bb;

  bb = acos(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_atan(const double* a, double* b)
{
  qd_real bb;

  bb = atan(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_atan2(const double* a, const double* b, double* c)
{
  qd_real cc;

  cc = atan2(qd_real(a), qd_real(b) );

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];
}

void c_qd_sinh(const double* a, double* b)
{
  qd_real bb;

  bb = sinh(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_cosh(const double* a, double* b)
{
  qd_real bb;

  bb = cosh(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_tanh(const double* a, double* b)
{
  qd_real bb;

  bb = tanh(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_asinh(const double* a, double* b)
{
  qd_real bb;

  bb = asinh(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_acosh(const double* a, double* b)
{
  qd_real bb;

  bb = acosh(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_atanh(const double* a, double* b)
{
  qd_real bb;

  bb = atanh(qd_real(a) );

  b[0] = bb.x[0];

  b[1] = bb.x[1];

  b[2] = bb.x[2];

  b[3] = bb.x[3];
}

void c_qd_sincos(const double* a, double* s, double* c)
{
  qd_real ss;

  qd_real cc;

  sincos(qd_real(a), ss, cc);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];

  s[0] = ss.x[0];

  s[1] = ss.x[1];

  s[2] = ss.x[2];

  s[3] = ss.x[3];
}

void c_qd_sincosh(const double* a, double* s, double* c)
{
  qd_real ss;

  qd_real cc;

  sincosh(qd_real(a), ss, cc);

  c[0] = cc.x[0];

  c[1] = cc.x[1];

  c[2] = cc.x[2];

  c[3] = cc.x[3];

  s[0] = ss.x[0];

  s[1] = ss.x[1];

  s[2] = ss.x[2];

  s[3] = ss.x[3];
}

void c_qd_read(const char* s, double* a)
{
  qd_real aa(s);

  a[0] = aa.x[0];

  a[1] = aa.x[1];

  a[2] = aa.x[2];

  a[3] = aa.x[3];
}

void c_qd_swrite(const double* a, int precision, char* s, int len)
{
  qd_real(a).write(s, len, precision);
}

void c_qd_write(const double* a)
{
  std::cout << qd_real(a).to_string(qd_real::_ndigits) << std::endl;
}

void c_qd_neg(const double* a, double* b)
{
  b[0] = -a[0];

  b[1] = -a[1];

  b[2] = -a[2];

  b[3] = -a[3];
}

void c_qd_rand(double* a)
{
  qd_real aa;

  aa = qdrand();

  a[0] = aa.x[0];

  a[1] = aa.x[1];

  a[2] = aa.x[2];

  a[3] = aa.x[3];
}

void c_qd_comp(const double* a, const double* b, int* result)
{
  qd_real aa(a);

  qd_real bb(b);

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

void c_qd_comp_qd_d(const double* a, double b, int* result)
{
  qd_real aa(a);

  if(aa < b)
  {
    *result = -1;
  }
  else if(aa > b)
  {
    *result = 1;
  }
  else
  {
    *result = 0;
  }
}

void c_qd_comp_d_qd(double a, const double* b, int* result)
{
  qd_real bb(b);

  if(a < bb)
  {
    *result = -1;
  }
  else if(a > bb)
  {
    *result = 1;
  }
  else
  {
    *result = 0;
  }
}

void c_qd_pi(double* a)
{
  a[0] = qd_real::_pi.x[0];

  a[1] = qd_real::_pi.x[1];

  a[2] = qd_real::_pi.x[2];

  a[3] = qd_real::_pi.x[3];
}

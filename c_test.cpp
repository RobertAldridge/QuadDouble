
// c_test.cpp

#include "include.h"

int test_1()
{
  double a[4], b[4], s[4], p[4], t[4], t2[4];

  double a_new[4], b_new[4], p_old[4];

  double m;

  int r, i;

  const int max_iter = 20;

  puts("Test 1.  (Salamin-Brent quadratically convergent formula for pi)");

  c_qd_copy_d(1.0, a);

  c_qd_copy_d(0.5, t);

  c_qd_sqrt(t, b);

  c_qd_copy_d(0.5, s);

  m = 1.0;

  c_qd_sqr(a, p);

  c_qd_selfmul_d(2.0, p);

  c_qd_selfdiv(s, p);

  printf("  iteration 0: ");

  c_qd_write(p);

  for(i = 1; i <= max_iter; i++)
  {
    m *= 2.0;

    c_qd_add(a, b, a_new);

    c_qd_selfmul_d(0.5, a_new);

    c_qd_mul(a, b, b_new);

    c_qd_sqr(a_new, t);

    c_qd_selfsub(b_new, t);

    c_qd_selfmul_d(m, t);

    c_qd_selfsub(t, s);

    c_qd_copy(a_new, a);

    c_qd_sqrt(b_new, b);

    c_qd_copy(p, p_old);

    c_qd_sqr(a, p);

    c_qd_selfmul_d(2.0, p);

    c_qd_selfdiv(s, p);

    c_qd_sub(p, p_old, t);

    c_qd_abs(t, t2);

    c_qd_comp_qd_d(t2, 1e-60, &r);

    if(r < 0)
    {
      break;
    }

    printf("  iteration %1d: ", i);

    c_qd_write(p);
  }

  c_qd_pi(p);

  printf(" _pi: ");

  c_qd_write(p);

  printf(" error: %.5e = %g eps\r\n", t2[0], t2[0] / ldexp(1.0, -209) );

  return 0;
}

int main1(int argc, const char* argv[] )
{
  int result = 0;

//unsigned int old_cw;

//fpu_fix_start(&old_cw);

  result = test_1();

//fpu_fix_end(&old_cw);

  return result;
}

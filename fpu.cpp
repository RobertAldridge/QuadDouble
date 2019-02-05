
// fpu.cpp

#include "include.h"

#ifdef _WIN32
#endif

#ifdef HAVE_FPU_CONTROL_H
#endif

#ifndef _FPU_GETCW
/*#define _FPU_GETCW(x) asm volatile ("fnstcw %0":"=m" (x) )*/
#endif

#ifndef _FPU_SETCW
/*#define _FPU_SETCW(x) asm volatile ("fldcw %0": :"m" (x) )*/
#endif

#ifndef _FPU_EXTENDED
/*#define _FPU_EXTENDED 0x0300*/
#endif

#ifndef _FPU_DOUBLE
/*#define _FPU_DOUBLE 0x0200*/
#endif

void fpu_fix_start(unsigned int* old_cw)
{
#ifdef X86

#ifdef _WIN32

#ifdef __BORLANDC__
  unsigned short cw = _control87(0, 0);
  _control87(0x0200, 0x0300);
  if(old_cw)
  {
    *old_cw = cw;
  }
#else
  unsigned int cw = _control87(0, 0);
  _control87(0x00010000, 0x00030000);
  if(old_cw)
  {
    *old_cw = cw;
  }
#endif

#else
  volatile unsigned short cw, new_cw;
  _FPU_GETCW(cw);

  new_cw = (cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw);

  if(old_cw)
  {
    *old_cw = cw;
  }
#endif

#endif
}

void fpu_fix_end(unsigned int* old_cw)
{
#ifdef X86

#ifdef _WIN32

#ifdef __BORLANDC__
  if(old_cw)
  {
    unsigned short cw = (unsigned short)*old_cw;
    _control87(cw, 0xFFFF);
  }
#else
  if(old_cw)
  {
    _control87(*old_cw, 0xFFFFFFFF);
  }
#endif

#else
  if(old_cw)
  {
    int cw;
    cw = *old_cw;
    _FPU_SETCW(cw);
  }
#endif

#endif
}

#ifdef HAVE_FORTRAN

#define f_fpu_fix_start FC_FUNC_(f_fpu_fix_start, F_FPU_FIX_START)

#define f_fpu_fix_end   FC_FUNC_(f_fpu_fix_end,   F_FPU_FIX_END)

void f_fpu_fix_start(unsigned int* old_cw)
{
  fpu_fix_start(old_cw);
}

void f_fpu_fix_end(unsigned int* old_cw)
{
  fpu_fix_end(old_cw);
}

#endif

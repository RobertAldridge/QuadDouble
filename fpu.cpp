
// fpu.cpp

#include "include.h"

//_MCW_DN - _DN_SAVE _DN_FLUSH

//_MCW_EM - _EM_INVALID _EM_DENORMAL _EM_ZERODIVIDE _EM_OVERFLOW _EM_UNDERFLOW _EM_INEXACT

//_MCW_IC - _IC_AFFINE _IC_PROJECTIVE

//_MCW_PC - _PC_24 _PC_53 _PC_64

//_MCW_RC - _RC_CHOP _RC_UP _RC_DOWN _RC_NEAR

bool fpu_fix_start(unsigned int* currentControl)
{
  errno_t result = 0;

  currentControl = 0;

  result = _controlfp_s(&currentControl, _DN_SAVE|_IC_PROJECTIVE|_PC_53|_RC_NEAR, _MCW_DN|_MCW_EM|_MCW_IC|_MCW_PC|_MCW_RC);

  return !result;
}

void fpu_fix_end(unsigned int* currentControl)
{
}

#ifdef HAVE_FORTRAN

#define f_fpu_fix_start FC_FUNC_(f_fpu_fix_start, F_FPU_FIX_START)

#define f_fpu_fix_end FC_FUNC_(f_fpu_fix_end, F_FPU_FIX_END)

void f_fpu_fix_start(unsigned int* currentControl)
{
  fpu_fix_start(currentControl);
}

void f_fpu_fix_end(unsigned int* currentControl)
{
  fpu_fix_end(currentControl);
}

#endif

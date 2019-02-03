
// fpu.h

/*
 * Contains functions to set and restore the round-to-double flag in the
 * control word of a x86 FPU.  The algorithms in the double-double and
 * quad-double package does not function with the extended mode found in
 * these FPU.
 */

/*
 * Set the round-to-double flag, and save the old control word in old_cw.
 * If old_cw is NULL, the old control word is not saved.
 */
void fpu_fix_start(unsigned int* old_cw);

/*
 * Restore the control word.
 */
void fpu_fix_end(unsigned int* old_cw);

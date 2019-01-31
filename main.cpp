
// file name: main.cpp

#include <cstdio>
#include <cstdlib>

#include "fpu.h"
#include "config.h"

//#define f_main FC_FUNC_(f_main, F_MAIN)

//extern "C" void f_main();

extern int main1(int argc, const char* argv[] );
extern int main2(int argc, const char* argv[] );
extern int main3(int argc, const char* argv[] );
extern int main4(int argc, const char* argv[] );
extern int main5(int argc, const char* argv[] );
extern int main6(int argc, const char* argv[] );
extern int main7(int argc, const char* argv[] );

int main(int argc, const char* argv[] )
{
  const char* argv1[] = {"c_test"};

  const char* argv2[] = {"compileExample"};

  const char* argv3[] = {"huge", "--all", "--verbose"};

  const char* argv4[] = {"pslq_test", "--all", "--verbose"};

  const char* argv5[] = {"qd_test", "--all", "--verbose"};

  const char* argv6[] = {"qd_timer", "--all", "--verbose"};

  const char* argv7[] = {"quadt_test", "--all", "--verbose"};

  int argc1 = sizeof(argv1) / sizeof(argv1[0] );

  int argc2 = sizeof(argv2) / sizeof(argv2[0] );

  int argc3 = sizeof(argv3) / sizeof(argv3[0] );

  int argc4 = sizeof(argv4) / sizeof(argv4[0] );

  int argc5 = sizeof(argv5) / sizeof(argv5[0] );

  int argc6 = sizeof(argv6) / sizeof(argv6[0] );

  int argc7 = sizeof(argv7) / sizeof(argv7[0] );

//unsigned int old_cw;
//fpu_fix_start(&old_cw);

printf("blah1\n");

//f_main();

main1(argc1, argv1);
main2(argc2, argv2);
main3(argc3, argv3);
main4(argc4, argv4);
main5(argc5, argv5);
main6(argc6, argv6);
main7(argc7, argv7);

printf("blah2\n");

//fpu_fix_end(&old_cw);

  return 0;
}


// include.h

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstring>
#include <limits>
#include <climits>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>

//#include <float.h>

//#include <math.h>

//#include <stdio.h>

//#include <string.h>

//#include <time.h>

//#include <builtins.h>

//#include <fpu_control.h>

#include <ieeefp.h>

//#include <sys/time.h>

//#include <windows.h>

//#if !defined(QD_INLINE)
#define inline
//#endif

namespace qd
{
}

#include "config.h"

#include "qd_config.h"

#include "inline.h"

#include "qd_real.h"
#include "dd_real.h"

#include "fpu.h"

//#if !defined(QD_INLINE)
#include "qd_inline.h"
#include "dd_inline.h"
//#endif

#include "util.h"

#include "bits.h"

#include "tictoc.h"

#include "pslq.h"

#include "c_qd.h"
#include "c_dd.h"

#include "quadt.h"

#ifndef __SSWFINDSIM__
#define __SSWFINDSIM__

#include "../lib/ssw/src/ssw.h"
#include "../lib/ssw/src/kseq.h"

#if (defined(__linux__) || defined(__APPLE__))
#include <getopt.h>
#include <unistd.h>
#else
#include "../lib/getopt9/include/getopt.h"
#include <io.h>
#include <process.h>
#endif

KSEQ_INIT(int, read)

#include "version.h"

#endif

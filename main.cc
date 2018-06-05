#include <cassert>
#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/mman.h>
#include <fcntl.h>
#include "variant.h"
#include "io_plink.h"
#include "methods.h"
#include "doctest.h"

int
main(int argc, char* argv[]) {
	doctest::Context context;
	context.run();
    return 0;
}

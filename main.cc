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

	assert(argc == 2);
	std::string bfile = argv[1];
	PLINKReader pr {bfile};
	double cr = call_rate(pr);
	printf("call rate is %f\n.", cr);
    return 0;
}

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"
#include "io_plink.h"
#include "methods.h"

PLINKReader pr {"data/example1"};

TEST_CASE("call rate") {
	CHECK(call_rate(pr) == 1.0);
}

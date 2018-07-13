#include <vector>
#include <functional>
#include <math.h>
#include "io_plink.h"

double
call_rate(PLINKReader &pr) {
	long total = 0;
	long n_missing = 0;
	for (auto i = 0; i < pr.n_variants; ++i) {
		for (auto j = 0; j < pr.n_samples; ++j) {
			auto gt = pr(i, j);
			if (gt == -1) {
				++n_missing;
			}
			++total;
		}
	}
	if (total != 0)
		return (double) (total - n_missing) / total;
	else {
        printf("divide by zero.\n");
        exit(1);
	}
}

#! /usr/bin/python

import sys
import hail as hl

n_samples = int(sys.argv[1])
n_variants = int(sys.argv[2])
path = sys.argv[3]

mt = hl.balding_nichols_model(1, n_samples, n_variants)
mt = mt.key_cols_by(s = hl.str(mt.sample_idx))
mt = mt.annotate_entries(GT = hl.unphased_diploid_gt_index_call(hl.rand_bool(0.5) * 2))
hl.export_vcf(mt, path + ".vcf")
hl.export_plink(mt, path)

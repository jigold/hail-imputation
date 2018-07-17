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

chimera0 = mt.filter_rows(mt.locus.position < n_variants / 2)
chimera0 = chimera0.filter_cols(chimera0.s == "0")

chimera1 = mt.filter_rows(mt.locus.position >= n_variants / 2)
chimera1 = chimera1.filter_cols(chimera1.s == "1")
chimera1 = chimera1.key_cols_by(s = "0")

mt2 = chimera0.union_rows(chimera1)
hl.export_vcf(mt2, path + "-chimera.vcf")
hl.export_plink(mt2, path + "-chimera")

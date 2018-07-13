#ifndef __IO_PLINK_H__
#define __IO_PLINK_H__

#include <vector>
#include <string>
#include <cassert>
#include <cmath>
#include <sys/mman.h>
#include "variant.h"
#include "site.h"

class PLINKReader {
	public:
		PLINKReader(const std::string &bfile): PLINKReader(bfile + ".bed", bfile + ".bim", bfile + ".fam") {}

		PLINKReader(const std::string &bed_file,
					const std::string &bim_file,
					const std::string &fam_file) {
						n_samples = read_fam(fam_file, samples);
						n_variants = read_bim(bim_file, variants, sites);
						row_sz_ = (n_samples + 3) / 4;
						bed_sz_ = 3 + (n_variants * row_sz_);
						data_sz = n_variants * row_sz_;
						data = read_bed(bed_file, bed_sz_);
					}

		~PLINKReader() { munmap((void *) data, bed_sz_); }

		std::vector<Variant> variants;
		std::vector<std::string> samples;
		std::vector<Site> sites;
		std::size_t n_samples;
		std::size_t n_variants;
		std::size_t data_sz;
		const char *data;
		inline int operator() (std::size_t v_idx, std::size_t s_idx) const;
		inline std::size_t distance(std::size_t i, std::size_t j) const;

	private:
		std::size_t row_sz_;
		std::size_t bed_sz_;

        std::size_t read_bim(const std::string &bim_file, std::vector<Variant> &variants, std::vector<Site> &sites);
        std::size_t read_fam(const std::string &fam_file, std::vector<std::string> &samples);
        char *read_bed(const std::string &bed_file, const std::size_t &sz);
};


inline
int
PLINKReader::operator() (std::size_t v_idx, std::size_t s_idx) const {
	assert(v_idx >= 0 && v_idx < n_variants && s_idx >= 0 && s_idx < n_samples);
	int gt;
	switch(data[(v_idx * row_sz_) + (s_idx >> 2)] >> ((s_idx & 3) << 1) & 3) {
		case 0:
			gt = 2;
			break;
		case 1:
			gt = -1;
			break;
		case 2:
			gt = 1;
			break;
		case 3:
			gt = 0;
			break;
		default:
			gt = -1;
	}
	return gt;
}


inline
std::size_t
PLINKReader::distance(std::size_t i, std::size_t j) const {
	return abs(variants[i].pos() - variants[j].pos());
}

ZippedResult
zip_sites(const PLINKReader &pr1, const PLINKReader &pr2);

#endif
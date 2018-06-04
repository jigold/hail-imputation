#ifndef __IO_PLINK_H__
#define __IO_PLINK_H__

#include <vector>
#include <string>
#include <cassert>
#include <sys/mman.h>
#include "variant.h"


class PLINKReader {
	public:
		PLINKReader(const std::string &bfile): PLINKReader(bfile + ".bed", bfile + ".bim", bfile + ".fam") {}

		PLINKReader(const std::string &bed_file,
					const std::string &bim_file,
					const std::string &fam_file) {
						n_samples = read_fam(fam_file, samples);
						n_variants = read_bim(bim_file, variants);
						row_sz = (n_samples + 3) / 4;
						bed_sz = 3 + (n_variants * row_sz);
						data = read_bed(bed_file, bed_sz);
					}

		~PLINKReader() { munmap((void *) data, bed_sz); }

		std::vector<Variant> variants;
		std::vector<std::string> samples;
		std::size_t n_samples;
		std::size_t n_variants;
		const char *data;
		inline const int operator() (const std::size_t &v_idx, const std::size_t &s_idx);

	private:
		std::size_t row_sz;
		std::size_t bed_sz;
        const std::size_t read_bim(const std::string &bim_file, std::vector<Variant> &variants);
        const std::size_t read_fam(const std::string &fam_file, std::vector<std::string> &samples);
        const char *read_bed(const std::string &bed_file, const std::size_t &sz);
};


inline
const int
PLINKReader::operator() (const std::size_t &v_idx, const std::size_t &s_idx) {
	assert(v_idx >= 0 && v_idx < n_variants && s_idx >= 0 && s_idx < n_samples);
	return data[(v_idx * row_sz) + (s_idx >> 2)] >> ((s_idx & 3) << 1);
}

#endif
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/mman.h>
#include <fcntl.h>
#include "variant.h"
#include "io_plink.h"
#include "site.h"


std::size_t
PLINKReader::read_bim(const std::string &bim_file, std::vector<Variant> &variants, std::vector<Site> &sites) {
	printf("Reading BIM file: '%s'.\n", bim_file.c_str());
    std::ifstream ifs (bim_file, std::ifstream::in);
    std::string line;
	long ct = 0;

    if (ifs.is_open()) {
        while (std::getline(ifs, line)) {
            std::istringstream iss(line);
            std::vector<std::string> v(std::istream_iterator<std::string>{iss},
                                       std::istream_iterator<std::string>());
			variants.push_back(Variant(v[0], std::stoi(v[3]), v[4], v[5], std::stold(v[2])));
			sites.push_back(Site(Variant(v[0], std::stoi(v[3]), v[4], v[5], std::stold(v[2])), ct));
            ct += 1;
        }
        ifs.close();
    } else {
        printf("Unable to open file: '%s'.\n", bim_file.c_str());
		exit(1);
    }

    printf("Read %lu variants.\n", ct);

    std::sort(sites.begin(), sites.end());

    return variants.size();
}

std::size_t
PLINKReader::read_fam(const std::string &fam_file, std::vector<std::string> &samples) {
	printf("Reading FAM file: '%s'.\n", fam_file.c_str());
    std::ifstream ifs (fam_file, std::ifstream::in);
    std::string line;
	long ct = 0;

    if (ifs.is_open()) {
        while (std::getline(ifs, line)) {
            std::istringstream iss(line);
            std::vector<std::string> s(std::istream_iterator<std::string>{iss},
                                       std::istream_iterator<std::string>());
			samples.push_back(s[1]);
            ct += 1;
        }
        ifs.close();
    } else {
        printf("Unable to open file: '%s'.\n", fam_file.c_str());
		exit(1);
    }

    printf("Read %lu samples.\n", ct);
    return samples.size();
}

char *
PLINKReader::read_bed(const std::string &bed_file, const std::size_t &sz) {
	int fd = open(bed_file.c_str(), O_RDONLY);
	if (fd == -1) {
        printf("Could not open file: %s -- %i\n", strerror(errno), errno);
        exit(EXIT_FAILURE);
    }

    void *data = mmap(NULL, sz, PROT_READ, MAP_SHARED, fd, 0);
    if (data  == MAP_FAILED) {
        printf("Could not map memory: %s -- %i\n", strerror(errno), errno);
        exit(EXIT_FAILURE);
    }

    assert(((char *) data)[0] == 108 && ((char *) data)[1] == 27 && ((char *) data)[2] == 1);
    return ((char *) data) + 3;
}

void
check_contig_boundary(std::set<std::size_t> &boundaries, std::string &last_contig, std::string contig, std::size_t idx) {
	if (last_contig == "" || contig != last_contig) {
		last_contig = contig;
		boundaries.insert(idx);
	}
}

ZippedResult
zip_sites(const PLINKReader &pr1, const PLINKReader &pr2) {
	std::size_t n_both = 0;
	std::size_t n_only1 = 0;
	std::size_t n_only2 = 0;
	std::size_t idx = 0;

	std::string last_contig = "";

	auto it1 = pr1.sites.begin();
	auto it2 = pr2.sites.begin();
	std::vector<ZippedSite> zipped_sites;
	std::set<std::size_t> contig_boundaries;

	while (it1 != pr1.sites.end() || it2 != pr2.sites.end()) {
		if (it1 != pr1.sites.end() && it2 != pr2.sites.end()) {
			if (it1->v() == it2->v()) {
				zipped_sites.push_back(ZippedSite {&*it1, &*it2});
                check_contig_boundary(contig_boundaries, last_contig, it1->v().contig(), idx);
				++it1;
				++it2;
				++n_both;
			} else if (it1->v() < it2->v()) {
				zipped_sites.push_back(ZippedSite {&*it1, nullptr});
                check_contig_boundary(contig_boundaries, last_contig, it1->v().contig(), idx);
				++it1;
				++n_only1;
			} else {
				zipped_sites.push_back(ZippedSite {nullptr, &*it2});
                check_contig_boundary(contig_boundaries, last_contig, it2->v().contig(), idx);
				++it2;
				++n_only2;
			}
		} else if (it1 != pr1.sites.end()) {
			zipped_sites.push_back(ZippedSite {&*it1, nullptr});
            check_contig_boundary(contig_boundaries, last_contig, it1->v().contig(), idx);
            ++it1;
            ++n_only1;
		} else {
			zipped_sites.push_back(ZippedSite {nullptr, &*it2});
            check_contig_boundary(contig_boundaries, last_contig, it2->v().contig(), idx);
            ++it2;
            ++n_only2;
		}
		++idx;
	}
	return ZippedResult {zipped_sites, contig_boundaries, n_both, n_only1, n_only2};
}
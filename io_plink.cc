#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/mman.h>
#include <fcntl.h>
#include "variant.h"
#include "io_plink.h"


const std::size_t
PLINKReader::read_bim(const std::string &bim_file, std::vector<Variant> &variants) {
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
            ct += 1;
        }
        ifs.close();
    } else {
        printf("Unable to open file: '%s'.\n", bim_file.c_str());
		exit(1);
    }

    printf("Read %lu variants.\n", ct);
    return variants.size();
}

const std::size_t
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


const char *
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

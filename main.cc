#include <cassert>
#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/mman.h>
#include <fcntl.h>

typedef struct {
  std::string contig;
  int pos;
  std::string ref;
  std::string alt;
} Variant;

typedef std::string Sample;

typedef struct {
  std::vector<Sample> samples;  
  std::vector<Variant> variants;
  std::size_t n_samples;
  std::size_t n_variants;
  std::size_t row_size;
  char* data;
} Dataset;

std::vector<std::string>
split(const std::string &s) {
  std::vector<std::string> result;
  std::istringstream iss(s);
  for (std::string s; iss >> s;) {
    result.push_back(s);
  }
  return result;
}

std::size_t
rowsize(int n_samples) {
  return (n_samples + 3) / 4;
}

std::size_t
bedfile_size(int n_samples, int n_variants) {
  return (std::size_t) 3 + (n_variants * rowsize(n_samples));
}

std::vector<Variant>
parse_bim(std::string file_path) {
  std::ifstream in(file_path);
  assert(in.is_open());
    
  std::vector<Variant> variants;
  std::string line;
  
  while(std::getline(in, line)) {
    std::vector<std::string> tokens = split(line);
    Variant v = {tokens[0], std::stoi(tokens[3]), tokens[4], tokens[5]};
    variants.push_back(v);
  }
  
  in.close();
  return variants;
}

std::vector<Sample>
parse_fam(std::string file_path) {
  std::ifstream in(file_path);
  assert(in.is_open());

  std::vector<Sample> samples;
  std::string line;

  while(std::getline(in, line)) {
    std::vector<std::string> tokens = split(line);
    Sample s = tokens[1];
    samples.push_back(s);
  }

  in.close();
  return samples;  
}

char*
parse_bed(std::string file_path, std::size_t sz) {
  int fd = open(file_path.c_str(), O_RDONLY);
  if (fd == -1) {
    printf("Could not open file: %s -- %i\n", strerror(errno), errno);
    exit(EXIT_FAILURE);
  }

  void *data = mmap(0, sz, PROT_READ, MAP_SHARED, fd, 0);
  if (data  == MAP_FAILED) {
    printf("Could not map memory: %s -- %i\n", strerror(errno), errno);
    exit(EXIT_FAILURE);
  }

  assert(((char *) data)[0] == 108 && ((char *) data)[1] == 27 && ((char *) data)[2] == 1);
  
  return ((char *) data) + 3;
}

double
call_rate(Dataset *d) {
  std::size_t n_genotypes = 0;
  std::size_t n_called = 0;
  std::size_t row_size = d->row_size;
  
  for (int vi = 0; vi < d->n_variants; vi++) {
    for (int si = 0; si < d->n_samples; si++) {
      int gt = (d->data[(vi * row_size) + (si >> 2)] >> ((si & 3) << 1)) & 3;
      if (gt != 1) {
	++n_called;
      }
      ++n_genotypes;
    }
  }

  return (double) n_called / n_genotypes;
}

int
main(int argc, char* argv[]) {
  assert(argc == 4);
  std::string bed_file = argv[1];
  std::string bim_file = argv[2];
  std::string fam_file = argv[3];

  std::vector<Sample> samples = parse_fam(fam_file);  
  std::vector<Variant> variants = parse_bim(bim_file);
  
  int n_samples = samples.size();  
  int n_variants = variants.size();
  
  char *data = parse_bed(bed_file, bedfile_size(n_samples, n_variants));
  
  Dataset d = {samples, variants, (std::size_t) n_samples, (std::size_t) n_variants, rowsize(n_samples), data};
  std::cout << "n_variants: " << d.n_variants << "\n";
  std::cout << "n_samples: " << d.n_samples << "\n";

  double cr = call_rate(&d);
  std::cout << "call rate: " << cr << "\n";
  
  return 0;
}

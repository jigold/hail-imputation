#ifndef __VARIANT_H__
#define __VARIANT_H__

#include <string>

class Variant {
    public:
        Variant(
        const std::string &contig_,
        const long &pos_,
        const std::string &ref_,
        const std::string &alt_,
        const double &pos_cm_) : contig(contig_), pos(pos_), ref(ref_), alt(alt_), pos_cm(pos_cm_) {}

        const std::string contig;
        const long pos;
        const std::string ref;
        const std::string alt;
        const double pos_cm;
        const std::string to_string();
};

inline
const std::string Variant::to_string() {
    return contig + ":" + std::to_string(pos) + ":" + ref + ":" + alt;
}

#endif

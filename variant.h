#ifndef __VARIANT_H__
#define __VARIANT_H__

#include <string>
#include <assert.h>

class Variant {
    public:
        Variant(
        const std::string &contig,
        const int &pos,
        const std::string &ref,
        const std::string &alt,
        const double &pos_cm) : contig_(contig), pos_(pos), ref_(ref), alt_(alt), pos_cm_(pos_cm) {
            assert(contig_ != "");
        }

        std::string contig() const { return contig_; }
        int pos() const { return pos_; }
        std::string ref() const { return ref_; }
        std::string alt() const { return alt_; }
        double pos_cm() const { return pos_cm_; }

        std::string to_string() const;

        bool operator==(const Variant &other) const;
		bool operator!=(const Variant &other) const;
		bool operator<(const Variant &other) const;
		bool operator>(const Variant &other) const;

	private:
		std::string contig_;
		int pos_;
		std::string ref_;
		std::string alt_;
		double pos_cm_;
};

inline
std::string Variant::to_string() const {
    return contig_ + ":" + std::to_string(pos_) + ":" + ref_ + ":" + alt_;
}

inline
bool Variant::operator==(const Variant& other) const {
    if (typeid(*this) != typeid(other))
        return false;

    return contig_ == other.contig() &&
        pos_ == other.pos() &&
        ref_ == other.ref() &&
        alt_ == other.alt() &&
        pos_cm_ == other.pos_cm();
}

inline
bool Variant::operator!=(const Variant& other) const {
    return !(*this == other);
}

inline
bool Variant::operator<(const Variant& other) const {
    if (contig_ != other.contig()) {
        return (contig_ < other.contig());
    } else if (pos_ != other.pos()) {
        return (pos_ < other.pos());
    } else if (ref_ != other.ref()) {
        return (ref_ < other.ref());
    } else {
        return (alt_ < other.alt());
    };
}

inline
bool Variant::operator>(const Variant& other) const {
	return !(*this < other);
}

#endif

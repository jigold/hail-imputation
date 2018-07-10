#ifndef __SITE_H__
#define __SITE_H__

#include <vector>
#include "variant.h"

class Site {
	public:
		Site(const Variant &v, const std::size_t &idx): v_(v), idx_(idx) {}

		Variant v() const { return v_; }
		std::size_t idx() const { return idx_; }

        bool operator==(const Site &other) const;
		bool operator!=(const Site &other) const;
		bool operator<(const Site &other) const;

		std::string to_string() const;

	private:
		Variant v_;
		std::size_t idx_;
};

inline
bool Site::operator==(const Site& other) const {
    if (typeid(*this) != typeid(other))
        return false;

    return v_ == other.v() && idx_ == other.idx();
}

inline
bool Site::operator!=(const Site& other) const {
    return !(*this == other);
}


inline
bool Site::operator<(const Site& other) const {
    if (v_ != other.v()) {
        return v_ < other.v();
    } else {
        return idx_ < other.idx();
    }
}

inline
std::string Site::to_string() const {
    return v_.to_string() + " " + std::to_string(idx_);
}

struct ZippedSite {
	ZippedSite(const Site* s1_, const Site* s2_): s1(s1_), s2(s2_) {}

	const Site* s1;
	const Site* s2;
};

#endif
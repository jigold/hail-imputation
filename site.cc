#include <vector>

std::vector<ZippedSite>
zip_sites(const PLINKReader &pr1, const PLINKReader &pr2) {
	auto it1 = pr1.sites.begin();
	auto it2 = pr2.sites.begin();
	std::vector<ZippedSite> zipped_sites;
	while (it1 != pr1.sites.end() || it2 != pr2.sites.end()) {
		if (it1 != pr1.sites.end() && it2 != pr2.sites.end()) {
			if (it1->v() == it2->v()) {
				zipped_sites.push_back(ZippedSite {&*it1, &*it2});
				++it1;
				++it2;
			} else if (it1->v() < it2->v()) {
				zipped_sites.push_back(ZippedSite {&*it1, nullptr});
				++it1;
			} else {
				zipped_sites.push_back(ZippedSite {nullptr, &*it2});
				++it2;
			}
		} else if (it1 != pr1.sites.end()) {
			zipped_sites.push_back(ZippedSite {&*it1, nullptr});
            ++it1;
		} else {
			zipped_sites.push_back(ZippedSite {nullptr, &*it2});
            ++it2;
		}
	}
	return zipped_sites;
}
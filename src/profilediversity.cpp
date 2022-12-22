#include "profilediversity.hpp"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "diversity.hpp"
#include "tools.hpp"

using namespace std;

ProfileDiversity::ProfileDiversity(string finput) : input(finput) {}

void ProfileDiversity::parseCsv() {
    ifstream stream(input);
    string line("");
    vector<string> vline;
    int loop(0);
    string identifiant;
    shared_ptr<Diversity> diversity = nullptr;

    while (getline(stream, line)) {
        if (loop == 0) {  // Version line
            vline = Tools::split(line, " ");

            if (vline.size() != 2 || vline[0].rfind("##", 0) != 0) {
                throw logic_error("Profile Diversity version malformated");
            }
            version = Tools::stringToFloat(vline[1]);
        } else if (loop == 1) {  // Header line
            vline = Tools::split(line, ",");
            if (vline.size() != 5 || vline[0].rfind("#", 0) != 0) {
                throw logic_error("Profile Diversity header malformated");
            }
        } else {
            vline = Tools::split(line, ",");

            identifiant = vline[0];
            diversity = make_shared<Diversity>();

            diversity->setMutRate(Tools::stringToDouble(vline[1]));
            diversity->setIndelFrac(Tools::stringToDouble(vline[2]));
            diversity->setIndelExtend(Tools::stringToDouble(vline[3]));
            diversity->setMaxInsSize(Tools::stringToInt(vline[4]));

            if (data.count(identifiant) > 0) {
                throw logic_error("Identifiant find " + identifiant +
                                  " isn't unique");
            }
            data.emplace(identifiant, diversity);
            diversity.reset();
        }
        loop++;
    }
}

int ProfileDiversity::count(string n) const noexcept { return data.count(n); }

shared_ptr<Diversity> ProfileDiversity::getDiversity(string id) noexcept {
    return data[id];
}


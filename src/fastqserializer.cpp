#include "fastqserializer.hpp"

#include <cstring>
#include <iostream>
#include <stdexcept>

#include "bgzf.h"
#include "fastq.hpp"

using namespace std;

FastqSerializer::FastqSerializer(string fid)
    : finput(fid), fp(nullptr), compression_level(1) {}
FastqSerializer::~FastqSerializer() {}
FastqSerializer& FastqSerializer::operator=(const FastqSerializer& other) {
    finput = other.finput;
    *fp = *(other.fp);
    compression_level = other.compression_level;
    return *this;
}

void FastqSerializer::open() {
    char mode[4] = "w";
    size_t len = finput.length();

    mode[2] = 0;
    mode[3] = 0;
    if (len > 3 && strstr(finput.c_str() + (len - 3), ".gz")) {
        mode[1] = 'g';
        mode[2] = compression_level + '0';
    } else if ((len > 4 && strstr(finput.c_str() + (len - 4), ".bgz")) ||
               (len > 5 && strstr(finput.c_str() + (len - 5), ".bgzf"))) {
        mode[1] = compression_level + '0';
    } else {
        mode[1] = 'u';
    }

    fp = bgzf_open(finput.c_str(), mode);
}
void FastqSerializer::close() { bgzf_close(fp); }
void FastqSerializer::writeRecord(const Fastq& record) {
    string str = record.toString();
    const char* ctr = str.c_str();
    if (bgzf_write(fp, ctr, str.length()) < 0) {
        throw logic_error("Unable to write record to fastq file");
    }
}

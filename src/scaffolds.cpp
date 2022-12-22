#include "scaffolds.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "index.hpp"
#include "scaffold.hpp"
#include "tools.hpp"

using namespace std;

Scaffolds::Scaffolds(string finput, string foutput)
    : input(finput), output(foutput) {}
Scaffolds::Scaffolds(string finput)
    : input(finput), output(finput + ".scaff") {}
Scaffolds::Scaffolds() {}
Scaffolds::~Scaffolds() { bgzf_close(bgzf); }

std::vector<Scaffold> Scaffolds::getVector() { return scaffolds; }

int Scaffolds::insertIndex(const string name, long long int start,
                           long long int end) {
    if (name.length() < 1) {
        cout << "Malformed line" << endl;
        return -1;
    }

    Scaffold scaffold(name, start, end);

    scaffolds.push_back(scaffold);

    return 0;
}

int Scaffolds::init(int min_len_to_keep) {
    int res;
    bgzf = bgzf_open(input.c_str(), "r");
    if (!Tools::isFileExist(output)) {
        if (buildCore() != 0) {
            throw logic_error("Error when read file");
        }
        res = save();
    } else {
        res = load(min_len_to_keep);
    }
    initMembers();
    return res;
}

int Scaffolds::buildCore() {
    stringstream bufname;
    int c(0), read_done(0), line_num(1), loop_name(0);
    // Faidx_t *idx;
    int seq_len(0), qual_len(0);
    int char_len(0), cl(0), line_len(0), ll(0);
    long long int start(1), stop(1);  // samtools fetch 1-based
    bool inN = true;

    enum read_state { OUT_READ, IN_NAME, IN_SEQ, SEQ_END, IN_QUAL } state;

    fai_format_options format = FAI_NONE;

    state = OUT_READ;

    while ((c = bgzf_getc(bgzf)) >= 0) {
        try {
            switch (state) {
                case OUT_READ:
                    switch (c) {
                        case '>':
                            if (format == FAI_FASTQ) {
                                cout << "Found '>' in a FASTQ file, error at "
                                        "line "
                                     << line_num << endl;
                                throw logic_error("");
                            }

                            format = FAI_FASTA;
                            state = IN_NAME;
                            break;

                        case '@':
                            if (format == FAI_FASTA) {
                                cout << "Found '@' in a FASTA file, error at "
                                        "line "
                                     << line_num << endl;
                                throw logic_error("");
                            }

                            format = FAI_FASTQ;
                            state = IN_NAME;
                            break;

                        case '\r':
                            // Blank line with cr-lf ending?
                            if ((c = bgzf_getc(bgzf)) == '\n') {
                                line_num++;
                            } else {
                                cout << "Format error, carriage return not "
                                        "followed by new line at line "
                                     << line_num << endl;
                                throw logic_error("");
                            }
                            break;

                        case '\n':
                            // just move onto the next line
                            line_num++;
                            break;

                        default: {
                            cout << "Format error, unexpected \"" << c
                                 << "\"\0 at line " << line_num << endl;
                            throw logic_error("");
                            break;
                        }
                    }
                    break;

                case IN_NAME:
                    if (read_done) {
                        if (!inN &&
                            insertIndex(bufname.str(), start, stop - 1) != 0)
                            throw logic_error("Create region index failded");
                        bufname.str(string());
                        read_done = 0;
                    }
                    loop_name = 0;
                    do {
                        if (!isspace(c)) {
                            bufname << (char)c;
                            loop_name++;
                        } else if (loop_name > 0 || c == '\n') {
                            break;
                        }
                    } while ((c = bgzf_getc(bgzf)) >= 0);

                    if (c < 0) {
                        cout << "The last entry has no sequence "
                             << bufname.str() << endl;
                        throw logic_error("");
                    }

                    // read the rest of the line if necessary
                    if (c != '\n')
                        while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n')
                            ;

                    state = IN_SEQ;
                    seq_len = qual_len = char_len = line_len = 0;
                    start = stop = 1;
                    line_num++;
                    break;

                case IN_SEQ:
                    if (format == FAI_FASTA) {
                        if (c == '\n') {
                            state = OUT_READ;
                            line_num++;
                            continue;
                        } else if (c == '>') {
                            state = IN_NAME;
                            continue;
                        }
                    } else if (format == FAI_FASTQ) {
                        if (c == '+') {
                            state = IN_QUAL;
                            if (c != '\n')
                                while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n')
                                    ;
                            line_num++;
                            continue;
                        } else if (c == '\n') {
                            cout << "Inlined empty line is not allowed in "
                                    "sequence "
                                 << bufname.str() << " at line " << line_num
                                 << endl;
                            throw logic_error("");
                        }
                    }

                    ll = cl = 0;

                    if (format == FAI_FASTA) read_done = 1;

                    do {
                        ll++;
                        cl++;
                        if (isgraph(c)) {
                            if (c == 'N' || c == 'n') {
                                if (!inN) {
                                    if (insertIndex(bufname.str(), start,
                                                    stop - 1) != 0) {
                                        throw logic_error(
                                            "Create region index failded");
                                    }
                                    start = stop;
                                    stop++;

                                } else {
                                    start++;
                                    stop++;
                                }
                                inN = true;
                            } else {
                                if (inN) {
                                    start = stop;
                                    stop++;
                                    inN = false;
                                } else {
                                    stop++;
                                }
                            }
                        }
                    } while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n');

                    ll++;
                    seq_len += cl;

                    if (line_len == 0) {  // first loop
                        line_len = ll;
                        char_len = cl;
                    } else if (line_len > ll) {  // end sequence

                        if (format == FAI_FASTA) {
                            state = OUT_READ;
                        } else {
                            state = SEQ_END;
                        }
                    } else if (line_len < ll) {  // line too long int
                        cout << "Different line length in sequence "
                             << bufname.str() << endl;
                        throw logic_error("");
                    }

                    line_num++;
                    break;

                case SEQ_END:
                    if (c == '+') {
                        state = IN_QUAL;
                        if (c != '\n')
                            while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n')
                                ;
                        line_num++;
                        continue;
                    } else {
                        cout << "Format error, expecting '+', got " << c
                             << " at line " << line_num << endl;
                        throw logic_error("");
                    }
                    break;

                case IN_QUAL:
                    if (c == '\n') {
                        if (!read_done) {
                            cout << "Inlined empty line is not allowed in "
                                    "quality of sequence "
                                 << bufname.str() << endl;
                            throw logic_error("");
                        }

                        state = OUT_READ;
                        line_num++;
                        continue;
                    } else if (c == '@' && read_done) {
                        state = IN_NAME;
                        continue;
                    }

                    ll = cl = 0;

                    do {
                        ll++;
                        if (isgraph(c)) cl++;
                    } while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n');

                    ll++;
                    qual_len += cl;

                    if (line_len < ll) {
                        cout << "Quality line length too long int in "
                             << bufname.str() << " at line " << line_num
                             << endl;
                        throw logic_error("");
                    } else if (qual_len == seq_len) {
                        read_done = 1;
                    } else if (qual_len > seq_len) {
                        cout << "Quality length long inter than sequence in "
                             << bufname.str() << " at line " << line_num
                             << endl;
                        throw logic_error("");
                    } else if (line_len > ll) {
                        cout << "Quality line length too short in "
                             << bufname.str() << " at line " << line_num
                             << endl;
                        throw logic_error("");
                    }

                    line_num++;
                    break;
            }

        }

        catch (const std::exception& e) {
            std::cerr << e.what();
        }
    }

    if (read_done) {
        if (!inN && insertIndex(bufname.str(), start, stop) != 0) {
            throw logic_error("Read done, error throw");
        }
    } else {
        throw logic_error("");
    }

    return 0;
}

string Scaffolds::toString() {
    stringstream buf;
    for (Scaffold f : scaffolds) {
        buf << f.getName() << '\t' << f.getStart() << '\t' << f.getStop()
            << endl;
    }
    return buf.str();
}

int Scaffolds::save() {
    ofstream flux(output.c_str());
    if (!flux) {
        throw logic_error("Output isn't writable");
    }
    flux << toString();
    return 0;
}

int Scaffolds::load(int gt) {
    string info;
    vector<string> tokens;

    // Read file
    ifstream flux(output.c_str());
    if (!flux) {
        throw logic_error("Error when read index file");
    }
    while (getline(flux, info)) {
        tokens = Tools::split(info.c_str(), "\t");
        if (tokens.size() != 3) {
            throw logic_error("Index file is malformated ?");
        }
        if (stol(tokens[2]) - stol(tokens[1]) >= gt) {
            insertIndex(tokens[0], stol(tokens[1]), stol(tokens[2]));
        }
    }

    return 0;
}

void Scaffolds::initMembers() {
    long long int length = 0, length_sum = 0;
    for (auto& f : scaffolds) {
        length = f.getStop() - f.getStart();
        length_sum += length;
        length_scaffolds += length;
        intervals.push_back(length_sum);
    }
}
long long int Scaffolds::getTotalLength() const noexcept {
    return length_scaffolds;
}

vector<long long int> Scaffolds::getIntervals() const noexcept {
    return intervals;
}


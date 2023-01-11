// Copyright 2022 guillaume-gricourt
#include "sequence.hpp"

#include <ctype.h>

#include <exception>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "diversity.hpp"
#include "randomgenerator.hpp"
#include "tools.hpp"

/*****************
 * Static member *
 * **************/
const char comp_tab[256] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'T', 'V', 'G', 'H', 'N', 'N', 'C', 'D', 'N', 'N',
    'M', 'N', 'K', 'N', 'N', 'N', 'N', 'Y', 'S', 'A', 'B', 'W', 'N', 'N', 'R',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T', 'v', 'G', 'v', 'N', 'N', 'C', 'd',
    'N', 'N', 'm', 'N', 'k', 'n', 'N', 'N', 'N', 'y', 's', 'A', 'N', 'b', 'w',
    'N', 'r', 'N', 'N', 'N', 'N', 'N', 'N',

    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'

    /* Code reverse complement
    ACGTRYKMSWBDHVN
    |||||||||||||||
    TGCAYRMKSWVHDBN
    */
};

const std::string nucleotides = "ATCG";

Sequence::Sequence(std::string &seq) : sequence(seq) {}

std::string Sequence::toString() const noexcept { return sequence; }
size_t Sequence::getLength() const noexcept { return sequence.size(); }
std::string Sequence::getSequence(int length, bool start, bool reverse_compl) {
  std::string sub("");
  int ix_start = 0;
  if (start) {
    sub = sequence.substr(0, length);
  } else {
    if (length < static_cast<int>(sequence.length())) {
      ix_start = sequence.length() - length;
    }
    sub = sequence.substr(ix_start, std::string::npos);
  }
  if (reverse_compl) {
    reverseComplement(&sub);
  }
  return sub;
}

void Sequence::chooseBase(RandomGenerator *random_generator, char *nucl,
                          char *remove) {
  std::string combinaison = "";
  if (remove) {
    if (islower(*remove)) {
      combinaison = Tools::toLower(&nucleotides);
    } else {
      combinaison = nucleotides;
    }
    size_t pos = combinaison.find_first_of(*remove);
    // Check if alphabet degenerated
    if (pos != std::string::npos) {
      combinaison.erase(combinaison.begin() + pos);
    }
  } else {
    combinaison = nucleotides;
  }

  *nucl = combinaison[random_generator->randomRange(
      0, static_cast<int>(combinaison.size()) - 1)];
}

void Sequence::makeMutation(RandomGenerator *random_generator,
                            std::shared_ptr<Diversity> diversity) {
  // Create data struct
  auto len_seq = sequence.size();
  auto ilen_seq = static_cast<int>(len_seq);
  std::vector<MutationType> mutations(len_seq, MutationType::NONE);
  std::vector<std::string> insertions;
  bool deleting = false;
  char nucl;
  int len_ins(0); //, total_len_ins(0), total_len_del(0);
  std::stringstream buf;

  // cout << "create vector" << endl;
  for (int i(0); i < ilen_seq; ++i) {
    if (deleting) {
      if (random_generator->random48() < diversity->getIndelExtend()) {
        mutations[i] = MutationType::DELETE;
        // total_len_del++;
        continue;
      } else {
        deleting = false;
      }
    }
    if (random_generator->random48() < diversity->getMutRate()) { // mutation
      if (random_generator->random48() >=
          diversity->getIndelFrac()) { // substitution
        mutations[i] = MutationType::SUBSTITUTE;
      } else {
        if (random_generator->random48() < 0.5 && i > 1) { // deletion
          mutations[i] = MutationType::DELETE;
          deleting = true;
          // total_len_del++;
        } else { // insertion
          len_ins = 1;
          auto p = random_generator->random48();
          while (len_ins <= diversity->getMaxInsSize() &&
                 p < diversity->getIndelExtend()) {
            len_ins++;
          }
          // total_len_ins += len_ins;
          mutations[i] = MutationType::INSERT;

          for (int j(0); j < len_ins; ++j) {
            Sequence::chooseBase(random_generator, &nucl);
            buf << (char)nucl;
          }
          insertions.push_back(buf.str());
          buf.str(std::string());
        }
      }
    }
  }
  // cout << "apply mutations" << endl;
  // Apply Mutations
  int pos_seq(0);
  //+ total_len_ins - total_len_del
  for (int i(0); i < ilen_seq; ++i) {
    if (mutations[i] == MutationType::INSERT) {
      sequence.insert(pos_seq, insertions[0]);
      pos_seq += insertions[0].length() - 1;
      insertions.erase(insertions.begin());
    }

    if (mutations[i] == MutationType::SUBSTITUTE) {
      Sequence::chooseBase(random_generator, &nucl, &sequence[pos_seq]);
      sequence[pos_seq] = nucl;
    }

    if (mutations[i] == MutationType::DELETE) {
      sequence.erase(sequence.begin() + pos_seq);
      --pos_seq;
    }

    pos_seq++;
  }
}

void Sequence::reverseComplement(std::string *str) noexcept {
  char c;
  int i = 0, j = str->length() - 1;

  while (i <= j) {
    c = (*str)[i];
    (*str)[i] = comp_tab[(unsigned char)(*str)[j]];
    (*str)[j] = comp_tab[(unsigned char)c];
    i++;
    j--;
  }
}
void Sequence::reverse(std::string *str) noexcept {
  char c;
  int i = 0, j = str->length() - 1;

  while (i < j) {
    c = (*str)[i];
    (*str)[i] = (*str)[j];
    (*str)[j] = c;
    i++;
    j--;
  }
}

char &Sequence::operator[](int i) { return sequence[i]; }
char Sequence::operator[](int i) const { return sequence[i]; }

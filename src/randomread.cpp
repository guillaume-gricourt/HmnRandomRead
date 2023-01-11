// Copyright 2022 guillaume-gricourt
#include "randomread.hpp"

#include <algorithm>
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

#include "args.hpp"
#include "faidxh.hpp"
#include "faidxr.hpp"
#include "randomgenerator.hpp"
#include "reference.hpp"
#include "scaffold.hpp"
#include "scaffolds.hpp"
#include "sequence.hpp"

// Helpers
RandomRead::RandomRead(Args &arg)
    : args(arg), random_generator(arg.getSeed()),
      read_ser_1(std::move(args.getRead1())),
      read_ser_2(std::move(args.getRead2())) {}

int RandomRead::init() {
  // Make index
  std::cout << "Make index for references" << std::endl;
  makeIndex();
  // Profile Diversity
  std::cout << "Parse - profil diversity" << std::endl;
  if (args.isProfileDiversity()) {
    args.getProfileDiversity()->parseCsv();
    for (auto &ref : *args.getReferences()) {
      if (args.getProfileDiversity()->count(ref->getIdDiversity()) > 1) {
        throw std::logic_error("Id profile diversity " + ref->getIdDiversity() +
                               " for reference " + ref->getFilePath() +
                               " does not match");
      }
    }
  }
  // Profile Error
  std::cout << "Parse - profil error" << std::endl;
  if (args.isProfileError()) {
    args.getProfileError()->parseCsv(args.getProfileErrorId(),
                                     args.isOutputPaired());
  }
  // Make reads
  std::cout << "Make reads" << std::endl;
  makeReads();
  return 0;
}
int RandomRead::makeIndex() {
  int min_length_scaffold = static_cast<int>(args.getLengthReads() * 2 / 3);

  for (auto &ref : *args.getReferences()) {
    // Fai
    ref->getFaidxr().init();
    // Scaff
    ref->getScaffolds().init(min_length_scaffold);
  }
  return 0;
}

int RandomRead::makeReads() {
  // Open files
  read_ser_1.open();
  read_ser_2.open();

  int64_t nbReads(0);
  int nb_ref(0);
  int ix_org(0);
  // Calculate total nb reads
  for (auto &ref : *args.getReferences()) {
    nbReads += ref->getNbReads();
  }

  // Run
  int count(0);
  std::vector<Scaffold> scaffolds;
  // vector<Scaffold> scaffolds_keep;
  std::vector<int64_t> intervals;
  int64_t length(0);
  int64_t ix_scaff(0);

  int64_t locationStart(0), locationStop(0);
  int64_t sizeInsert(0);
  std::string seq = "";
  bool scaffold_find = false;

  for (int64_t i = 0; i < nbReads; ++i) {
    Scaffold scaffold;
    // Find Org - Get ref name unique
    nb_ref = static_cast<int>(args.getReferences()->size());
    ix_org = 0;
    if (nb_ref > 1) {
      ix_org = random_generator.randomRange(0, nb_ref - 1);
    }
    auto &choiceRef = args.getReferences()->at(ix_org);
    while (count < 10) {
      // Find Scaffold
      // From:
      // https://softwareengineering.stackexchange.com/questions/150616/get-weighted-random-item
      scaffolds = choiceRef->getScaffolds().getVector();

      length = choiceRef->getScaffolds().getTotalLength();
      intervals = choiceRef->getScaffolds().getIntervals();

      ix_scaff = random_generator.randomRangeLong(0L, length);

      for (size_t j = 0; j < intervals.size(); j++) {
        if (intervals.at(j) >= ix_scaff) {
          scaffold = scaffolds.at(j);
          scaffold_find = true;
          break;
        }
      }
      if (!scaffold_find) {
        count++;
        continue;
      }

      // Defined Insert length/Loc
      sizeInsert = std::round(random_generator.randomNormal(
          (double)args.getMeanInsertSize(), (double)args.getStdInsertSize()));
      // TODO(guillaume-gricourt): Loop on choice
      // Check if size in scaffolds

      // Defined loc
      if (sizeInsert > scaffold.getLength()) {
        locationStart = scaffold.getStart();
        locationStop = scaffold.getStop();
      } else {
        // normalize
        locationStart = random_generator.randomRangeLong(
            scaffold.getStart(),
            scaffold.getStop() - static_cast<int64_t>(sizeInsert));
        locationStop = locationStart + sizeInsert;
      }

      // Retrieve Insert-Sequence
      std::stringstream bufRoi;
      bufRoi << scaffold.getName() << ":" << locationStart << "-"
             << locationStop;
      seq = choiceRef->getFaidxr().fetch(bufRoi.str());

      if (seq == "") {
        count++;
      } else {
        break;
      }

      intervals.clear();
      length = 0;
    }
    // Add mutations : deletion/insertion/sub/N
    Sequence sequence(seq);
    if (args.getProfileDiversity() && choiceRef->getIdDiversity() != "") {
      sequence.makeMutation(&random_generator,
                            args.getProfileDiversity()->getDiversity(
                                choiceRef->getIdDiversity()));
    }
    // Strand : F or R
    std::string read1_seq(""), read2_seq("");

    double strand = random_generator.randomRange(0.0, 1.0);

    if (strand >= 0.5) {
      read1_seq = sequence.getSequence(args.getLengthReads(), true, false);
      read2_seq = sequence.getSequence(args.getLengthReads(), false, true);
    } else {
      read1_seq = sequence.getSequence(args.getLengthReads(), false, true);
      read2_seq = sequence.getSequence(args.getLengthReads(), true, false);
    }
    // Create record
    Fastq read1(read1_seq, 33, i, (strand >= 0.5) ? true : false,
                choiceRef->getFileName(), scaffold.getName(), locationStart,
                locationStart + static_cast<int>(read1_seq.length()));
    Fastq read2(read2_seq, 33, i, (strand < 0.5) ? true : false,
                choiceRef->getFileName(), scaffold.getName(),
                locationStop - static_cast<int>(read2_seq.length()),
                locationStop);
    // Quality
    read1.initQual(&random_generator);
    read2.initQual(&random_generator);
    // Errors sequencing
    if (args.isProfileError()) {
      read1.makeErrors(&random_generator, args.getProfileError());
      read2.makeErrors(&random_generator, args.getProfileError());
    }
    // Write output
    if (read1.getStrand() == 0) {
      read_ser_1.writeRecord(read1);
      read_ser_2.writeRecord(read2);
    } else {
      read_ser_2.writeRecord(read1);
      read_ser_1.writeRecord(read2);
    }

    choiceRef->minus();
    if (choiceRef->getNbReadsRemaining() == 0) {
      args.removeReferences(ix_org);
    }

    count = 0;
    scaffolds.clear();
    intervals.clear();
    length = 0;
    ix_scaff = 0;
    locationStart = 0;
    locationStop = 0;
    sizeInsert = 0;
    seq = "";
    scaffold_find = false;
  }
  read_ser_1.close();
  read_ser_2.close();
  return 0;
}

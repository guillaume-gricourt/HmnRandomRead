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

using namespace std;

// Helpers
RandomRead::RandomRead(Args& arg)
    : args(arg),
      random_generator(arg.getSeed()),
      read_ser_1(std::move(args.getRead1())),
      read_ser_2(std::move(args.getRead2())) {}

int RandomRead::init() {
    // Make index
    makeIndex();
    cout << "Make index" << endl;
    // Profile Diversity
    if (args.isProfileDiversity()) {
        args.getProfileDiversity()->parseCsv();
        for (auto& ref : *args.getReferences()) {
            if (args.getProfileDiversity()->count(ref->getIdDiversity()) > 1) {
                throw logic_error("Id profile diversity " +
                                  ref->getIdDiversity() + " for reference " +
                                  ref->getFilePath() + " doesn't exist");
            }
        }
    }
    // Profile Error
    if (args.isProfileError()) {
        args.getProfileError()->parseCsv(args.getProfileErrorId(),
                                         args.isOutputPaired());
    }
    // Make reads
    makeReads();
    return 0;
}
int RandomRead::makeIndex() {
    int min_length_scaffold = static_cast<int>(args.getLengthReads() * 2 / 3);

    for (auto& ref : *args.getReferences()) {
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

    long long int nbReads(0);
    int nb_ref(0);
    int ix_org(0);
    // Calculate total nb reads
    for (auto& ref : *args.getReferences()) {
        nbReads += ref->getNbReads();
    }

    // cout << "total reads " << nbReads << endl;
    // Run
    int count(0);
    vector<Scaffold> scaffolds;
    // vector<Scaffold> scaffolds_keep;
    vector<long long int> intervals;
    long long int length(0);
    long long int ix_scaff(0);

    long long int locationStart(0), locationStop(0);
    long long int sizeInsert(0);
    string seq = "";
    bool scaffold_find = false;

    for (long int i = 0; i < nbReads; ++i) {
        Scaffold scaffold;

        // cout << "Start for loop" << endl;
        // Find Org - Get ref name unique
        nb_ref = static_cast<int>(args.getReferences()->size());

        // cout << "Nb ref " << nb_ref << endl;
        ix_org = 0;

        if (nb_ref > 1) {
            ix_org = random_generator.randomRange(0, nb_ref - 1);
        }
        // cout << "ix org " << ix_org << endl;

        auto& choiceRef = args.getReferences()->at(ix_org);

        // cout << "Choice ref " << choiceRef->getFileName() << endl;

        // cout << "min len scaff " << min_length_scaffold << endl;

        while (count < 10) {
            // Find Scaffold

            // From :
            // https://softwareengineering.stackexchange.com/questions/150616/get-weighted-random-item
            // cout << "choice scaff" << endl;
            scaffolds = choiceRef->getScaffolds().getVector();

            length = choiceRef->getScaffolds().getTotalLength();
            intervals = choiceRef->getScaffolds().getIntervals();

            // cout << "Length " << length << endl;
            // cout << "intervals " << intervals.size() << endl;
            /*
            for(auto& scaff : scaffolds){

                    scaffolds_keep.push_back(scaff);
                    length += scaff.getLength();
                    intervals.push_back(length);
            }
            */
            // cout << "Length scaff " << length << endl;

            ix_scaff = random_generator.randomRangeLong(0L, length);

            // cout << "Ix scaff " << ix_scaff << endl;

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
            sizeInsert = round(
                random_generator.randomNormal((double)args.getMeanInsertSize(),
                                              (double)args.getStdInsertSize()));
            // cout << "Size insert " << sizeInsert << " " <<
            // scaffold.getLength() << " " << scaffold.getName() << endl;
            // TODO Loop on choice
            // Check if size in scaffolds

            // Defined loc
            if (sizeInsert > scaffold.getLength()) {
                locationStart = scaffold.getStart();
                locationStop = scaffold.getStop();
            } else {
                // normalize
                // cout << "Scaff start " << scaffold.getStart() << " " <<
                // scaffold.getStop() - static_cast<long long int>(sizeInsert)
                // << endl;

                locationStart = random_generator.randomRangeLong(
                    scaffold.getStart(),
                    scaffold.getStop() -
                        static_cast<long long int>(sizeInsert));
                locationStop = locationStart + sizeInsert;
            }

            // Retrieve Insert-Sequence
            stringstream bufRoi;
            bufRoi << scaffold.getName() << ":" << locationStart << "-"
                   << locationStop;

            // cout << "Before seq " << bufRoi.str() <<  endl;

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
        // cout << "Before sequence" << endl;
        Sequence sequence(seq);
        // cout << "Before make mutation" << endl;
        // cout << "Sequence " << sequence.toString() << endl;
        if (args.getProfileDiversity() && choiceRef->getIdDiversity() != "") {
            sequence.makeMutation(&random_generator,
                                  args.getProfileDiversity()->getDiversity(
                                      choiceRef->getIdDiversity()));
        }
        // Strand : F or R
        // cout << "Before create string read" << endl;
        string read1_seq(""), read2_seq("");

        double strand = random_generator.randomRange(0.0, 1.0);

        // cout << "Strand " << strand << endl;
        // cout << "Sequence " << sequence.toString() << endl;
        if (strand >= 0.5) {
            // cout << "Before r1" << endl;
            read1_seq =
                sequence.getSequence(args.getLengthReads(), true, false);
            // cout << "Before r2" << endl;
            read2_seq =
                sequence.getSequence(args.getLengthReads(), false, true);
            // cout << "after r2" << endl;
        } else {
            read1_seq =
                sequence.getSequence(args.getLengthReads(), false, true);
            read2_seq =
                sequence.getSequence(args.getLengthReads(), true, false);
        }
        // cout << "R1 " << read1_seq << endl;
        // cout << "R2 " << read2_seq << endl;

        // Create record
        // cout << "Before create Fastq" << endl;
        Fastq read1(read1_seq, 33, i, (strand >= 0.5) ? true : false,
                    choiceRef->getFileName(), scaffold.getName(), locationStart,
                    locationStart + static_cast<int>(read1_seq.length()));
        Fastq read2(read2_seq, 33, i, (strand < 0.5) ? true : false,
                    choiceRef->getFileName(), scaffold.getName(),
                    locationStop - static_cast<int>(read2_seq.length()),
                    locationStop);

        // Quality
        // cout << "Before create quality" << endl;
        read1.initQual(&random_generator);
        read2.initQual(&random_generator);

        // Errors sequencing
        // cout << "Before make error" << endl;
        if (args.isProfileError()) {
            read1.makeErrors(&random_generator, args.getProfileError());
            read2.makeErrors(&random_generator, args.getProfileError());
        }

        // Write output
        // cout << "Before write output" << endl;
        if (read1.getStrand() == 0) {
            read_ser_1.writeRecord(read1);
            read_ser_2.writeRecord(read2);
        } else {
            read_ser_2.writeRecord(read1);
            read_ser_1.writeRecord(read2);
        }

        /********
        ** MAJ **
        ********/
        // cout << "Before MAJ" << endl;
        // cout << "choice ref " << choiceRef->getNbReadsRemaining()<< endl;
        choiceRef->minus();
        if (choiceRef->getNbReadsRemaining() == 0) {
            cout << "Remove " << choiceRef->getFilePath() << " " << ix_org
                 << endl;
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
    // cout << "Before close" << endl;
    read_ser_1.close();
    read_ser_2.close();
    return 0;
}


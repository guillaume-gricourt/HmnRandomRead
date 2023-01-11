// Copyright 2022 guillaume-gricourt
#include <getopt.h>

#include <string>

#include "args.hpp"
#include "randomread.hpp"
#include "version.hpp"

int main(int argc, char *argv[]) {
  Args args(argv[0]);

  // Command line
  std::string arg;

  if (argc < 2) {
    args.showUsage("");
    return 0;
  } else {
    for (int i = 1; i < argc; ++i) {
      arg = argv[i];
      if (arg == "help" || arg == "--help" || arg == "-h") {
        args.showUsage("");
        return 0;
      } else if (arg == "-v" || arg == "--version") {
        args.showVersion(HMNRANDOMREAD_VERSION);
        return 0;
      } else if (arg == "-r" || arg == "--reference") {
        if (i + 1 < argc) {
          args.addReference(argv[++i]);
        } else {
          args.showUsage("Reference must be filled");
        }
      } else if (arg == "-lengthReads" || arg == "--length-reads") {
        if (i + 1 < argc) {
          args.setLengthReads(argv[++i]);
        } else {
          args.showUsage("Length Reads must be indicated");
        }
      } else if (arg == "-meanInsert" || arg == "--mean-insert-size") {
        if (i + 1 < argc) {
          args.setMeanInsertSize(argv[++i]);
        } else {
          args.showUsage("Mean Insert Size is not validated");
        }
      } else if (arg == "-stdInsert" || arg == "--std-insert-size") {
        if (i + 1 < argc) {
          args.setStdInsertSize(argv[++i]);
        } else {
          args.showUsage("Standard Deviation Insert Size is not valid");
        }
      } else if (arg == "-r1" || arg == "--read-forward") {
        if (i + 1 < argc) {
          args.setRead1(argv[++i]);
        } else {
          args.showUsage("Format Read forward is not valid");
        }
      } else if (arg == "-r2" || arg == "--read-reverse") {
        if (i + 1 < argc) {
          args.setRead2(argv[++i]);
        } else {
          args.showUsage("Format Read reverse is not valid");
        }
      } else if (arg == "-profileDiversity" || arg == "--profile-diversity") {
        if (i + 1 < argc) {
          args.setProfileDiversity(argv[++i]);
        } else {
          args.showUsage("File profile diversity must be indicated");
        }
      } else if (arg == "-profileError" || arg == "--profile-error") {
        if (i + 1 < argc) {
          args.setProfileError(argv[++i]);
        } else {
          args.showUsage("File profile error must be indicated");
        }
      } else if (arg == "-profileErrorId" || arg == "--profile-error-id") {
        if (i + 1 < argc) {
          args.setProfileErrorId(argv[++i]);
        } else {
          args.showUsage("Error profile ID must be indicated");
        }
      } else if (arg == "-s" || arg == "--seed") {
        if (i + 1 < argc) {
          args.setSeed(argv[++i]);
        } else {
          args.showUsage("Seed number is not valid");
        }
      } else {
        args.showUsage("Unrecognized argument: " + arg);
      }
    }
  }

  // Check if ok
  if (!args.isValid()) {
    args.showUsage("Arguments are not valid");
  } else {
    if (!args.getUsageVisible()) {
      args.showStart(HMNRANDOMREAD_VERSION);
      RandomRead randomRead(args);
      int res = randomRead.init();
      args.showEnd();
      return res;
    }
  }
  return 0;
}

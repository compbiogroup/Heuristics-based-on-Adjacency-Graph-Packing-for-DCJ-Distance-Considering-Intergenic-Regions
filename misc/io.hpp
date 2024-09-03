#pragma once
#include "../cycle/cycles.hpp"
#include "../heur/ga.hpp"
#include "genome.hpp"
#include <iostream>

template <class T> void print_vec(ostream &os, vector<T> vec) {
  for (auto a : vec) {
    os << a << " ";
  }
}

vector<string> *read_lines(istream &is);
GenomePair input(string &line1, string &line2, bool extend);
GenomePair input(string &line1, string &line2, string &line3, string &line4,
                 bool extend);
MultiGenomePair input_mg(string &line1, string &line2, bool extend);
MultiGenomePair input_mg(string &line1, string &line2, string &line3,
                         string &line4, bool extend);
void output(ostream &os, int dist, double time);
void output(ostream &os, string &dist_etc, double time);
void output(ostream &os, const CycleGraph &, double time);
void output(ostream &os, GenomePair &data);
void output(ostream &os, PermsIrs);
void output(ostream &os, CommonPartition);

const string input_description =
    "input file (if not provided stdin is "
    "used). Each 4 lines of the input file correspond to a instance, "
    "each line has a list of space separated values, and represent in "
    "order the origin strings, the origin intergenic regions lists, the "
    "target strings, and the target intergenic regions lists. Mutliple "
    "chromossomes must be separeted by the charactere ; both in the strings "
    "and intergenic regions lists. An lenear chromossome may be proceeded by "
    "the charactere L and an circular chromossome must be proceeded by an "
    "charactere C. The intergenic region between the first and last gene of a "
    "circular chromossome must be the last on the list. If the -z"
    "option is used, each 2 lines correspond to a instance.";

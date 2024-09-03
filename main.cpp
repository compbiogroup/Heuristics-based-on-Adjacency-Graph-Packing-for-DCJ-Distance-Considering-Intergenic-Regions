#include <getopt.h>

#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>

#include "cycle/cycles.hpp"
#include "heur/ga.hpp"
#include "misc/io.hpp"
#include "misc/timer.hpp"
namespace fs = experimental::filesystem;
using namespace std;

#define N_POS_ARGS 1

struct Args {
  string heuristic;
  string input_file;
  string output_folder;
  int iterations = 100;
  bool extend = false;
  bool all = false;
  int tournament_size = 2;
  int mutation_rate = 50;
  int crossover_rate = 50;
  bool fill_zero = false;
};

void help(char *name) {
  cout
      << "usage: Find a cycle packing for the generalized adjacency graph "
         "formed "
         "from the origin and target signed genomes."
      << endl
      << "\t" << name << " HEUR [OPTIONS]" << endl
      << endl
      << "positional arguments:" << endl
      << "\tHEUR                    the heuristic to use (ga|rand)" << endl
      << endl
      << "optional arguments:" << endl
      << "\t-h, --help              show this help message and exit" << endl
      << "\t-i, --input INPUT       " << input_description
      << " Multiple Chromossome are separated dy semi-colon. An L in the "
         "beginning of the Chromossome indicates that it is linear and an C "
         "indicates that it is circular (If there is no letter it is assumed "
         "linear)."
      << endl
      << "\t-o, --output OUTPUT     output folder (if not provided stdout is "
         "used)"
      << endl
      << "\t-k, --iterations ITER   number of iterations (default 100)" << endl
      << "\t-m, --mutation MUT      mutation rate in percentage for GA "
         "(default 50%)"
      << endl
      << "\t-c, --crossover CROS    crossover rate in percentage for GA "
         "(default 50%)"
      << endl
      << "\t-t, --tournament TOUR   number of elements in each tournament for "
         "GA (default 2)"
      << endl
      << "\t-e, --extend            whether to extend (cap) the genomes before "
         "apply the algorithm. Should always be used with multiple "
         "chromossomes and is ignore in circular chromossomes. In "
         "unichromossomal linear genomes the genomes may "
         "already be extended."
      << endl
      << "\t-a, --all               print all cycles" << endl
      << "\t-z, --zero              fill intergenic regions with zero instead "
         "of read them"
      << endl
      << endl;

  exit(EXIT_SUCCESS);
}

void get_args(Args &args, int argc, char *argv[]) {
  extern char *optarg;
  extern int optind;
  int n_pos_args = 0;

  struct option longopts[] = {
      {"input", 1, NULL, 'i'},      {"output", 1, NULL, 'o'},
      {"iterations", 1, NULL, 'k'}, {"mutation", 1, NULL, 'm'},
      {"crossover", 1, NULL, 'c'},  {"tournament", 1, NULL, 't'},
      {"extend", 0, NULL, 'e'},     {"all", 0, NULL, 'a'},
      {"help", 0, NULL, 'h'},       {"zero", 0, NULL, 'z'}};

  char op;
  while ((op = getopt_long(argc, argv, "i:o:k:m:c:t:heaz", longopts, NULL)) !=
         -1) {
    switch (op) {
    case 'i':
      args.input_file = optarg;
      break;
    case 'o':
      args.output_folder = optarg;
      break;
    case 'k':
      args.iterations = atoi(optarg);
      break;
    case 'm':
      args.mutation_rate = atoi(optarg);
      break;
    case 'c':
      args.crossover_rate = atoi(optarg);
      break;
    case 't':
      args.tournament_size = atoi(optarg);
      break;
    case 'e':
      args.extend = true;
      break;
    case 'a':
      args.all = true;
      break;
    case 'z':
      args.fill_zero = true;
      break;
    default:
      help(argv[0]);
    }
  }
  for (int i = optind; i < argc; i++) {
    args.heuristic = argv[i];
    n_pos_args++;
  }

  if (n_pos_args != N_POS_ARGS) {
    help(argv[0]);
  }
}

CycleGraph *get_best_cg(CycleGraph *cg1, CycleGraph *cg2) {
  if (cg1->gain() > cg2->gain()) {
    delete cg2;
    return cg1;
  } else {
    delete cg1;
    return cg2;
  }
}

int main(int argc, char *argv[]) {
  Args args;
  ifstream is;
  unique_ptr<vector<string>> input_lines;

  get_args(args, argc, argv);

  // set generator
  random_device rd;
  mt19937 rand_gen(rd());
  /* mt19937 rand_gen(0); */

  if (args.input_file != "") {
    is.open(args.input_file);
    input_lines.reset(read_lines(is));
    is.close();
  } else {
    input_lines.reset(read_lines(cin));
  }

  // div controls how to read lines, either an instance is 2 lines (if
  // fill_zero) or an instance is 4 lines.
  int div = (args.fill_zero) ? 2 : 4;
  try {
    if (input_lines->size() % div == 1) {
      throw invalid_argument(
          "Number of lines is not multiple of 4 (or 2 if fill_zero).");
    }

#pragma omp parallel for // Process each instance in parallel
    for (size_t i = 0; i < input_lines->size(); i += div) {
      Timer timer;
      ofstream os;
      unique_ptr<CycleGraph> cg, cg_aux, cg_best;
      int name_idx = i / div;

      MultiGenomePair data;
      if (args.fill_zero) {
        data = input_mg((*input_lines)[i], (*input_lines)[i + 1], args.extend);
      } else {
        data =
            input_mg((*input_lines)[i], (*input_lines)[i + 1],
                     (*input_lines)[i + 2], (*input_lines)[i + 3], args.extend);
      }

      cg = make_unique<CycleGraph>(*data.g, *data.h, rand_gen);

      if (args.output_folder != "" && args.all) {
        os.open((args.output_folder / fs::path(args.input_file).filename())
                    .string() +
                string(5 - to_string(name_idx).size(), '0') +
                to_string(name_idx) + "-all");
      }

      if (args.heuristic == "rand") {
        cg_aux.reset(new CycleGraph(*cg));
        cg_aux->decompose(false);
        cg_best = move(cg_aux);
        if (args.all) {
          output((args.output_folder != "") ? os : cout, *cg_aux,
                 timer.since_last_mark());
        }

        for (int i = 1; i < args.iterations; ++i) {
          timer.mark_time();
          cg_aux.reset(new CycleGraph(*cg));
          cg_aux->decompose(true);
          if (args.all) {
            output((args.output_folder != "") ? os : cout, *cg_aux,
                   timer.since_last_mark());
          }
          if (cg_aux->gain() > cg_best->gain()) {
            cg_best = move(cg_aux);
          }
        }

      } else if (args.heuristic == "ga") {
        ostream *ga_os;
        if (args.all) {
          ga_os = (args.output_folder != "") ? &os : &cout;
        } else {
          ga_os = nullptr;
        }
        int start = args.iterations / 10;
        if (start % 2 == 1) {
          start += 1;
        }
        GA ga =
            GA(make_unique<Chromossome>(*cg), args.mutation_rate / 100.0,
               args.crossover_rate / 100.0, args.tournament_size, start, start,
               (args.iterations - start) / start, ga_os, timer, rand_gen);
        ga.solve(timer);
        cg_best = ga.get_best_chr();
      } else {
        help(argv[0]);
      }

      if (args.output_folder != "") {
        os.close();
        os.open((args.output_folder / fs::path(args.input_file).filename())
                    .string() +
                string(5 - to_string(name_idx).size(), '0') +
                to_string(name_idx) + "-best");
      }

      output((args.output_folder != "") ? os : cout, *cg_best,
             timer.elapsed_time());
      timer.mark_time();

      // Add output for dcj distance
      if (args.output_folder != "") {
        os.close();
        os.open((args.output_folder / fs::path(args.input_file).filename())
                    .string() +
                string(5 - to_string(name_idx).size(), '0') +
                to_string(name_idx) + "-dcj");
      }
      output((args.output_folder != "") ? os : cout, cg_best->dcj_ir_distance(),
             timer.since_last_mark());
    }

  } catch (const invalid_argument &e) {
    cerr << "Something went wrong!!!" << endl;
    cerr << e.what() << endl;
  }

  return 0;
}

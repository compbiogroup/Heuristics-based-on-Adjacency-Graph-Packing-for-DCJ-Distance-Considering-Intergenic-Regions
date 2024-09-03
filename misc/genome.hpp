#pragma once

#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>
using namespace std;

typedef int Gene;
typedef pair<Gene, bool> Genea; // Gene with an indicator whether it is alpha or
                                // not (alpha genes must be deleted/inserted)
typedef int IR;

/* String representing a unichromossomal and linear or circular genome (first
 * gane has index 1). Circular genomes ignore the extend argument. In circular
 * chromosomes the intergenic regions between the extremities should be the last
 * of the list.
 */

class Genome {
protected:
  bool circular;
  unique_ptr<vector<Genea>> genes;
  unique_ptr<vector<IR>> intergenic_regions;
  unique_ptr<vector<vector<Gene>>> positions;
  vector<Gene> empty_vec; // empty_vector to return when a gene is not found
  Gene op_max;
  Genome(){};
  void record_positions();

public:
  Genome(string str_g, bool extend, int op_max_init = 1, bool circular = false);
  Genome(string str_g, string str_i, bool extend, int op_max_init = 1,
         bool circular = false);
  Genome(vector<Genea> gs, vector<IR> irs, int op_max_init = 1,
         bool circular = false);
  Genome(int *g_vals, int *irs, int n, int op_max_init = 1,
         bool circular = false);
  Genome(const Genome&);
  Genome& operator=(const Genome&) = delete;
  virtual ~Genome() = default;
  size_t size() const { return genes->size(); }
  /* Only read access to the elements. */
  Gene operator[](int i) const { return (*genes)[i - 1].first; }
  /* Only read access to intergenic regions. */
  IR get_ir(int i) const { return (*intergenic_regions)[i - 1]; }
  /* Only read intergenic regions, but returns 0 for invalid positions. */
  IR safe_get_ir(int i) const {
      if (i <= 0 || i > size()) return 0;
      else return (*intergenic_regions)[i - 1];
  }
  /* Check if the gene in position i is alpha (must be deleted) */
  IR check_alpha(int i) const { return (*genes)[i - 1].second; }
  /* Get maximum occurrence of a label */
  int occ(Gene label) const;
  /* Get maximum occurrence */
  int occ_max() const;
  /* Get positions of a given label (read only). */
  const vector<int> &pos(Gene label) const;
  /* Insertion operation
  (insert each genes of the gene vector followed by each ir from the ir vector)
*/
  void insertion(int, vector<Gene> &);
  void insertion(int, vector<Gene> &, vector<IR> &);
  /* Deletion operation */
  void deletion(int, int, IR);
  /* Reversal operation */
  void reversal(int, int, IR, IR);
  /* Transposition operation */
  void transposition(int, int, int, IR, IR, IR);
  /* Change the label of a gene */
  void replace_label(int, Gene);
  /* Pretty print genome */
  virtual void serialize(ostream &) const;
  /* Get and set a value bigger then all the labels */
  Gene get_op_max() const { return op_max; }
  /* Get the genes alphabet */
  void alphabet(unordered_map<int, int> &alf, bool negative) const;
  /* Check if two genomes hava the same genes */
  bool balanced(Genome &) const;
  /* Check if genome is linear */
  bool is_linear() const { return !circular; }
  /* Check if genome is circular */
  bool is_circular() const { return circular; }
  /* Transforme linear in a circular genome */
  void make_circular();
  /* Get the permutation as a c array */
  int *get_perm() const;
  /* Get list of genes */
  vector<Genea> get_genes() const { return *genes; }
  /* Get list of intergenic regions */
  vector<IR> get_irs() const { return *intergenic_regions; }
};

ostream &operator<<(ostream &os, const Genome &g);

struct GenomePair {
  unique_ptr<Genome> g;
  unique_ptr<Genome> h;
};

/* String representing a multichromossomal genome (first chromossome has index
 * 1) */
class MultiGenome {
private:
  vector<unique_ptr<Genome>> chromosomes;
  Gene op_max;

public:
  MultiGenome(unique_ptr<Genome> &&g) {
    op_max = g->get_op_max();
    chromosomes.push_back(std::move(g));
  }
  MultiGenome(const MultiGenome&);
  MultiGenome(string str_g, bool extend);
  MultiGenome(string str_g, string str_i, bool extend);
  /* The bool in the genome pair indicates if the genome is circular or not */
  MultiGenome(vector<pair<bool, vector<Genea>>> gss, vector<vector<IR>> irss);
  /* Full access to the elements. */
  Genome &get_genome(int i) { return *(chromosomes[i - 1]); }
  /* Only read access to the elements. */
  Genome &operator[](int i) const { return *(chromosomes[i - 1]); }
  /* Break chromosome c at position i in two genomes (add extensions if the chromosome is linear) */
  void break_genome(int c, int i);
  /* Get and set a value bigger then all the labels */
  Gene get_op_max() const { return op_max; }
  size_t size() const;
  size_t n_chrom() const { return chromosomes.size(); }
  /* Pretty print genomes */
  void serialize(ostream &) const;
  /* Get maximum occurrence of a label */
  int occ(Gene label) const {
    int sum = 0;
    for (auto &chrom : chromosomes)
      sum += chrom->occ(label);
    return sum;
  }
  /* Count the number of linear chromosomes in the genome */
  int num_lin() const {
    int sum = 0;
    for (size_t i = 1; i <= (*this).n_chrom(); ++i) {
      if ((*this)[i].is_linear()) {
        sum++;
      }
    }
    return sum;
  }
};

ostream &operator<<(ostream &os, const MultiGenome &g);

struct MultiGenomePair {
  unique_ptr<MultiGenome> g;
  unique_ptr<MultiGenome> h;
};

#pragma once

#include "../misc/genome.hpp"
#include <bitset>
#include <map>
#include <queue>
#include <random>
#include <set>

typedef int Vtx_id;
#define NO_EDGE -1

typedef enum { Head, Tail, Cap } VtxType;

struct Vertex {
  Gene gene_val;   // value of the gene correspondent to this vertex
  bool sign_positive = false; // record the sign of the correspondent gene
  VtxType vtx_type;     // if it is the head or tail vertex of the correspondent
                        // gene, or if it is a gene added during capping
  Vtx_id black = NO_EDGE;
  int weight = 0; // weight of the incident black edge, negative values for
                  // origin genome and positive for target
  Vtx_id fix_gray = NO_EDGE;
  Vtx_id red = NO_EDGE;
  vector<Vtx_id> grays;
  bool in_cycle = false;      // indicate if the vertex is in a cycle
  bool is_indel = false;
  int indel_update = 0; // how to update the indel count if adding the red edge
                        // incident to this vertex. Can be used to indicate if
                        // the correspondent gene is from G (-1) or H (1).
  Vertex(Gene gene_val, bool sign_positive, VtxType vtx_type, Vtx_id black, int weight) : gene_val(gene_val), sign_positive(sign_positive), vtx_type(vtx_type), black(black), weight(weight), grays() {}
};

/* Raw representation of permutations to interact with other programs */
struct PermsIrs {
  size_t s_n, p_n;
  int *s;
  int *s_ir;
  int *p;
  int *p_ir;
};

struct CommonPartition {
  vector<int> breaks_s;
  vector<int> breaks_p;
};

struct QEntry {
  Vtx_id vtx;
  map<Vtx_id, Vtx_id> fixed;
  set<Vtx_id> vizited;
  map<Gene, int> indel_count;

  // Constructor for the first entry
  QEntry(Vtx_id start, map<Gene, int> indel_count)
      : fixed(), vizited(), indel_count(indel_count) {
    vtx = start;
  }

  // Constructor for entry with head-tail correspondence
  QEntry(Vtx_id v, Vtx_id u, Vtx_id alter_v, Vtx_id alter_u,
         map<Vtx_id, Vtx_id> fixed_old, set<Vtx_id> vizited_old,
         map<Gene, int> indel_count_old)
      : fixed(fixed_old), vizited(vizited_old), indel_count(indel_count_old) {
    vtx = u;
    fixed[v] = u;
    fixed[u] = v;
    fixed[alter_v] = alter_u;
    fixed[alter_u] = alter_v;
    vizited.insert(u);
    vizited.insert(v);
  }

  // Constructor for entry without head-tail correspondence
  QEntry(Vtx_id v, Vtx_id u, map<Vtx_id, Vtx_id> fixed_old,
         set<Vtx_id> vizited_old, map<Gene, int> indel_count_old)
      : fixed(fixed_old), vizited(vizited_old), indel_count(indel_count_old) {
    vtx = u;
    fixed[v] = u;
    fixed[u] = v;
    vizited.insert(u);
    vizited.insert(v);
  }

  bool not_yet_visited(const set<Vtx_id> &vizited_in_level, Vtx_id v) {
    return vizited_in_level.find(v) == vizited_in_level.end() &&
           vizited.find(v) == vizited.end();
  }
};

struct Run {
  Gene fst_gene; // first gene of the run
  char genome;   // G or H
  vector<Gene> genes_to_add;
};

class CycleGraph {

private:
  // Counters for types of cycles
  int balanced_cycles = 0;
  int divergent_cycles = 0;
  int divergent_balanced_or_neg_cycles = 0;
  int unitary_balanced_cycles = 0;

  int indel_potation = 0;
  /* id of the last vertex correspondent to a gene in
   * the origin/target genome */
  int last_origin_vertex;
  int last_target_vertex;

  /* Whether to consider only unitary cycles in the gain and decomposition
   * functions */
  bool only_unitary = false;
  /* Value bigger then all the labels */
  int op_max;
  /* Map from position on the genomes to the vertex id */
  vector<Vertex> vertices;
  vector<tuple<size_t, int, Vtx_id>>
      cycles; // we indentify cycles by one of their vertices the first value is
              // used to sort the cycles (may be their size of may involve the
              // potation for instance, see the definiton in make_cycle)
  /* Count the number of genes with positive sign for origin genome and negative
   * sign for target. The final map will indicate the number of each gene that
   * must be inserted or deleted. */
  map<Gene, int> indel_count;
  mt19937 rand_gen;

  void follow_indel_edge(set<Vtx_id> &vizited_in_level,
                         vector<unique_ptr<QEntry>> &q1);
  bool is_chromosome_end(Vtx_id v) const;
  /* Select a cycle with a bfs
   * Arguments:
   *     start - initial vertex
   *     is_random - whether to use random approach
   */
  /* Get string with the cycles from the decomposition */
  void bfs(Vtx_id start, bool is_random);
  void get_trivial_cycle(Vtx_id start, bool is_random);

public:
  /* Initial Constructor. */
  CycleGraph(const MultiGenome &origin, const MultiGenome &target,
             mt19937 rand_gen, bool only_unitary = false);
  /* Copy Constructor. */
  CycleGraph(const CycleGraph &that)
      : balanced_cycles(that.balanced_cycles),
        divergent_cycles(that.divergent_cycles),
        divergent_balanced_or_neg_cycles(that.divergent_balanced_or_neg_cycles),
        unitary_balanced_cycles(that.unitary_balanced_cycles),
        indel_potation(that.indel_potation),
        last_origin_vertex(that.last_origin_vertex),
        last_target_vertex(that.last_target_vertex),
        only_unitary(that.only_unitary), op_max(that.op_max),
        vertices(that.vertices), cycles(that.cycles),
        indel_count(that.indel_count), rand_gen(that.rand_gen) {}
  size_t size() const { return last_target_vertex + 1; };
  void decompose_random_map();
  void decompose(bool is_random);
  string show_cycles() const;
  /* Recover cycles from string */
  void read_cycles(string);
  /* Make assigment correspondent to the current decomposition by changing the
   * values of each vertex */
  /* indels are represented with a zero. */
  void make_assignment();
  /* Get genomes correspondent to the graph */
  MultiGenomePair get_genomes() const;
  /* Recover permutations from decomposition (only if we are using unichromosome
   * genomes) */
  PermsIrs get_perms() const;
  /* Recover common partition from decomposition */
  CommonPartition get_part(bool allow_rev) const;
  /* Get number of cycles from the decomposition */
  int dec_size() const { return cycles.size(); }
  void rem_cycle(tuple<size_t, int, Vtx_id> c);
  /* Get indel potation */
  int potation() const { return indel_potation; }
  /* Gain function, represents the value of the cycle packing */
  vector<int> gain() const {
    vector<int> val(1);
    /* val[1] = this->dec_size() - this->potation(); */
    /* val[0] = this->dec_size() - this->potation(); */
    // val[0] = this->dec_size() + this->balanced_cycles;
    /* val[0] = this->dec_size() + this->divergent_balanced_or_neg_cycles; */
    if (only_unitary) {
      val[0] = this->unitary_balanced_cycles;
    } else {
      val[0] = this->balanced_cycles;
    }
    return val;
  }
  /* Make representation of a cycle with the first value used for sorting */
  tuple<size_t, int, Vtx_id> make_cycle(vector<Vtx_id> cycle, int weight, int pot) {
    /* return pair<size_t, Vtx_id>(cycle.size() + pot, weight, cycle[0]); */
    return tuple<size_t, int, Vtx_id>(cycle.size(), weight, cycle[0]);
  }

  /* Verify if a cycle can be added */
  bool check_cycle(vector<Vtx_id> cycle);
  void add_cycle(vector<Vtx_id> cycle);
  void check_and_add_cycle(vector<Vtx_id> cycle);

  vector<tuple<size_t, int, Vtx_id>> cycle_list() const { return cycles; }
  vector<Vtx_id> get_cycle(Vtx_id i) const;
  int is_divergent(Vtx_id) const;
  int cycle_weight(Vtx_id i) const;
  int cycle_potation(Vtx_id i) const;
  Run cycle_run(int i) const;
  void serialize(ostream &) const;
  // Get vertex info by index (starting at zero)
  int get_gray_degree(Vtx_id i) const;
  int get_black_degree(Vtx_id i) const;
  int get_red_degree(Vtx_id i) const;
  int is_indel(Vtx_id i) const { return vertices[i].is_indel; }
  VtxType get_type(Vtx_id i) const { return vertices[i].vtx_type; }
  int get_value(Vtx_id i) const { return vertices[i].gene_val; }
  Vtx_id get_correspondent(Vtx_id i) const { return vertices[i].fix_gray; }
  bool in_G(Vtx_id i) const { return vertices[i].indel_update == -1; }

  // Distance algorithms based on the cycle graph
  // 2017-fertin-etal algorithm
  int dcj_ir_distance() const;
};

ostream &operator<<(std::ostream &os, const CycleGraph &cg);
ostream &operator<<(std::ostream &os, const Run &r);

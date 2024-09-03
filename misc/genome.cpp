#include "genome.hpp"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <numeric>
#include <sstream>
#include <unordered_map>

bool hasOnlySpaces(const string &str) {
  return all_of(str.begin(), str.end(), ::isspace);
}

Genome::Genome(string str_g, bool extend, int op_max_init, bool circular)
    : Genome(str_g, "", extend, op_max_init, circular) {}

Genome::Genome(string str_g, string str_i, bool extend, int op_max_init,
               bool circular_init)
    : circular(circular_init), empty_vec(), op_max(op_max_init) {
  string token;
  genes.reset(new vector<Genea>());
  intergenic_regions.reset(new vector<IR>());

  if (circular && extend) {
    throw invalid_argument("Circular genomes cannot be extended.");
  }

  /* Read each element of the string. */
  if (extend)
    genes->push_back(Genea(0, false));
  stringstream ssg(str_g);
  while (getline(ssg, token, ' ')) {
    if (!hasOnlySpaces(token)) {
      Genea a = Genea(stoi(token), false);
      genes->push_back(a);
      if (abs(genes->back().first) > op_max)
        op_max = abs(genes->back().first);
    }
  }
  for (int i = 0; i < int(genes->size()); i++) {
    Genea a = (*genes)[i];
    if (a.first == 0) { // Alpha (gene to be deleted)
      op_max++;
      (*genes)[i] = Genea(op_max, true);
    } else {
      /* Make sure the absolute value of each gene is at least 2 so we can use 1
       * in the extremity. */
      (*genes)[i] = Genea((a.first < 0) ? a.first - 1 : a.first + 1, false);
    }
  }
  if (extend) {
    (*genes)[0] = Genea(0, false);
    genes->push_back(Genea(1, false));
  }
  if (!extend && !circular) {
    (*genes)[0] = Genea(0, false);
    genes->back() = Genea(1, false);
  }
  op_max += 2; // Sum 1 to compensate for the dislocation and 1 to be above the
               // current position

  /* Read each intergenic regions. */
  if (str_i == "") {
    for (int i = 0; i < (circular) ? genes->size() : int(genes->size()) - 1;
         i++) {
      intergenic_regions->push_back(0);
    }
  } else {
    stringstream ssi(str_i);
    while (getline(ssi, token, ' ')) {
      if (!hasOnlySpaces(token)) {
        IR r = stoi(token);
        intergenic_regions->push_back(r);
      }
    }
  }

  record_positions();

  if (circular && intergenic_regions->size() != genes->size()) {
    throw invalid_argument("Circular genomes must have equal quantities of "
                           "intergenic regions and genes.");
  } else if (!circular && intergenic_regions->size() != genes->size() - 1) {
    throw invalid_argument("Linear genomes must have one less intergenic "
                           "region than the number of genes.");
  }
}

Genome::Genome(const vector<Genea> gs, const vector<IR> irs, int op_max_init,
               bool circular_init)
    : circular(circular_init), empty_vec(), op_max(op_max_init) {
  for (Genea gene : gs) {
    if (abs(gene.first) > op_max)
      op_max = abs(gene.first);
  }
  op_max++;

  genes = make_unique<vector<Genea>>(gs);
  intergenic_regions = make_unique<vector<IR>>(irs);

  record_positions();

  if (circular && intergenic_regions->size() == genes->size()) {
    throw invalid_argument("Circular genomes must have equal quantities of "
                           "intergenic regions and genes.");
  } else if (!circular && intergenic_regions->size() != genes->size() - 1) {
    throw invalid_argument("Linear genomes must have one less intergenic "
                           "region than the number of genes.");
  }
}

vector<Genea> build_Genea_vector(int *g_vals, int n) {
  vector<Genea> gs(n);
  for (int i = 0; i < n; i++) {
    gs[i] = Genea(g_vals[i], false);
  }
  return gs;
}

Genome::Genome(int *g_vals, int *irs, int n, int op_max, bool circular)
    : Genome(build_Genea_vector(g_vals, n), vector<IR>(irs, irs + n - 1),
             op_max, circular) {}

Genome::Genome(const Genome &g) : circular(g.circular), op_max(g.op_max) {
  genes = make_unique<vector<Genea>>(*g.genes);
  intergenic_regions = make_unique<vector<IR>>(*g.intergenic_regions);
  positions = make_unique<vector<vector<int>>>(*g.positions);
}

void Genome::record_positions() {
  /* Record list of positions for each label. */
  positions.reset(new vector<vector<Gene>>(op_max + 1));
  for (size_t i = 0; i < genes->size(); ++i) {
    (*positions)[abs((*genes)[i].first)].push_back(i + 1);
  }
}

int Genome::occ(Gene label) const { return this->pos(label).size(); }

int Genome::occ_max() const {
  vector<int> count(this->positions->size(), 0);
  for (size_t i = 1; i <= this->size(); ++i) {
    count[abs((*this)[i])]++;
  }
  return *max_element(count.begin(), count.end());
}

const vector<int> &Genome::pos(Gene label) const {
  if (abs(label) >= positions->size())
    return empty_vec;
  return (*positions)[abs(label)];
}

void Genome::insertion(int i, vector<Gene> &new_genes) {
  vector<IR> new_irs(new_genes.size());
  insertion(i, new_genes, new_irs);
}

void Genome::insertion(int i, vector<Gene> &new_genes_, vector<IR> &new_irs) {
  assert(1 <= i);
  assert(i <= int(size()));
  vector<Genea> new_genes;

  for (Gene a : new_genes_) {
    if (abs(a) >= op_max)
      op_max = abs(a) + 1;
    new_genes.push_back(Genea(a, false));
  }

  genes->insert(genes->begin() + i, new_genes.begin(), new_genes.end());
  intergenic_regions->insert(intergenic_regions->begin() + i - 1,
                             new_irs.begin(), new_irs.end());
  record_positions();
}

void Genome::deletion(int i, int j, IR x) {
  assert(2 <= i);
  assert(i < j);
  assert(j <= int(size()));
  assert(0 <= x &&
         x <= (*intergenic_regions)[i - 2] + (*intergenic_regions)[j - 2]);

  int k = i - 1, l = j - 1;
  (*intergenic_regions)[i - 2] = x;
  for (; l < int(size()) - 1; k++, l++) {
    (*genes)[k] = (*genes)[l];
    (*intergenic_regions)[k] = (*intergenic_regions)[l];
  }
  (*genes)[k] = (*genes)[l];
  genes->resize(k + 1);
  intergenic_regions->resize(k);
  record_positions();
}

void Genome::reversal(int i, int j, IR x, IR y) {
  assert(2 <= i);
  assert(i <= j);
  assert(j <= int(size()));
  assert(0 <= x && x <= (*intergenic_regions)[i - 2]);
  assert(0 <= y && y <= (*intergenic_regions)[j - 1]);

  int x_rest = (*intergenic_regions)[i - 2] - x;
  int y_rest = (*intergenic_regions)[j - 1] - y;

  for (int k = 0; k < (j - i + 1) / 2; k++) {
    swap((*genes)[i + k - 1], (*genes)[j - k - 1]);
    swap((*intergenic_regions)[i + k - 1], (*intergenic_regions)[j - k - 2]);
  }

  (*intergenic_regions)[i - 2] = x + y;
  (*intergenic_regions)[j - 1] = x_rest + y_rest;
  record_positions();
}

void Genome::transposition(int i, int j, int k, IR x, IR y, IR z) {
  assert(2 <= i);
  assert(i < j);
  assert(j < k);
  assert(k <= int(size()));
  assert(0 <= x && x <= (*intergenic_regions)[i - 2]);
  assert(0 <= y && y <= (*intergenic_regions)[j - 2]);
  assert(0 <= z && z <= (*intergenic_regions)[k - 2]);

  int x_rest = (*intergenic_regions)[i - 2] - x;
  int y_rest = (*intergenic_regions)[j - 2] - y;
  int z_rest = (*intergenic_regions)[k - 2] - z;

  vector<Genea> aux1;
  vector<IR> aux2;
  for (int t = i - 1; t <= j - 2; t++) {
    aux1.push_back((*genes)[t]);
    aux2.push_back((*intergenic_regions)[t]);
  }
  int t = i - 1;
  for (int r = j - 1; r <= k - 2; r++) {
    (*genes)[t] = (*genes)[r];
    (*intergenic_regions)[t] = (*intergenic_regions)[r];
    t++;
  }
  for (size_t r = 0; r < aux1.size(); r++) {
    (*genes)[t] = aux1[r];
    (*intergenic_regions)[t] = aux2[r];
    t++;
  }

  (*intergenic_regions)[i - 2] = x + y_rest;
  (*intergenic_regions)[i + k - j - 2] = z + x_rest;
  (*intergenic_regions)[k - 2] = y + z_rest;
  record_positions();
}

void Genome::replace_label(int idx, Gene label) {
  Genea old = (*genes)[idx - 1];
  (*genes)[idx - 1] = Genea((old.first < 0) ? -label : label, old.second);

  vector<int> &pos_vec = (*positions)[abs(old.first)];

  pos_vec.erase(find(pos_vec.begin(), pos_vec.end(), idx));
  if (int(positions->size()) <= label)
    positions->resize(label + 1);
  (*positions)[label].push_back(idx);

  if (abs(label) >= op_max)
    op_max = abs(label) + 1;
}

void Genome::serialize(ostream &os) const {
  for (size_t i = 1; i <= size() - 1; ++i) {
    os << "(" << (*this)[i] << ((check_alpha(i)) ? "*" : "") << ") - "
       << get_ir(i) << " - ";
  }
  os << "(" << (*this)[size()] << ")";
  if (circular) {
    os << " - " << get_ir(size()) << " - ";
  }
}

ostream &operator<<(ostream &os, const Genome &g) {
  g.serialize(os);
  return os;
}

void Genome::alphabet(unordered_map<int, int> &alf, bool negative) const {
  for (Genea a : *(this->genes)) {
    if (negative) {
      alf[abs(a.first)]--;
    } else {
      alf[abs(a.first)]++;
    }
  }
}

bool Genome::balanced(Genome &h) const {
  unordered_map<int, int> alf;
  this->alphabet(alf, false);
  h.alphabet(alf, true);
  bool ok = true;
  for (auto i : alf) {
    if (i.second != 0)
      ok = false;
  }
  int sum_ir_g = accumulate(this->intergenic_regions->begin(),
                            this->intergenic_regions->end(), 0);
  int sum_ir_h =
      accumulate(h.intergenic_regions->begin(), h.intergenic_regions->end(), 0);
  return ok and sum_ir_g == sum_ir_h;
}

void Genome::make_circular() {
  if (circular)
    return;
  circular = true;
  genes->pop_back();
  genes->erase(genes->begin());
  intergenic_regions->back() += intergenic_regions->front();
  intergenic_regions->erase(intergenic_regions->begin());
}

int *Genome::get_perm() const {
  int *perm = new int[size()];
  for (int i = 0; i < int(genes->size()); i++) {
    perm[i] = (*genes)[i].first;
    if (perm[i] >= 0)
      perm[i] += 1;
    else
      perm[i] -= 1;
  }
  return perm;
}

MultiGenome::MultiGenome(const MultiGenome &that) {
  for (size_t i = 0; i < that.chromosomes.size(); i++) {
    chromosomes.push_back(make_unique<Genome>(*that.chromosomes[i+1]));
  }
  op_max = that.op_max;
}

MultiGenome::MultiGenome(string str_g, bool extend)
    : MultiGenome(str_g, "", extend) {}

MultiGenome::MultiGenome(string str_g, string str_i, bool extend) {
  string token;
  vector<string> tokens_g;
  vector<string> tokens_i;

  stringstream ssg(str_g);
  while (getline(ssg, token, ';')) {
    if (!hasOnlySpaces(token)) {
      tokens_g.push_back(token);
    }
  }

  if (str_i != "") {
    stringstream ssi(str_i);
    while (getline(ssi, token, ';')) {
      if (!hasOnlySpaces(token)) {
        tokens_i.push_back(token);
      }
    }
    assert(tokens_g.size() == tokens_i.size());
  }

  bool circular;
  int beg;
  op_max = 1;
  for (size_t i = 0; i < tokens_g.size(); i++) {

    circular = false;
    beg = 0;
    int pos = tokens_g[i].find('C');
    if (pos != string::npos) {
      circular = true;
      beg = pos + 1;
    }
    pos = tokens_g[i].find('L');
    if (pos != string::npos) {
      beg = pos + 1;
    }

    if (str_i == "") {
      chromosomes.push_back(make_unique<Genome>(
          tokens_g[i].substr(beg), extend && !circular, op_max, circular));
    } else {
      chromosomes.push_back(unique_ptr<Genome>(
          new Genome(tokens_g[i].substr(beg), tokens_i[i], extend && !circular,
                     op_max, circular)));
    }
    op_max = chromosomes.back()->get_op_max();
  }
}

MultiGenome::MultiGenome(vector<pair<bool, vector<Genea>>> gss,
                         vector<vector<IR>> irss) {
  op_max = 1;
  for (size_t i = 0; i < gss.size(); i++) {
    chromosomes.push_back(unique_ptr<Genome>(
        new Genome(gss[i].second, irss[i], op_max, gss[i].first)));
    op_max = chromosomes.back()->get_op_max();
  }
}

void MultiGenome::break_genome(int c, int i) {
  Genome &g = get_genome(c);
  if (i == 0) {
    return;
  }
  vector<Genea> new_genes_1;
  vector<IR> new_irs_1;
  vector<Genea> new_genes_2;
  vector<IR> new_irs_2;

  if (g.is_linear()) {
    new_genes_2.push_back(Genea(0, false));
    new_irs_2.push_back(0);
  }
  vector<Genea> genes = g.get_genes();
  new_genes_1.insert(new_genes_1.end(), genes.begin(), genes.begin() + i);
  new_genes_2.insert(new_genes_2.end(), genes.begin() + i, genes.end());
  vector<IR> irs = g.get_irs();
  new_irs_1.insert(new_irs_1.end(), irs.begin(), irs.begin() + i);
  new_irs_2.insert(new_irs_2.end(), irs.begin() + i, irs.end());
  if (g.is_linear()) {
    new_genes_1.push_back(Genea(1, false));
  }

  bool circular = g.is_circular();
  chromosomes.erase(chromosomes.begin() + c - 1);
  chromosomes.insert(
      chromosomes.begin() + c - 1,
      make_unique<Genome>(new_genes_1, new_irs_1, get_op_max(), circular));
  chromosomes.insert(
      chromosomes.begin() + c,
      make_unique<Genome>(new_genes_2, new_irs_2, get_op_max(), circular));
}

size_t MultiGenome::size() const {
  int sum = 0;
  for (int i = 0; i < int(chromosomes.size()); i++) {
    sum += chromosomes[i]->size();
  }
  return sum;
}

void MultiGenome::serialize(ostream &os) const {
  for (size_t i = 0; i < chromosomes.size(); i++) {
    cout << *chromosomes[i] << " | ";
  }
  cout << endl;
}

ostream &operator<<(ostream &os, const MultiGenome &g) {
  g.serialize(os);
  return os;
}

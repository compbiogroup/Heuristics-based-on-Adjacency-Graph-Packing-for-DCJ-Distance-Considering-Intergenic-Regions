#include "cycles.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <sstream>

#include "../misc/genome.hpp"

CycleGraph::CycleGraph(const MultiGenome &origin, const MultiGenome &target,
                       mt19937 rand_gen_init, bool only_unitary_init)
    : only_unitary(only_unitary_init), vertices(), cycles(), indel_count(),
      rand_gen(rand_gen_init) {
  int a = max(origin.size() + 2 * origin.n_chrom(),
              target.size() + 2 * target.n_chrom());
  int b = max(origin.get_op_max(), target.get_op_max());
  op_max = max(a, b) + 1;
  int cap = op_max;
  op_max++;
  int v_id = 0;
  size_t beg, end;

  /* Count the number of linear chromosomes in the origin and target genomes */
  int n_lin_origin = origin.num_lin();
  int n_lin_target = target.num_lin();

  /* Using gene values to label the vertices, also indicate how to update indels
   * of each vertice. Furthermore, include the black and indel edges. */
  // vertices and edges of origin chromosomes
  for (size_t i = 1; i <= origin.n_chrom(); ++i) {
    if (origin[i].is_linear()) {
      vertices.push_back(
          Vertex(cap, true, Cap, v_id + 1, -origin[i].get_ir(1)));
      vertices[v_id].indel_update = -1;
      indel_count[cap]++;
      v_id++;
      beg = 2;
      end = origin[i].size() - 1;
    } else {
      beg = 1;
      end = origin[i].size();
    }
    int fst_v_id = v_id;
    for (size_t j = beg; j <= end; ++j) {
      if (origin[i][j] >= 0) {
        vertices.push_back(Vertex(abs(origin[i][j]), true, Tail, v_id - 1,
                                  -origin[i].safe_get_ir(j - 1)));
        vertices.push_back(Vertex(abs(origin[i][j]), true, Head, v_id + 2,
                                  -origin[i].get_ir(j)));
      } else {
        vertices.push_back(Vertex(abs(origin[i][j]), false, Head, v_id - 1,
                                  -origin[i].safe_get_ir(j - 1)));
        vertices.push_back(Vertex(abs(origin[i][j]), false, Tail, v_id + 2,
                                  -origin[i].get_ir(j)));
      }
      vertices[v_id].indel_update = -1;
      vertices[v_id].red = v_id + 1;
      vertices[v_id + 1].indel_update = -1;
      vertices[v_id + 1].red = v_id;
      indel_count[abs(origin[i][j])]++;
      v_id += 2;
    }
    if (origin[i].is_linear()) {
      vertices.push_back(Vertex(cap, true, Cap, v_id - 1,
                                -origin[i].get_ir(origin[i].size() - 1)));
      vertices[v_id].indel_update = -1;
      indel_count[cap]++;
      v_id++;
    } else {
      /* In circular chromosomes the first and last vertices are connected by a
       * black edge */
      vertices[v_id - 1].black = fst_v_id;
      vertices[fst_v_id].black = v_id - 1;
      vertices[v_id - 1].weight = origin[i].get_ir(origin[i].size());
      vertices[fst_v_id].weight = origin[i].get_ir(origin[i].size());
    }
  }
  // Add empty origin chromosomes
  for (int i = 0; i < n_lin_target - n_lin_origin; ++i) {
    vertices.push_back(Vertex(cap, true, Cap, v_id + 1, 0));
    vertices[v_id].indel_update = -1;
    indel_count[cap]++;
    v_id++;
    vertices.push_back(Vertex(cap, true, Cap, v_id - 1, 0));
    vertices[v_id].indel_update = -1;
    indel_count[cap]++;
    v_id++;
  }
  last_origin_vertex = v_id - 1;
  // vertices and edges of target chromosomes
  for (size_t i = 1; i <= target.n_chrom(); ++i) {
    if (target[i].is_linear()) {
      vertices.push_back(Vertex(cap, true, Cap, v_id + 1, target[i].get_ir(1)));
      vertices[v_id].indel_update = 1;
      indel_count[cap]--;
      v_id++;
      beg = 2;
      end = target[i].size() - 1;
    } else {
      beg = 1;
      end = target[i].size();
    }
    int fst_v_id = v_id;
    for (size_t j = beg; j <= end; ++j) {
      if (target[i][j] >= 0) {
        vertices.push_back(Vertex(abs(target[i][j]), true, Tail, v_id - 1,
                                  target[i].safe_get_ir(j - 1)));
        vertices.push_back(Vertex(abs(target[i][j]), true, Head, v_id + 2,
                                  target[i].get_ir(j)));
      } else {
        vertices.push_back(Vertex(abs(target[i][j]), false, Head, v_id - 1,
                                  target[i].safe_get_ir(j - 1)));
        vertices.push_back(Vertex(abs(target[i][j]), false, Tail, v_id + 2,
                                  target[i].get_ir(j)));
      }
      vertices[v_id].indel_update = 1;
      vertices[v_id].red = v_id + 1;
      vertices[v_id + 1].indel_update = 1;
      vertices[v_id + 1].red = v_id;
      indel_count[abs(target[i][j])]--;
      v_id += 2;
    }
    if (target[i].is_linear()) {
      vertices.push_back(Vertex(cap, true, Cap, v_id - 1,
                                target[i].get_ir(target[i].size() - 1)));
      vertices[v_id].indel_update = 1;
      indel_count[cap]--;
      v_id++;
    } else {
      /* In circular chromosomes the first and last vertices are connected by a
       * black edge */
      vertices[v_id - 1].black = fst_v_id;
      vertices[fst_v_id].black = v_id - 1;
      vertices[v_id - 1].weight = target[i].get_ir(target[i].size());
      vertices[fst_v_id].weight = target[i].get_ir(target[i].size());
    }
  }
  // Add empty target chromosomes
  for (int i = 0; i < n_lin_origin - n_lin_target; ++i) {
    vertices.push_back(Vertex(cap, true, Cap, v_id + 1, 0));
    vertices[v_id].indel_update = 1;
    indel_count[cap]--;
    v_id++;
    vertices.push_back(Vertex(cap, true, Cap, v_id - 1, 0));
    vertices[v_id].indel_update = 1;
    indel_count[cap]--;
    v_id++;
  }
  last_target_vertex = v_id - 1;

  /* Gray Edges */
  for (int i = 0; i <= last_origin_vertex; i++) {
    for (int j = last_origin_vertex + 1; j <= last_target_vertex; j++) {
      if (vertices[i].gene_val == vertices[j].gene_val &&
          vertices[i].vtx_type == vertices[j].vtx_type) {
        vertices[i].grays.push_back(j);
        vertices[j].grays.push_back(i);
      }
    }
  }
}

/* Decompose the graph by using a random assignment of the genes */
void CycleGraph::decompose_random_map() {
  Vtx_id headtailcorresp_v, headtailcorresp_u;

  for (int v = 0; v <= last_origin_vertex; v++) {
    vector<Vtx_id> neigs;
    for (Vtx_id u : vertices[v].grays) {
      if (vertices[u].fix_gray == NO_EDGE) {
        neigs.push_back(u);
      }
    }
    shuffle(neigs.begin(), neigs.end(), rand_gen);
    if (neigs.empty()) {
      throw runtime_error("No neighbor found in decompose_random_map.");
    }
    Vtx_id u = neigs.front();
    vertices[v].fix_gray = u;
    vertices[u].fix_gray = v;
    if (vertices[v].vtx_type != Cap) {
      headtailcorresp_v = vertices[v].red;
      headtailcorresp_u = vertices[u].red;
      assert(vertices[headtailcorresp_u].fix_gray == NO_EDGE);
      assert(vertices[headtailcorresp_v].fix_gray == NO_EDGE);
      vertices[headtailcorresp_v].fix_gray = headtailcorresp_u;
      vertices[headtailcorresp_u].fix_gray = headtailcorresp_v;
      v++;
    }
  }
  this->decompose(true);
}

/* Decompose the remaning graph */
void CycleGraph::decompose(bool is_random) {
  vector<Vtx_id> idxs(size());
  iota(idxs.begin(), idxs.end(), 0);
  if (is_random) {
    shuffle(idxs.begin(), idxs.end(), rand_gen);
  }
  if (only_unitary) {
    for (Vtx_id i : idxs) {
      get_trivial_cycle(i, is_random);
    }
  }
  for (Vtx_id i : idxs) {
    bfs(i, is_random);
  }
}

void CycleGraph::get_trivial_cycle(Vtx_id start, bool is_random) {
  vector<Vtx_id> red1, red2, cycle;
  Vtx_id b1 = vertices[start].black;
  Vtx_id b1_ = b1;
  map<Gene, int> indel_count = this->indel_count;

  if (vertices[start].in_cycle)
    return;

  /* Loop to maybe follow red edges */
  while (true) {
    vector<Vtx_id> neigs;
    if (vertices[b1_].fix_gray != NO_EDGE) {
      neigs.push_back(vertices[b1_].fix_gray);
    } else if (!vertices[b1_].is_indel) {
      for (Vtx_id u : vertices[b1_].grays) {
        if (vertices[u].fix_gray == NO_EDGE) {
          neigs.push_back(u);
        }
      }
    }

    for (Vtx_id g1 : neigs) {
      if (vertices[g1].in_cycle) {
        continue;
      }
      Vtx_id b2 = vertices[g1].black;
      Vtx_id b2_ = b2;

      /* Loop to maybe follow red edges */
      while (true) {
        vector<Vtx_id> neigs2;
        if (vertices[b2_].fix_gray != NO_EDGE) {
          neigs2.push_back(vertices[b2_].fix_gray);
        } else if (!vertices[b2_].is_indel) {
          for (Vtx_id u : vertices[b2_].grays) {
            if (vertices[u].fix_gray == NO_EDGE) {
              neigs2.push_back(u);
            }
          }
        }
        for (Vtx_id g2 : neigs2) {
          if (g2 == start && vertices[b1].weight + vertices[b2].weight == 0) {
            cycle.resize(0);
            cycle.push_back(start);
            cycle.push_back(b1);
            cycle.insert(cycle.end(), red1.begin(), red1.end());
            cycle.push_back(g1);
            cycle.push_back(b2);
            cycle.insert(cycle.end(), red2.begin(), red2.end());
          }
        }
        if (vertices[b2_].fix_gray != NO_EDGE ||
            vertices[vertices[b2_].red].in_cycle) {
          break;
        }
        Vtx_id r2 = vertices[b2_].red;
        if (indel_count[vertices[r2].gene_val] * vertices[r2].indel_update <
            0) {
          red2.push_back(r2);
          b2_ = vertices[r2].black;
          red2.push_back(b2_);
          indel_count[vertices[r2].gene_val] += vertices[r2].indel_update;
        } else {
          break;
        }
      }
    }

    if (vertices[b1_].fix_gray != NO_EDGE ||
        vertices[vertices[b1_].red].in_cycle) {
      break;
    }
    Vtx_id r1 = vertices[b1_].red;
    if (indel_count[vertices[r1].gene_val] * vertices[r1].indel_update < 0) {
      red1.push_back(r1);
      b1_ = vertices[r1].black;
      red1.push_back(b1_);
      indel_count[vertices[r1].gene_val] += vertices[r1].indel_update;
    } else {
      break;
    }
  }

  if (cycle.size() > 0) {
    add_cycle(cycle);
  }
}

bool CycleGraph::is_chromosome_end(Vtx_id v) const {
  return vertices[v].vtx_type == Cap ||
         (vertices[v].black < v == vertices[v].red < v);
}

/* Follow indel edges from the last entry of queue and add all vertices reached
 * to the queue */
void CycleGraph::follow_indel_edge(set<Vtx_id> &vizited_in_level,
                                   vector<unique_ptr<QEntry>> &q1) {
  Vtx_id v = vertices[q1.back()->vtx].black;
  while (vertices[v].fix_gray == NO_EDGE && !is_chromosome_end(v) &&
         !vertices[vertices[v].red].in_cycle) {
    Vtx_id u = vertices[v].red;
    if (q1.back()->indel_count[vertices[u].gene_val] *
                vertices[u].indel_update <
            0 &&
        q1.back()->not_yet_visited(vizited_in_level, u)) {
      q1.push_back(make_unique<QEntry>(
          v, u, q1.back()->fixed, q1.back()->vizited, q1.back()->indel_count));
      q1.back()->indel_count[vertices[u].gene_val] += vertices[u].indel_update;
    } else {
      break;
    }
    v = vertices[u].black;
  }
}

void CycleGraph::bfs(Vtx_id start, bool is_random) {
  unique_ptr<QEntry> entry;
  vector<unique_ptr<QEntry>> q1;
  vector<unique_ptr<QEntry>> q2;
  set<Vtx_id> vizited_in_level;
  Vtx_id headtailcorresp_v, headtailcorresp_u;
  size_t pos = 0;

  if (vertices[start].in_cycle)
    return;

  q1.push_back(make_unique<QEntry>(start, this->indel_count));
  vizited_in_level.insert(start);
  follow_indel_edge(vizited_in_level, q1);

  do {
    entry = std::move(q1[pos]);
    pos++;

    Vtx_id v = vertices[entry->vtx].black;

    /* Follow each gray edge */
    vector<Vtx_id> neigs;
    if (vertices[v].fix_gray != NO_EDGE) {
      neigs.push_back(vertices[v].fix_gray);
    } else if (!vertices[v].is_indel) {
      auto vu = entry->fixed.find(v);
      if (vu != entry->fixed.end()) {
        neigs.push_back(vu->second);
      } else {
        for (auto u : vertices[v].grays) {
          if (entry->fixed.find(u) == entry->fixed.end() &&
              vertices[u].fix_gray == NO_EDGE) {
            neigs.push_back(u);
          }
        }
      }
    }

    vector<Vtx_id> notDone;
    for (Vtx_id u : neigs) {
      if (!vertices[u].in_cycle &&
          entry->not_yet_visited(vizited_in_level, u)) {
        notDone.push_back(u);
      }
    }

    // if v doest have a fixed gray edge and is not a cap vertex
    if (vertices[v].fix_gray == NO_EDGE && vertices[v].vtx_type != Cap) {
      for (Vtx_id u : notDone) {
        headtailcorresp_v = vertices[v].red;
        headtailcorresp_u = vertices[u].red;
        q2.push_back(make_unique<QEntry>(v, u, headtailcorresp_v,
                                         headtailcorresp_u, entry->fixed,
                                         entry->vizited, entry->indel_count));
        vizited_in_level.insert(u);
        follow_indel_edge(vizited_in_level, q2);
      }
    } else {
      for (Vtx_id u : notDone) {
        q2.push_back(make_unique<QEntry>(v, u, entry->fixed, entry->vizited,
                                         entry->indel_count));
        vizited_in_level.insert(u);
        follow_indel_edge(vizited_in_level, q2);
      }
    }

    if (pos == q1.size()) {
      q1.clear();
      for (unique_ptr<QEntry> &e : q2) {
        q1.push_back(std::move(e));
      }
      if (is_random) {
        shuffle(q1.begin(), q1.end(), rand_gen);
      }
      q2.clear();
      vizited_in_level.clear();
      pos = 0;
    }

  } while (q1[pos]->vtx != start);

  /* Find balanced cycle in level if it exists */
  vector<Vtx_id> cycle;
  for (; pos < q1.size(); pos++) {
    if (q1[pos]->vtx == start) {
      entry = std::move(q1[pos]);
      cycle.clear();
      Vtx_id u, v;
      int weight = 0;
      v = start;
      do {
        u = vertices[v].black;
        cycle.push_back(v);
        cycle.push_back(u);
        v = entry->fixed[u];
        if (v != vertices[u].red) {
          weight += vertices[v].weight;
        }
      } while (v != start);
      if (weight == 0) {
        break;
      }
    }
  }

  /* Add the balanced cycle, if there is none take any cycle in level */
  add_cycle(cycle);
}

string CycleGraph::show_cycles() const {
  ostringstream ss;

  ss << '[';
  for (size_t i = 0; i < cycles.size(); ++i) {
    Vtx_id v = get<2>(cycles[i]);
    ss << '[' << v;
    v = vertices[v].black;
    ss << "," << v;
    v = vertices[v].is_indel ? vertices[v].red : vertices[v].fix_gray;
    while (v != get<2>(cycles[i])) {
      ss << "," << v;
      v = vertices[v].black;
      ss << "," << v;
      v = vertices[v].is_indel ? vertices[v].red : vertices[v].fix_gray;
    }
    ss << ']';
    if (i != cycles.size() - 1)
      ss << ',';
  }
  ss << ']';

  return ss.str();
}

void CycleGraph::read_cycles(string str) {
  string str_cycle, str_el;
  vector<Vtx_id> cycle;

  /* Remove '[[' and ']]'. */
  str = str.substr(2, str.size() - 4);

  /* Read each cycle. */
  string delimiter = "],[";
  while (not str.empty()) {
    size_t idx = str.find(delimiter);
    if (idx != str.npos) {
      str_cycle = str.substr(0, idx);
      str.erase(0, idx + delimiter.length());
    } else {
      str_cycle = str;
      str = "";
    }

    /* Read each element. */
    cycle.clear();
    stringstream ss(str_cycle);
    while (getline(ss, str_el, ',')) {
      cycle.push_back(stoi(str_el));
    }

    add_cycle(cycle);
  }
}

void CycleGraph::make_assignment() {
  int val = 1;
  for (size_t i = last_origin_vertex + 1; i <= last_target_vertex; ++i) {
    if (vertices[i].is_indel) {
      vertices[i].gene_val = 0;
    } else {
      vertices[i].gene_val = val;
      if (vertices[i].vtx_type == Cap || vertices[i].red < i) {
        val++;
      }
    }
  }
  for (size_t i = 0; i <= last_origin_vertex; ++i) {
    if (vertices[i].is_indel) {
      vertices[i].gene_val = 0;
      vertices[i].sign_positive = true;
    } else {
      vertices[i].gene_val = vertices[vertices[i].fix_gray].gene_val;
      vertices[i].sign_positive = (vertices[vertices[i].fix_gray].sign_positive)
                                      ? vertices[i].sign_positive
                                      : !vertices[i].sign_positive;
    }
  }
  for (size_t i = last_origin_vertex + 1; i <= last_target_vertex; ++i) {
    vertices[i].sign_positive = true;
  }
}

Genea vertex_to_gene(Vertex v) {
  return Genea(((v.sign_positive) ? 1 : -1) * v.gene_val, v.is_indel);
}

MultiGenomePair CycleGraph::get_genomes() const {
  MultiGenomePair data;
  vector<pair<bool, vector<Genea>>> gss;
  vector<vector<IR>> irss;

  for (size_t i = 0; i <= last_target_vertex;) {
    if (vertices[i].vtx_type == Cap) {
      gss.push_back(pair<bool, vector<Genea>>(false, vector<Genea>()));
      irss.push_back(vector<IR>());
      assert(!vertices[i].is_indel);
      gss.back().second.push_back(vertex_to_gene(vertices[i]));
      i++;
      while (!is_chromosome_end(i)) {
        gss.back().second.push_back(vertex_to_gene(vertices[i]));
        irss.back().push_back(vertices[i].weight);
        if (vertices[i].is_indel) {
          while (vertices[i].is_indel) {
            i += 2;
          }
        } else {
          i += 2;
        }
      }
      assert(!vertices[i].is_indel);
      gss.back().second.push_back(vertex_to_gene(vertices[i]));
      irss.back().push_back(vertices[i].weight);
      i++;
    } else {
      gss.push_back(pair<bool, vector<Genea>>(true, vector<Genea>()));
      irss.push_back(vector<IR>());
      while (i <= last_target_vertex && vertices[i + 1].black >= i) {
        gss.back().second.push_back(vertex_to_gene(vertices[i]));
        irss.back().push_back(vertices[i].weight);
        if (vertices[i].is_indel) {
          while (i <= last_target_vertex && vertices[i].is_indel &&
                 vertices[i + 1].black >= i) {
            i += 2;
          }
        } else {
          i += 2;
        }
      }
    }

    if (i == last_origin_vertex + 1) {
      data.g = make_unique<MultiGenome>(gss, irss);
      gss.clear();
      irss.clear();
    }
    if (i == last_target_vertex + 1) {
      data.h = make_unique<MultiGenome>(gss, irss);
    }
  }

  return data;
}

PermsIrs CycleGraph::get_perms() const {
  PermsIrs perms_irs;

  MultiGenomePair data = get_genomes();
  if (data.g->n_chrom() != 1 || data.h->n_chrom() != 1) {
    throw runtime_error("get_perms can only be used if the cycle was produce "
                        "from single chromosomes.");
  }

  perms_irs.s_n = (*data.g)[1].size();
  perms_irs.s_n = (*data.h)[1].size();
  perms_irs.s = (int *)malloc(perms_irs.s_n * sizeof(int));
  perms_irs.s_ir = (int *)malloc((perms_irs.s_n - 1) * sizeof(int));
  perms_irs.p = (int *)malloc(perms_irs.p_n * sizeof(int));
  perms_irs.p_ir = (int *)malloc((perms_irs.p_n - 1) * sizeof(int));

  for (int i = 1; i <= perms_irs.s_n; i++) {
    perms_irs.s[i] = (*data.g)[1][i];
  }
  for (int i = 1; i < perms_irs.s_n; i++) {
    perms_irs.s[i] = (*data.g)[1].get_ir(i);
  }
  for (int i = 1; i <= perms_irs.p_n; i++) {
    perms_irs.p[i] = (*data.h)[1][i];
  }
  for (int i = 1; i < perms_irs.p_n; i++) {
    perms_irs.p[i] = (*data.h)[1].get_ir(i);
  }

  return perms_irs;
}

PermsIrs flip_perms_irs(PermsIrs perms_irs) {
  PermsIrs perms_irs_flip;
  perms_irs_flip.s = perms_irs.p;
  perms_irs_flip.p = perms_irs.s;
  perms_irs_flip.s_ir = perms_irs.p_ir;
  perms_irs_flip.p_ir = perms_irs.s_ir;
  perms_irs_flip.s_n = perms_irs.p_n;
  perms_irs_flip.p_n = perms_irs.s_n;
  return perms_irs_flip;
}

void CycleGraph::serialize(ostream &os) const {
  for (size_t i = 0; i <= last_target_vertex; ++i) {
    os << i << "(" << (vertices[i].sign_positive ? "+" : "-")
       << vertices[i].gene_val << "," << vertices[i].black << ","
       << vertices[i].weight << "," << (vertices[i].in_cycle ? "o" : "_")
       << ",";
    switch (vertices[i].vtx_type) {
    case Head:
      os << "H";
      break;
    case Tail:
      os << "T";
      break;
    case Cap:
      os << "C";
      break;
    }
    os << "|[";
    if (vertices[i].is_indel) {
      os << "!" << vertices[i].red;
    } else if (vertices[i].fix_gray == NO_EDGE) {
      for (auto gray : vertices[i].grays) {
        os << gray << ",";
      }
      os << "!" << vertices[i].red;
    } else {
      os << "*" << vertices[i].fix_gray;
    }
    os << "])"
       << " ";
    if (int(i) == last_origin_vertex)
      os << endl;
  }
}

ostream &operator<<(std::ostream &os, const CycleGraph &cg) {
  cg.serialize(os);
  return os;
}

ostream &operator<<(ostream &os, const Run &r) {
  os << r.genome << r.fst_gene << ":";
  for (Gene a : r.genes_to_add) {
    os << " " << a;
  }
  cout << endl;
  return os;
}

void CycleGraph::rem_cycle(tuple<size_t, int, Vtx_id> c) {
  Vtx_id headtailcorresp_v, headtailcorresp_u;
  Vtx_id i = get<2>(c);

  indel_potation -= cycle_potation(i);
  vector<Vtx_id> cycle = get_cycle(i);

  if (cycle_weight(i) == 0) {
    balanced_cycles--;
    if (cycle.size() == 4) {
      unitary_balanced_cycles--;
    }
  }

  if (is_divergent(i)) {
    divergent_cycles--;
    if (cycle_weight(i) <= 0) {
      divergent_balanced_or_neg_cycles--;
    }
  }

  assert(cycle.size() % 2 == 0);
  for (size_t i = 0; i < cycle.size(); ++i) {
    Vtx_id v = cycle[i];
    assert(vertices[v].in_cycle);
    vertices[v].in_cycle = false;
    if (i % 2 == 1) {
      Vtx_id u = (i == cycle.size() - 1) ? cycle[0] : cycle[i + 1];
      if (vertices[v].is_indel) {
        vertices[v].is_indel = false;
        vertices[u].is_indel = false;
        this->indel_count[vertices[v].gene_val] -= vertices[v].indel_update;
      } else if (vertices[v].vtx_type == Cap) {
        assert(vertices[v].fix_gray != NO_EDGE);
        assert(vertices[u].fix_gray != NO_EDGE);
        vertices[v].fix_gray = NO_EDGE;
        vertices[u].fix_gray = NO_EDGE;
      } else {
        headtailcorresp_v = vertices[v].red;
        headtailcorresp_u = vertices[u].red;
        if (not vertices[headtailcorresp_v].in_cycle &&
            not vertices[headtailcorresp_u].in_cycle) {
          assert(vertices[v].fix_gray != NO_EDGE);
          assert(vertices[u].fix_gray != NO_EDGE);
          assert(vertices[headtailcorresp_v].fix_gray != NO_EDGE);
          assert(vertices[headtailcorresp_u].fix_gray != NO_EDGE);
          vertices[v].fix_gray = NO_EDGE;
          vertices[u].fix_gray = NO_EDGE;
          vertices[headtailcorresp_v].fix_gray = NO_EDGE;
          vertices[headtailcorresp_u].fix_gray = NO_EDGE;
        }
      }
    }
  }
  cycles.erase(find(cycles.begin(), cycles.end(), c));
}

bool CycleGraph::check_cycle(vector<Vtx_id> cycle) {
  for (size_t i = 0; i < cycle.size(); ++i) {
    Vtx_id v = cycle[i];
    if (vertices[v].in_cycle)
      return false;
    if (i % 2 == 1) {
      Vtx_id u = (i == cycle.size() - 1) ? cycle[0] : cycle[i + 1];
      if (vertices[v].is_indel) {
        if (vertices[v].red != u)
          return false;
      } else if (vertices[v].fix_gray != NO_EDGE) {
        if (vertices[v].fix_gray != u)
          return false;
      } else if (vertices[u].fix_gray != NO_EDGE) {
        return false;
      }
    }
  }
  return true;
}

void CycleGraph::add_cycle(vector<Vtx_id> cycle) {
  Vtx_id headtailcorresp_v, headtailcorresp_u;

  assert(cycle.size() % 2 == 0);
  for (size_t i = 0; i < cycle.size(); ++i) {
    Vtx_id v = cycle[i];
    assert(not vertices[v].in_cycle);
    vertices[v].in_cycle = true;
    if (i % 2 == 1 && vertices[v].fix_gray == NO_EDGE) {
      Vtx_id u = (i == cycle.size() - 1) ? cycle[0] : cycle[i + 1];
      if (u == vertices[v].red) {
        vertices[v].is_indel = true;
        vertices[u].is_indel = true;
        this->indel_count[vertices[v].gene_val] += vertices[v].indel_update;
      } else {
        if (vertices[u].fix_gray == NO_EDGE) {
          vertices[v].fix_gray = u;
          vertices[u].fix_gray = v;
          if (vertices[v].vtx_type != Cap) {
            headtailcorresp_v = vertices[v].red;
            headtailcorresp_u = vertices[u].red;
            assert(vertices[headtailcorresp_u].fix_gray == NO_EDGE);
            assert(vertices[headtailcorresp_v].fix_gray == NO_EDGE);
            vertices[headtailcorresp_v].fix_gray = headtailcorresp_u;
            vertices[headtailcorresp_u].fix_gray = headtailcorresp_v;
          }
        } else {
          if (vertices[u].fix_gray != v || vertices[v].fix_gray != u) {
            throw runtime_error("Atempting to add a cycle with gray edge "
                                "distinct from the fixed one.");
          }
        }
      }
    }
  }
  int pot = cycle_potation(cycle[0]);
  int weight = cycle_weight(cycle[0]);
  cycles.push_back(make_cycle(cycle, weight, pot));

  indel_potation += pot;
  if (weight == 0) {
    balanced_cycles++;
    if (cycle.size() == 4) {
      unitary_balanced_cycles++;
    }
  }

  if (is_divergent(cycle[0])) {
    divergent_cycles++;
    if (weight <= 0) {
      divergent_balanced_or_neg_cycles++;
    }
  }
}

void CycleGraph::check_and_add_cycle(vector<Vtx_id> cycle) {
  if (check_cycle(cycle)) {
    add_cycle(cycle);
  }
}

vector<Vtx_id> CycleGraph::get_cycle(Vtx_id i) const {
  vector<Vtx_id> cycle;

  Vtx_id v = i;
  cycle.push_back(v);
  v = vertices[v].black;
  cycle.push_back(v);
  v = vertices[v].is_indel ? vertices[v].red : vertices[v].fix_gray;
  while (v != cycle[0]) {
    cycle.push_back(v);
    v = vertices[v].black;
    cycle.push_back(v);
    v = vertices[v].is_indel ? vertices[v].red : vertices[v].fix_gray;
  }

  return cycle;
}

int CycleGraph::is_divergent(Vtx_id i) const {
  bool ans = false;
  int flips = 0;
  int count = 0;
  bool dir = i < vertices[i].black;

  Vtx_id v = i;
  int end = v;
  do {
    if (dir != (v < vertices[v].black)) {
      flips += 1;
      if (flips == 1) {
        end = v;
      }
      if (flips >= 3) {
        ans = true;
      }
    }
    if (flips == 1) {
      count += 1;
    }
    if (flips == 2) {
      count -= 1;
    }

    v = vertices[v].black;
    if (vertices[v].is_indel) {
      v = vertices[v].red;
    } else {
      v = vertices[v].fix_gray;
    }
  } while (v != end && !ans);

  if (ans || count != 0)
    return true;
  else
    return false;
}

int CycleGraph::cycle_weight(Vtx_id i) const {
  int weight = 0;

  Vtx_id v = i;
  do {
    v = vertices[v].black;
    if (vertices[v].is_indel) {
      v = vertices[v].red;
    } else {
      v = vertices[v].fix_gray;
      weight += vertices[v].weight;
    }
  } while (v != i);

  return weight;
}

int CycleGraph::cycle_potation(int i) const {
  int runs = 0;
  int first_state = 0, last_state = 0, state = 0;

  Vtx_id v = i;
  do {
    v = vertices[v].black;
    if (vertices[v].is_indel) {
      state = vertices[v].indel_update;
      if (last_state != state) {
        if (last_state != 0) {
          runs++;
        } else {
          first_state = state;
        }
        last_state = state;
      }
      v = vertices[v].red;
    } else {
      v = vertices[v].fix_gray;
    }
  } while (v != i);
  if (last_state != 0 && (last_state != first_state || runs == 0)) {
    runs++;
  }

  if (runs == 0) {
    return 0;
  } else {
    return ceil((runs + 1) / 2.0);
  }
}

Run CycleGraph::cycle_run(int i) const {
  Run run;
  int run_state = 0,
      state = 0; // States indicate if a run is a insertion or a deletion run
  Vtx_id v = i;

  // Find one end of the run
  int first_end = -1;
  do {
    v = vertices[v].black;
    if (vertices[v].is_indel) {
      state = vertices[v].indel_update;
      if (run_state != state) {
        if (run_state == 0) {
          run_state = state;
        } else {
          first_end = v;
          break;
        }
      }
      v = vertices[v].red;
    } else {
      v = vertices[v].fix_gray;
    }
  } while (v != i);

  if (run_state == 0)
    return run; // No run on the cycle

  bool revert = false;
  if (first_end == -1) { // single run
    while (vertices[v].indel_update == run_state) {
      v = vertices[v].black;
      if (vertices[v].is_indel) {
        v = vertices[v].red;
      } else {
        v = vertices[v].fix_gray;
      }
    }
    first_end = v;
  }

  // Get insertion position and check if the list of insertions must be inverted
  assert(vertices[first_end].indel_update != run_state);
  if (first_end > vertices[first_end].black) {
    run.fst_gene = vertices[vertices[v].black].gene_val;
  } else {
    run.fst_gene = vertices[v].gene_val;
    revert = true;
  }

  // Get genome where the insertion will occur
  if (run_state == -1) {
    run.genome = 'H';
  } else {
    run.genome = 'G';
  }

  // Get genes to insert
  v = first_end;
  do {
    v = vertices[v].black;
    if (vertices[v].is_indel) {
      state = vertices[v].indel_update;
      if (run_state != state) {
        break;
      }

      if (vertices[v].sign_positive == (v < vertices[v].red)) {
        run.genes_to_add.push_back(vertices[v].gene_val);
      } else {
        run.genes_to_add.push_back(-vertices[v].gene_val);
      }

      v = vertices[v].red;
    } else {
      v = vertices[v].fix_gray;
    }
  } while (v != first_end);

  if (revert) {
    reverse(run.genes_to_add.begin(), run.genes_to_add.end());
    for (int i = 0; i < int(run.genes_to_add.size()); i++) {
      run.genes_to_add[i] = -run.genes_to_add[i];
    }
  }

  return run;
}

int CycleGraph::get_gray_degree(Vtx_id i) const {
  if (vertices[i].fix_gray != NO_EDGE) {
    return 1;
  } else {
    return vertices[i].grays.size();
  }
}

int CycleGraph::get_black_degree(Vtx_id i) const {
  return (vertices[i].black == NO_EDGE) ? 0 : 1;
}

int CycleGraph::get_red_degree(Vtx_id i) const {
  return (vertices[i].red == NO_EDGE) ? 0 : 1;
}

int CycleGraph::dcj_ir_distance() const {
  int num_bal_cycles = this->balanced_cycles;
  int merges = 0;
  int splits = 0;

  vector<int> unbalanced_weights;
  for (auto c : cycles) {
    if (get<1>(c) != 0) {
      unbalanced_weights.push_back(get<1>(c));
    }
  }
  int num_unbal_cycles = unbalanced_weights.size();

  // Merge possible pairs of unbalanced cycles
  for (int i = 0; i < unbalanced_weights.size(); i++) {
    for (int j = i + 1; j < unbalanced_weights.size(); j++) {
      if (unbalanced_weights[i] != 0 && unbalanced_weights[j] != 0 &&
          unbalanced_weights[i] + unbalanced_weights[j] == 0) {
        unbalanced_weights[i] = 0;
        unbalanced_weights[j] = 0;
        merges += 1;
        num_bal_cycles += 1;
        num_unbal_cycles -= 2;
      }
    }
  }

  // 1/3-approximation to merge possible triples of unbalanced cycles
  for (int i = 0; i < unbalanced_weights.size(); i++) {
    for (int j = 0; j < unbalanced_weights.size(); j++) {
      if (j == i)
        continue;
      for (int k = j + 1; k < unbalanced_weights.size(); k++) {
        if (k == i)
          continue;
        if (unbalanced_weights[i] != 0 && unbalanced_weights[j] != 0 &&
            unbalanced_weights[k] != 0 &&
            -unbalanced_weights[i] ==
                unbalanced_weights[j] + unbalanced_weights[k]) {
          unbalanced_weights[i] = 0;
          unbalanced_weights[j] = 0;
          unbalanced_weights[k] = 0;
          merges += 2;
          num_bal_cycles += 1;
          num_unbal_cycles -= 3;
        }
      }
    }
  }

  // All remaining unbalanced cycles are merged into one
  merges += num_unbal_cycles - 1;

  splits = (last_origin_vertex / 2) - num_bal_cycles;
  return merges + splits;
}

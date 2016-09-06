#ifndef HAPLOTYPE_ENUMERATOR_H
#define HAPLOTYPE_ENUMERATOR_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "vg.pb.h"
#include "xg.hpp"

using namespace std;

//  RRMemo functions
//  Created by Jordan Eizenga on 6/21/16.
struct RRMemo {
private:

  std::vector<double> S_multipliers;
  double T_multiplier;

  std::vector< std::vector<double> > S;
  std::vector<double> T;

  double rho;
  double exp_rho;

  double S_value(int height, int width);
  double T_value(int width);

public:
  RRMemo(double recombination_penalty);
  ~RRMemo(void);

  double recombination_penalty();

  double rr_diff(int height, int width);
  double rr_same(int height, int width);
  double rr_adj(int width);
  double rr_all(int height, int width);
};

class rectangle {
private:
  xg::XG::ThreadSearchState state;
  // We don't use these yet (we're using relative indices instead) but they will
  // be used in edit-propsal
  int a_index;
public:
  ~rectangle(void) {};
  // Switched to indices from pointers as quick fix for messy bug caused by
  // reindexing on removal of empty rectangles... TODO? fix this
  int prev = -1;
  int next = -1;
  int J = 0;
  int I = 0;
  double R = 0;
  // Computes J at next_id for the strip corresponding to state
  // NB that this also calls rectangle::extend
  int get_next_J(xg::XG::ThreadMapping next_node, xg::XG& graph);
  // Extends state by node next_id
  void extend(xg::XG::ThreadMapping next_node, xg::XG& graph);
};

// A cross-section is a column of rectangles S^a_b, a <= b. Each "rectangle" in
// the sense of recomb-rectangle functions is a whole cross_section
struct cross_section {
private:
  xg::XG::ThreadMapping node;
  int b_index;
public:
  cross_section(int64_t new_height,int b,xg::XG::ThreadMapping new_node);
  ~cross_section(void) {};
  vector<rectangle> S;
  int height; // height (in consistent thread_ts)
  int width = 1; // width (in base pairs)
  inline xg::XG::ThreadMapping get_node();
};

using thread_t = vector<xg::XG::ThreadMapping>;

// A haplo_d indexes |A| + 1 columns of rectangles S^*_b according in A-order
class haplo_d {
public:
  rectangle empty_rect;
  vector<cross_section> cs;
  haplo_d(const thread_t& t, xg::XG& graph);
  ~haplo_d(void) {};
  // calculate_Is() needs to be called for the cross_sections have I values in
  // their rectangles. The haplo_d constructor only builds the most recent (in
  // terms of node history) 1 or 2 rectangles at each node
  void calculate_Is(xg::XG& graph);
  double probability(double recombination_penalty);
  // Prints csv containing: number of threads joining, number leaving, number of
  // consistent threads, length in base pairs, |A_curr| for each a->a+1 interval
  // in the set A of "active" nodes. Returns |A_curr|^max and total length in
  // base pairs
  pair<int,int> print_decomposition_stats(string output_path);
  // Prints csv containing: J-values, prev indices, next indices--ie smallest
  // human-interpretable representation of the haplo_d
  void print_decomposition(string output_path);
  // We want to make a big array where rows are strips and cells are rectangles'
  // I-values. This is for making pretty pictures of haplo_ds

  // This array will be relatively sparse since it will hold all-zero entries
  // prior to strips joining our query haplotype and all-one entries after. This
  // is for output to R or a similar program for visualization; haplo_d entities
  // as implemented are more compact but not easily human-readable
  void unfold_rectangles(string output_path);
  inline double prev_R(int b, int a);
  inline int prev_I(int b, int a);
};

thread_t path_to_thread_t(vg::Path& path);

// A function to take an xg index and extract all of its embedded threads into
// haplo_ds. This is SLOW to run to completion since there are in general tons
// of embededded threads in a haplo_d
void extract_threads_into_haplo_ds(xg::XG& index, string output_path, int64_t start_node, int64_t end_node, bool make_graph);
void decompose_and_print(const thread_t& t, xg::XG& graph, string output_path);
bool check_for_edges(int64_t old_node_id, bool old_node_is_reverse, int64_t new_node_id, bool new_node_is_reverse, xg::XG& index);
bool check_if_thread_t_broken(const thread_t& t, XG& graph);
#endif

#ifndef HAPLOTYPE_EXTRACTER_H
#define HAPLOTYPE_EXTRACTER_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "vg.pb.h"
#include "xg.hpp"
#include "haplotypes.hpp"

using namespace std;
using namespace vg;
using thread_t = vector<xg::XG::ThreadMapping>;

void output_weighted_haplotype_list(string output_path, vector<pair<thread_t,int>>& haplotype_list, xg::XG& index);
void thread_to_graph_spanned(thread_t& t, Graph& graph, xg::XG& index);
Path path_from_thread_t(thread_t& t);
vector<pair<thread_t,int> > list_haplotypes(xg::XG& index, xg::XG::ThreadMapping start_node, int extend_distance);

#endif

#include <iostream>
#include "haplotypes.hpp"
#include "path.hpp"
#include "position.hpp"

using namespace std;
using namespace xg;

RRMemo::RRMemo(double recombination_penalty)  {
  rho = recombination_penalty;
  exp_rho = exp(-rho);
  assert(exp_rho < 1);
  continue_probability = 1 - population_size * exp_rho;
  exp_rho = exp_rho/continue_probability;
  rho = -log(exp_rho);
  S.push_back(std::vector<double>(1, 1.0));
  S_multipliers.push_back(1.0);
  T.push_back(1.0);
  T_multiplier = 1.0 - exp_rho;

  // log versions
  logT_base = log1p(-exp_rho);
  for(int i = 0; i < population_size; i++) {
    logS_bases.push_back(log1p(i*exp_rho));
  }
}

RRMemo::~RRMemo(void) {

}

double logdiff(double a, double b) {
  if(b > a) {
    double c = a;
    a = b;
    b = c;
  }
  return a + log1p(-exp(b - a));
}

double logsum(double a, double b) {
  if(b > a) {
    double c = a;
    a = b;
    b = c;
  }
  return a + log1p(exp(b - a));
}

double RRMemo::logSN(vector<double> logRs, vector<int> Is) {
  if(logRs.size() == 1) {
    return logRs[0] + log(Is[0]);
  } else {
    double max_summand = logRs[0] + log(Is[0]);
    int max_index = 0;
    vector<double> summands;
    for(int i = 0; i < logRs.size(); i++){
      summands.push_back(logRs[i] + log(Is[i]));
      if(summands.back() > max_summand) {
        max_summand = summands.back();
        max_index = i;
      }
    }
    double sum = 0;
    for(int i = 0; i < summands.size(); i++) {
      if(i != max_index) {
        sum += exp(summands[i]-max_summand);
      }
    }
    return max_summand + log1p(sum);
  }
}

double RRMemo::logT(int width) {
  return (width-1)*logT_base; //logT_base = log(1 - exp_rho)
}

double RRMemo::logS(int height, int width) {
  return (width-1)*logS_bases[height-1]; //logS_base = log(1 + i*exp_rho)
}

double RRMemo::logRRDiff(int height, int width) {

  return logdiff(logS(height,width),logT(width)) - log(height);
}

double RRMemo::log_continue_factor(int64_t totwidth) {
  return totwidth * log(continue_probability);
}

double RRMemo::continue_factor(int64_t totwidth) {
  return exp(log_continue_factor(totwidth));
}

double RRMemo::recombination_penalty() {
  return exp_rho;
}

double RRMemo::log_recombination_penalty() {
  return -rho;
}

double RRMemo::cont_probability() {
  return continue_probability;
}

double RRMemo::S_value(int height, int width) {

  while (S.size() < height) {
    S_multipliers.push_back(S_multipliers[S_multipliers.size() - 1] + exp_rho);
    S.push_back(std::vector<double>(1, 1.0));
  }
  std::vector<double>& S_row = S[height - 1];
  double S_multiplier = S_multipliers[height - 1];

  while (S_row.size() < width) {
    S_row.push_back(S_row[S_row.size() - 1] * S_multiplier);
  }

  return S_row[width - 1];
}

double RRMemo::T_value(int width) {

  while (T.size() < width) {
    T.push_back(T[T.size() - 1] * T_multiplier);
  }

  return T[width - 1];
}

double RRMemo::rr_diff(int height, int width) {

  if (height < 1 || width < 1) {
    cerr << "error:[RRMemo] height and width of recombination rectangle must be >= 1" << endl;
  }

  return (S_value(height, width) - T_value(width)) / height;
}

double RRMemo::rr_same(int height, int width) {

  if (height < 1 || width < 1) {
    cerr << "error:[RRMemo] height and width of recombination rectangle must be >= 1" << endl;
  }

  double T_val = T_value(width);
  return (S_value(height, width) - T_val) / height + T_val;
}

double RRMemo::rr_adj(int width) {

  if (width < 1) {
    cerr << "error:[RRMemo] height and width of recombination rectangle must be >= 1" << endl;
  }

  return T_value(width);
}

double RRMemo::rr_all(int height, int width) {

  if (height < 1 || width < 1) {
    cerr << "error:[RRMemo] height and width of recombination rectangle must be >= 1" << endl;
  }

  return S_value(height, width);
}

// Unmemoized implementations of same polynomials

double rr_diff(int height, int width, double recombination_penalty) {
  double exp_rho = exp(-recombination_penalty);
  return (pow(1.0 + (height - 1.0) * exp_rho, width - 1.0) - pow(1.0 - exp_rho, width - 1.0)) / height;
}

double rr_same(int height, int width, double recombination_penalty) {
  double exp_rho = exp(-recombination_penalty);
  double T_val = pow(1.0 - exp_rho, width - 1.0);
  return (pow(1.0 + (height - 1.0) * exp_rho, width - 1.0) - T_val) / height + T_val;
}

double rr_adj(int width, double recombination_penalty) {
  return pow(1.0 - exp(-recombination_penalty), width - 1.0);
}

double rr_all(int height, int width, double recombination_penalty) {
  double exp_rho = exp(-recombination_penalty);
  return exp_rho * pow(1.0 + (height - 1.0) * exp_rho, width - 1.0);
}

// Rectangular decomposition building:

cross_section::cross_section(int64_t new_height, int b, XG::ThreadMapping new_node) {
  b_index = b;
  height = new_height;
  node = new_node;
}

void rectangle::extend(XG::ThreadMapping next_node, XG& graph) {
  int64_t next_side = graph.id_to_rank(next_node.node_id) * 2 + next_node.is_reverse;
  if(state.current_side == 0) {
    // We're extending an empty state
    state.range_start = 0;
    state.range_end = graph.node_height(next_node);
  } else {
    // bool edge_exists = check_for_edges(graph.rank_to_id(state.current_side / 2),state.current_side % 2,
    //        next_node.node_id, next_node.is_reverse, graph);
    // Else, look at where the path goes to and apply the where_to function to
    // shrink the range down.
    //if(edge_exists) {
      state.range_start = graph.where_to(state.current_side, state.range_start, next_side);
      state.range_end = graph.where_to(state.current_side, state.range_end, next_side);
    //} else {
    //  state.range_end = state.range_start;
    //}
  }
  state.current_side = next_side;
}

void rectangle::simple_extend(thread_t& extension, xg::XG& graph) {
  xg::XG::ThreadMapping next_node = extension.back();
  int64_t next_side = graph.id_to_rank(next_node.node_id) * 2 + next_node.is_reverse;
  state.current_side = next_side;
}

void rectangle::simple_extend(xg::XG::ThreadMapping next_node, xg::XG& graph) {
  int64_t next_side = graph.id_to_rank(next_node.node_id) * 2 + next_node.is_reverse;
  state.current_side = next_side;
}

int rectangle::get_next_J(XG::ThreadMapping next_node, XG& graph) {
  extend(next_node, graph);
  return state.count();
}

int rectangle::get_next_J(thread_t& extension, XG& graph) {
  if(extension.size() == 1) {
    return get_next_J(extension.back(), graph);
  } else {
    xg::XG::ThreadMapping second_last_node = extension.end()[-2];
    state.current_side = graph.id_to_rank(second_last_node.node_id) * 2 + second_last_node.is_reverse;
    extend(extension.back(), graph);
    return state.count();
  }
}

inline XG::ThreadMapping cross_section::get_node() {
  return node;
}

bool check_for_edges(int64_t old_node_id, bool old_node_is_reverse, int64_t new_node_id,
          bool new_node_is_reverse, xg::XG& index) {
  // What edge are we following
  Edge edge_taken = make_edge(old_node_id, old_node_is_reverse, new_node_id, new_node_is_reverse);

  // Make sure we find it
  bool edge_found = false;

  vector<Edge> edges = new_node_is_reverse ? index.edges_on_end(new_node_id) : index.edges_on_start(new_node_id);

  for(auto& edge : edges) {
    // Look at every edge in order.
    if(edges_equivalent(edge, edge_taken)) {
      // If we found the edge we're taking, break.
      edge_found = true;
      break;
    }
  }
  if(edge_found == false) {  cerr << "did not find edge between" << old_node_id << " and " << new_node_id << endl;}
  return edge_found;
}

bool check_if_thread_t_broken(const thread_t& t, XG& graph) {
  bool broken = false;
  XG::ThreadMapping current_node = t[0];
  for(int i = 1; i < t.size(); i++) {
    XG::ThreadMapping next_node = t[i];
    bool edge_exists = check_for_edges(current_node.node_id, current_node.is_reverse,
            next_node.node_id, next_node.is_reverse, graph);
    if(!edge_exists) {
      broken = true;
      break;
    }
    current_node = next_node;
  }
  return broken;
}

pair<vector<pair<double,int>>,vector<pair<double,int>>> haplo_d::identify_structures(double leave_threshold,
        double return_threshold, int timeout, xg::XG& graph) {
  vector<pair<double,int>> recombinations;
  vector<pair<double,int>> private_polys;
  // timer ensures that you don't identify the "same" complex variant twice
  int timer = 0;
  for(int b = 1; b < cs.size()-1; b++) {
    if(timer < timeout) {
      timer += cs[b].width;
    } else {
      // does a lot of thread leave?
      int leavers = cs[b-1].S[0].J - cs[b].S[cs[b-1].S[0].next].J;
      double leave_ratio = leavers / cs[b-1].height;
      int starters = 0;
      if(cs[b].S[0].prev == -1) {
        starters = cs[b].S[0].I;
      }
      double return_ratio = 0;
      if(leavers > 0) {
        return_ratio = starters / leavers;
      }
      if(leave_ratio > leave_threshold) {
        if(return_ratio > return_threshold) {
          // returning "defines" recombination structurally--therefore we count it instead of leaving
          recombinations.push_back(make_pair(return_ratio,leavers));
        } else {
          private_polys.push_back(make_pair(leave_ratio,leavers));
        }
        // reset the timer
        timer = cs[b].width;
      }
    }
  }
  return make_pair(recombinations,private_polys);
}

haplo_d::haplo_d() {

}

haplo_d::haplo_d(const thread_t& t, XG& graph) {
  rectangle rect;
  rect.J = rect.get_next_J(t[0],graph);
  // At the leftmost node there is only one strip, so I = J
  rect.I = rect.J;
  int last_height = rect.J;
  cs.push_back(cross_section(rect.J,0,t[0]));
  cs.back().S.push_back(rect);
  int width = 0;
  int new_height;
  bool add_rectangle;
  bool add_A;
  for(int i = 1; i < t.size(); i++) {
    // Count the number of base pairs since the last entry or exit node
    width += graph.node_length(t[i-1].node_id);
    new_height = graph.node_height(t[i]);
    if(cs.back().S.size() != 0) {
      rect = cs.back().S[0];
      rect.J = rect.get_next_J(t[i],graph); // step this strip forward
      // Did any threads leave?
      if(last_height > rect.J) {
        add_A = 1;
      }
      // Are there any threads here which didn't come from the previous node?
      if(rect.J < new_height) {
        add_rectangle = 1;
        add_A = 1;
      }
      // This is an entry or exit node, add a cross-section to the vector of
      // cross-sections (which corresponds to the "A" set in the theory doc)
      if(add_A) {
        cs.back().width = width;
        width = 0;
        cs.push_back(cross_section(new_height,i,t[i]));
      } else {
        // This isn't a node where anything leaves or joins, let's skip over it
        cs.back().bridge.push_back(t[i]);
        for (size_t a = 0; a < cs.back().S.size(); a++) {
          cs.back().S[a].extend(t[i],graph);
        }
      }
      // This is an entry node; we also need a new rectangle corresponding to the
      // new strip. We need to do this *before* we populate since cross_sections
      // arrange rectangles newest -> oldest
      // NB that add_rectangle implies add_A
      if(add_rectangle) {
        rectangle new_rect;
        new_rect.extend(t[i],graph);
        new_rect.J = new_height;
        cs.back().height = new_rect.J;
        cs.back().S.push_back(new_rect);
        cs.back().S.back().I = new_rect.J - rect.J;
      }
      if(add_A) {
        int b = cs.size()-1;
        if(rect.J > 0) {
          cs[b].S.push_back(rect);
          cs[b].S.back().prev = 0;
          cs[b-1].S[0].next = cs[b].S.size()-1;
        }
      }
      last_height = new_height;
      add_A = 0;
      add_rectangle = 0;
    } else {
      cs.back().width = width;
      width = 0;
      cs.push_back(cross_section(new_height,i,t[i]));
      if(new_height > 0) {
        rectangle new_rect;
        new_rect.extend(t[i],graph);
        new_rect.J = new_height;
        cs.back().height = new_rect.J;
        cs.back().S.push_back(new_rect);
        cs.back().S.back().I = new_rect.J - rect.J;
      }
    }
  }
  if(cs.size() == 1) {
    cs.back().width = width;
  }
  cs.back().width += graph.node_length(t.back().node_id) - 1;
  for(int i = 0; i < cs.size(); i++) {
    tot_width += cs[i].width;
  }
}

void haplo_d::calculate_Is(XG& graph) {
  // node 0 was done in the haplo_d constructor; start at node 1
  for(int b = 1; b < cs.size(); b++) {
    // make sure that there is at least one rectangle here
    if(cs[b].S.size() != 0) {
      // get side and orientation of the next element in our query thread_t
      XG::ThreadMapping next_node = cs[b].get_node();
      // if J = 0 for a rectangle, then J must be 0 for all older rectangles
      bool nonempty_J = (cs[b].S.back().J > 0);
      if(nonempty_J) {
        int new_J;
        // start at a = 1 since the haplo_d initializer handles a = 0
        for(int a = 1; a < cs[b-1].S.size(); a++) {
          rectangle new_rect = cs[b-1].S[a];
          if (cs[b-1].bridge.size() != 0 &&
                new_rect.state.current_side != graph.id_to_rank(cs[b-1].bridge.back().node_id) * 2 + cs[b-1].bridge.back().is_reverse) {
            for(int i = 0; i < cs[b-1].bridge.size(); i++) {
              new_rect.extend(cs[b-1].bridge[i],graph);
            }
          }
          new_J = new_rect.get_next_J(next_node,graph);
          new_rect.J = new_J;
          if(new_J != 0) {
            if(new_J == cs[b].S.back().J) {
              // There was nothing in the last rectangle, so let's remove it
              cs[b].S.pop_back();
            } else {
              // There were threads in the last rectangle, so let's count them
              cs[b].S.back().I = cs[b].S.back().J - new_J;
              cs[b-1].S[a].next = cs[b].S.size();
            }
            // Let's assume for now that there is something in this rectangle;
            // if it turns out at the next step that we're wrong, then we'll
            // remove it
            cs[b].S.push_back(new_rect);
            cs[b].S.back().prev = a;
          } else {
            // don't try to add more rectangles, we know that J = 0 from here on
            break;
          }
        }
      } else {
        // this shouldn't be here
        cs[b].S.pop_back();
      }
      cs[b].S.back().I = cs[b].S.back().J;
    }
  }
}

// Querying haplo_d structure:

thread_t path_to_thread_t(Path& path) {
  thread_t t;
  for(size_t i = 0; i < path.mapping_size(); i++) {
    Mapping mapping = path.mapping(i);
    auto pos = mapping.position();
    XG::ThreadMapping m = {pos.node_id(), pos.is_reverse()};
    t.push_back(m);
  }
  return t;
}

void decompose_and_print(const thread_t& t, XG& graph, string output_path) {
  haplo_d decomposition = haplo_d(t, graph);
  decomposition.calculate_Is(graph);
  decomposition.print_decomposition_stats(output_path);
}

void haplo_d::unfold_rectangles(string output_path) {
  ofstream haplo_d_out (output_path);
  int norm_A = cs.size();
  // Scan through the haplo_d and record where initial rectangles of strips lie
  vector<pair<int,int>> start_nodes;
  for(int i = 0; i < cs.size(); i++) {
    for(int j = 0; j < cs[i].S.size(); j++) {
      if(cs[i].S[j].prev == -1) {
        start_nodes.push_back(make_pair(i,j));
        break; // You can't have more than one "starting" rectangle
      }
    }
  }
  cerr << start_nodes.size() << " start nodes" << endl;
  // Iterate through all initial nodes and build the strips in a row of the
  // "rectangles" array
  sort(start_nodes.begin(),start_nodes.end());
  // Need to use a vector of vectors then expand it at the time of printing
  // since you'd be trying to allocate a huge 2D array otherwise. Be careful of
  // indexing!!
  vector<vector<int>> rectangles(start_nodes.size(),vector<int>());
  for(int i = 0; i < start_nodes.size(); i++) {
    int current_rect = start_nodes[i].second;
    int j = 0;
    while(current_rect != -1) {
      rectangles[i].push_back(cs[start_nodes[i].first + j].S[current_rect].I);
      current_rect = cs[start_nodes[i].first + j].S[current_rect].next;
      j++;
    }
  }
  // Print the rectangles array as tab-delimited plaintext
  for(int i = 0; i < start_nodes.size(); i++) {
    for(int j = 0; j < start_nodes[i].first; j++) {
      haplo_d_out << 0 << "\t";
    }
    for(int j = 0; j < rectangles[i].size(); j++) {
      haplo_d_out << rectangles[i][j] << "\t";
    }
    for(int j = start_nodes[i].first + rectangles[i].size(); j < norm_A; j++) {
      haplo_d_out << 0 << "\t";
    }
    haplo_d_out << endl;
  }
  haplo_d_out.close();
}

void haplo_d::print_decomposition(string output_path) {
  ofstream haplo_d_out (output_path);
  for(int i = 0; i < cs.size(); i++) {
    haplo_d_out << cs[i].get_node().node_id<< ", J's: \t";
    for(int j = 0; j < cs[i].S.size(); j++) {
      haplo_d_out << cs[i].S[j].J << "\t\t";
    }
    haplo_d_out << endl;
    haplo_d_out << cs[i].get_node().node_id<< ", prev: \t";
    for(int j = 0; j < cs[i].S.size(); j++) {
      haplo_d_out << cs[i].S[j].prev << "\t";
    }
    haplo_d_out << endl;
    haplo_d_out << cs[i].get_node().node_id<< ", next: \t";
    for(int j = 0; j < cs[i].S.size(); j++) {
      haplo_d_out << cs[i].S[j].next << "\t";
    }
    haplo_d_out << endl;
  }
  haplo_d_out.close();
}

pair<int,int> haplo_d::print_decomposition_stats(string output_path) {
  ofstream haplo_d_out (output_path);
  int Acurrmax = 0;
  int tot_length = cs[0].width;
  haplo_d_out << cs[0].S.size() << "\t" << cs[0].height << "\t" << cs[0].height << "\t" << 0 << "\t"
        << cs[0].width << "\t" << (cs[0].S[0].prev == -1) << "\n";
  for(int i = 1; i < cs.size(); i++) {
    int joiners = (cs[i].S[0].prev == -1) ? cs[i].S[0].I : 0;
    haplo_d_out << cs[i].S.size() << "\t" << cs[i].height << "\t" << joiners << "\t"
          << cs[i-1].height + joiners - cs[i].height << "\t" << cs[i].width << "\t"
          << (cs[i].S[0].prev == -1) << "\n";
    if(cs[i].S.size() > Acurrmax) {Acurrmax = cs[i].S.size();}
    tot_length += cs[i].width;
  }
  cerr << ": |A_curr|^max = " << Acurrmax << "; length = " << cs.size() << " nodes, " << tot_length << " bp"<< endl;
  haplo_d_out.close();
  pair<int,int> return_val;
  return_val.first = Acurrmax;
  return_val.second = tot_length;
  return return_val;
}

// Calculating probabilitiees

inline double haplo_d::prev_R(int b, int a) {
  if(cs[b].S[a].prev == -1) {
    return 0;
  } else {
    return cs[b-1].S[cs[b].S[a].prev].R;
  }
}

inline double haplo_d::prev_logR(int b, int a) {
  return cs[b-1].S[cs[b].S[a].prev].logR;
}

inline int haplo_d::prev_I(int b, int a) {
  if(cs[b].S[a].prev == -1) {
    return 0;
  } else {
    return cs[b-1].S[cs[b].S[a].prev].I;
  }
}

vector<double> haplo_d::prev_logRs(int b) {
  vector<double> returnRs;
  for(int i = 0; i < cs[b].S.size(); i++) {
    if(cs[b].S[i].prev != -1) {
      returnRs.push_back(prev_logR(b,i));
    }
  }
  return returnRs;
}

vector<int> haplo_d::prev_Is(int b) {
  vector<int> returnIs;
  for(int i = 0; i < cs[b].S.size(); i++) {
    if(cs[b].S[i].prev != -1) {
      returnIs.push_back(prev_I(b,i));
    }
  }
  return returnIs;
}

vector<double> haplo_d::current_logRs(int b) {
  vector<double> returnRs;
  for(int i = 0; i < cs[b].S.size(); i++) {
    returnRs.push_back(cs[b].S[i].logR);
  }
  return returnRs;
}

vector<int> haplo_d::current_Is(int b) {
  vector<int> returnIs;
  for(int i = 0; i < cs[b].S.size(); i++) {
    returnIs.push_back(cs[b].S[i].I);
  }
  return returnIs;
}


double haplo_d::probability(RRMemo& memo) {
  // Make sure that we have nonempty rectangles through which to propagate our
  // calculation
  for (int b = 0; b < cs.size(); b++) {
    if(cs[b].height == 0) {
      return 0;
    }
  }
  // defined same as in writeup
  double S1 = 0;
  double S1S2 = 0;
  // compute R for the first interval (which has no predecessor)
  // NB: we are always working at the left edge of a cross_section
  cerr << "starting probability calculation. Tot length " << tot_width << ", |A| "<< cs.size() << endl;
  cs[0].S[0].R = memo.rr_all(cs[0].height,cs[0].width);
  for (int b = 1; b < cs.size(); b++) {
    S1 = 0;
    S1S2 = 0;
    if(cs[b].width != 1) {
      for(int a = 0; a < cs[b].S.size(); a++) {
        // N.B. that R's are r^a_b's rather that R^a_b's. Thus the I factor
        S1 += (prev_R(b,a)) * (prev_I(b,a));
      }
    }
    for(int a = 0; a < cs[b-1].S.size(); a++) {
      S1S2 += cs[b-1].S[a].R * cs[b-1].S[a].I;
    }
    // calculate contributions from all continuing strips
    for(int a = 0; a < cs[b].S.size(); a++) {
      if(cs[b].S[a].prev != -1) {
        cs[b].S[a].R =
        (1 - memo.recombination_penalty()) * ((S1 * memo.rr_diff(cs[b].height, cs[b].width)) + (prev_R(b,a) * memo.rr_adj(cs[b].width)) +
        (memo.recombination_penalty() * S1S2 * memo.rr_all(cs[b].height,cs[b].width)));
        cs[b].S[a].R = cs[b].S[a].R * memo.continue_factor(cs[b].width);
      } else {
        cs[b].S[a].R = (memo.recombination_penalty() * S1S2 * memo.rr_all(cs[b].height,cs[b].width));
        cs[b].S[a].R = cs[b].S[a].R * memo.continue_factor(cs[b].width);
      }
    }
  }
  double total_probability_haplotype = 0;
  for(int a = 0; a < cs.back().S.size(); a++) {
    total_probability_haplotype += cs.back().S[a].R * cs.back().S[a].I;
  }
  return total_probability_haplotype;
}

double haplo_d::log_probability(RRMemo& memo) {
  // Make sure that we have nonempty rectangles through which to propagate our
  // calculation
  for (int b = 0; b < cs.size(); b++) {
    if(cs[b].height == 0 || cs[b].S.size() == 0) {
      return nan("");
    }
  }
  // defined same as in writeup
  double logS1 = 0;
  double logS1S2 = 0;
  double logpS1S2RRS = 0;
  double logS1RRD = 0;
  double logLHS = 0;
  // compute R for the first interval (which has no predecessor)
  // we are always working at the left edge of a cross_section
  cerr << "starting log-probability calculation. Tot length " << tot_width << ", |A| " << cs.size() << endl;
  cs[0].S[0].logR = memo.logS(cs[0].height,cs[0].width);
  for (int b = 1; b < cs.size(); b++) {
    logS1S2 = memo.logSN(current_logRs(b-1),current_Is(b-1));
    logpS1S2RRS = logS1S2 + memo.log_recombination_penalty() + memo.logS(cs[b].height,cs[b].width);
    if(prev_logRs(b).size() == 0) {
      return nan("");
    }
    logS1 = memo.logSN(prev_logRs(b),prev_Is(b));
    logS1RRD = logS1 + memo.logRRDiff(cs[b].height,cs[b].width);
    // calculate contributions from all continuing strips
    for(int a = 0; a < cs[b].S.size(); a++) {
      if(cs[b].S[a].prev != -1) {
        if(cs[b].width == 1) {
          logLHS = memo.logT_base + prev_logR(b,a) + memo.logT(cs[b].width);
        } else {
          logLHS = memo.logT_base + logsum(logS1RRD,prev_logR(b,a) + memo.logT(cs[b].width));
        }
        cs[b].S[a].logR = logsum(logLHS,logpS1S2RRS);
      } else {
        cs[b].S[a].logR = logpS1S2RRS;
      }
    }
  }
  double total_probability_haplotype = memo.logSN(current_logRs(cs.size()-1),current_Is(cs.size()-1));
  return total_probability_haplotype;
}

void logRR_tests(double recombination_penalty) {
  RRMemo memo = RRMemo(recombination_penalty);
  vector<int> Is = {4000,1000,500,10,5};
  vector<int> lengths = {1000,20,5,3,1};
  for(int i = 0; i < lengths.size(); i++) {

    cerr << "logsum(" << Is[i] << "," << lengths[i] << ") = " << logsum(log(Is[i]),log(lengths[i])) << " | ";
    cerr << "log(" << Is[i] << "+" << lengths[i] << ") = " << log(Is[i] + lengths[i]) << endl;
  }
  for(int i = 0; i < lengths.size(); i++) {
    cerr << "logdiff(" << Is[i] << "," << lengths[i] << ") = " << logdiff(log(Is[i]),log(lengths[i])) << " | ";
    cerr << "log(" << Is[i] << "-" << lengths[i] << ") = " << log(Is[i] - lengths[i]) << endl;
  }

  for(int i = 0; i < lengths.size(); i++) {
    for(int j = 0; j < Is.size(); j++) {
      cerr << "logDiff(" << Is[j] <<","<< lengths[i] << ") = " << memo.logRRDiff(Is[j],lengths[i]) << " | ";
      cerr << "log[Diff(" << Is[j] <<","<< lengths[i] << ")] = " << log(memo.rr_diff(Is[j],lengths[i])) << endl;
    }
  }
  vector<double> Rs = {0.00001,0.000002,0.000000003,0.0000000000004,0.000000000000000000007};
  vector<double> logRs;
  for(int i = 0; i < Rs.size(); i++) {
    logRs.push_back(log(Rs[i]));
    cerr << "exp(logR) = " << logRs[i] << " | " << "R = " << Rs[i] << endl;
  }
  cerr << "logSN = " << memo.logSN(logRs,Is) << " | " << "log[SN] = ";
  double result = 0;
  double intermediate = 0;
  for(int i = 0; i < Rs.size(); i++) {
    intermediate = Rs[i]*Is[i];
    result += intermediate;
  }
  cerr << log(result) << endl;
  cerr << "------------" << endl;
}

bool RR_tests(void) {
  // RRMemo tests
  double epsilon = 0.0000001;

  double memo_val;
  double direct_val;

  for (double rho = 1.0; rho < 5.0; rho += 1.0) {

    RRMemo memo(rho);

    for (int c = 1; c < 10; c++) {
      for (int n = 1; n < 10; n++) {

        memo_val = memo.rr_diff(n, c);
        direct_val = rr_diff(n, c, rho);

        if (fabs(memo_val - direct_val) > epsilon) {
          cerr << "FAIL: rr_diff, n = " << n << ", c = " << c << ", rho = " << rho
          << ", memo value = " << memo_val << ", direct value = " << direct_val << endl;
          exit(1);
        }

        memo_val = memo.rr_same(n, c);
        direct_val = rr_same(n, c, rho);

        if (fabs(memo_val - direct_val) > epsilon) {
          cerr << "FAIL: rr_same, n = " << n << ", c = " << c << ", rho = " << rho
          << ", memo value = " << memo_val << ", direct value = " << direct_val << endl;
          exit(1);
        }

        memo_val = memo.rr_all(n, c);
        direct_val = rr_all(n, c, rho);

        if (fabs(memo_val - direct_val) > epsilon) {
          cerr << "FAIL: rr_all, n = " << n << ", c = " << c << ", rho = " << rho
          << ", memo value = " << memo_val << ", direct value = " << direct_val << endl;
          exit(1);
        }
      }

      memo_val = memo.rr_adj(c);
      direct_val = rr_adj(c, rho);

      if (fabs(memo_val - direct_val) > epsilon) {
        cerr << "FAIL: rr_adj, c = " << c << ", rho = " << rho
        << ", memo value = " << memo_val << ", direct value = " << direct_val << endl;
        exit(1);
      }
    }
  }

  cerr << "RR tests passed!" << endl;
  return true;
}


void probabilities_of_all_theads_in_index(xg::XG& index, int64_t start_node, int64_t end_node, int64_t internal_index, double recombination_penalty) {
  RRMemo memo(recombination_penalty);
  end_node = (end_node == -1) ? index.ts_iv.size() : end_node + 1;
  for(int64_t i = 1; i < end_node; i++) {
    // Skip it if no threads start at it
    if(index.ts_iv[i] == 0) {
        continue;
    }
    // If there are threads starting,
    internal_index = internal_index % (index.ts_iv[i] - 1);
    for(int64_t j = internal_index; j < index.ts_iv[i]; j++) {
      // For every thread starting there
      thread_t path;
      int64_t side = i;
      int64_t offset = j;
      cerr << "computing probability for thread starting at \t" << i << " " << j << endl;
      while(true) {
          // Unpack the side into a node traversal
          xg::XG::ThreadMapping m = {index.rank_to_id(side / 2), (bool) (side % 2)};
          // Add the mapping to the thread
          path.push_back(m);
          // Work out where we go:
          // What edge of the available edges do we take?
          int64_t edge_index = index.bs_get(side, offset);
          // If we find a separator, we're very broken.
          assert(edge_index != index.BS_SEPARATOR);
          // Does the path end?
          if(edge_index == index.BS_NULL) {
              // Path ends here.
              break;
          } else {
              // If we've got an edge, convert to an actual edge index
              edge_index -= 2;
          }
          // We also should not have negative edges.
          assert(edge_index >= 0);
          // Look at the edges we could have taken next
          vector<Edge> edges_out = side % 2 ? index.edges_on_start(index.rank_to_id(side / 2)) :
                index.edges_on_end(index.rank_to_id(side / 2));
          assert(edge_index < edges_out.size());
          Edge& taken = edges_out[edge_index];
          // Follow the edge
          int64_t other_node = taken.from() == index.rank_to_id(side / 2) ? taken.to() :
                taken.from();
          bool other_orientation = (side % 2) != taken.from_start() != taken.to_end();
          // Get the side
          int64_t other_side = index.id_to_rank(other_node) * 2 + other_orientation;
          // Go there with where_to
          offset = index.where_to(side, offset, other_side);
          side = other_side;
          // Continue the process from this new side
      }
    // We have a thread to follow, take it
    haplo_d h = haplo_d(path, index);
    h.calculate_Is(index);
    cout << h.log_probability(memo) << endl;
    }
  }
}

void report_threads_and_breaks(xg::XG& index, int64_t start_node, int64_t end_node) {
  end_node = (end_node == -1) ? index.ts_iv.size() : end_node + 1;
  for(int64_t i = 1; i < end_node; i++) {
    // Skip it if no threads start at it
    if(index.ts_iv[i] == 0) {
        continue;
    }
    // If there are threads starting,
    for(int64_t j = 0; j < index.ts_iv[i]; j++) {
      // For every thread starting there
      thread_t path;
      int64_t side = i;
      int64_t offset = j;
      while(true) {
          // Unpack the side into a node traversal
          xg::XG::ThreadMapping m = {index.rank_to_id(side / 2), (bool) (side % 2)};
          // Add the mapping to the thread
          path.push_back(m);
          // Work out where we go:
          // What edge of the available edges do we take?
          int64_t edge_index = index.bs_get(side, offset);
          // If we find a separator, we're very broken.
          assert(edge_index != index.BS_SEPARATOR);
          // Does the path end?
          if(edge_index == index.BS_NULL) {
              // Path ends here.
              break;
          } else {
              // If we've got an edge, convert to an actual edge index
              edge_index -= 2;
          }
          // We also should not have negative edges.
          assert(edge_index >= 0);
          // Look at the edges we could have taken next
          vector<Edge> edges_out = side % 2 ? index.edges_on_start(index.rank_to_id(side / 2)) :
                index.edges_on_end(index.rank_to_id(side / 2));
          assert(edge_index < edges_out.size());
          Edge& taken = edges_out[edge_index];
          // Follow the edge
          int64_t other_node = taken.from() == index.rank_to_id(side / 2) ? taken.to() :
                taken.from();
          bool other_orientation = (side % 2) != taken.from_start() != taken.to_end();
          // Get the side
          int64_t other_side = index.id_to_rank(other_node) * 2 + other_orientation;
          // Go there with where_to
          offset = index.where_to(side, offset, other_side);
          side = other_side;
          // Continue the process from this new side
      }
    // We have a thread to follow, take it

    cout << i << " " << j << "\t broken?: \t" << check_if_thread_t_broken(path, index) << endl;
    }
  }
}

void extract_threads_into_threads(xg::XG& index, string output_path,
        int64_t start_node, int64_t end_node) {
  end_node = (end_node == -1) ? index.ts_iv.size() : end_node + 1;
  for(int64_t i = start_node; i < end_node; i++) {
    // Skip it if no threads start at it
    if(index.ts_iv[i] == 0) {
        continue;
    }
    // If there are threads starting,
    for(int64_t j = 0; j < index.ts_iv[i]; j++) {
      ofstream threadout(output_path+to_string(i)+"_"+to_string(j)+".thread.csv");
      // For every thread starting there
      thread_t path;
      int64_t side = i;
      int64_t offset = j;
      while(true) {
          // Unpack the side into a node traversal
          xg::XG::ThreadMapping m = {index.rank_to_id(side / 2), (bool) (side % 2)};
          // Add the mapping to the thread
          path.push_back(m);
          // Work out where we go:
          // What edge of the available edges do we take?
          int64_t edge_index = index.bs_get(side, offset);
          // If we find a separator, we're very broken.
          assert(edge_index != index.BS_SEPARATOR);
          // Does the path end?
          if(edge_index == index.BS_NULL) {
              // Path ends here.
              break;
          } else {
              // If we've got an edge, convert to an actual edge index
              edge_index -= 2;
          }
          // We also should not have negative edges.
          assert(edge_index >= 0);
          // Look at the edges we could have taken next
          vector<Edge> edges_out = side % 2 ? index.edges_on_start(index.rank_to_id(side / 2)) :
                index.edges_on_end(index.rank_to_id(side / 2));
          assert(edge_index < edges_out.size());
          Edge& taken = edges_out[edge_index];
          // Follow the edge
          int64_t other_node = taken.from() == index.rank_to_id(side / 2) ? taken.to() :
                taken.from();
          bool other_orientation = (side % 2) != taken.from_start() != taken.to_end();
          // Get the side
          int64_t other_side = index.id_to_rank(other_node) * 2 + other_orientation;
          // Go there with where_to
          offset = index.where_to(side, offset, other_side);
          side = other_side;
          // Continue the process from this new side
      }
      int64_t total_width = 0;
      threadout << "tot length \t length \t ID" << endl;
      for(int64_t k = 0; k < path.size(); k++) {
        total_width += index.node_length(path[k].node_id);
        threadout << total_width << "\t" << index.node_length(path[k].node_id) << "\t" << path[k].node_id << endl;
      }
      threadout.close();
    }
  }
}

void B_experiment(xg::XG& index, int64_t internal_index, double recombination_penalty) {
  vector<int> permuted5008 = {4861,1611,905,4649,2367,3058,3301,103,1729,718,3651,1845,1506,70,2005,4059,4206,3974,3241,3956,1070,4248,4797,3539,4983,184,1525,411,909,3202,4986,2941,2226,3924,4623,4722,4752,1985,3660,68,4155,4270,1948,3669,881,4600,2035,3023,3401,3076,3297,2520,2486,1708,1477,1362,349,1674,2586,2268,2826,1794,2570,4212,3763,260,1104,3845,3321,4202,1968,1510,4825,2093,627,3256,1865,4557,1879,2984,396,4530,2911,863,1097,1058,1239,776,385,4221,3175,3836,216,3536,4659,616,2417,2193,1254,3365,2411,2398,3441,1366,4764,1818,2694,2407,4302,3621,3715,4421,4362,932,1388,1647,1726,3351,3125,1997,3132,3730,3210,2003,4590,2621,4833,3152,173,2573,4867,2327,1831,2914,1050,237,570,2169,3567,4943,3059,3992,1650,3211,840,2155,2129,470,1312,1393,2775,368,2121,3546,3563,2745,1534,4817,3744,2795,4081,2266,1813,3990,3603,3803,4584,3307,917,4340,4631,1943,4343,3649,469,631,2236,3483,417,653,1670,971,2372,2041,4400,1203,2466,2819,1392,2654,4068,2243,2898,2689,3473,923,4157,953,4761,2162,2592,3666,526,3533,4341,3917,2316,448,4278,4087,2502,2583,3162,2824,535,4509,238,3813,2605,1093,3290,4120,477,1586,3063,3228,4978,1952,226,1530,232,1667,752,1479,4318,3849,1325,1435,242,4121,331,1374,4553,4187,3398,831,3421,3406,1942,4668,1167,220,4874,1617,2613,1649,1193,4812,3959,2918,871,3423,4828,1103,3678,1503,3553,434,207,3637,4790,4385,2763,3517,792,1624,3525,1343,3869,4805,5,4745,3901,2143,878,1003,2040,324,1588,1207,4650,3575,3551,1485,452,621,2410,3516,4185,980,2543,3016,1979,2464,1078,4045,4931,337,1893,808,3843,665,3706,4743,3594,3915,3775,2614,1522,2103,946,416,741,224,1024,2432,4555,3288,3896,4208,890,1815,4625,593,1499,338,4878,2724,5002,494,1359,1556,310,2828,3863,3212,1846,3166,2020,4879,1402,708,2298,970,3449,97,2094,2909,1286,3895,1155,4244,4533,2860,3518,1955,3543,743,4733,4243,2234,356,2663,2498,4258,2288,2390,2087,4314,3781,3201,4644,4756,3757,1644,1596,4052,1810,4670,4502,4560,8,2132,914,1944,3809,1959,1480,3477,573,1626,1069,4562,431,2031,1785,872,3820,447,2056,4123,4147,4487,891,4474,109,4150,1230,2674,3719,3299,2871,4429,2539,4731,1921,2991,1299,4495,2071,894,2246,3116,1528,3850,4379,2835,2913,2348,3489,294,2695,1213,4231,479,2986,2537,4786,1211,4271,3647,459,2665,4174,141,1702,1509,2135,3071,1966,630,4974,2753,620,658,95,503,3664,3750,1907,45,1296,467,4472,4637,309,3395,2831,3326,2680,841,3120,2044,2947,2713,818,4768,4973,2552,4443,4464,2712,2054,4383,2856,3587,1087,4488,533,1880,3189,3639,3903,3494,536,3865,2208,3972,3034,9,4995,4655,1151,1306,709,2163,738,2346,4788,1894,1373,4900,3573,4626,4681,1587,2461,2997,1037,1864,2953,244,4371,4413,1987,1051,3571,2116,1541,1323,122,2962,3557,3653,3938,3002,4130,4053,290,3045,4133,4911,1672,4217,476,790,3812,1146,3472,4005,268,3625,2604,3870,595,1773,137,4425,4992,2287,1866,193,1106,960,816,964,4207,4057,868,1209,398,1628,851,2286,4388,1314,858,4457,988,4635,2386,1380,3981,4585,3889,4807,2066,4143,2253,2876,1620,1915,3431,2730,3876,712,3953,1196,2211,3402,2350,4370,1884,3964,710,2459,4914,4384,2703,703,1935,4607,3740,488,4213,4706,1074,553,3382,734,1016,287,4404,853,2472,4862,2099,453,199,4090,392,1625,3137,274,3902,2827,4897,332,4107,2780,4859,4489,3048,4398,3355,1294,2973,1927,4895,3413,3476,1938,1859,4292,3844,1780,1909,4188,2748,679,615,4066,2667,3887,1621,4969,2152,3875,2115,1924,4432,1424,1782,2364,4442,3165,4824,1236,3114,2258,1095,999,3276,2946,3081,731,2837,1266,1720,1476,2252,1903,2563,4159,450,4042,2050,471,4769,2320,2940,1229,4729,2339,2972,83,2467,1933,2423,828,2749,463,320,4918,2983,3259,772,1006,92,4134,3360,4638,2323,4027,1249,58,23,3315,2091,4500,3793,4420,2633,591,386,4679,1124,1189,4426,3658,926,1597,4849,3822,2476,3541,3012,866,4054,1946,2218,191,3963,3258,688,1129,3692,3227,1886,380,1032,1690,3097,3936,1900,3975,4964,3075,4257,4526,2340,2128,3263,266,4539,4799,2859,2949,3108,2293,1018,1222,3690,497,3537,412,422,2166,4407,2971,4740,969,4279,1953,4793,4480,826,2900,3601,3458,4930,486,3453,3286,2095,2772,1134,880,3420,2684,2676,3050,2391,3735,825,3128,1411,3339,4282,3912,699,2083,2782,3239,2274,4505,1793,212,1996,864,656,3617,2739,2981,4801,4186,4139,3347,3052,1434,435,1015,4658,1645,1518,2922,4392,163,2812,5006,2493,1163,1231,1825,472,4098,2912,31,1700,2064,116,4677,1132,3198,4576,4572,3316,601,317,1875,2928,2230,3633,4490,99,2247,3324,2815,4988,198,4030,1263,2295,3631,2403,3265,351,4858,1652,650,177,2282,751,3761,3226,2186,418,1012,4775,546,3652,3283,2113,3554,2671,4112,2863,1839,921,4506,2239,4069,4194,4842,2825,4062,4844,1526,62,91,2761,2223,629,2433,2519,4273,3407,748,1361,2963,4800,1445,3149,4285,3507,4071,2052,3616,128,4711,3511,2988,3646,3756,1961,3067,2059,4713,208,1415,1717,4755,2501,793,4397,1899,1805,2741,611,2178,4760,506,2369,2478,1532,4510,4813,311,2475,325,4734,2936,61,1932,3147,3759,2470,3612,3547,3344,1668,4077,282,1116,4427,1404,2487,2025,4599,4558,
  1075,3578,4226,262,2776,1114,1705,73,2349,538,3090,2608,2515,2203,4129,626,3099,2207,585,1327,4070,4055,931,957,3470,363,2778,4393,892,1264,178,1876,3861,3331,2889,4611,2560,182,3774,3303,4923,1444,4762,3806,1228,3135,1368,4276,684,3697,1489,292,4972,2175,3827,3703,1800,1118,4883,2529,2832,2736,1439,721,2996,2630,2039,3526,4535,4703,4532,3245,2522,4154,3929,3946,2877,2474,414,942,1468,2471,833,3514,3856,394,1740,2841,1165,3899,847,4742,1265,330,1660,22,3934,1406,3455,1338,3169,2636,4961,608,3392,2027,1888,4689,2732,1320,285,992,1017,3104,4446,2140,3529,1967,4571,4263,897,633,3379,4598,3583,1096,2345,1746,2652,3885,170,814,378,1988,313,4201,1627,391,918,3388,3950,4349,194,1770,4665,2381,4220,784,848,2409,2257,704,2596,2861,4806,263,3383,1766,1914,3270,2750,1046,3312,966,4044,1081,2304,374,2495,3440,529,1497,1040,635,4380,16,2915,1085,4467,4135,4477,1549,144,1348,2104,2752,4727,3790,1395,4164,3010,520,2894,3696,4105,2188,1131,1716,2437,1956,1159,2707,3508,3359,2738,2396,1492,3540,499,3223,930,4091,1637,3997,1577,4305,2773,254,4757,487,1014,4680,204,1562,1687,4675,1162,1322,1454,4189,2994,10,3275,2958,2755,3466,2799,464,4732,3555,53,3415,4873,1809,3354,496,460,769,636,4078,298,1527,2431,2931,2233,4826,1778,3656,2078,1308,3888,4643,3582,2711,3005,3606,4419,2842,302,3746,672,3465,250,3599,2507,3574,2131,1842,2453,4015,1777,4329,1272,1110,778,4816,4953,794,3598,2976,3911,1183,1796,377,369,3717,457,716,3491,346,1334,3036,3568,3001,4687,4238,1204,2238,3852,4511,1566,40,4886,1817,2285,3797,4554,272,2852,2294,458,3641,783,1315,4304,4335,4223,4253,3920,2924,4360,1646,4766,4851,2598,3017,654,599,555,1593,1615,1555,3609,4857,4470,247,2662,774,2120,86,3833,4326,192,4431,3126,4043,1446,4628,1575,1975,468,44,2053,3397,2525,4579,1355,421,2764,4915,4934,767,1053,4968,256,3179,861,1052,3139,906,4926,1680,4482,2292,1589,4991,842,4430,3766,2006,1502,1635,4250,1482,2769,139,559,4763,3243,4433,2205,152,76,165,2235,820,1648,278,1139,4496,4109,695,1195,1862,3792,1853,3436,4483,1135,1750,4190,2130,3714,2271,3308,339,4449,681,4088,3786,1472,2100,1284,669,1565,4104,5007,3987,3122,3588,4701,1689,2624,3057,4798,2666,3289,3433,2183,1601,4320,304,2765,3754,893,2905,2506,4439,2046,1000,4381,3253,4126,2418,3531,4888,4158,2377,3368,4298,174,2618,2187,2321,1570,4466,785,1925,4662,1995,624,2623,1430,2314,2848,3163,2254,2942,1694,4512,997,4609,3093,3839,1963,2698,4471,1631,2927,1923,1111,3429,3335,1677,4375,1285,2932,671,483,4058,911,579,2561,359,1937,131,2021,1849,1623,889,4693,4264,3006,3009,502,4887,4352,2419,2114,2944,4622,43,1950,1341,500,2048,293,3862,2118,3035,1177,3485,2541,2643,2134,3770,4096,2111,1375,2599,3173,415,4182,3976,291,3313,157,93,2261,1252,2569,1067,1977,1047,2157,3084,2421,4454,3913,2546,4272,3512,41,4863,1784,3236,2307,4013,3700,3682,1618,3387,3768,2329,2720,228,3513,4685,2625,2733,1803,1713,884,2195,2504,2886,580,3872,2603,57,395,1744,1904,4315,1391,1537,1679,249,1656,4219,4676,2967,2851,4405,257,3131,2158,542,3607,3464,4073,2679,341,1610,4249,2347,4696,2336,4440,1086,4593,1043,899,2180,3142,2906,706,958,952,941,4889,4082,1936,3475,4100,2650,2933,2255,3831,3796,1227,713,1429,1964,3958,571,3743,4854,2033,2060,2930,3500,2547,2670,82,4445,2653,603,582,693,428,4901,348,2489,1727,2838,1153,4436,475,760,998,2542,243,1332,2627,4804,443,1661,3281,836,4699,1023,2793,1757,1801,610,4789,929,4191,111,2975,1170,1344,3073,4620,3148,215,1186,4086,4447,241,1309,550,1636,214,286,4291,4815,3320,2509,4721,2399,740,1364,2473,1798,3340,59,3734,4196,233,164,3435,4424,3866,3167,3502,3877,2262,2303,3255,1372,3523,4725,4589,4333,301,3828,564,1695,3600,651,4602,1926,2675,787,3880,835,189,4939,2839,3954,739,3538,1013,2378,50,4690,2305,4956,1984,4485,4307,4256,745,3041,560,742,2875,1401,2172,252,3068,3091,2127,3904,1523,4976,3906,1045,3086,4712,972,1814,1084,2989,510,1042,990,649,55,547,3261,2538,3024,4993,4795,3699,4493,4402,3367,1812,623,2937,3390,1571,1501,3881,3414,2374,950,4461,1837,3444,3832,3988,750,694,4355,4832,612,1038,517,2845,2935,1438,5003,1475,2081,2817,4898,4835,1721,2766,549,934,2149,261,2393,2481,2810,4779,2395,1304,3025,1261,148,4621,2007,1260,4141,4567,229,3372,3930,3751,2890,1199,4410,1682,1137,991,4151,4728,3498,3384,1531,4630,676,4601,4142,4947,4925,1598,5005,3925,3306,3688,3088,1642,1370,2225,1516,4313,3859,4136,1887,1068,4970,2656,4673,4020,1166,4092,1957,2276,2979,2380,1592,4739,4652,922,2678,3562,2867,642,1579,2445,1191,2792,807,2960,4784,1947,90,561,3107,2209,4274,511,1483,4083,4508,2558,2012,2803,3998,138,2043,4587,1769,1504,3353,2062,2184,1678,35,295,1371,2224,4499,652,1986,1569,4076,2584,1890,89,4034,3143,4075,667,824,441,3238,4909,1790,4038,1182,4596,3252,2768,4559,3893,4546,4037,677,2198,2231,4160,1414,3693,1970,1302,2759,3520,1816,2204,2999,2688,3648,3,4441,3478,4880,832,3695,1514,2868,2206,
  3753,3891,305,4927,4366,3300,1326,2469,2426,637,2306,429,2664,2737,2897,3463,3207,874,4707,2096,4286,3782,3271,3229,2562,4168,4516,910,1683,2516,1913,1554,3342,3952,4000,154,2805,3620,3657,3622,2352,1863,156,2494,4200,632,3848,1169,267,690,2951,4080,114,4267,3769,3282,210,4765,1346,4700,3799,3685,3418,3480,819,3687,2993,1572,1607,4089,1161,3867,782,3155,2256,3882,2872,1405,2325,4039,4198,4654,360,437,117,333,3701,4564,1843,1233,3105,3853,2297,3302,2119,2492,1619,1697,754,3213,3965,3640,508,4917,4575,4336,551,735,3446,4406,3650,1664,1428,862,2794,3451,3860,327,4955,4297,2880,3077,442,4178,258,3522,4145,2197,444,732,1885,4810,3947,1447,1735,2216,1062,3264,1083,4697,1004,523,74,3986,3442,3983,1709,161,1609,3140,3835,2404,859,4919,2760,3404,4582,607,686,1466,3443,3065,1469,3287,119,4347,345,223,993,541,4773,2221,2449,1467,3608,3403,584,726,3118,1073,1901,2388,3565,1298,2047,4952,2042,12,2080,2575,3018,150,2416,4006,3749,1951,105,2865,3409,2565,4299,3710,4547,1407,1835,1319,2057,1929,786,936,2553,379,401,4811,4416,505,4162,675,393,2908,4369,1174,1889,4242,2995,575,115,1140,554,3493,52,3136,2428,4580,1363,945,2578,1441,3185,3939,1009,4476,4912,2069,2601,3328,1019,2658,2702,1872,1548,3760,1091,1545,1246,668,765,3805,1423,1739,449,432,2708,697,2086,123,3311,2572,678,1663,4672,4774,928,4780,1671,4577,167,4549,4708,130,4479,545,4565,4933,2462,4551,1156,3061,1715,1928,108,3894,5001,438,2384,3062,1747,4018,1507,1258,439,2693,4855,1300,2125,30,1115,1666,3434,279,1574,1457,2714,531,3356,1206,1178,3430,4847,979,1478,1714,1123,4717,2500,4056,1079,1804,1685,2526,4834,2434,4591,1384,2588,4262,867,2788,2309,1698,1871,4351,3357,1529,2798,895,3309,3437,2444,4656,3008,1347,4653,259,2334,1100,2324,2146,1653,3842,1908,2451,4613,2330,2550,1808,2036,2063,3868,2723,2559,1255,3053,4990,771,838,995,1063,2219,4684,1399,3371,4177,2356,81,3821,221,1452,3771,2503,3234,647,4131,1860,927,297,2530,2065,2106,4338,3180,1789,3659,2026,2201,2610,136,3244,4165,281,1243,2791,3618,4205,3295,2365,3991,2232,2150,4688,3235,3826,3115,1173,2957,3333,2722,1982,4061,1722,4758,1712,2137,387,587,3221,1126,4636,2862,3231,2001,1775,2655,3391,548,3694,1242,3215,2729,1749,2754,4140,2301,2338,1099,3741,1991,1060,746,3033,2430,3447,590,3109,2632,2904,4821,4486,4715,3919,2635,3183,1293,3314,4515,3362,2227,3222,4246,592,2248,3378,1599,2406,3994,4183,4819,3817,440,4295,3671,4184,532,1488,1638,4232,370,855,4913,3462,3948,4032,3361,1989,4840,1867,350,3847,1838,4498,3003,1253,4588,646,2771,4522,2801,1524,3310,4581,1278,1686,2508,248,4875,1420,1791,2370,4962,1271,1919,4666,4024,2734,1389,405,1767,4885,1259,2482,1147,985,4563,2954,733,1762,4280,3298,2333,1873,3519,2884,4929,723,79,885,994,4437,1154,4592,3720,273,4632,2296,924,4225,1175,4365,1273,756,2640,2429,3171,406,4403,577,1753,1776,3762,4782,3247,3683,3293,1641,4339,1241,3772,1426,2161,245,3878,4064,552,2921,3655,3857,3134,2505,4114,2319,2250,2070,1494,1465,3798,4881,1408,125,4615,1844,2685,2058,445,2217,3830,4783,2523,3072,2661,645,124,2145,2122,3200,3129,3000,4337,606,4289,4293,4294,4944,1168,3589,2945,2725,4719,717,1883,2600,719,4345,3412,3490,1719,1202,2626,1329,3279,308,870,4016,788,2245,3718,4640,4785,66,3585,2580,1869,4019,4751,4528,4269,1458,522,3079,326,3046,1357,4959,775,493,700,2156,1349,2645,4456,2284,419,2797,3804,3898,2107,4951,253,673,4830,4144,4869,1544,4011,195,3580,1270,4079,4868,574,3779,3156,3322,4767,3426,796,3481,3334,3051,4236,4967,2686,2709,1533,777,691,4958,3445,3905,3619,2846,3224,2022,4357,21,4239,4152,619,2556,2313,919,2443,1897,1550,3944,2170,3044,1337,1041,528,682,101,3935,4386,2342,485,655,648,4,2606,4907,1351,3123,2789,4050,4452,1681,813,2181,2800,4455,1806,2452,3192,4310,2977,518,219,4227,4021,3794,4211,534,1906,1317,657,1011,2067,180,3724,1941,3272,1036,670,3951,1035,355,299,581,1379,2485,18,1,2998,1836,2840,303,1277,4169,1462,4606,2992,4072,3336,2355,3672,2691,1377,662,811,2153,1310,773,829,3945,54,4910,803,1912,3488,2045,799,1728,1830,3343,4928,1214,1827,515,4281,3967,4605,3274,4214,3712,2289,4041,614,1755,896,504,1896,1535,2705,3858,2544,4171,3510,3654,981,817,2344,2631,4544,251,1410,2177,2888,2807,1811,2318,321,2196,104,3542,2540,3257,340,4527,1194,312,4266,2796,4667,3277,3610,3550,235,202,200,3778,4678,2015,2456,1394,3352,2574,3043,2110,1590,3174,4822,2512,3748,4460,4412,3417,3187,2200,2277,1133,3280,2415,1760,4468,4845,1834,3854,1313,3758,2581,1007,3188,4237,2034,605,4642,1148,628,4209,288,3054,2142,3824,2275,3521,977,1090,2394,2174,2533,155,3196,4288,424,4179,4067,1998,1217,353,976,3363,602,1662,986,1451,4542,3424,597,2521,3499,959,3632,3908,1376,2641,3082,4960,1634,382,1205,4945,1684,3254,2450,3698,2866,2049,3890,33,4074,576,1898,72,4548,2382,3883,3840,3837,4009,277,4778,2072,947,1692,4998,3020,4317,1318,1157,1232,800,85,2970,3738,413,3742,4989,2785,1064,2535,
  2439,1561,2075,625,3205,2424,1221,1736,381,4503,4772,2402,3569,1881,1823,634,2622,2124,2594,3979,2955,727,4128,1703,3549,2966,318,1676,938,4890,2463,1490,3823,4113,3292,4040,179,1002,1993,1731,1745,3214,2192,1324,225,1238,1560,3100,1180,407,1188,408,300,63,334,3232,2013,4569,618,4245,1696,4669,3846,1235,3918,1333,3602,3393,3816,206,3630,4648,2549,3267,188,2978,3492,2089,2647,2213,383,118,1216,2701,598,3926,3593,1848,1336,4709,4448,3940,4414,1553,1149,761,2536,3661,240,4529,1220,322,2881,3096,4010,4116,1520,1723,1117,3345,1076,3181,1629,933,3907,1763,1655,544,3524,1419,2629,3584,352,4633,1515,2965,265,2528,4698,3897,516,1330,1200,2595,3191,1756,1508,4594,3195,46,3532,873,4417,19,168,2593,2943,3923,2692,2566,565,2683,3022,3150,898,768,316,846,639,902,1852,4809,2082,0,4148,400,2853,2568,3765,1874,1138,2874,2210,3728,2854,4724,3564,3592,2527,375,2165,2762,3190,4012,4691,1730,2907,3993,1799,2779,1613,643,3157,4122,264,3636,4193,1127,1275,4531,4997,4197,4550,3102,1752,1150,4312,3250,4163,711,3739,875,4975,763,4791,2436,3119,2809,4823,1910,4843,1917,4634,1643,720,2777,2425,939,3242,2585,2108,3380,1437,1030,143,1181,3810,1916,1121,3348,3980,4368,25,4692,617,255,3172,3338,3572,2389,166,4309,1754,1125,3130,1978,4364,1594,4501,1585,2833,984,102,1733,3999,2770,2638,426,1847,4570,361,2743,2873,2283,315,4938,1102,4356,3330,3675,96,1807,1197,539,4453,4514,1279,367,4216,1345,4195,3825,935,3591,3623,3909,328,1417,2272,2878,2609,3932,3015,4866,1055,1781,4475,3094,961,3733,3220,296,3559,67,2616,1365,1774,2806,2076,3176,456,925,965,1031,3419,2513,4520,2646,749,1765,4741,3197,3141,3204,489,2577,436,2774,3829,2484,2483,883,2109,1707,4753,1418,2164,3037,1008,4376,482,2879,779,1795,3127,4736,4180,60,2077,1659,805,1215,42,2029,1962,203,2822,3038,2811,3505,3168,1854,14,821,1665,3721,3159,4491,1094,4099,2138,3708,1779,3673,4513,2212,1331,1413,4716,2756,687,1250,1761,1071,4814,3785,2244,329,4604,2051,3106,4478,3251,3970,3160,2435,4170,13,3117,3262,1737,509,2136,3535,1604,4330,4215,937,3624,4950,1276,3874,1699,513,908,88,2602,879,2787,1706,4517,3386,4358,1109,1891,4595,1421,213,2844,3266,3456,4574,2758,1471,3233,674,1383,1824,4838,4265,4321,4428,11,3400,3795,1287,146,1369,1026,4856,4894,1176,913,3004,3626,1974,2455,1710,1797,1512,2816,3070,1612,1639,1552,4902,1295,2028,466,2496,3087,962,4481,989,38,4261,4348,1339,4418,2160,762,1367,3732,3085,4230,3029,519,3663,2849,2229,987,27,2882,3428,3013,1305,4796,4469,354,4451,3121,4746,568,1119,4573,4118,1980,1725,1930,49,1931,916,3237,3031,4996,3791,3319,4540,4695,3679,3145,276,3039,1464,2144,2673,2299,2648,3377,683,2264,364,4777,3203,2335,1940,4494,4389,140,3219,3452,2929,4268,572,1741,3808,478,4173,996,4260,698,3030,239,1080,65,2024,827,3783,3548,2016,2273,2176,4802,4332,3369,4137,2084,480,4382,3153,2901,1122,2392,4537,4608,1616,3375,4748,1403,4860,2923,3745,1311,1861,1918,2098,4334,1291,659,1971,2263,3729,409,112,4322,2189,1981,1390,2969,3984,75,1658,2985,3552,1461,3707,1500,887,2950,1614,4111,705,725,1511,4882,2744,3851,1251,514,1922,2154,4492,3604,2836,129,3922,2781,2829,948,1436,1828,954,1257,4283,3381,4614,4046,844,2214,3056,2746,4093,1342,3784,1711,2376,2491,1820,3503,1855,2642,1840,4981,3978,69,3527,2337,1455,4850,3614,830,4965,3955,1029,4218,4870,967,537,4204,1198,4871,4401,4199,3635,3028,3170,4084,4434,2167,2088,1742,2373,2249,222,3780,2740,4942,1972,3815,1990,3937,4730,4839,1210,121,3208,3206,186,1992,4001,809,2037,3366,1301,2576,2032,4316,4770,1945,1340,2171,2018,4241,3556,1328,4387,1010,3586,3370,3776,2510,147,4411,3871,2038,404,4646,4235,4004,3217,1751,724,1772,1453,4893,2308,3723,4377,764,1201,4138,4616,4921,4541,4106,1969,3482,3566,2401,2690,907,3273,2290,1563,1584,4735,2607,974,3186,1540,2715,3349,4538,1179,153,1764,1142,4119,558,1583,2699,4710,3471,3879,2813,3787,3838,2717,1460,4222,798,4458,1101,465,1573,1303,1112,3422,3095,4534,1473,2422,1851,4181,36,4578,797,1386,187,2202,1290,1141,4932,491,1378,4660,1521,2820,1732,2802,4233,900,4507,4311,2704,4396,4827,4586,230,2517,3921,3629,2353,190,4033,2447,2446,4350,1487,1911,3681,865,107,2420,289,4008,4275,4957,2332,474,2858,1321,770,4342,801,604,4462,2055,1288,1089,4629,666,955,4319,2079,2126,3892,781,843,461,2728,4395,1269,29,1895,1005,956,2011,4747,3777,1976,4247,2959,1651,3611,4146,4029,583,1397,3755,3138,2190,4047,3230,4603,3931,3479,3725,3112,3676,4612,492,3942,596,3576,2315,2269,2847,4325,3595,1381,2892,132,1144,342,4115,4744,1400,1223,1442,4994,3886,3385,2068,3133,1748,3291,4754,1297,1578,2074,521,3689,3680,886,2267,4224,384,4394,660,644,2023,903,2617,2554,2397,3284,755,4035,849,4545,4175,4714,3113,557,2441,397,1219,2830,4940,3722,1983,17,3506,1724,3066,856,481,2783,569,757,2682,2362,2123,2358,973,1832,2242,2987,3914,2823,1841,2117,2916,3515,3083,4980,3026,2644,4346,4948,98,1870,2251,4372,877,358,3613,3570,1470,
  753,28,2002,1316,4568,280,3487,171,512,1994,3240,2651,2363,2857,4327,2731,795,2910,839,3711,3146,4726,2310,766,78,4303,1001,1450,4831,3364,1771,4051,1743,3686,4523,2442,2343,48,4624,1120,4161,1358,789,2311,4172,3962,2280,3927,3736,2460,158,3432,169,3884,3394,1130,4971,3643,373,197,2564,1657,3709,1958,2518,420,4837,4829,4892,3249,3496,4408,84,1513,2557,376,578,4022,707,4750,1021,2677,2571,4852,3662,963,3968,4117,728,307,4125,347,1164,1688,6,2465,3332,943,1999,4808,362,876,1234,3484,2408,1675,4985,2148,4017,1595,3713,246,4094,37,1542,2073,430,4864,366,3638,4705,3802,2133,904,4543,7,1020,3329,1564,399,4324,3177,715,4374,4049,912,2185,3376,2619,722,1640,2719,4818,4963,1082,1608,3460,1443,1240,4300,860,3049,176,1281,3534,2322,501,2742,4254,3425,2260,3032,3408,4720,507,2887,1934,4999,1538,2920,2767,3040,944,882,3427,2948,2488,4884,3801,4663,2405,3704,4484,2302,56,47,2085,451,1066,484,1190,2956,3642,3305,3199,4936,1048,1136,1654,3560,2009,2413,2934,71,4803,4284,120,2855,2895,3285,2326,3074,1225,4399,284,4097,4036,2968,4949,3731,4977,1519,1274,3269,4048,4627,1551,1581,3811,901,4354,1856,4954,806,729,2628,490,594,3973,3154,371,2357,4228,2587,2222,1822,2341,1734,3727,2090,4674,127,2279,3504,4103,3098,1226,3943,3528,162,4014,275,4641,4556,149,1352,3246,1416,77,4920,1622,1877,3995,1307,1558,978,495,172,804,4344,1098,689,3864,4617,5000,780,2620,34,2030,3497,2925,39,3665,3834,4891,661,1145,3103,3218,1580,3151,4450,15,4619,3060,1788,3677,2312,4853,1431,80,4941,1768,181,1160,638,3374,4702,641,2710,3457,3439,3078,1353,2361,390,4922,269,1248,525,1493,744,4916,1065,4378,4287,3080,2939,3971,1481,1158,562,822,524,4664,4903,3716,812,2,4435,4987,446,692,2697,2590,4060,2727,4023,4229,4525,1783,462,2747,3933,4583,142,1960,3773,3789,145,4738,4149,2440,664,3577,314,1787,2454,1396,2902,1965,2551,126,344,701,3318,1539,271,4561,4610,4937,175,1360,1718,4657,2964,2179,4110,1557,3969,3628,1486,1517,3341,4124,2004,1172,3055,2368,1802,3461,2700,588,1850,888,2790,1433,425,1382,1463,2438,1280,423,1398,4524,4277,2191,218,1829,1758,1025,2241,2784,4841,2903,4536,183,3900,201,4085,3358,4423,3579,3800,3346,2477,4373,4234,1474,4007,236,3014,527,4463,1459,4003,1868,4519,2657,4063,3225,3767,1224,455,1054,4792,134,20,834,2611,1107,4028,3089,4794,2615,4749,1061,3590,852,4979,4904,2706,2982,815,2883,1559,227,2270,2173,1088,4651,1354,4908,4836,4737,1826,1546,1484,2597,758,2387,1022,217,410,185,3764,4683,3019,3182,113,498,2534,3350,32,1245,4661,357,1701,2354,4359,389,1892,4877,1882,4872,1143,1268,4671,4308,3317,1244,4176,1505,2010,1432,4924,2834,4363,4459,1759,4166,3405,1044,3841,3645,1128,837,2843,433,4776,133,2019,920,3910,1567,4820,135,2961,3674,1630,1184,2669,87,3158,3110,3327,1496,3644,196,3448,4566,2278,2524,4301,2366,696,4618,2151,2360,2097,3605,205,2427,160,3325,4982,730,3468,2240,2220,1498,1536,1858,3007,3989,2215,810,2938,3818,3021,1056,663,586,1693,3509,3194,2105,3416,3396,3389,4095,3411,2589,3047,1108,4694,3530,3184,3691,2468,1283,2726,3916,1605,2696,1491,1237,1821,2000,1902,1669,3982,2497,4132,403,4192,473,3819,563,2637,1412,4848,4504,1973,4025,2659,106,1954,3069,2919,4865,2499,3193,4723,3807,1606,2893,2168,4391,857,3399,915,951,2818,3042,1456,640,4597,2300,530,4896,3454,335,4259,3684,3855,3737,2102,2672,3960,2757,372,983,323,3977,4026,2480,2974,336,2147,3495,1356,2385,1033,4966,234,64,2545,4306,4065,1335,1208,2870,3966,2721,4240,1633,3949,3668,4645,1582,1289,2228,100,4409,2808,1028,1673,823,3544,4296,2639,1105,1576,4328,3064,3985,3957,2383,2751,1171,1282,714,2379,609,4251,2112,343,4290,3702,1495,365,94,1039,3670,2328,3459,1262,3450,4946,702,2681,4252,2014,949,1256,2869,4153,4906,1387,270,2821,1049,543,2281,2804,1591,982,2359,2917,1059,2412,51,2896,1632,3178,2926,2017,2555,680,2490,3788,2814,4876,151,3144,2182,1185,5008,2400,737,3941,3561,759,2885,2850,3597,1034,4781,4323,2718,4935,1192,3501,3304,4984,2634,2735,209,3248,2514,2864,4444,1425,2579,1792,4203,1786,4102,2891,4361,26,4647,1547,4101,3101,869,427,231,2716,306,3486,3581,1603,1057,4438,2582,2591,4422,5004,2511,388,4899,1187,3323,1440,2980,2371,4367,2259,2061,3545,3961,975,4108,1602,2375,1427,4255,1218,747,3027,1878,4031,2649,3438,3260,4473,3558,1568,3278,968,3011,2567,1212,3092,1833,4156,2199,2414,4465,2237,1350,736,283,2899,2952,940,3469,1949,3124,4905,1247,3667,1448,540,3209,3216,4771,3615,3752,3111,791,3410,802,2458,1819,4718,613,2159,1152,4704,2660,1072,4497,845,3337,4846,4331,1267,4639,3726,3164,1543,1857,2668,159,3873,3705,4127,2786,1704,3467,4415,2139,4518,2101,4787,3161,2531,2351,2457,3928,600,2092,1409,567,2194,4682,2008,4390,566,2141,1738,1449,1905,1385,3627,3294,2265,2331,2612,1920,4167,850,2990,685,854,4210,1422,211,1113,2548,622,3996,1027,1691,3296,4759,3373,4521,2532,3474,1292,3268,4353,1077,1092,2479,2687,589,3634,4002,1939,556,3814,3747,4552,4686,3596,
  2317,110,24,1600,402,454,2291,319,2448};
  RRMemo memo(recombination_penalty);
  for(int64_t i = 1; i < index.ts_iv.size(); i++) {
    // Skip it if no threads start at it
    if(index.ts_iv[i] == 0) {
        continue;
    }
    // If there are threads starting,
    for(int64_t k = internal_index; k < permuted5008.size(); k++) {
      int64_t j = permuted5008[k];
      if (j >= index.ts_iv[i]) {
        continue;
      }
      // For every thread starting there
      thread_t path;
      int64_t side = i;
      int64_t offset = j;
      cerr << "computing probability for thread starting at \t" << i << " " << j << endl;
      while(true) {
          // Unpack the side into a node traversal
          xg::XG::ThreadMapping m = {index.rank_to_id(side / 2), (bool) (side % 2)};
          // Add the mapping to the thread
          path.push_back(m);
          // Work out where we go:
          // What edge of the available edges do we take?
          int64_t edge_index = index.bs_get(side, offset);
          // If we find a separator, we're very broken.
          assert(edge_index != index.BS_SEPARATOR);
          // Does the path end?
          if(edge_index == index.BS_NULL) {
              // Path ends here.
              break;
          } else {
              // If we've got an edge, convert to an actual edge index
              edge_index -= 2;
          }
          // We also should not have negative edges.
          assert(edge_index >= 0);
          // Look at the edges we could have taken next
          vector<Edge> edges_out = side % 2 ? index.edges_on_start(index.rank_to_id(side / 2)) :
                index.edges_on_end(index.rank_to_id(side / 2));
          assert(edge_index < edges_out.size());
          Edge& taken = edges_out[edge_index];
          // Follow the edge
          int64_t other_node = taken.from() == index.rank_to_id(side / 2) ? taken.to() :
                taken.from();
          bool other_orientation = (side % 2) != taken.from_start() != taken.to_end();
          // Get the side
          int64_t other_side = index.id_to_rank(other_node) * 2 + other_orientation;
          // Go there with where_to
          offset = index.where_to(side, offset, other_side);
          side = other_side;
          // Continue the process from this new side
      }
    // We have a thread to follow, take it
    haplo_d h = haplo_d(path, index);
    h.calculate_Is(index);
    double probability = h.log_probability(memo);
    cout << "[" << i << "," << j << "]\t" << h.tot_width << "\t" << probability << "\t" <<  probability - memo.logS(memo.population_size,h.tot_width) - log(memo.population_size) << endl;
    }
  }
}


void A_experiment(xg::XG& index,int seed) {
  vector<int> cuts = {10,100,1000,10000,100000,1000000,10,100,1000,10000,100000,1000000};

  for(int64_t i = 1; i < index.ts_iv.size(); i++) {
    // Skip it if no threads start at it
    if(index.ts_iv[i] <10) {
        continue;
    }
    for(int cut_i = 0; cut_i < cuts.size(); cut_i++) {
      srand(cuts[cut_i]+seed+cut_i);
      vector<int64_t> randomstarts;
      for(int j = 0; j < 10; j++) {
        randomstarts.push_back(rand() % index.ts_iv[i]);
      }
      for(int k = 0; k < 10; k++) {
        int64_t j = randomstarts[k];
        thread_t path;
        int64_t side = i;
        int64_t offset = j;
        while(true) {
            // Unpack the side into a node traversal
            xg::XG::ThreadMapping m = {index.rank_to_id(side / 2), (bool) (side % 2)};
            // Add the mapping to the thread
            path.push_back(m);
            // Work out where we go:
            // What edge of the available edges do we take?
            int64_t edge_index = index.bs_get(side, offset);
            // If we find a separator, we're very broken.
            assert(edge_index != index.BS_SEPARATOR);
            // Does the path end?
            if(edge_index == index.BS_NULL) {
                // Path ends here.
                break;
            } else {
                // If we've got an edge, convert to an actual edge index
                edge_index -= 2;
            }
            // We also should not have negative edges.
            assert(edge_index >= 0);
            // Look at the edges we could have taken next
            vector<Edge> edges_out = side % 2 ? index.edges_on_start(index.rank_to_id(side / 2)) :
                  index.edges_on_end(index.rank_to_id(side / 2));
            assert(edge_index < edges_out.size());
            Edge& taken = edges_out[edge_index];
            // Follow the edge
            int64_t other_node = taken.from() == index.rank_to_id(side / 2) ? taken.to() :
                  taken.from();
            bool other_orientation = (side % 2) != taken.from_start() != taken.to_end();
            // Get the side
            int64_t other_side = index.id_to_rank(other_node) * 2 + other_orientation;
            // Go there with where_to
            offset = index.where_to(side, offset, other_side);
            side = other_side;
            // Continue the process from this new side
        }
        int64_t last_possible_start = path.size() - cuts[cut_i];
        int64_t start_i = rand() % last_possible_start;
        thread_t path_subset;
        for(int64_t l = 0; l < cuts[cut_i]; l++) {
          path_subset.push_back(path[start_i + l]);
        }
        haplo_d h = haplo_d(path_subset, index);
        h.calculate_Is(index);
        int Acurrmax = 1;
        int64_t avg_A = 1;
        for(int l = 1; l < h.cs.size(); l++) {
          if(h.cs[l].S.size() > Acurrmax) {Acurrmax = h.cs[l].S.size();}
          avg_A += h.cs[l].S.size();
        }
        cout << "[" << cuts[cut_i] << "]\t" << Acurrmax << "\t" << avg_A/h.cs.size() << "\t" << h.tot_width << "\t" << h.cs.size() << endl;
      }
    }
  }
}


void extract_threads_into_haplo_ds(xg::XG& index, string output_path,
        int64_t start_node, int64_t end_node, int64_t inner_index, bool make_graph) {
  //ts_iv is a vector of the # of threads starting at each side
  end_node = (end_node == -1) ? index.ts_iv.size() : end_node + 1;
  ofstream all_thread_stats (output_path+"summary."+to_string(start_node)+"to"+to_string(end_node)+".csv");
  for(int64_t i = start_node; i < end_node; i++) {
    // Skip it if no threads start at it
    if(index.ts_iv[i] == 0) {
        continue;
    }
    // If there are threads starting,
    for(int64_t j = inner_index; j < index.ts_iv[i]; j++) {
      // For every thread starting there
      thread_t path;
      int64_t side = i;
      int64_t offset = j;
      while(true) {
          // Unpack the side into a node traversal
          xg::XG::ThreadMapping m = {index.rank_to_id(side / 2), (bool) (side % 2)};
          // Add the mapping to the thread
          path.push_back(m);
          // Work out where we go:
          // What edge of the available edges do we take?
          int64_t edge_index = index.bs_get(side, offset);
          // If we find a separator, we're very broken.
          assert(edge_index != index.BS_SEPARATOR);
          // Does the path end?
          if(edge_index == index.BS_NULL) {
              // Path ends here.
              break;
          } else {
              // If we've got an edge, convert to an actual edge index
              edge_index -= 2;
          }
          // We also should not have negative edges.
          assert(edge_index >= 0);
          // Look at the edges we could have taken next
          vector<Edge> edges_out = side % 2 ? index.edges_on_start(index.rank_to_id(side / 2)) :
                index.edges_on_end(index.rank_to_id(side / 2));
          assert(edge_index < edges_out.size());
          Edge& taken = edges_out[edge_index];
          // Follow the edge
          int64_t other_node = taken.from() == index.rank_to_id(side / 2) ? taken.to() :
                taken.from();
          bool other_orientation = (side % 2) != taken.from_start() != taken.to_end();
          // Get the side
          int64_t other_side = index.id_to_rank(other_node) * 2 + other_orientation;
          // Go there with where_to
          offset = index.where_to(side, offset, other_side);
          side = other_side;
          // Continue the process from this new side
      }
      // We have a thread_t to follow; let's make a haplo_d
      haplo_d h = haplo_d(path, index);
      h.calculate_Is(index);
      cerr << i << " / " << j;
      pair<int,int> hd_output =
            h.print_decomposition_stats(output_path+to_string(i)+"_"+to_string(j)+".stats.csv");
      h.print_decomposition(output_path+to_string(i)+"_"+to_string(j)+".structure.csv");
      if(make_graph) {
        h.unfold_rectangles(output_path+to_string(i)+"_"+to_string(j)+".plot.csv");
      }
      all_thread_stats << i << ", " << j << "\t" << hd_output.first << "\t" << hd_output.second << endl;
    }
  }
  all_thread_stats.close();
}

void haplo_d::log_calculate_Is(xg::XG& graph) {
  // node 0 was done in the haplo_d constructor; start at node 1
  for(int b = 1; b < cs.size(); b++) {
    // make sure that there is at least one rectangle here
    if(cs[b].S.size() != 0) {
      // get side and orientation of the next element in our query thread_t
      XG::ThreadMapping next_node = cs[b].get_node();
      // if J = 0 for a rectangle, then J must be 0 for all older rectangles
      bool nonempty_J = (cs[b].S.back().J > 0);
      if(nonempty_J) {
        int new_J;
        // start at a = 1 since the haplo_d initializer handles a = 0
        int total_deltaJ = cs[b-1].S[1].J - cs[b-1].S[cs[b-1].S[1].next].J;
        if(total_deltaJ == 0) {
          for(int a = 1; a < cs[b-1].S.size(); a++) {
            rectangle new_rect = cs[b-1].S[a];
            if (cs[b-1].bridge.size() != 0 &&
                new_rect.state.current_side != graph.id_to_rank(cs[b-1].bridge.back().node_id) * 2 + cs[b-1].bridge.back().is_reverse) {
              new_rect.simple_extend(cs[b-1].bridge,graph);
              new_rect.simple_extend(next_node,graph);
              new_rect.prev = a;
              cs[b].S.push_back(new_rect);
              cs[b-1].S[a].next = cs[b].S.size()-1;
            }
          }
        } else {
          thread_t extension = cs[b-1].bridge;
          extension.push_back(next_node);
          binaryI(graph, extension, b, cs[b].S.size()-1, cs[b-1].S[0].next + cs[b-1].S.size(), total_deltaJ, 0, cs[b].S.back().J, 0);
          for(int a = 1; a < cs[b].S.size() - 1; a++) {
            cs[b].S[a].I = cs[b].S[a].J - cs[b].S[a+1].J;
          }
          cs[b].S.back().I = cs[b].S.back().J;
        }
      } else {
        // this shouldn't be here
        cs[b].S.pop_back();
      }
    }
  }
}

void haplo_d::binaryI(xg::XG& graph, thread_t extension, int b, int atop, int abottom, int deltaItop, int deltaIbottom, int Jtop, int Jbottom) {
  if(atop == abottom + 1) {
    // You've reached max recursion depth!
    return;
  } else if(deltaItop == deltaIbottom) {
    // The nexting behavior of the J-counts ensures that there are no changes to
    // thread membership in the interval we're evaluating here. We can safely
    // extend everything
    for(int i = atop + 1; i < abottom; i++) {
      rectangle rect = cs[b-1].S[i];
      rect.simple_extend(extension,graph);
      rect.prev = i;
      cs[b].S.push_back(rect);
      cs[b-1].S[i].next = cs[b].S.size()-1;
    }
  } else {
    // At least one thread diverges within this interval; let's find it
    int mid = atop + (abottom-atop)/2;
    rectangle rect_mid = cs[b-1].S[mid];
    int Jmid = rect_mid.get_next_J(extension,graph);
    int deltaImid = cs[b-1].S[mid].J - Jmid;
    if(Jmid == 0) {
      // all smaller rectangles disappear
      // handled by boundary behavior
    }
    if(Jmid == Jtop) {
      // all rectangles between disappear
    } else {
      binaryI(graph, extension, b, atop, mid, deltaItop, deltaImid, Jtop, Jmid);
    }
    if(Jmid != 0 && Jmid != Jbottom) {
      rect_mid.prev = mid;
      cs[b].S.push_back(rect_mid);
      cs[b-1].S[mid].next = cs[b].S.size()-1;
    }
    if(Jbottom == Jmid) {
      // all rectangles between disappear
    } else {
      binaryI(graph, extension, b, mid, abottom, deltaImid, deltaIbottom, Jmid, Jbottom);
    }
  }
}

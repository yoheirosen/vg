#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <vector>
#include <sstream>
#include <omp.h>
#include <cstring>
#include <algorithm>
#include <unistd.h>
#include "vg.pb.h"
#include "sha1.hpp"
#include "Variant.h"

namespace vg {

using namespace std;

char reverse_complement(const char& c);
string reverse_complement(const string& seq);
int get_thread_count(void);
string wrap_text(const string& str, size_t width);
bool is_number(const string& s);

// split a string on any character found in the string of delimiters (delims)
std::vector<std::string>& split_delims(const std::string &s, const std::string& delims, std::vector<std::string> &elems);
std::vector<std::string> split_delims(const std::string &s, const std::string& delims);

const std::string sha1sum(const std::string& data);
const std::string sha1head(const std::string& data, size_t head);

bool allATGC(const string& s);
string nonATGCNtoN(const string& s);
void mapping_cigar(const Mapping& mapping, vector<pair<int, char> >& cigar);
string cigar_string(vector<pair<int, char> >& cigar);
string mapping_string(const string& source, const Mapping& mapping);
double median(std::vector<int> &v);

// Convert a probability to a natural log probability.
inline double prob_to_logprob(double prob) {
    return log(prob);
}
// Convert natural log probability to a probability
inline double logprob_to_prob(double logprob) {
    return exp(logprob);
}
// Add two probabilities (expressed as logprobs) together and return the result
// as a logprob.
inline double logprob_add(double logprob1, double logprob2) {
    // Pull out the larger one to avoid underflows
    double pulled_out = max(logprob1, logprob2);
    return pulled_out + prob_to_logprob(logprob_to_prob(logprob1 - pulled_out) + logprob_to_prob(logprob2 - pulled_out));
}
// Invert a logprob, and get the probability of its opposite.
inline double logprob_invert(double logprob) {
    return prob_to_logprob(1.0 - logprob_to_prob(logprob));
}

// Convert integer Phred quality score to probability of wrongness.
inline double phred_to_prob(int phred) {
    return pow(10, -((double)phred) / 10);
}

// Convert probability of wrongness to integer Phred quality score.
inline int prob_to_phred(double prob) {
    return round(-10.0 * log10(prob));
}

// Convert a Phred quality score directly to a natural log probability of wrongness.
inline double phred_to_logprob(int phred) {
    return (-((double)phred) / 10) / log10(exp(1.0));
}

// Convert a natural log probability of wrongness directly to a Phred quality score.
inline int logprob_to_phred(double logprob ) {
    return round(-10.0 * logprob * log10(exp(1.0)));
}

template<typename T, typename V>
set<T> map_keys_to_set(const map<T, V>& m) {
    set<T> r;
    for (auto p : m) r.insert(p.first);
    return r;
}

// pairwise maximum
template<typename T>
vector<T> pmax(const std::vector<T>& a, const std::vector<T>& b) {
    std::vector<T> c;
    assert(a.size() == b.size());
    c.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(c),
                   [](T a, T b) { return std::max<T>(a, b); });
    return c;
}

// maximum of all vectors
template<typename T>
vector<T> vpmax(const std::vector<std::vector<T>>& vv) {
    std::vector<T> c;
    if (vv.empty()) return c;
    c = vv.front();
    typename std::vector<std::vector<T> >::const_iterator v = vv.begin();
    ++v; // skip the first element
    for ( ; v != vv.end(); ++v) {
        c = pmax(c, *v);
    }
    return c;
}

/**
 * Compute the sum of the values in a collection. Values must be default-
 * constructable (like numbers are).
 */
template<typename Collection>
typename Collection::value_type sum(const Collection& collection) {

    // Set up an alias
    using Item = typename Collection::value_type;

    // Make a new zero-valued item to hold the sum
    auto total = Item();
    for(auto& to_sum : collection) {
        total += to_sum;
    }

    return total;

}

/**
 * Compute the sum of the values in a collection, where the values are log
 * probabilities and the result is the log of the total probability. Items must
 * be convertible to/from doubles for math.
 */
template<typename Collection>
typename Collection::value_type logprob_sum(const Collection& collection) {

    // Set up an alias
    using Item = typename Collection::value_type;

    // Pull out the minimum value
    auto min_iterator = min_element(begin(collection), end(collection));

    if(min_iterator == end(collection)) {
        // Nothing there, p = 0
        return Item(prob_to_logprob(0));
    }

    auto check_iterator = begin(collection);
    ++check_iterator;
    if(check_iterator == end(collection)) {
        // We only have a single element anyway. We don't want to subtract it
        // out because we'll get 0s.
        return *min_iterator;
    }

    // Pull this much out of every logprob.
    Item pulled_out = *min_iterator;

    if(logprob_to_prob(pulled_out) == 0) {
        // Can't divide by 0!
        // TODO: fix this in selection
        pulled_out = prob_to_logprob(1);
    }

    Item total(0);
    for(auto& to_add : collection) {
        // Sum up all the scaled probabilities.
        total += logprob_to_prob(to_add - pulled_out);
    }

    // Re-log and re-scale
    return pulled_out + prob_to_logprob(total);
}

string tmpfilename(const string& base);

// Code to detect if a variant lacks an ID and give it a unique but repeatable
// one.
string get_or_make_variant_id(vcflib::Variant variant);
string make_variant_id(vcflib::Variant variant);

// Simple little tree
template<typename T>
struct TreeNode {
    T v;
    vector<TreeNode<T>*> children;
    ~TreeNode() { for (auto c : children) { delete c; } }
    void for_each_preorder(function<void(TreeNode<T>*)> lambda) {
        lambda(this);
        for (auto c : children) {
            c->for_each_preorder(lambda);
        }
    }
    void for_each_postorder(function<void(TreeNode<T>*)> lambda) {
        for (auto c : children) {
            c->for_each_preorder(lambda);
        }
        lambda(this);
    }
};

template<typename T>
struct Tree {
    typedef TreeNode<T> Node;
    Node* root;
    Tree(Node* r = 0) : root(r) {}
    ~Tree() { delete root; }
    void for_each_preorder(function<void(Node*)> lambda) {
        if (root) root->for_each_preorder(lambda);
    }
    void for_each_postorder(function<void(Node*)> lambda) {
       if (root) root->for_each_postorder(lambda);
    }

};


}

#endif

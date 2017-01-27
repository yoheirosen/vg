//
//  suffix_tree.hpp
//  
// Suffix tree implementation including Ukkonen's linear time construction algorithm
//


#ifndef suffix_tree_hpp
#define suffix_tree_hpp

#include <stdio.h>
#include <unordered_map>
#include <list>
#include <vector>
#include <sstream>
#include <iostream>

using namespace std;

namespace vg {
    
    /**
     * An implementation of a suffix tree with linear time and space complexity for construction.
     *
     */
    class SuffixTree {
    
    public:
        /// Linear time constructor.
        ///
        /// Note: string cannot have null characters, but this is not checked.
        SuffixTree(string::const_iterator begin, string::const_iterator end);
        ~SuffixTree() = default;
        
        /// Returns the length of the longest prefix of str that exactly matches
        /// a suffix of the string used to construct the suffix tree.
        ///
        /// Note: string cannot have null characters, but this is not checked.
        size_t longest_overlap(const string& str);
        
        /// Retuns a vector of all of the indices where a string occurs as a substring
        /// of the string used to construct the suffix tree. Indices are ordered arbitrarily.
        ///
        /// Note: string cannot have null characters, but this is not checked.
        vector<size_t> substring_locations(const string& str);
        
        /// String representation for debugging.
        string to_string();
        
    private:
        struct STNode;
        
        // keep nodes in list to avoid difficulties with pointers and reallocations
        list<STNode> nodes;
        
        unordered_map<char, STNode*> root;
        
        string::const_iterator begin;
        string::const_iterator end;
        
        inline char get_char(size_t i) {
            return i == end - begin ? '\0' : *(begin + i);
        }
        
        // debugging functions for constructor
        
        string partial_tree_to_string(int64_t phase);
        string label_string(unordered_map<char, STNode*>* branch_point);
        string node_string(STNode* node, int64_t phase);
        string active_point_string(unordered_map<char, STNode*>* active_branch_point,
                                   STNode* active_node,
                                   int64_t active_length,
                                   int64_t phase);
        string suffix_links_string(unordered_map<unordered_map<char, STNode*>*,
                                                 unordered_map<char, STNode*>*>& suffix_links);
    };
    /**
     * A node of a suffix tree corresponding a substring of the string
     *
     */
    struct SuffixTree::STNode {
        
        STNode(int64_t first, int64_t last);
        
        unordered_map<char, STNode*> children;
        
        // inclusive range of string on this node, -1 indicates end sentinel
        int64_t first;
        int64_t last;
        
        inline int64_t length(int64_t phase) {
            return last >= 0 ? last - first + 1 : phase - first + 1;
        }
        
        inline int64_t final_index(int64_t phase) {
            return last >= 0 ? last - first : phase - first;
        }
    };
    
}


#endif /* suffix_tree_hpp */

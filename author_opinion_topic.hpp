#include "common.hpp"
#include "progress.hpp"
#include <cassert>
#include <unordered_map>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>

struct Word {
  int word_id;
  int topic;
};

struct Document {
  std::vector<Word> neutral_words;
  std::vector<Word> opinion_words;
  int author_id;
};

class AuthorTopic{

public:
  AuthorTopic(){
  }

  AuthorTopic(double a, double b, int t, int loop, int show_limit);

  // Add a new document to the corpus, update all of the necessary data structures
  void add_document(string author, vector<string> neutral_words, vector<string> opinion_words);

  void sample_all();

  void output(char* filename);

  // Print the current state of the topic model. Use for debugging purposes
  void print_state();

private:

  // get the corresponding id from the author or word
  int get_author_id(string author);
  int get_neutral_word_id(string word);
  int get_opinion_word_id(string word);

  // add a topic assignment for a word, increment the private counts
  void add_neutral_topic_assignment(int word_id, int topic_id, int author_id);
  void add_opinion_topic_assignment(int word_id, int topic_id, int author_id, int doc_id);

  // remove a topic assignment for a word, decrement the private counts
  void decrease_neutral_topic_assignment(int word_id, int topic_id, int author_id);
  void decrease_opinion_topic_assignment(int word_id, int topic_id, int author_id, int doc_id);

  // return the probability that this word is assigned to this author and topic
  // Requires: c_wt, c_at, sum_c_wt, and sum_c_at all do not include the assignment for this word
  double neutral_sampling_prob(int author_id, int word_id, int topic_id);

  // return the probability that this opinion word is assigned to this author and topic
  // Requires: c_owt, c_at, sum_c_owt, and sum_c_at all do not include the assignment for this word
  double opinion_sampling_prob(int author_id, int word_id, int topic_id, int doc_id);

  // sample the distributions and update the assignment for this opinion word
  void sample_word(Word &w, int author_id, int doc_id, bool neutral);

  // THETA Data Structures

  // key: <author_id, topic_id>
  // value: # of assign
  unordered_map<key, int, myhash, myeq> c_at;  // count for author-topic assignments
  // key: author_id
  unordered_map<int, int> sum_c_at;  // sum of total author mentions


  // PHI Data Structures

  // key: <word_id, topic_id>
  // value: # of assignments
  unordered_map<key, int, myhash, myeq> c_wt;  // count for word-topic assignments
  // key: topic_id
  unordered_map<int, int> sum_c_wt;  // sum of total neutral topic mentions


  // RHO Data Structures

  // key: [author]<word_id, topic_id>
  // value: # of assignments
  unordered_map<int, unordered_map<key, int, myhash, myeq>> c_owt; // count for opinion word-topic assignments per author
  // key: [author]<topic_id>
  unordered_map<int, unordered_map<int, int>> sum_c_owt; // sum of total opinion topic mentions per author
  

  // DEBUG 

  // key: <doc_id, topic_id>
  // value: # of assignments
  unordered_map<key, int, myhash, myeq> c_dt;
  // key: doc_id
  unordered_map<int, int> sum_c_dt;

  // TODO: turn these vectors into maps

  // author names (index is used as author_id)
  vector<string> authors;
  // neutral words (index is word_id)
  vector<string> neutral_words;
  // opinion words
  vector<string> opinion_words;

  // vector of Documents, which contain the author info and Words of the document
  vector<Document> documents;  

  // model parameters
  double alpha;
  double beta;
  int T;
  int loop_count;
  int show_limit;
};

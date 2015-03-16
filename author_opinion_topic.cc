#include "common.hpp"
#include "progress.hpp"
#include <cassert>
#include <unordered_map>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

struct Word {
  int word_id;
  int topic;
};

struct Document {
  vector<Word> neutral_words;
  vector<Word> opinion_words;
  int author_id;
};

class AuthorTopic{
public:
  AuthorTopic(){
  }

  AuthorTopic(double a, double b, int t, int loop, int show_limit){
    this->alpha = a;
    this->beta = b;
    this->T = t;
    this->loop_count = loop;
    this->show_limit = show_limit;

    // set the random seed
    // DEBUG
    //srand(time(0));
    srand(0);
  }

  // Print the current state of the topic model. Use for debugging purposes
  void print_state() {

    /*
    // print the word-topic assignments in a table format
    for (int doc_id = 0; doc_id < documents.size(); ++doc_id) {

      Document &doc = documents[doc_id];
      cout << "Document " << doc_id << endl;

      cout << "Word:\t";
      for (int neutral_word_index = 0; neutral_word_index < doc.neutral_words.size(); ++neutral_word_index) {
        cout << neutral_words.at(doc.neutral_words[neutral_word_index].word_id) << "\t";
      }
      for (int opinion_word_index = 0; opinion_word_index < doc.opinion_words.size(); ++opinion_word_index) {
        cout << opinion_words.at(doc.opinion_words[opinion_word_index].word_id) << "\t";
      }

      cout << endl << "Topic:\t";
      for (int neutral_word_index = 0; neutral_word_index < doc.neutral_words.size(); ++neutral_word_index) {
        cout << doc.neutral_words[neutral_word_index].topic << "\t";
      }
      for (int opinion_word_index = 0; opinion_word_index < doc.opinion_words.size(); ++opinion_word_index) {
        cout << doc.opinion_words[opinion_word_index].topic << "\t";
      }
      cout << endl;
    }
    cout << endl;
    */

    // print out the word-topic counts (c_wt)
    cout << "c_wt" << endl;
    for (int neutral_id = 0; neutral_id < neutral_words.size(); ++neutral_id) {
      cout << "Word " << neutral_words[neutral_id];
      for (int topic_id = 0; topic_id < T; ++topic_id) {
        cout << ": topic " << topic_id << ", count " << c_wt[make_pair(neutral_id, topic_id)];
      }
      cout << endl;
    }
    cout << endl;
    // print sum_c_wt
    cout << "sum_c_wt" << endl;
    for (int topic_id = 0; topic_id < T; ++topic_id) {
      cout << "Topic " << topic_id << " assigned " << sum_c_wt[topic_id] << " times" << endl;
    }
    cout << endl;

    // print out the author-topic counts
    cout << "c_at" << endl;
    for (int author_id = 0; author_id < authors.size(); ++author_id) {
      cout << "Author " << authors.at(author_id);
      for (int topic_id = 0; topic_id < T; ++topic_id) {
        cout << ": topic " << topic_id << ", count " << c_at[make_pair(author_id, topic_id)];
      }
      cout << endl;
    }
    // print out total author mentions (sum_c_at)
    cout << "sum_c_at" << endl;
    for (int author_id = 0; author_id < authors.size(); ++author_id) {
      cout << "Author: " << authors[author_id] << " total mentions: " << sum_c_at[author_id] << endl;
    }
    cout << endl;

    // print out the author-opinion word topic counts
    cout << "c_owt" << endl;
    for (int author_id = 0; author_id < authors.size(); ++author_id) {
      cout << "Author " << authors[author_id] << endl;
      for (int opinion_word_id = 0; opinion_word_id < opinion_words.size(); ++opinion_word_id) {
        cout << "Opinion Word " << opinion_words.at(opinion_word_id);
        for (int topic_id = 0; topic_id < T; ++topic_id) { 
          cout << " topic " << topic_id << " count: " << c_owt[author_id][make_pair(opinion_word_id, topic_id)];
        }
        cout << endl;
      }
      cout << endl;
    }
    cout << "sum_c_owt" << endl;
    for (int author_id = 0; author_id < authors.size(); ++author_id) {
      cout << "Author " << authors[author_id] << " =>";
      for (int topic_id = 0; topic_id < T; ++topic_id) { 
        cout << " topic " << topic_id << ", opinion sum " << sum_c_owt[author_id][topic_id] << ":";
      }
      cout << endl;
    }
    cout << endl << endl;


  }

  int get_author_id(string author){
    for(int i = 0; i < authors.size(); ++i){
      if(authors.at(i) == author){
       return i;
      }
    }
    authors.push_back(author);
    return authors.size() - 1;
  }

  int get_neutral_word_id(string word){
    // find the index of the word in the vector
    for(int i = 0; i < neutral_words.size(); ++i){
      if(neutral_words.at(i) == word){
       return i;
      }
    }
    // if the word is not in the vector, then add it to the words vector
    neutral_words.push_back(word);
    return neutral_words.size() - 1;
  }

  int get_opinion_word_id(string word){
    // find the index of the word in the vector
    for(int i = 0; i < opinion_words.size(); ++i){
      if(opinion_words.at(i) == word){
       return i;
      }
    }
    // if the word is not in the vector, then add it to the words vector
    opinion_words.push_back(word);
    return opinion_words.size() - 1;
  }

  void add_neutral_topic_assignment(int word_id, int topic_id, int author_id) {
    // increment the count for this word_id/topic and author/topic pair
    ++c_wt[make_pair(word_id, topic_id)];
    ++c_at[make_pair(author_id, topic_id)];

    // increment the word count for this topic
    ++sum_c_wt[topic_id];
    // increment the word count for this author
    ++sum_c_at[author_id];
  }

  void add_opinion_topic_assignment(int word_id, int topic_id, int author_id, int doc_id) {
    // increment the count for this word_id/topic and author/topic pair
    ++c_owt[author_id][make_pair(word_id, topic_id)];
    ++c_at[make_pair(author_id, topic_id)];

    // increment the word count for this topic
    ++sum_c_owt[author_id][topic_id];
    // increment the word count for this author
    ++sum_c_at[author_id];

    //++sum_c_owt[1-author_id][topic_id];
    //++c_owt[1-author_id][make_pair(word_id, topic_id)];

    ++c_dt[make_pair(doc_id, topic_id)];
    ++sum_c_dt[doc_id];
  }

  void decrease_neutral_topic_assignment(int word_id, int topic_id, int author_id) {
    // increment the count for this word_id/topic and author/topic pair
    --c_wt[make_pair(word_id, topic_id)];
    --c_at[make_pair(author_id, topic_id)];

    // increment the word count for this topic
    --sum_c_wt[topic_id];
    // increment the word count for this author
    --sum_c_at[author_id];
  }

  void decrease_opinion_topic_assignment(int word_id, int topic_id, int author_id, int doc_id) {
    // increment the count for this word_id/topic and author/topic pair
    --c_owt[author_id][make_pair(word_id, topic_id)];
    --c_at[make_pair(author_id, topic_id)];

    // increment the word count for this topic
    --sum_c_owt[author_id][topic_id];
    // increment the word count for this author
    --sum_c_at[author_id];

    //--sum_c_owt[1-author_id][topic_id];
    //--c_owt[1-author_id][make_pair(word_id, topic_id)];
    --c_dt[make_pair(doc_id, topic_id)];
    --sum_c_dt[doc_id];
  }

  // Add a new document to the corpus, update all of the necessary data structures
  void add_document(string author, vector<string> neutral_words, vector<string> opinion_words){
    
    Document doc;
    int doc_id = documents.size();

    // set the author
    doc.author_id = get_author_id(author);

    // process the neutral word ids
    for(vector<string>::iterator i = neutral_words.begin(); i != neutral_words.end(); ++i){
      // add the word to the document, assign a random topic [0, t)
      Word w;
      w.word_id = get_neutral_word_id(*i);
      w.topic = rand() % T;

      add_neutral_topic_assignment(w.word_id, w.topic, doc.author_id);

      doc.neutral_words.push_back(w);
    }

    // process the opinion word ids
    for(vector<string>::iterator i = opinion_words.begin(); i != opinion_words.end(); ++i){
      // add the word to the document
      Word w;
      w.word_id = get_opinion_word_id(*i);

      // assign a random topic [0, t)
      w.topic = rand() % T;

      add_opinion_topic_assignment(w.word_id, w.topic, doc.author_id, doc_id);

      doc.opinion_words.push_back(w);
    }

    documents.push_back(doc);
  }

  // find the probability that this word is assigned to this author and topic
  // Requires: _c_wt, _c_at, _sum_c_wt, and _sum_c_at all do not include the assignment for this word
  double neutral_sampling_prob(int author_id, int word_id, int topic_id){

    int W = neutral_words.size();
    unordered_map<key, int, myhash, myeq>::iterator i;

    // find the number of assignments for this word-topic pair
    int c_wt_count = 0;
    if(c_wt.find(make_pair(word_id, topic_id)) != c_wt.end()){
      c_wt_count = c_wt[make_pair(word_id, topic_id)];
    }
    
    // find the number of assignments for this author-topic pair
    int c_at_count = 0;
    if(c_at.find(make_pair(author_id, topic_id)) != c_at.end()){
      c_at_count = c_at[make_pair(author_id, topic_id)];
    }
    // (wt_count + beta)/(t_sum + W*beta) * (at_count + alpha)/(author_sum + T*alpha)
    double prob = (c_wt_count + beta) / (sum_c_wt[topic_id] + (W * beta));
    prob *= (c_at_count + alpha) / (sum_c_at[author_id] + (T * alpha));

    return prob;
  }

  // find the probability that this opinion word is assigned to this author and topic
  // Requires: _c_wt, _c_at, _sum_c_wt, and _sum_c_at all do not include the assignment for this word
  double opinion_sampling_prob(int author_id, int word_id, int topic_id, int doc_id){
    
    int O = opinion_words.size();

    // find the number of assignments for this word-topic pair
    int c_owt_count = 0;
    if(c_owt[author_id].find(make_pair(word_id, topic_id)) != c_owt[author_id].end()){
      c_owt_count = c_owt[author_id][make_pair(word_id, topic_id)];
    }
    
    // find the number of assignments for this author-topic pair
    int c_at_count = 0;
    if(c_at.find(make_pair(author_id, topic_id)) != c_at.end()){
      c_at_count = c_at.at(make_pair(author_id, topic_id));
    }
    // (wt_count + beta)/(t_sum + W*beta) * (at_count + alpha)/(author_sum + T*alpha)
    double prob = (c_owt_count + beta) / (sum_c_owt.at(author_id).at(topic_id) + (O * beta));
    // DEBUG
    prob *= (c_at_count + alpha) / (sum_c_at.at(author_id) + (T * alpha));

    int c_dt_count = 0;
    if(c_dt.find(make_pair(doc_id, topic_id)) != c_dt.end()){
      c_dt_count = c_dt.at(make_pair(doc_id, topic_id));
    }

    //prob *= (c_dt_count + alpha) / (sum_c_dt.at(doc_id) + (T * alpha));
    return prob;
  }

  // sample the distributions and update the assignment for this opinion word
  void sample_word(Word &w, int author_id, int doc_id, bool neutral){

    // store the previous assignment for the topic
    int prev_topic = w.topic;

    // vector contains prob density over topics
    // image
    //  |------t_1-----|---t_2---|-t_3-|------t_4-----|
    // 0.0                ^^                         1.0
    vector<double> prob;

    // decrement the word-topic/topic assignment for the previous topic assignment
    if (neutral){
      decrease_neutral_topic_assignment(w.word_id, prev_topic, author_id);
    } else {
      decrease_opinion_topic_assignment(w.word_id, prev_topic, author_id, doc_id);
    } 
    // get the prob density over all topic combinations
    for(int topic_id = 0; topic_id < T; ++topic_id){
      double topic_prob = 0;
      if (neutral){
        topic_prob = neutral_sampling_prob(author_id, w.word_id, topic_id);
      } else {
        topic_prob = opinion_sampling_prob(author_id, w.word_id, topic_id, doc_id);
      }
      prob.push_back(topic_prob);

      // keep a cumulative distribution of probabilities
      if(prob.size() > 1){
        prob.at(prob.size() - 1) += prob.at(prob.size() - 2);
      }
    }
    
    /*
    if (!neutral) {
      if (opinion_words.at(w.word_id) == "0") {
        for(int topic_id = 0; topic_id < T; ++topic_id){
          cout << "Opinion " << opinion_words.at(w.word_id) << " topic " << topic_id << " probability " << prob.at(topic_id) << endl;
        }
        cout << endl;
      }
    }
    */
    
    // scaling [0,  1]
    double sum = prob.at(prob.size() - 1);
    for(int i = 0; i < prob.size(); ++i){
      prob.at(i) /= sum;
    }
    /*
    if (!neutral) {
      if (opinion_words.at(w.word_id) == "5" || opinion_words.at(w.word_id) == "3") {
        for(int topic_id = 0; topic_id < T; ++topic_id){
          cout << "Opinion " << opinion_words.at(w.word_id) << " topic " << topic_id << " probability " << prob.at(topic_id) << endl;
        }
        cout << endl;
      }
    }
    */
    
    // get a random uniform sample from the author/topic distribution
    double pos_prob = uniform_rand();
    /*
    if (!neutral) {
      if (opinion_words.at(w.word_id) == "5" || opinion_words.at(w.word_id) == "3") {
        cout << pos_prob << endl;
      }
    }
    */
    int new_topic = 0;
    if(pos_prob > prob.at(0)){
      for(int i = 1; i < prob.size(); ++i){
        if((pos_prob <= prob.at(i)) && (pos_prob > prob.at(i - 1))){
          new_topic = i;
          break;
        }
      }
    }

    // update the topic and author assignment
    w.topic = new_topic;

    if (!neutral) {
      if (w.word_id == 0 && new_topic == 0) {
        cout << "New topic: " << new_topic << ", prob: " << pos_prob << endl;
        for(int topic_id = 0; topic_id < T; ++topic_id){
          cout << "Opinion " << opinion_words.at(w.word_id) << " topic " << topic_id << " probability " << opinion_sampling_prob(author_id, w.word_id, topic_id, doc_id) << endl;
        }
        cout << "c_owt[0][0,0] = " << c_owt.at(author_id).at(make_pair(w.word_id, w.topic)) << endl;
        cout << "sum_c_owt[0][0] = " << sum_c_owt.at(author_id).at(w.topic) << endl;
        cout << "c_at[0, 0] = " << c_at.at(make_pair(author_id, w.topic)) << endl;
        cout << "sum_c_at[0] = " << sum_c_at.at(author_id) << endl;
        cout << endl;
      } 
    }
  
    // increment the word-topic and author-topic assignment counts
    if (neutral) {
      add_neutral_topic_assignment(w.word_id, w.topic, author_id);
    } else {
      add_opinion_topic_assignment(w.word_id, w.topic, author_id, doc_id);
    }
  }

  void sample_all(){  

    // initialize the progress bar
    boost::progress_display progress( loop_count * documents.size() );

    // do the prescribed number of iterations
    for(int i = 0; i < loop_count; ++i){

      // DEBUG
      //print_state();

      // loop through each document
      for(int pos_doc = 0; pos_doc < documents.size(); ++pos_doc){

        // loop through neutral each word and sample the topic assignment
        for(int pos_word = 0; pos_word < (documents.at(pos_doc)).neutral_words.size(); ++pos_word){
          // sample the word
          sample_word(documents.at(pos_doc).neutral_words.at(pos_word), documents.at(pos_doc).author_id, pos_doc, true);
        }
        // loop through each opinion word and sample the topic assignment
        for(int pos_word = 0; pos_word < (documents.at(pos_doc)).opinion_words.size(); ++pos_word){
          sample_word(documents.at(pos_doc).opinion_words.at(pos_word), documents.at(pos_doc).author_id, pos_doc, false);
        }
        ++progress;
      }
    }
  }

  void output(char* filename){

    // DEBUG - topic assignments
    /*
    for (Document d: documents) {
      for (Word neutral: d.neutral_words) {
        cout << neutral_words[neutral.word_id] << " topic " << neutral.topic << endl;
      }
      for (Word opinion: d.opinion_words) {
        cout << opinion_words[opinion.word_id] << " topic " << opinion.topic << endl;
      }
      cout << endl;
    }
    */

    // output theta
    ostringstream oss_theta;
    oss_theta << filename << "_theta.txt" ;
    ofstream ofs_theta;
    ofs_theta.open((oss_theta.str()).c_str());

    ofs_theta << "Author" << "\t" << "Topic" << "\t" << "Score" << endl;
    for(int author_id = 0; author_id < authors.size(); ++author_id){
      // ソート

      // key: score
      // value: topic_id
      vector<pair<double, int> > theta;

      int author_count_sum = sum_c_at[author_id];

      for(int topic_id = 0; topic_id < T; ++topic_id){

        int c_at_count = 0;
        if(c_at.find(make_pair(author_id, topic_id)) != c_at.end()){
          c_at_count = c_at[make_pair(author_id, topic_id)];
        }
        // score = (at_count + alpha)/(author_sum + T*alpha)
        double score = (c_at_count + alpha)/(author_count_sum + (T * alpha));
        theta.push_back(make_pair(score, topic_id));
        //all_theta[make_pair(author_id, topic_id)] = score;
      }

      // sort theta by topic score and output the show_limit topics for this author
      sort(theta.begin(), theta.end(), std::greater<pair<double, int> >());
      vector<pair<double, int> >::iterator j;
      int count = 0;
      for(j = theta.begin(); j != theta.end() && count < show_limit; ++j, ++count){
        // output: Author topic_id score \n
        ofs_theta << authors.at(author_id) << "\t" << (*j).second << "\t" << (*j).first << endl;
      }
    }
    ofs_theta.close();

    // output phi
    int W = neutral_words.size();
    ostringstream oss_phi;
    oss_phi << filename << "_phi.txt" ;
    ofstream ofs_phi;
    ofs_phi.open((oss_phi.str()).c_str());

    ofs_phi << "Topic" << "\t" << "Neutral Word" << "\t" << "Score" << endl;
    for(int topic_id = 0; topic_id < T; ++topic_id){
      
      // key: score
      // value: word
      vector<pair<double, string> > phi;

      for(int word_id = 0; word_id < W; ++word_id){
        int c_wt_count = 0;
        if(c_wt.find(make_pair(word_id, topic_id)) != c_wt.end()){
          c_wt_count = c_wt[make_pair(word_id, topic_id)];
        }
        // score = (wt_count + beta)/(topic_sum + W*beta)
        double score = (c_wt_count + beta)/(sum_c_wt[topic_id] + (W * beta));
        phi.push_back(make_pair(score, neutral_words.at(word_id)));
      }

      // sort phi from high score to low
      sort(phi.begin(), phi.end(), greater<pair<double, string> >());
      vector<pair<double, string> >::iterator j;
      int count = 0;
      for(j = phi.begin(); j != phi.end() && count < show_limit; ++j, ++count){
        // Output: Topic_id neutral_word score
        ofs_phi << topic_id << "\t" << (*j).second << "\t" << (*j).first << endl;
      }

      /*
      // sort all theta by topic/author
      vector<pair<double, string> > topic_given_theta;
      for(int author_id = 0; author_id < _authors.size(); ++author_id){
        double score = all_theta[make_pair(author_id, topic_id)];
        topic_given_theta.push_back(make_pair(score, _authors.at(author_id)));
      }

      sort(topic_given_theta.begin(), topic_given_theta.end());
      count = 0;
      for(j = topic_given_theta.rbegin(); j != topic_given_theta.rend(); ++j){
        if(count >= this->show_limit){
          break;
        }else{
          count++;
        }
      }
      */
    }
    // output rho
    ostringstream oss_rho;
    oss_rho << filename << "_rho.txt" ;
    ofstream ofs_rho;
    ofs_rho.open((oss_rho.str()).c_str());
    int opinion_size = opinion_words.size();

    // compute rho
    ofs_rho << "Author" << "\t" << "Topic" << "\t" << "Opinion Word" << "\t" << "Score" << endl;

    for (int author_id = 0; author_id < authors.size(); ++author_id) {
      for (int topic_id = 0; topic_id < T; ++topic_id) {

        int author_opinion_topic_count = sum_c_owt[author_id][topic_id];

        // key: score, value: opinion word
        vector<pair<double, string> > rho;
        for(int word_id = 0; word_id < opinion_size; ++word_id){
          int c_oat_count = 0;
          if(c_owt[author_id].find(make_pair(word_id, topic_id)) != c_owt[author_id].end()){
            c_oat_count = c_owt[author_id][make_pair(word_id, topic_id)];
          }
          // score = (owt_count + beta)/(author_topic_sum + O*beta)
          double score = (c_oat_count + beta)/(author_opinion_topic_count + (opinion_size * beta));
          rho.push_back(make_pair(score, opinion_words.at(word_id)));
        }

        // output rho (per author)
        sort(rho.begin(), rho.end());
        vector<pair<double, string> >::reverse_iterator j;
        int count = 0;
        for(j = rho.rbegin(); j != rho.rend() && count < show_limit; ++j, ++count){
          ofs_rho << author_id << "\t" << topic_id << "\t" << (*j).second << "\t" << (*j).first << endl;
        }
      }
    }

    ofs_phi.close();
  }
  
private:

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

// Usage: ./opinion_topic input.tsv t iterations output_lim output_path
// input.tsv: input file of format: author:..:author \t word:word...:word \t opinion:opinion:.. \n
// t: # topics
// iter: # of gibbs sampling iterations
// lim: limit of output
// output_path: path to the directory to 
int main(int argc, char** argv){

  assert(argc == 6);
  char *input_filename = argv[1];
  int topic_size = atoi(argv[2]);
  int num_iterations = atoi(argv[3]);
  int output_limit = atoi(argv[4]);
  char *output_filename = argv[5];

  // set the hyperparameters according to current standards
  double alpha = 50/topic_size;
  double beta = 0.01;

  AuthorTopic t(alpha, beta, topic_size, num_iterations, output_limit);
  
  // load the input documents into vectors
  ifstream ifs;
  ifs.open(input_filename, ios::in);
  string line;
  while(getline(ifs, line)){
    // file format
    // author \t neutral_1:neutral_2:... \t opinion_1:opinion_2:...
    vector<string> elem = split_string(line, "\t");
    vector<string> authors = split_string(elem.at(0), ":");
    assert(authors.size() == 1);
    string author = authors[0];
    vector<string> words = split_string(elem.at(1), ":");
    vector<string> opinion_words = split_string(elem.at(2), ":");

    t.add_document(author, words, opinion_words);

  }
  ifs.close();
  t.sample_all();
  t.print_state();
  t.output(output_filename);
}

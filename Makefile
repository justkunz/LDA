# Makefile

CXX=g++
CXX_FLAGS=-g -std=c++11 -I ~/Library/boost/

opinion_topic.o : author_opinion_topic.cc
	${CXX} ${CXX_FLAGS} -c author_opinion_topic.cc

common.o: common.hpp common.cc
	${CXX} ${CXX_FLAGS} -c common.cc

all: author_topic cpt labeled_lda

author_topic: 
	${CXX} ${CXX_FLAGS} author_topic.cc -o author_topic

opinion_topic: opinion_topic.o common.o
	${CXX} ${CXX_FLAGS} author_opinion_topic.o common.o -o opinion_topic

cpt: 
	${CXX} ${CXX_FLAGS} cpt.cc -o cpt

labeled_lda: 
	${CXX} ${CXX_FLAGS} labeled_lda.cc -o labeled_lda


clean: 
	rm author_topic cpt labeled_lda

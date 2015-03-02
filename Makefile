# Makefile

CXX=g++
CXX_FLAGS=-O3 -std=c++11 -I ~/Library/boost/

all: author_topic cpt labeled_lda

author_topic: 
	${CXX} ${CXX_FLAGS} ./src/author_topic.cc -o author_topic

cpt: 
	${CXX} ${CXX_FLAGS} ./src/cpt.cc -o cpt

labeled_lda: 
	${CXX} ${CXX_FLAGS} ./src/labeled_lda.cc -o labeled_lda


clean: 
	rm author_topic cpt labeled_lda

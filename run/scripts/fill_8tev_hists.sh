#!/bin/sh

rm ../results/latest.root

source COMPILE

./bin/AnalyzeHxx ../data/data/HZZ4lTree_DoubleEle.root ../results/latest.root ../results data 100
./bin/AnalyzeHxx ../data/data/HZZ4lTree_DoubleMu.root  ../results/latest.root ../results data 101
./bin/AnalyzeHxx ../data/data/HZZ4lTree_DoubleOr.root  ../results/latest.root ../results data 102
#./bin/AnalyzeHxx ../data/current_data.root ../results/latest.root ../results

# ls ../data/data/
# HZZ4lTree_DoubleEle.root  HZZ4lTree_DoubleMu.root  HZZ4lTree_DoubleOr.root

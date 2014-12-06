#!/bin/sh

rm ../results/latest.root

source COMPILE

# Data
./bin/AnalyzeHxx ../data/data/HZZ4lTree_DoubleEle.root ../results/latest.root ../results data 0
./bin/AnalyzeHxx ../data/data/HZZ4lTree_DoubleMu.root  ../results/latest.root ../results data 0
./bin/AnalyzeHxx ../data/data/HZZ4lTree_DoubleOr.root  ../results/latest.root ../results data 0

# Background
./bin/AnalyzeHxx ../data/2e2mu/HZZ4lTree_ZZJetsTo4L.root  ../results/latest.root ../results zzjets 20
./bin/AnalyzeHxx ../data/4e/HZZ4lTree_ZZJetsTo4L.root     ../results/latest.root ../results zzjets 20
./bin/AnalyzeHxx ../data/4mu/HZZ4lTree_ZZJetsTo4L.root    ../results/latest.root ../results zzjets 20
./bin/AnalyzeHxx ../data/2e2mu/HZZ4lTree_ZZTo2e2mu.root   ../results/latest.root ../results 2e2mu 21
./bin/AnalyzeHxx ../data/4e/HZZ4lTree_ZZTo2e2mu.root      ../results/latest.root ../results 2e2mu 21
./bin/AnalyzeHxx ../data/4mu/HZZ4lTree_ZZTo2e2mu.root     ../results/latest.root ../results 2e2mu 21
./bin/AnalyzeHxx ../data/2e2mu/HZZ4lTree_ZZTo2e2tau.root  ../results/latest.root ../results 2e2tau 22
./bin/AnalyzeHxx ../data/4e/HZZ4lTree_ZZTo2e2tau.root     ../results/latest.root ../results 2e2tau 22
./bin/AnalyzeHxx ../data/4mu/HZZ4lTree_ZZTo2e2tau.root    ../results/latest.root ../results 2e2tau 22
./bin/AnalyzeHxx ../data/2e2mu/HZZ4lTree_ZZTo2mu2tau.root ../results/latest.root ../results 2mu2tau 23
./bin/AnalyzeHxx ../data/4e/HZZ4lTree_ZZTo2mu2tau.root    ../results/latest.root ../results 2mu2tau 23
./bin/AnalyzeHxx ../data/4mu/HZZ4lTree_ZZTo2mu2tau.root   ../results/latest.root ../results 2mu2tau 23
./bin/AnalyzeHxx ../data/2e2mu/HZZ4lTree_ZZTo4e.root      ../results/latest.root ../results 4e 24
./bin/AnalyzeHxx ../data/4e/HZZ4lTree_ZZTo4e.root         ../results/latest.root ../results 4e 24
./bin/AnalyzeHxx ../data/4mu/HZZ4lTree_ZZTo4e.root        ../results/latest.root ../results 4e 24
./bin/AnalyzeHxx ../data/2e2mu/HZZ4lTree_ZZTo4mu.root     ../results/latest.root ../results 4mu 25
./bin/AnalyzeHxx ../data/4e/HZZ4lTree_ZZTo4mu.root        ../results/latest.root ../results 4mu 25
./bin/AnalyzeHxx ../data/4mu/HZZ4lTree_ZZTo4mu.root       ../results/latest.root ../results 4mu 25
./bin/AnalyzeHxx ../data/2e2mu/HZZ4lTree_ZZTo4tau.root    ../results/latest.root ../results 4tau 26
./bin/AnalyzeHxx ../data/4e/HZZ4lTree_ZZTo4tau.root       ../results/latest.root ../results 4tau 26
./bin/AnalyzeHxx ../data/4mu/HZZ4lTree_ZZTo4tau.root      ../results/latest.root ../results 4tau 26

# Signal It's probably easier to add these in from hxx at plotting stage
#./bin/AnalyzeHxx ../data/signal/hxx_8TeV_1GeV.root    ../results/latest.root ../results hxx1 100
#./bin/AnalyzeHxx ../data/signal/hxx_8TeV_10GeV.root   ../results/latest.root ../results hxx10 101
#./bin/AnalyzeHxx ../data/signal/hxx_8TeV_100GeV.root  ../results/latest.root ../results hxx100 102
#./bin/AnalyzeHxx ../data/signal/hxx_8TeV_500GeV.root  ../results/latest.root ../results hxx500 103
#./bin/AnalyzeHxx ../data/signal/hxx_8TeV_1000GeV.root ../results/latest.root ../results hxx1000 104


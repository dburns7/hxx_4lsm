#!/bin/sh

source COMPILE

#./bin/AnalyzeHxx --8tev --fake_rate=1E-5 --num_smear=100 --met_smear=20 ../data/Analysis/CMSSW/DAS/8TeV/all.root ../results/8TeV/latest.root ../results/8TeV/

#./bin/AnalyzeHxx --8tev --fake_rate=1E-5 --num_smear=100 --met_smear=20 ../data/Analysis/CMSSW/DAS_100k/8TeV/all.root ../results/8TeV/latest_100k.root ../results/8TeV/

#./bin/AnalyzeHxx --8tev --fake_rate=1E-5 --num_smear=100 --met_smear=20 ../data/Analysis/CMSSW/DAS_Full/8TeV/all.root ../results/8TeV/latest_Full.root ../results/8TeV/

#./bin/AnalyzeHxx --8tev --fake_rate=1E-5 --num_smear=100 --met_smear=20 ../data/Analysis/CMSSW/DAS_Full/8TeV/all.root ../results/8TeV/latest_Full_Combo.root ../results/8TeV/

./bin/AnalyzeHxx --8tev --fake_rate=1E-5 --num_smear=100 --met_smear=20 ../data/Analysis/CMSSW/MCValidation/all.root ../results/8TeV/latest_MCValidation_4emuta.root ../results/8TeV/

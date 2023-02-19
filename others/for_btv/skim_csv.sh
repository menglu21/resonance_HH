#!/bin/bash

jessource=("AbsoluteMPFBias" "AbsoluteScale" "AbsoluteStat" "RelativeBal" "RelativeFSR" "RelativeJEREC1" "RelativeJEREC2" "RelativeJERHF" "RelativePtBB" "RelativePtEC1" "RelativePtEC2" "RelativePtHF" "RelativeStatEC" "RelativeStatFSR" "RelativeStatHF" "PileUpDataMC" "PileUpPtBB" "PileUpPtEC1" "PileUpPtEC2" "PileUpPtHF" "PileUpPtRef" "FlavorQCD" "Fragmentation" "SinglePionECAL" "SinglePionHCAL" "TimePtEta")

cp DeepJet_106XUL16postVFPSF_v3.csv DeepJet_106XUL16postVFPSF_v3_skimmed.csv
cp DeepJet_106XUL16preVFPSF_v2.csv DeepJet_106XUL16preVFPSF_v2_skimmed.csv

for ic in ${jessource[@]}
do
#  printf "%s\n" $ic
  sed -i "/$ic/d" DeepJet_106XUL16postVFPSF_v3_skimmed.csv
  sed -i "/$ic/d" DeepJet_106XUL16preVFPSF_v2_skimmed.csv
done

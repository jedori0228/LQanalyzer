#!/bin/bash

###########################
# run trilepton_mumumu.sh #
###########################
#runMC=$1
#runSignal=$2
#runDoubleMuon=$3
#runFRMethod=$4
#whichFR=$5
###########################

#1) Signal
#source ./trilepton_mumumu.sh true true false false '' 'runCutOptNtuple'
#2) Bkg
source ./trilepton_mumumu.sh true false false false '' 'runCutOptNtuple'
#3) Data
#source ./trilepton_mumumu.sh false false true false '' 'runCutOptNtuple'
#4) FR_weighted
source ./trilepton_mumumu.sh false false true true 'dijet_topology' 'runCutOptNtuple'
source ./trilepton_mumumu.sh false false true true 'HighdXY' 'runCutOptNtuple'
source ./trilepton_mumumu.sh false false true true 'DiMuon_HighdXY' 'runCutOptNtuple'
source ./trilepton_mumumu.sh false false true true 'DiMuon_HighdXY_n_jets' 'runCutOptNtuple'





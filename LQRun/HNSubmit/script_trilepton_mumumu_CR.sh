#!/bin/bash

##############################
# run trilepton_mumumu_CR.sh #
##############################
# runMC=$1
# runSignal=$2
# runDoubleMuon=$3
# runFRMethod=$4
# whichFR=$5
# whichCR=$6
##############################

## CR1 ##
#1) Signal
#source ./trilepton_mumumu_CR.sh true true false false '' 'CR1'
#2) Bkg
source ./trilepton_mumumu_CR.sh true false false false '' 'CR1'
#3) Data
#source ./trilepton_mumumu_CR.sh false false true false '' 'CR1'
#4) FR_weighted
source ./trilepton_mumumu_CR.sh false false true true 'dijet_topology' 'CR1'
source ./trilepton_mumumu_CR.sh false false true true 'HighdXY' 'CR1'
source ./trilepton_mumumu_CR.sh false false true true 'DiMuon_HighdXY' 'CR1'
source ./trilepton_mumumu_CR.sh false false true true 'DiMuon_HighdXY_n_jets' 'CR1'

## CR2 ##
#1) Signal
#source ./trilepton_mumumu_CR.sh true true false false '' 'CR2'
#2) Bkg
source ./trilepton_mumumu_CR.sh true false false false '' 'CR2'
#3) Data
#source ./trilepton_mumumu_CR.sh false false true false '' 'CR2'
#4) FR_weighted
source ./trilepton_mumumu_CR.sh false false true true 'dijet_topology' 'CR2'
source ./trilepton_mumumu_CR.sh false false true true 'HighdXY' 'CR2'
source ./trilepton_mumumu_CR.sh false false true true 'DiMuon_HighdXY' 'CR2'
source ./trilepton_mumumu_CR.sh false false true true 'DiMuon_HighdXY_n_jets' 'CR2'




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
source ./trilepton_mumumu.sh true true false false ''
#2) Bkg
source ./trilepton_mumumu.sh true false false false ''
#3) Data
source ./trilepton_mumumu.sh false false true false ''
#4) FR_weighted
source ./trilepton_mumumu.sh false false true true 'Dijet'
source ./trilepton_mumumu.sh false false true true 'HighdXY'
source ./trilepton_mumumu.sh false false true true 'DiMuonHighdXY'
source ./trilepton_mumumu.sh false false true true 'DiMuonHighdXYnjets'

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
source ./trilepton_mumumu_CR.sh true true false false '' 'CR1'
#2) Bkg
source ./trilepton_mumumu_CR.sh true false false false '' 'CR1'
#3) Data
source ./trilepton_mumumu_CR.sh false false true false '' 'CR1'
#4) FR_weighted
source ./trilepton_mumumu_CR.sh false false true true 'Dijet' 'CR1'
source ./trilepton_mumumu_CR.sh false false true true 'HighdXY' 'CR1'
source ./trilepton_mumumu_CR.sh false false true true 'DiMuonHighdXY' 'CR1'
source ./trilepton_mumumu_CR.sh false false true true 'DiMuonHighdXYnjets' 'CR1'

## CR2 ##
#1) Signal
source ./trilepton_mumumu_CR.sh true true false false '' 'CR2'
#2) Bkg
source ./trilepton_mumumu_CR.sh true false false false '' 'CR2'
#3) Data
source ./trilepton_mumumu_CR.sh false false true false '' 'CR2'
#4) FR_weighted
source ./trilepton_mumumu_CR.sh false false true true 'Dijet' 'CR2'
source ./trilepton_mumumu_CR.sh false false true true 'HighdXY' 'CR2'
source ./trilepton_mumumu_CR.sh false false true true 'DiMuonHighdXY' 'CR2'
source ./trilepton_mumumu_CR.sh false false true true 'DiMuonHighdXYnjets' 'CR2'

######################################
# run FakeRateCalculator_Musubmit.sh #
######################################
# runMC=$1
# runData=$2
# runDataLowPt=$3
######################################

#1) MC
source ./trilepton_mumumu_CR.sh true false false
#2) Data
source ./trilepton_mumumu_CR.sh false true false
#3) Data_low_pt
source ./trilepton_mumumu_CR.sh false true true










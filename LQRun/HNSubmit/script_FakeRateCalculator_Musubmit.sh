#!/bin/bash

######################################
# run FakeRateCalculator_Musubmit.sh #
######################################
# runMC=$1
# runData=$2
# runDataLowPt=$3
######################################

#1) MC
source ./FakeRateCalculator_Musubmit.sh true false false
#2) Data
source ./FakeRateCalculator_Musubmit.sh false true false
#3) Data_low_pt
source ./FakeRateCalculator_Musubmit.sh false true true










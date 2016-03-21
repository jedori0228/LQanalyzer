#!/bin/sh

###### SET WHAT JOBS TO RUN
runMC=false
runDoubleMuon=true

if [[ $runMC  == "true" ]]; 
then
    source functions.sh
    cycle="FakeRateCalculator_Mu"
    skinput="True"
#    useskim="NoCut"
    outputdir=$LQANALYZER_DIR"/data/output/Muon/"
    #### JOB CONFIGURATION
    njobs=1
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    nevents=100000
 #   declare -a input_samples=("DY10to50")
# "DY50plus" "ttbar" "Wjets" "WZ" "ZZ" "WW" "QCD_mumu")
#		declare -a input_samples=("HN100_new")
    #declare -a input_samples=("ttbar_central")
    source submit.sh $1
fi
    

################ DOUBLEMUON DATA
if [[ $runDoubleMuon  == "true" ]];
then
    source functions.sh
    cycle="FakeRateCalculator_Mu"
    skinput="True"
    stream="muon"
    useskim="DiLep"
    outputdir=$LQANALYZER_DIR"/data/output/Muon/"
    #### JOB CONFIGURATION
    njobs=30
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    nevents=100000
    declare -a input_samples=("C" "D")
    source submit.sh $1
fi     





echo ""
echo "End of example_submit.sh script."

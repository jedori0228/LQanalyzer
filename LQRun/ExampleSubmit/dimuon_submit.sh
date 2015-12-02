#!/bin/sh

###### SET WHAT JOBS TO RUN
runMC=true
runDoubleMuon=false

if [[ $runMC  == "true" ]]; 
then
    source functions.sh
    cycle="ExampleAnalyzerDiMuon"
    skinput="True"
    useskim=""
    outputdir=$LQANALYZER_DIR"/data/output/Muon/"
    #### JOB CONFIGURATION
    njobs=1
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    nevents=-1
    #declare -a input_samples=("DY10to50" "DY50plus" "ttbar" "Wjets" "WZ" "ZZ" "WW" "QCD_mumu")
    declare -a input_samples=("HN200_mumumu_new")
    source submit.sh $1
fi
    

################ DOUBLEMUON DATA
if [[ $runDoubleMuon  == "true" ]];
then
    source functions.sh
    cycle="ExampleAnalyzerDiMuon"
    skinput="True"
    stream="muon"
    useskim="DiLep"
    outputdir=$LQANALYZER_DIR"/data/output/Muon/"
    #### JOB CONFIGURATION
    njobs=1
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    nevents=100000
    declare -a input_samples=("A")
    source submit.sh $1
fi     





echo ""
echo "End of example_submit.sh script."

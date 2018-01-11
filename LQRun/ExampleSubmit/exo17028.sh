#!/bin/sh

###### SET WHAT JOBS TO RUN
runMC=true
runDoubleMuon=false

if [[ $runMC  == "true" ]]; 
then
    source functions.sh
    cycle="ExampleAnalyzer"
    skinput="True"
    outputdir=/data4/LQAnalyzerCode/jskim/LQ8TeVBatch/
    #### JOB CONFIGURATION
    njobs=50
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    usebatch="True"
    queue="middleq"
    #declare -a input_samples=("DY10to50" "DY50plus" "ttbar" "Wjets" "WZ" "ZZ" "WW" "QCD_mumu")
    declare -a input_samples=("ttbar")
    source submit.sh $1
fi
    

################ DOUBLEMUON DATA
if [[ $runDoubleMuon  == "true" ]];
then
    source functions.sh
    cycle="ExampleAnalyzerDiMuon"
    skinput="True"
    stream="singlemuon"
#    useskim="DiLep"
    outputdir=$LQANALYZER_DIR"/data/output/Muon/"
    #### JOB CONFIGURATION
    njobs=30
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    declare -a input_samples=("A" "B" "C" "D")
    source submit.sh $1
fi     





echo ""
echo "End of example_submit.sh script."

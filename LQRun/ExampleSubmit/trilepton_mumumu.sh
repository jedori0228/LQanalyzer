#!/bin/sh

###### SET WHAT JOBS TO RUN
runMC=false
runDoubleMuon=true

if [[ $runMC  == "true" ]]; 
then
    source functions.sh
    cycle="trilepton_mumumu"
    skinput="True"
#    useskim="NoCut"
    outputdir=$LQANALYZER_DIR"/data/output/trilepton/mumumu/LLL/"
    #outputdir=$LQANALYZER_DIR"/"
    #### JOB CONFIGURATION
    njobs=30
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    nevents=-1
    #declare -a input_samples=("HN40_mumumu_new" "HN50_mumumu_new" "HN60_mumumu_new" "HN90_mumumu_new" "HN100_mumumu_new" "HN150_mumumu_new" "HN200_mumumu_new" "HN300_mumumu_new" "HN400_mumumu_new" "HN500_mumumu_new" "HN700_mumumu_new" "HN1000_mumumu_new")
    #declare -a input_samples=("HN40_mumumu_new")
    declare -a input_samples=("DY10to50" "DY50plus" "topDIL" "Wbb" "ttW" "ttZ" "WWW" "TTWW" "TTG" "ZZZ" "WZZ" "WWZ" "WWG" "WW_mg" "WZtollqq_mg" "WZtoqqln_mg" "WZtollln_mg" "ZZtollnn_mg" "ZZtollqq_mg" "ZZtollll_mg" "Zbb" "HtoWW" "ggHtoZZ" "Wtollln_new" "ttbar")
    source submit.sh $1
fi
    

################ DOUBLEMUON DATA
if [[ $runDoubleMuon  == "true" ]];
then
    source functions.sh
    cycle="trilepton_mumumu"
    skinput="True"
    stream="muon"
    useskim="DiLep"
    outputdir=$LQANALYZER_DIR"/data/output/trilepton/mumumu/TTL/weighted/"
    #outputdir=$LQANALYZER_DIR"/"
    #### JOB CONFIGURATION
    njobs=30
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    nevents=-1
    declare -a input_samples=("A" "B" "C" "D")
    #declare -a input_samples=("A")
    source submit.sh $1
fi     





echo ""
echo "End of example_submit.sh script."

#!/bin/sh

###### SET WHAT JOBS TO RUN
runMC=true
runTrileptonMuMuE=false

if [[ $runMC  == "true" ]]; 
then
    source functions.sh
    cycle="trilepton_mumue"
    skinput="True"
    useskim=""
    outputdir=$LQANALYZER_DIR"/data/output/trilepton/mumue/"
    #### JOB CONFIGURATION
    njobs=30
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    
    #declare -a input_samples=("HN40_lep" "Wtollln")
		declare -a input_samples=("Wtollln_new")    
    #declare -a input_samples=("WWZ")
		#declare -a input_samples=("DY10to50" "DY50plus" "topDIL" "Wbb" "ttW" "ttZ" "WWW" "TTWW" "TTG" "ZZZ" "WZZ" "WWZ" "WWG" "WW_mg" "WZtollqq_mg" "WZtoqqln_mg" "WZtollln_mg" "ZZtollnn_mg" "ZZtollqq_mg" "ZZtollll_mg" "Zbb" "HtoWW")
    source submit.sh $1
fi
    

################ DOUBLEMUON DATA
if [[ $runTrileptonMuMuE  == "true" ]];
then
    source functions.sh
    cycle="trilepton_mumue"
    skinput="True"
    stream="emu"
    useskim="DiLep"
    outputdir=$LQANALYZER_DIR"/data/output/trilepton/mumue/"
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

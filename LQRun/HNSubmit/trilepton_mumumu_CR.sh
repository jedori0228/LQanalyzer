#!/bin/sh

###### SET WHAT JOBS TO RUN
runMC=false
runDoubleMuon=true
runFRMethod=true

if [[ $runMC  == "true" ]]; 
then
    source functions.sh
    cycle="trilepton_mumumu_CR"
    skinput="True"
#    useskim="NoCut"
    outputdir=$LQANALYZER_DIR"/data/output/trilepton_mumumu_CR/"
    #outputdir=$LQANALYZER_DIR"/"
    #### JOB CONFIGURATION
    njobs=10
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    nevents=-1
    declare -a input_samples=("HN40_mumumu_new" "HN50_mumumu_new" "HN60_mumumu_new" "HN90_mumumu_new" "HN100_mumumu_new" "HN150_mumumu_new" "HN200_mumumu_new" "HN300_mumumu_new" "HN400_mumumu_new" "HN500_mumumu_new" "HN700_mumumu_new" "HN1000_mumumu_new")
    #declare -a input_samples=("DY10to50" "DY50plus" "Wjets" "WW_mg" "WZtollqq_mg" "WZtoqqln_mg" "WZtollln_mg" "ZZtollnn_mg" "ZZtollqq_mg" "ZZtollll_mg" "topDIL")
    #declare -a input_samples=("ggHtoZZ" "WZtollln_mg" "ZZtollll_mg" "Wtollln_new")
    #declare -a input_samples=("Wtollln_new")
    #declare -a input_samples=("DY10to50" "DY50plus" "topDIL" "Wbb" "ttW" "ttZ" "WWW" "TTWW" "TTG" "ZZZ" "WZZ" "WWZ" "WWG" "WW_mg" "WZtollqq_mg" "WZtoqqln_mg" "WZtollln_mg" "ZZtollnn_mg" "ZZtollqq_mg" "ZZtollll_mg" "Zbb" "HtoWW" "ggHtoZZ" "Wtollln_new" "ttbar" "Wjets")

    #declare -a input_samples=("WW_mg" "WZtollqq_mg" "WZtoqqln_mg" "WZtollln_mg" "ZZtollnn_mg" "ZZtollqq_mg" "ZZtollll_mg" "Zbb" "HtoWW" "ggHtoZZ" "Wtollln_new") #cms5
    #declare -a input_samples=("Zbb") #cms6

    source submit.sh $1
fi
    

################ DOUBLEMUON DATA
if [[ $runDoubleMuon  == "true" ]];
then
    source functions.sh
    cycle="trilepton_mumumu_CR"
    outputdir=$LQANALYZER_DIR"/data/output/trilepton_mumumu_CR/period/"
    if [[ $runFRMethod == "true" ]];
    then
        cycle="trilepton_mumumu_CR_FR_method"
        outputdir=$LQANALYZER_DIR"/data/output/trilepton_mumumu_CR/FR_weighted/HighdXY/"
        #outputdir=$LQANALYZER_DIR"/"
    fi
    skinput="True"
    stream="muon"
    useskim="DiLep"
    #useskim="Lepton"
    #### JOB CONFIGURATION
    njobs=100
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

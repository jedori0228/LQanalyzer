#!/bin/sh

###### SET WHAT JOBS TO RUN
runMC=$1
runSignal=$2
runDoubleMuon=$3
runFRMethod=$4
whichFR=$5
runCutOpt=$6

if [[ $runMC  == "true" ]]; 
then
    source functions.sh
    cycle="trilepton_mumumu"
    outdirname="trilepton_mumumu"
    skinput="True"
    #### JOB CONFIGURATION
    njobs=100
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    nevents=-1
    usebatch="True"
    declare -a input_samples=("DY10to50" "DY50plus" "topDIL" "Wbb" "ttW" "ttZ" "WWW" "TTWW" "TTG" "ZZZ" "WZZ" "WWZ" "WWG" "WW_mg" "WZtollqq_mg" "WZtoqqln_mg" "WZtollln_mg" "ZZtollnn_mg" "ZZtollqq_mg" "ZZtollll_mg" "Zbb" "HtoWW" "ggHtoZZ" "Wtollln_new" "ttbar" "Wjets")

    if [[ $runSignal == "true" ]];
    then
      declare -a input_samples=("HN40_mumumu_new" "HN50_mumumu_new" "HN60_mumumu_new" "HN90_mumumu_new" "HN100_mumumu_new" "HN150_mumumu_new" "HN200_mumumu_new" "HN300_mumumu_new" "HN400_mumumu_new" "HN500_mumumu_new" "HN700_mumumu_new" "HN1000_mumumu_new")
    fi

    if [[ $runCutOpt == "runCutOptNtuple" ]];
    then
      k_jskim_flag_2="runCutOptNtuple"
      outdirname="trilepton_mumumu_cutop"
    fi

    outputdir=$LQANALYZER_DIR"/data/output/"$outdirname"/"

    #### for debugging
    declare -a input_samples=("DY10to50")
    #outputdir=$LQANALYZER_DIR"/"

    source submit.sh
fi
    

################ DOUBLEMUON DATA
if [[ $runDoubleMuon  == "true" ]];
then
    source functions.sh
    cycle="trilepton_mumumu"
    outdirname="trilepton_mumumu/period"
    jskimflag1=$whichFR
    if [[ $runFRMethod == "true" ]];
    then
      cycle="trilepton_mumumu_FR_method"
      outdirname="trilepton_mumumu/FR_weighted/"$whichFR
      if [[ $runCutOpt == "runCutOptNtuple" ]];
      then
        k_jskim_flag_2="runCutOptNtuple"
        outdirname="trilepton_mumumu_cutop/FR_weighted/"$whichFR
      fi
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
    usebatch="True"
    declare -a input_samples=("A" "B" "C" "D")

    outputdir=$LQANALYZER_DIR"/data/output/"$outdirname"/"

    #### for debugging
    #declare -a input_samples=("A")
    #outputdir=$LQANALYZER_DIR"/"

    #### 18 GeV check
    #declare -a input_samples=("A" "B" "C" "D")
    #outputdir=$LQANALYZER_DIR"/data/output/18GeVmumumumu/"
    #skinput="False"

    source submit.sh
fi     





echo ""
echo "End of example_submit.sh script."

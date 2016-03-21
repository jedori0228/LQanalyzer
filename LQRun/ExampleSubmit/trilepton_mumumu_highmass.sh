#!/bin/sh

###### SET WHAT JOBS TO RUN
runMC=true
runDoubleMuon=true

if [[ $runMC  == "true" ]]; 
then
    source functions.sh
    cycle="trilepton_mumumu_high_mass"
    skinput="True"
    useskim=""
    outputdir=$LQANALYZER_DIR"/data/output/trilepton/mumumu/"
    #### JOB CONFIGURATION
    njobs=30
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    nevents=-1
    declare -a input_samples=("Wbb" "WWW" "ZZZ" "WZZ" "WWZ" "WW_mg" "WZtollqq_mg" "WZtoqqln_mg" "WZtollln_mg" "ZZtollnn_mg" "ZZtollqq_mg" "ZZtollll_mg" "Zbb" "Wtollln" "DY10to50" "DY50plus" "topDIL" "ttW" "ttZ" "TTWW" "TTG" "WWG" "HtoWW" "stbar_sch" "stbar_tch" "stbar_tW" "st_sch" "st_tch" "st_tW" "topHAD" "topLJ" "ggHtoZZ" "HtoZZ" "HtoTauTau" "QCD_mumu" "WgammaMu")
    #declare -a input_samples=("DY10to50" "DY50plus" "topDIL" "Wbb" "ttW" "ttZ" "WWW" "TTWW" "TTG" "ZZZ" "WZZ" "WWZ" "WWG" "WW_mg" "WZtollqq_mg" "WZtoqqln_mg" "WZtollln_mg" "ZZtollnn_mg" "ZZtollqq_mg" "ZZtollll_mg" "Zbb" "HtoWW" "Wtollln_new" "HN40_mumumu_new" "HN50_mumumu_new" "HN60_mumumu_new")
    source submit.sh $1
fi
    

################ DOUBLEMUON DATA
if [[ $runDoubleMuon  == "true" ]];
then
    source functions.sh
    cycle="trilepton_mumumu_high_mass"
    skinput="True"
    stream="muon"
    useskim=""
    outputdir=$LQANALYZER_DIR"/data/output/trilepton/mumumu/"
    #### JOB CONFIGURATION
    njobs=30
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    nevents=-1
    declare -a input_samples=("A" "B" "C" "D")
    source submit.sh $1
fi     





echo ""
echo "End of example_submit.sh script."

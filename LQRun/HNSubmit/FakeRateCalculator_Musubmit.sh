#!/bin/sh

runMC=$1
runData=$2
runDataLowPt=$3

if [[ $runMC  == "true" ]];
then
    source functions.sh
    cycle="FakeRateCalculator_Mu"
    skinput="True"
    usebatch="True"
    njobs=100
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    outputdir=$LQANALYZER_DIR"/data/output/MuonFakes/dXY_0p01_dZ_0p5/"
    #outputdir=$LQANALYZER_DIR"/"

    declare -a input_samples=("DY10to50" "DY50plus" "Wjets" "Wgamma" "stbar_sch" "stbar_tch" "stbar_tW" "st_sch" "st_tch" "st_tW" "ttbarMS" "Wbb" "ttW" "ttZ" "WWW" "TTWW" "TTG" "ZZZ" "WZZ" "WWZ" "WWG" "WW_mg" "WZtollqq_mg" "WZtoqqln_mg" "WZtollln_mg" "ZZtollnn_mg" "ZZtollqq_mg" "ZZtollll_mg" "Zbb" "HtoWW" "ggHtoZZ" "Wtollln_new" "QCD_mumu")
    #declare -a input_samples=("QCD_1000_mu" "QCD_15-20_mu" "QCD_20-30_mu" "QCD_30-50_mu" "QCD_50-80_mu" "QCD_800-1000_mu" "QCD_120-170_mu" "QCD_170-300_mu" "QCD_300-470_mu" "QCD_470-600_mu" "QCD_600-800_mu" "QCD_80-120_mu")
    #declare -a input_samples=("DY50plus")

    source submit.sh  

fi

if [[ $runData  == "true" ]];
then
    source functions.sh
    cycle="FakeRateCalculator_Mu"
    skinput="True"
    #skinput="False"

    njobs=100
    usebatch="True"
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    #usebatch="False"

    nevents=-1

    #### for pT > 15 GeV bin ####
    outputdir=$LQANALYZER_DIR"/data/output/MuonFakes/dXY_0p01_dZ_0p5/period/"
    #outputdir=$LQANALYZER_DIR"/"
    stream="muon"
    #stream="singlemuon"
    declare -a input_samples=("A" "B" "C" "D")

    #### for 10-15 GeV pT bin ####
    if [[ $runDataLowPt == "true" ]];
    then
        outputdir=$LQANALYZER_DIR"/data/output/MuonFakes/dXY_0p01_dZ_0p5/"
        stream="muon_lowpt"
        declare -a input_samples=("D")
    fi

    source submit.sh $1

fi



echo ""
echo "End of example_submit.sh script."

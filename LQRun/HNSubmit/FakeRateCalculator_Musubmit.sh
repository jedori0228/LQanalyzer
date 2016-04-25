#!/bin/sh

runData=true
runMC=false
runQCD=false
runQCD2=false


if [[ $1  == "MC"  ]];
then
    runMC=true
fi

if [[ $1  == "QCD"  ]];
then
    runQCD=true
fi


if [[ $1  == "QCD2"  ]];
then
    runQCD2=true
fi



if [[ $1  == "DATA"  ]];
then

    runData=true
fi


if [[ $runMC  == "true" ]];
then
    source functions.sh
    cycle="FakeRateCalculator_Mu"
    skinput="True"

    njobs=5
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    #outputdir=$LQANALYZER_DIR"/data/output/MuonFakes/MCTruth/QCD_single_mu_enriched/"
    #outputdir=$LQANALYZER_DIR"/data/output/MuonFakes/dijet_topology/dXY_0p01_dZ_0p5/"
    outputdir=$LQANALYZER_DIR"/data/output/MuonFakes/Large_dXY/"

    #declare -a input_samples=("DY10to50" "DY50plus" "Wjets" "Wgamma" "stbar_sch" "stbar_tch" "stbar_tW" "st_sch" "st_tch" "st_tW" "ttbarMS")

    #declare -a input_samples=("stbar_sch" "stbar_tch" "stbar_tW" "st_sch" "st_tch" "st_tW") #cms5
    declare -a input_samples=("ttbarMS") #cms6

    #declare -a input_samples=("ttbarMSpow" "ttbarMS" "ttbarMS_chs")
    #declare -a input_samples=("QCD_1000_mu" "QCD_15-20_mu" "QCD_20-30_mu" "QCD_30-50_mu" "QCD_50-80_mu" "QCD_800-1000_mu" "QCD_120-170_mu" "QCD_170-300_mu" "QCD_300-470_mu" "QCD_470-600_mu" "QCD_600-800_mu" "QCD_80-120_mu")
    source submit.sh  
    #source hadd.sh /home/chasejeon/LQanalyzer_Oct2015_8TeV/LQanalyzer/data/output/ElectronFakes/MC/  FakeRateCalculator_El_mc_5_3_14.root  FakeRateCalculator_El_SK*
    #mv /home/chasejeon/LQanalyzer_Oct2015_8TeV/LQanalyzer/data/output/ElectronFakes/MC/FakeRateCalculator_El_mc_5_3_14.root /home/jskim/LQanalyzer_Oct2015_8TeV/LQanalyzer/data/output/ElectronFakes/

fi

if [[ $runQCD  == "true" ]];
then
    source functions.sh
    cycle="FakeRateCalculator_Mu"
    skinput="True"

    njobs=30
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    outputdir=$LQANALYZER_DIR"/data/output/MuonFakes/QCD/"
    
    declare -a input_samples=("QCD_20_30_EM" "QCD_20_30_BCtoE" "QCD_30_80_EM" "QCD_30_80_BCtoE" "QCD_80_170_EM" "QCD_80_170_BCtoE" "QCD_170_250_EM" "QCD_170_250_BCtoE" "QCD_250_350_EM" "QCD_250_350_BCtoE" "QCD_350_EM" "QCD_350_BCtoE" )

   source submit.sh
   # source hadd.sh /home/chasejeon/LQanalyzer_Oct2015_8TeV/LQanalyzer/data/output/ElectronFakes/QCD/ FakeRateCalculator_El_SKQCD_5_3_14.root FakeRateCalculator_El_SKQCD*
   # mv /home/chasejeon/LQanalyzer_Oct2015_8TeV/LQanalyzer/data/output/ElectronFakes/QCD/FakeRateCalculator_El_SKQCD_5_3_14.root /home/jskim/LQanalyzer_Oct2015_8TeV/LQanalyzer/data/output/ElectronFakes/
fi

if [[ $runQCD2  == "true" ]];
then
    source functions.sh
    cycle="FakeRateCalculator_Mu"
    skinput="True"

    njobs=30
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    outputdir=$LQANALYZER_DIR"/data/output/MuonFakes/QCD2/"

    declare -a input_samples=( "QCD_30-40_EM2" "QCD_40_EM2")


    source submit.sh
    #source hadd.sh /home/chasejeon/LQanalyzer_Oct2015_8TeV/LQanalyzer/data/output/ElectronFakes/QCD2/ FakeRateCalculator_El_SKQCD_5_3_14.root FakeRateCalculator_El_SKQCD*

fi


if [[ $runData  == "true" ]];
then
    source functions.sh
    cycle="FakeRateCalculator_Mu"
    skinput="True"
    #skinput="False"

    njobs=30
    usebatch="True"
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    #stream="muon"
    stream="muon_lowpt"
    #stream="singlemuon"
    nevents=-1
    #outputdir=$LQANALYZER_DIR"/data/output/MuonFakes/dijet_topology/dXY_0p01_dZ_0p5/period/"
    outputdir="/data4/LQAnalyzerCode/jskim/LQanalyzer/data/output/MuonFakes/"
    
    #declare -a input_samples=("A" "B" "C" "D")
    declare -a input_samples=("D")

    source submit.sh $1

    #source hadd.sh /home/chasejeon/LQanalyzer_Oct2015_8TeV/LQanalyzer/data/output/ElectronFakes/Data/ FakeRateCalculator_El_data_5_3_14.root FakeRateCalculator_El_period*  
   # mv /home/chasejeon/LQanalyzer_Oct2015_8TeV/LQanalyzer/data/output/ElectronFakes/Data/FakeRateCalculator_El_data_5_3_14.root /home/jskim/LQanalyzer_Oct2015_8TeV/LQanalyzer/data/output/ElectronFakes/
fi



echo ""
echo "End of example_submit.sh script."

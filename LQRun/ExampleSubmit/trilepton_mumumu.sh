#!/bin/sh

###### SET WHAT JOBS TO RUN
runMC=true
runTrileptonMuMuMu=true

if [[ $runMC  == "true" ]]; 
then
    source functions.sh
    cycle="trilepton_mumumu"
    skinput="True"
    useskim=""
    outputdir=$LQANALYZER_DIR"/data/output/trilepton/mumumu/no_lep_assigning_here/"
    #### JOB CONFIGURATION
    njobs=30
  	data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    nevents=-1
    remove="False"
		#declare -a input_samples=("HN40_mumumu_new" "HN50_mumumu_new" "HN60_mumumu_new")
		#declare -a input_samples=("ttbar" "stbar_sch" "stbar_tch" "stbar_tW" "st_sch" "st_tch" "st_tW" "topHAD" "topLJ")
		#declare -a input_samples=("HN40_new" "HN50_new" "HN60_new")
		#declare -a input_samples=("ggHtoZZ" "HtoZZ")
		declare -a input_samples=("DY10to50" "DY50plus" "topDIL" "Wbb" "ttW" "ttZ" "WWW" "TTWW" "TTG" "ZZZ" "WZZ" "WWZ" "WWG" "WW_mg" "WZtollqq_mg" "WZtoqqln_mg" "WZtollln_mg" "ZZtollnn_mg" "ZZtollqq_mg" "ZZtollll_mg" "Zbb" "HtoWW" "Wtollln_new")
		source submit.sh $1
fi
    

################ DOUBLEMUON DATA
if [[ $runTrileptonMuMuMu  == "true" ]];
then
    source functions.sh
    cycle="trilepton_mumumu"
    skinput="True"
    stream="muon"
    useskim="DiLep"
    outputdir=$LQANALYZER_DIR"/data/output/trilepton/mumumu/no_lep_assigning_here/"
    #### JOB CONFIGURATION
    njobs=30
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
		nevents=-1
		declare -a input_samples=("A" "B" "C" "D")
		#declare -a input_samples=("B")
    source submit.sh $1
fi     





echo ""
echo "End of example_submit.sh script."

#!/bin/sh

###### SET WHAT JOBS TO RUN
runMC=true
runTrileptonMuMuMu=false

if [[ $runMC  == "true" ]]; 
then
    source functions.sh
    cycle="trilepton_eee"
    skinput="True"
    useskim=""
    outputdir=$LQANALYZER_DIR"/data/output/trilepton/eee/"
    #### JOB CONFIGURATION
    njobs=30
  	data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
    #nevents=-1
    remove="False"
		#declare -a input_samples=("ggHtoZZ" "HtoZZ" "HtoTauTau" "HtoWW" "QCD_mumu")
		#declare -a input_samples=("ttbar" "stbar_sch" "stbar_tch" "stbar_tW" "st_sch" "st_tch" "st_tW" "topHAD" "topLJ")
		declare -a input_samples=("HN40_eee_new" "HN50_eee_new" "HN60_eee_new")
		#declare -a input_samples=("Wgamma" "WgammaMu")
		#declare -a input_samples=("DY10to50" "DY50plus" "WZ_py" "ZZ_py" "WW_py" "topDIL" "Wbb" "ttW" "ttZ" "WWW" "TTWW" "TTG" "ZZZ" "WZZ" "WWZ" "WWG" "WW_mg" "WZtollqq_mg" "WZtoqqln_mg" "WZtollln_mg" "ZZtollnn_mg" "ZZtollqq_mg" "ZZtollll_mg" "Zbb")
		#declare -a input_samples=("DY10to50" "DY50plus" "topDIL" "Wbb" "ttW" "ttZ" "WWW" "TTWW" "TTG" "ZZZ" "WZZ" "WWZ" "WWG" "WW_mg" "WZtollqq_mg" "WZtoqqln_mg" "WZtollln_mg" "ZZtollnn_mg" "ZZtollqq_mg" "ZZtollll_mg" "Zbb" "HtoWW" "Wtollln_new" "ggHtoZZ" "HtoZZ")
		source submit.sh $1
fi
    

################ DOUBLEMUON DATA
if [[ $runTrileptonMuMuMu  == "true" ]];
then
    source functions.sh
    cycle="trilepton_eee"
    skinput="True"
    stream="egamma"
    useskim="DiLep"
    outputdir=$LQANALYZER_DIR"/data/output/trilepton/eee/"
    #### JOB CONFIGURATION
    njobs=30
    data_lumi="AtoD"
    loglevel="INFO"
    logstep=1000
		#nevents=100000
		declare -a input_samples=("A" "B" "C" "D")
		#declare -a input_samples=("B")
    source submit.sh $1
fi     





echo ""
echo "End of example_submit.sh script."

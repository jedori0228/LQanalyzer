#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################

### FR ###
declare -a trilep_fake_bkg=(
'DYJets_10to50' 'DYJets'
'SingleTop_s' 'SingleTop_t' 'SingleTbar_t' 'SingleTop_tW' 'SingleTbar_tW'
'WJets'
'WZ' 'ZZ' 'WW'
'WWW' 'WZZ' 'WWZ' 'ZZZ'
'TTJets_aMC'
)
declare -a QCD_FR=(
'QCD_Pt-1000toInf_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-15to20_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-80to120_MuEnriched'
)

### FR MC Closure ###
declare -a FR_MC_Closure=(
'DYJets_10to50' 'DYJets'
'WJets'
)

### SR ###
declare -a trilep_bkg=(
'WZTo3LNu_powheg' 'ZZTo4L_powheg'
'WWW' 'WZZ' 'WWZ' 'ZZZ'
'ttW' 'ttZ'
'ttH_nonbb'
'WGtoLNuG' 'WGtoLNuMM' 'ZGto2LG'
)
declare -a trilep_signal=('HN40_mumumu_VmuN_0p1' 'HN60_mumumu_VmuN_0p1' 'HN150_mumumu_VmuN_0p1' 'HN700_mumumu_VmuN_0p1')

### CR ###
declare -a trilep_CR_bkg=(
'WZTo3LNu_powheg' 'ZZTo4L_powheg'
'WWW' 'WZZ' 'WWZ' 'ZZZ'
'ttW' 'ttZ'
'ttH_nonbb'
'WGtoLNuG' 'WGtoLNuMM' 'ZGto2LG'
)
declare -a trilep_signal=('HN40_mumumu_VmuN_0p1' 'HN60_mumumu_VmuN_0p1' 'HN150_mumumu_VmuN_0p1' 'HN700_mumumu_VmuN_0p1')



declare -a trilep_bkg_miss=(
'DY10to50_MCatNLO' 'DY50plus_MCatNLO' 'ZZ_llll_powheg' 'WZ_lllnu_powheg' 'WG_lnuG_madgraph' 'WGstarToLNuMuMu' 'ZG_llG_MCatNLO'
)




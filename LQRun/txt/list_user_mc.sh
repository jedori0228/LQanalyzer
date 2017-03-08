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
'TTJets_aMC'
)
declare -a trilep_fake_bkg_short=(
'DYJets_10to50'
'SingleTop_s' 'SingleTop_t' 'SingleTbar_t' 'SingleTop_tW' 'SingleTbar_tW'
'WZ' 'ZZ' 'WW'
)
declare -a trilep_fake_bkg_long=(
'DYJets'
'WJets'
'TTJets_aMC'
)
declare -a QCD_FR=(
'QCD_Pt-1000toInf_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-15to20_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-80to120_MuEnriched'
)

### FR MC Closure ###
declare -a FR_MC_Closure=(
'TTJets_aMC'
'DYJets_10to50' 'DYJets'
'WJets'
'QCD_Pt-1000toInf_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-15to20_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-80to120_MuEnriched'
)

### SR ###
declare -a trilep_bkg=(
'WZTo3LNu_powheg' 'ZZTo4L_powheg'
'WWW' 'WZZ' 'WWZ' 'ZZZ'
'ttW' 'ttZ'
'ttH_nonbb'
)
declare -a trilep_nonprompt_bkg=(
'DYJets_10to50' 'DYJets' 'WJets' 'TTJets_aMC'
)
declare -a trilep_diboson_had=(
'WZ' 'ZZ' 'WW'
)
#declare -a trilep_signal=('HN40_mumumu_VmuN_0p1' 'HN60_mumumu_VmuN_0p1' 'HN150_mumumu_VmuN_0p1' 'HN700_mumumu_VmuN_0p1')
declare -a trilep_signal=(
'HN_MuMuMu_5'
'HN_MuMuMu_10'
'HN_MuMuMu_20'
'HN_MuMuMu_30'
'HN_MuMuMu_40'
'HN_MuMuMu_50'
'HN_MuMuMu_60'
'HN_MuMuMu_70'
'HN_MuMuMu_90'
'HN_MuMuMu_100'
'HN_MuMuMu_150'
'HN_MuMuMu_200'
'HN_MuMuMu_300'
'HN_MuMuMu_400'
'HN_MuMuMu_500'
'HN_MuMuMu_700'
'HN_MuMuMu_1000')
declare -a Official_trilep_signal=(
'Official_HN_MuMuMu_5'
'Official_HN_MuMuMu_10'
'Official_HN_MuMuMu_20'
'Official_HN_MuMuMu_30'
'Official_HN_MuMuMu_40'
'Official_HN_MuMuMu_50'
'Official_HN_MuMuMu_60'
'Official_HN_MuMuMu_70'
'Official_HN_MuMuMu_90'
'Official_HN_MuMuMu_100'
'Official_HN_MuMuMu_150'
'Official_HN_MuMuMu_200'
'Official_HN_MuMuMu_300'
'Official_HN_MuMuMu_400'
'Official_HN_MuMuMu_500'
'Official_HN_MuMuMu_700'
'Official_HN_MuMuMu_1000')

### CR ###
declare -a trilep_CR_bkg=(
'WZTo3LNu_powheg' 'ZZTo4L_powheg'
'WWW' 'WZZ' 'WWZ' 'ZZZ'
'ttW' 'ttZ'
'ttH_nonbb'
)
declare -a dilep_CR_bkg=(
'DYJets_10to50' 'DYJets'
)


declare -a trilep_bkg_miss=(
'TTJets_aMC' 'QCD_Pt-1000toInf_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-15to20_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-80to120_MuEnriched'
)

declare -a VGamma=(
'ZGto2LG')
# 'WGtoLNuG' 'WGtoLNuMM'
#)

declare -a FourLep=(
'WZZ' 'WWZ' 'ZZZ' 'ttZ' 'ttH_nonbb'
)

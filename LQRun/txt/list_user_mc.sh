#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################

declare -a example=("" "WW" "" "WZ")

declare -a tmplist=('WpWp_qcd_madgraph' 'ZG_llG_MCatNLO' 'ZZ_llnunu_powheg' 'ZZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'ZZ_llll_powheg' 'ZZ_pythia8' 'ttHnobb_Powheg' 'ttHtobb_Powheg')


declare -a hn=('HN_EE_M40' 'HN_EE_M100' 'HN_EE_M500' 'HN_EE_M1500')

declare -a trilep_bkg=(
'WZTo3LNu_powheg' 'ZZTo4L_powheg'
'WWW' 'WZZ' 'WWZ' 'ZZZ'
'ttW' 'ttZ'
'ttH_nonbb'
'WGtoLNuG' 'WGtoLNuMM' 'ZGto2LG'
)

declare -a trilep_leptonskim=(
'ZZ_llll_powheg' 'WZ_lllnu_powheg' 'WG_lnuG_madgraph' 'WGstarToLNuMuMu' 'ZG_llG_MCatNLO'
)

declare -a trilep_CR_bkg=(
'WZTo3LNu_powheg' 'ZZTo4L_powheg'
'WWW' 'WZZ' 'WWZ' 'ZZZ'
'ttW' 'ttZ'
'ttH_nonbb'
'WGtoLNuG' 'WGtoLNuMM' 'ZGto2LG'
)

declare -a trilep_bkg_miss=(
'DY10to50_MCatNLO' 'DY50plus_MCatNLO' 'ZZ_llll_powheg' 'WZ_lllnu_powheg' 'WG_lnuG_madgraph' 'WGstarToLNuMuMu' 'ZG_llG_MCatNLO'
)

declare -a trilep_fake_bkg=(
'DYJets_10to50' 'DYJets'
'SingleTop_s' 'SingleTop_t' 'SingleTbar_t' 'SingleTop_tW' 'SingleTbar_tW'
'ttW' 'ttZ'
'WJets'
'WZ' 'ZZ' 'WW'
'WWW' 'WZZ' 'WWZ' 'ZZZ'
'TTJets_aMC'
'WGtoLNuG'
'ttH_nonbb'
)

declare -a trilep_signal=('HN40_mumumu_VmuN_0p1' 'HN60_mumumu_VmuN_0p1' 'HN150_mumumu_VmuN_0p1' 'HN700_mumumu_VmuN_0p1')

declare -a QCD_FR=(
'QCD_Pt-1000toInf_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-15to20_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-80to120_MuEnriched'
)

## QCD em
## 'QCD_em15to20_pythia8' 'QCD_em20to30_pythia8' 'QCD_em30to50_pythia8' 'QCD_em50to80_pythia8' 'QCD_em80to120_pythia8' 'QCD_em120to170_pythia8' 'QCD_em170to300_pythia8' 'QCD_em300toINF_pythia8'
## QCD
## 'QCD_15to20_bcToE_pythia8' 'QCD_20to30_bcToE_pythia8' 'QCD_30to80_bcToE_pythia8' 'QCD_80to170_bcToE_pythia8' 'QCD_170to250_bcToE_pythia8' 'QCD_250toInf_bcToE_pythia8'


#declare -a QCD_FR=('QCD_mu1000toINF_pythia8' 'QCD_em120to170_pythia8' 'QCD_mu120to170_pythia8' 'QCD_em15to20_pythia8' 'QCD_mu15to20_pythia8' 'QCD_em170to300_pythia8' 'QCD_mu170to300_pythia8' 'QCD_em20to30_pythia8' 'QCD_mu20to30_pythia8' 'QCD_mu300to470_pythia8' 'QCD_em300toINF_pythia8' 'QCD_DoubleEM_30to40_pythia8' 'QCD_em30to50_pythia8' 'QCD_mu30to50_pythia8' 'QCD_DoubleEM_30toInf_pythia8' 'QCD_DoubleEM_40toInf_pythia8' 'QCD_mu470to600_pythia8' 'QCD_em50to80_pythia8' 'QCD_mu50to80_pythia8' 'QCD_mu600to800_pythia8' 'QCD_mu800to1000_pythia8' 'QCD_em80to120_pythia8' 'QCD_mu80to120_pythia8' 'QCD_15to20_bcToE_pythia8' 'QCD_170to250_bcToE_pythia8' 'QCD_20to30_bcToE_pythia8' 'QCD_250toInf_bcToE_pythia8' 'QCD_30to80_bcToE_pythia8' 'QCD_80to170_bcToE_pythia8')


#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################

declare -a example=('WZ_pythia8')

declare -a tmplist=('WpWp_qcd_madgraph' 'ZG_llG_MCatNLO' 'ZZ_llnunu_powheg' 'ZZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'ZZ_llll_powheg' 'ZZ_pythia8' 'ttHnobb_Powheg' 'ttHtobb_Powheg')


declare -a hn=('HN_EE_M40' 'HN_EE_M100' 'HN_EE_M500' 'HN_EE_M1500')

declare -a trilep_bkg=(
'WZ_pythia8' 'ZZ_pythia8'
'ttHnobb_Powheg'
'WZZ_MCatNLO' 'WWZ_MCatNLO' 'ZZZ_MCatNLO'
'ttWJets' 'ttZJets'
'ttHnobb_Powheg'
)
declare -a trilep_CR_bkg=(
'WW_pythia8' 'WZ_pythia8' 'ZZ_pythia8'
'WZZ_MCatNLO' 'WWZ_MCatNLO' 'ZZZ_MCatNLO'
'ttWJets' 'ttZJets'
'ttHnobb_Powheg'
)

declare -a trilep_bkg_miss=(
'ttWJets' 'ttZJets'
)

declare -a trilep_fake_bkg=(
'DY10to50_MCatNLO' 'DY50plus_MCatNLO'
'singletop_s_MCatNLO' 'singletop_tbar_Powheg' 'singletop_t_Powheg' 'singletop_tbarW_Powheg' 'singletop_tW_Powheg'
'ttWJets' 'ttZJets'
'WJets_MCatNLO'
'WZ_pythia8' 'ZZ_pythia8' 'WW_pythia8'
'WZZ_MCatNLO' 'WWZ_MCatNLO' 'ZZZ_MCatNLO'
'TT_MCatNLO'
'ttHnobb_Powheg'
)

declare -a trilep_signal=('HN40_mumumu' 'HN60_mumumu' 'HN150_mumumu' 'HN700_mumumu')

declare -a QCD_FR=(
'QCD_mu15to20_pythia8' 'QCD_mu20to30_pythia8' 'QCD_mu30to50_pythia8' 'QCD_mu50to80_pythia8' 'QCD_mu80to120_pythia8' 'QCD_mu120to170_pythia8' 'QCD_mu170to300_pythia8' 'QCD_mu300to470_pythia8' 'QCD_mu470to600_pythia8' 'QCD_mu600to800_pythia8' 'QCD_mu800to1000_pythia8' 'QCD_mu1000toINF_pythia8' 
'QCD_DoubleEM_30to40_pythia8' 'QCD_DoubleEM_40toInf_pythia8' 'QCD_DoubleEM_30toInf_pythia8'
)

## QCD em
## 'QCD_em15to20_pythia8' 'QCD_em20to30_pythia8' 'QCD_em30to50_pythia8' 'QCD_em50to80_pythia8' 'QCD_em80to120_pythia8' 'QCD_em120to170_pythia8' 'QCD_em170to300_pythia8' 'QCD_em300toINF_pythia8'
## QCD
## 'QCD_15to20_bcToE_pythia8' 'QCD_20to30_bcToE_pythia8' 'QCD_30to80_bcToE_pythia8' 'QCD_80to170_bcToE_pythia8' 'QCD_170to250_bcToE_pythia8' 'QCD_250toInf_bcToE_pythia8'


#declare -a QCD_FR=('QCD_mu1000toINF_pythia8' 'QCD_em120to170_pythia8' 'QCD_mu120to170_pythia8' 'QCD_em15to20_pythia8' 'QCD_mu15to20_pythia8' 'QCD_em170to300_pythia8' 'QCD_mu170to300_pythia8' 'QCD_em20to30_pythia8' 'QCD_mu20to30_pythia8' 'QCD_mu300to470_pythia8' 'QCD_em300toINF_pythia8' 'QCD_DoubleEM_30to40_pythia8' 'QCD_em30to50_pythia8' 'QCD_mu30to50_pythia8' 'QCD_DoubleEM_30toInf_pythia8' 'QCD_DoubleEM_40toInf_pythia8' 'QCD_mu470to600_pythia8' 'QCD_em50to80_pythia8' 'QCD_mu50to80_pythia8' 'QCD_mu600to800_pythia8' 'QCD_mu800to1000_pythia8' 'QCD_em80to120_pythia8' 'QCD_mu80to120_pythia8' 'QCD_15to20_bcToE_pythia8' 'QCD_170to250_bcToE_pythia8' 'QCD_20to30_bcToE_pythia8' 'QCD_250toInf_bcToE_pythia8' 'QCD_30to80_bcToE_pythia8' 'QCD_80to170_bcToE_pythia8')


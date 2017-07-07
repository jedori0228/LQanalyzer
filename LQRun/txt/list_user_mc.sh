#!/bin/sh

declare -a tmp=(
'TTJets_aMC' 'DYJets_MG' 'TTLL_powheg'
)

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
'ZGto2LG' 'WGtoLNuG'
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
'DYJets_10to50' 'DYJets' 'DYJets_MG_10to50' 'DYJets_MG'
'TTJets_aMC' 'TTLL_powheg' 'TTLJ_powheg'
)

### SR ###
declare -a trilep_bkg=(
'WZTo3LNu_powheg' 'ZZTo4L_powheg'
'WWW' 'WZZ' 'WWZ' 'ZZZ'
'ttW' 'ttZ'
'ttH_nonbb'
'ZGto2LG' 'WGtoLNuG'
'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau' 'ggZZto2mu2nu' 'ggZZto2mu2tau' 'ggZZto4e' 'ggZZto4mu' 'ggZZto4tau'
)
declare -a trilep_nonprompt_bkg=(
'DYJets_10to50' 'DYJets' 'WJets' 'TTJets_aMC' 'TTLL_powheg'
'ttW' 'ttZ'
'ttH_nonbb' 'ttH_bb'
)
declare -a trilep_diboson_had=(
'WZ' 'ZZ' 'WW'
)

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
'HN_MuMuMu_1000'
'HN_SSSF_MuMuE_5'
'HN_SSSF_MuMuE_10'
'HN_SSSF_MuMuE_20'
'HN_SSSF_MuMuE_30'
'HN_SSSF_MuMuE_40'
'HN_SSSF_MuMuE_50'
'HN_SSSF_MuMuE_60'
'HN_SSSF_MuMuE_70'
'HN_SSSF_MuMuE_90'
'HN_SSSF_MuMuE_100'
'HN_SSSF_MuMuE_150'
'HN_SSSF_MuMuE_200'
'HN_SSSF_MuMuE_300'
'HN_SSSF_MuMuE_400'
'HN_SSSF_MuMuE_500'
'HN_SSSF_MuMuE_700'
'HN_SSSF_MuMuE_1000'
'HN_SSSF_MuEMu_5'
'HN_SSSF_MuEMu_10'
'HN_SSSF_MuEMu_20'
'HN_SSSF_MuEMu_30'
'HN_SSSF_MuEMu_40'
'HN_SSSF_MuEMu_50'
'HN_SSSF_MuEMu_60'
'HN_SSSF_MuEMu_70'
'HN_SSSF_MuEMu_90'
'HN_SSSF_MuEMu_100'
'HN_SSSF_MuEMu_150'
'HN_SSSF_MuEMu_200'
'HN_SSSF_MuEMu_300'
'HN_SSSF_MuEMu_400'
'HN_SSSF_MuEMu_500'
'HN_SSSF_MuEMu_700'
'HN_SSSF_MuEMu_1000'
)

### CR ###
declare -a trilep_CR_bkg=(
'WZTo3LNu_powheg' 'ZZTo4L_powheg'
'WWW' 'WZZ' 'WWZ' 'ZZZ'
'ttW' 'ttZ'
'ttH_nonbb'
'ZGto2LG' 'WGtoLNuG'
'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau' 'ggZZto2mu2nu' 'ggZZto2mu2tau' 'ggZZto4e' 'ggZZto4mu' 'ggZZto4tau'
)
declare -a dilep_CR_bkg=(
'DYJets_10to50' 'DYJets'
)


declare -a trilep_bkg_miss=(
'TTJets_aMC' 'QCD_Pt-1000toInf_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-15to20_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-80to120_MuEnriched'
)

declare -a VGamma=(
'ZGto2LG'
'WGtoLNuG')
#'WGtoLNuMM')

declare -a FourLep=(
'WZZ' 'WWZ' 'ZZZ' 'ttZ' 'ttH_nonbb'
)

declare -a ghent=(
'Ghent_HN_MuMuMu_emu_mixing_1'
'Ghent_HN_MuMuMu_emu_mixing_2'
'Ghent_HN_MuMuMu_emu_mixing_5'
'Ghent_HN_MuMuMu_emu_mixing_10'
'Ghent_HN_MuMuMu_emu_mixing_20'
'Ghent_HN_MuMuMu_emu_mixing_30'
'Ghent_HN_MuMuMu_emu_mixing_40'
'Ghent_HN_MuMuMu_emu_mixing_50'
'Ghent_HN_MuMuMu_emu_mixing_60'
'Ghent_HN_MuMuMu_emu_mixing_80'
'Ghent_HN_MuMuMu_emu_mixing_100'
'Ghent_HN_MuMuMu_emu_mixing_130'
'Ghent_HN_MuMuMu_emu_mixing_150'
'Ghent_HN_MuMuMu_emu_mixing_200'
'Ghent_HN_MuMuMu_emu_mixing_400'
'Ghent_HN_MuMuMu_emu_mixing_600'
'Ghent_HN_MuMuMu_emu_mixing_800'
'Ghent_HN_MuMuMu_emu_mixing_1000'
'Ghent_HN_MuMuMu_mu_mixing_1'
'Ghent_HN_MuMuMu_mu_mixing_2'
'Ghent_HN_MuMuMu_mu_mixing_5'
'Ghent_HN_MuMuMu_mu_mixing_10'
'Ghent_HN_MuMuMu_mu_mixing_20'
'Ghent_HN_MuMuMu_mu_mixing_30'
'Ghent_HN_MuMuMu_mu_mixing_40'
'Ghent_HN_MuMuMu_mu_mixing_50'
'Ghent_HN_MuMuMu_mu_mixing_60'
'Ghent_HN_MuMuMu_mu_mixing_80'
'Ghent_HN_MuMuMu_mu_mixing_100'
'Ghent_HN_MuMuMu_mu_mixing_130'
'Ghent_HN_MuMuMu_mu_mixing_150'
'Ghent_HN_MuMuMu_mu_mixing_200'
'Ghent_HN_MuMuMu_mu_mixing_400'
'Ghent_HN_MuMuMu_mu_mixing_600'
'Ghent_HN_MuMuMu_mu_mixing_800'
'Ghent_HN_MuMuMu_mu_mixing_1000'
)


declare -a test=(
'DYJets' 'TTJets_aMC' 'ZGto2LG'
'DYJets_MG' 'TTLL_powheg'
'WJets' 'WGtoLNuG'
'WZTo3LNu_powheg' 'ZZTo4L_powheg'
)

declare -a test2=(
'HNMoriondLLMumMum_200'
'HNMoriondLLMupMup_200'
'HNMoriondLLEmEm_200'
'HNMoriondLLEpEp_200'
'HNMoriondLL_Tchannel_MumMum_200'
'HNMoriondLL_Tchannel_MupMup_200'
'HNMoriondLL_Tchannel_EmEm_200'
'HNMoriondLL_Tchannel_EpEp_200'
)

declare dilepbkg_use=(
'WWTo2L2Nu_DS' 'ww_ds'
'WWTo2L2Nu' 'ggWWto2L2Nu'
'WZTo3LNu_powheg' 'WZto2L2Q_amcatnlo' 'WZTo3LNu_amcatnlo'
'ZZTo4L_powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau' 'ggZZto2mu2nu' 'ggZZto2mu2tau' 'ggZZto4e' 'ggZZto4mu' 'ggZZto4tau'
'ZZTo2L2Nu_Powheg'
'ZZTo2L2Q_Powheg'
'TTLL_powheg'
'ttH_nonbb' 'ttZ' 'ttW'
'WWW' 'WWZ' 'WZZ' 'ZZZ'
)

declare dilepbkg_missing=(
'DYJets_10to50' 'DYJets'
'vbhHtoZZ' 'ggHtoZZ'
'ZGto2LG'
'WGtoLNuG' 'WgstarToLNuEE' 'WgstarToLNuMuMu'
'ttWToLNu'
'ttZToLL_M-1to10'
)

declare dilepbkg_all=(
'DYJets_10to50' 'DYJets'
'WWTo2L2Nu_DS' 'ww_ds'
'WWTo2L2Nu' 'ggWWto2L2Nu'
'WZTo3LNu_powheg' 'WZto2L2Q_amcatnlo' 'WZTo3LNu_amcatnlo'
'ZZTo4L_powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau' 'ggZZto2mu2nu' 'ggZZto2mu2tau' 'ggZZto4e' 'ggZZto4mu' 'ggZZto4tau'
'ZZTo2L2Nu_Powheg'
'ZZTo2L2Q_Powheg'
'vbhHtoZZ' 'ggHtoZZ'
'ZGto2LG'
'WGtoLNuG' 'WgstarToLNuEE' 'WgstarToLNuMuMu'
'ttWToLNu'
'TTLL_powheg'
'ttH_nonbb' 'ttZ' 'ttW' 'ttZToLL_M-1to10'
'WWW' 'WWZ' 'WZZ' 'ZZZ'
)

## HNEmEm_250
declare -a DiMuSig=(
'HNMoriondLLMumMum_200'
'HNMoriondLLMupMup_200'
'HNMoriondLLMumMum_100'
'HNMoriondLLMupMup_100'
'HNMoriondLLMumMum_50'
'HNMoriondLLMupMup_50'
'HNMoriondLL_Tchannel_MumMum_100'
'HNMoriondLL_Tchannel_MupMup_100'
'HNMoriondLL_Tchannel_MumMum_200'
'HNMoriondLL_Tchannel_MupMup_200'
'HNMoriondLL_Tchannel_MumMum_500'
'HNMoriondLL_Tchannel_MupMup_500'
'HNMoriondLLEmEm_200'
'HNMoriondLLEpEp_200'
'HNMoriondLLEmEm_100'
'HNMoriondLLEpEp_100'
'HNMoriondLLEmEm_50'
'HNMoriondLLEpEp_50'
'HNMoriondLL_Tchannel_EmEm_100'
'HNMoriondLL_Tchannel_EpEp_100'
'HNMoriondLL_Tchannel_EmEm_200'
'HNMoriondLL_Tchannel_EpEp_200'
'HNMoriondLL_Tchannel_EmEm_500'
'HNMoriondLL_Tchannel_EpEp_500'
)
declare -a DiMuSig_all=(
'HNMupMup_40'  'HNMupMup_50'  'HNMupMup_60'  'HNMupMup_70'  'HNMupMup_80'  'HNMupMup_90'  'HNMupMup_100'  'HNMupMup_125'  'HNMupMup_150'  'HNMupMup_200'  'HNMupMup_250'  'HNMupMup_300'  'HNMupMup_400'  'HNMupMup_500'  'HNMupMup_600'  'HNMupMup_700'  'HNMupMup_800'  'HNMupMup_900'  'HNMupMup_1000'  'HNMupMup_1100'  'HNMupMup_1200'  'HNMupMup_1300'  'HNMupMup_1400'  'HNMupMup_1500'
'HNMumMum_40'  'HNMumMum_50'  'HNMumMum_60'  'HNMumMum_70'  'HNMumMum_80'  'HNMumMum_90'  'HNMumMum_100'  'HNMumMum_125'  'HNMumMum_150'  'HNMumMum_200'  'HNMumMum_250'  'HNMumMum_300'  'HNMumMum_400'  'HNMumMum_500'  'HNMumMum_600'  'HNMumMum_700'  'HNMumMum_800'  'HNMumMum_900'  'HNMumMum_1000'  'HNMumMum_1100'  'HNMumMum_1200'  'HNMumMum_1300'  'HNMumMum_1400'  'HNMumMum_1500'
'HNEpEp_40'  'HNEpEp_50'  'HNEpEp_60'  'HNEpEp_70'  'HNEpEp_80'  'HNEpEp_90'  'HNEpEp_100'  'HNEpEp_125'  'HNEpEp_150'  'HNEpEp_200'  'HNEpEp_250'  'HNEpEp_300'  'HNEpEp_400'  'HNEpEp_500'  'HNEpEp_600'  'HNEpEp_700'  'HNEpEp_800'  'HNEpEp_900'  'HNEpEp_1000'  'HNEpEp_1100'  'HNEpEp_1200'  'HNEpEp_1300'  'HNEpEp_1400'  'HNEpEp_1500'
'HNEmEm_40'  'HNEmEm_50'  'HNEmEm_60'  'HNEmEm_70'  'HNEmEm_80'  'HNEmEm_90'  'HNEmEm_100'  'HNEmEm_125'  'HNEmEm_150'  'HNEmEm_200' 'HNEmEm_250' 'HNEmEm_300'  'HNEmEm_400'  'HNEmEm_500'  'HNEmEm_600'  'HNEmEm_700'  'HNEmEm_800'  'HNEmEm_900'  'HNEmEm_1000'  'HNEmEm_1100'  'HNEmEm_1200'  'HNEmEm_1300'  'HNEmEm_1400'  'HNEmEm_1500'
'HNMupMup_Tchannel_300'  'HNMupMup_Tchannel_400'  'HNMupMup_Tchannel_600'  'HNMupMup_Tchannel_700'  'HNMupMup_Tchannel_800'  'HNMupMup_Tchannel_900'  'HNMupMup_Tchannel_1000'  'HNMupMup_Tchannel_1200'  'HNMupMup_Tchannel_1300'  'HNMupMup_Tchannel_1400'  'HNMupMup_Tchannel_1500' 'HNMoriondLL_Tchannel_MupMup_100'  'HNMoriondLL_Tchannel_MupMup_200'  'HNMoriondLL_Tchannel_MupMup_500'  'HNMoriondLL_Tchannel_MupMup_1100'
'HNMumMum_Tchannel_300'  'HNMumMum_Tchannel_400'  'HNMumMum_Tchannel_600'  'HNMumMum_Tchannel_700'  'HNMumMum_Tchannel_800'  'HNMumMum_Tchannel_900'  'HNMumMum_Tchannel_1000'  'HNMumMum_Tchannel_1200'  'HNMumMum_Tchannel_1300'  'HNMumMum_Tchannel_1400'  'HNMumMum_Tchannel_1500' 'HNMoriondLL_Tchannel_MumMum_100'  'HNMoriondLL_Tchannel_MumMum_200'  'HNMoriondLL_Tchannel_MumMum_500'  'HNMoriondLL_Tchannel_MumMum_1100'
'HNEpEp_Tchannel_300'  'HNEpEp_Tchannel_400'  'HNEpEp_Tchannel_600'  'HNEpEp_Tchannel_700'  'HNEpEp_Tchannel_800'  'HNEpEp_Tchannel_900'  'HNEpEp_Tchannel_1000'  'HNEpEp_Tchannel_1200'  'HNEpEp_Tchannel_1300'  'HNEpEp_Tchannel_1400'  'HNEpEp_Tchannel_1500' 'HNMoriondLL_Tchannel_EpEp_100'  'HNMoriondLL_Tchannel_EpEp_200'  'HNMoriondLL_Tchannel_EpEp_500'  'HNMoriondLL_Tchannel_EpEp_1100'
'HNEmEm_Tchannel_300'  'HNEmEm_Tchannel_400'  'HNEmEm_Tchannel_600'  'HNEmEm_Tchannel_700'  'HNEmEm_Tchannel_800'  'HNEmEm_Tchannel_900'  'HNEmEm_Tchannel_1000'  'HNEmEm_Tchannel_1200'  'HNEmEm_Tchannel_1300'  'HNEmEm_Tchannel_1400'  'HNEmEm_Tchannel_1500' 'HNMoriondLL_Tchannel_EmEm_100'  'HNMoriondLL_Tchannel_EmEm_200'  'HNMoriondLL_Tchannel_EmEm_500'  'HNMoriondLL_Tchannel_EmEm_1100'
)

declare -a trilep_signal_mme=(
'HN_SSSF_MuEMu_5'
'HN_SSSF_MuEMu_10'
'HN_SSSF_MuEMu_20'
'HN_SSSF_MuEMu_30'
'HN_SSSF_MuEMu_40'
'HN_SSSF_MuEMu_50'
'HN_SSSF_MuEMu_60'
'HN_SSSF_MuEMu_70'
'HN_SSSF_MuEMu_90'
'HN_SSSF_MuEMu_100'
'HN_SSSF_MuEMu_150'
'HN_SSSF_MuEMu_200'
'HN_SSSF_MuEMu_300'
'HN_SSSF_MuEMu_400'
'HN_SSSF_MuEMu_500'
'HN_SSSF_MuEMu_700'
'HN_SSSF_MuEMu_1000'
'HNMoriondLLMumMum_200'
'HNMoriondLLMupMup_200'
'HNMoriondLLEmEm_200'
'HNMoriondLLEpEp_200'
'HNMoriondLL_Tchannel_MumMum_200'
'HNMoriondLL_Tchannel_MupMup_200'
'HNMoriondLL_Tchannel_EmEm_200'
'HNMoriondLL_Tchannel_EpEp_200'
)



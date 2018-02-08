#!/bin/sh

declare -a samples8TeV=(
'WW' 'ZZ' 'WZ'
'WWW' 'WWZ' 'WZZ' 'ZZZ'
'ttW' 'ttZ' 'ttH_nonbb'
)

declare -a privTch=(
'HeavyNeutrinoToEpEp_Tchannel_M300' 'HeavyNeutrinoToEpEp_Tchannel_M400'  
'HeavyNeutrinoToEpEp_Tchannel_M500'  'HeavyNeutrinoToEpEp_Tchannel_M600'  'HeavyNeutrinoToEpEp_Tchannel_M700'  'HeavyNeutrinoToEpEp_Tchannel_M800'  
'HeavyNeutrinoToEpEp_Tchannel_M900'  'HeavyNeutrinoToEpEp_Tchannel_M1000'  'HeavyNeutrinoToEpEp_Tchannel_M1100'  'HeavyNeutrinoToEpEp_Tchannel_M1200'  
'HeavyNeutrinoToEpEp_Tchannel_M1300'  'HeavyNeutrinoToEpEp_Tchannel_M1400'  'HeavyNeutrinoToEpEp_Tchannel_M1500'  'HeavyNeutrinoToEpEp_Tchannel_M1700'  
'HeavyNeutrinoToEpEp_Tchannel_M2000'  'HeavyNeutrinoToEmEm_Tchannel_M300'  'HeavyNeutrinoToEmEm_Tchannel_M400'  'HeavyNeutrinoToEmEm_Tchannel_M500'  
'HeavyNeutrinoToEmEm_Tchannel_M600'  'HeavyNeutrinoToEmEm_Tchannel_M700'  'HeavyNeutrinoToEmEm_Tchannel_M800'  'HeavyNeutrinoToEmEm_Tchannel_M900'  
'HeavyNeutrinoToEmEm_Tchannel_M1000'  'HeavyNeutrinoToEmEm_Tchannel_M1100'  'HeavyNeutrinoToEmEm_Tchannel_M1200'  'HeavyNeutrinoToEmEm_Tchannel_M1300'  
'HeavyNeutrinoToEmEm_Tchannel_M1400'  'HeavyNeutrinoToEmEm_Tchannel_M1500'  'HeavyNeutrinoToEmEm_Tchannel_M1700'  'HeavyNeutrinoToEmEm_Tchannel_M2000'  
'HeavyNeutrinoToEpMup_Tchannel_M300'  'HeavyNeutrinoToEpMup_Tchannel_M400'  'HeavyNeutrinoToEpMup_Tchannel_M500'  'HeavyNeutrinoToEpMup_Tchannel_M600'  
'HeavyNeutrinoToEpMup_Tchannel_M700'  'HeavyNeutrinoToEpMup_Tchannel_M800'  'HeavyNeutrinoToEpMup_Tchannel_M900'  'HeavyNeutrinoToEpMup_Tchannel_M1000'  
'HeavyNeutrinoToEpMup_Tchannel_M1100'  'HeavyNeutrinoToEpMup_Tchannel_M1200'  'HeavyNeutrinoToEpMup_Tchannel_M1300'  'HeavyNeutrinoToEpMup_Tchannel_M1400'  
'HeavyNeutrinoToEpMup_Tchannel_M1500'  'HeavyNeutrinoToEpMup_Tchannel_M1700'  'HeavyNeutrinoToEpMup_Tchannel_M2000'  'HeavyNeutrinoToEmMum_Tchannel_M300'  
'HeavyNeutrinoToEmMum_Tchannel_M400'  'HeavyNeutrinoToEmMum_Tchannel_M500'  'HeavyNeutrinoToEmMum_Tchannel_M600'  'HeavyNeutrinoToEmMum_Tchannel_M700'  
'HeavyNeutrinoToEmMum_Tchannel_M800'  'HeavyNeutrinoToEmMum_Tchannel_M900'  'HeavyNeutrinoToEmMum_Tchannel_M1000'  'HeavyNeutrinoToEmMum_Tchannel_M1100'  
'HeavyNeutrinoToEmMum_Tchannel_M1200'  'HeavyNeutrinoToEmMum_Tchannel_M1300'  'HeavyNeutrinoToEmMum_Tchannel_M1400'  'HeavyNeutrinoToEmMum_Tchannel_M1500'  
'HeavyNeutrinoToEmMum_Tchannel_M1700'  'HeavyNeutrinoToEmMum_Tchannel_M2000'  'HeavyNeutrinoToMupMup_Tchannel_M300'  'HeavyNeutrinoToMupMup_Tchannel_M400'  
'HeavyNeutrinoToMupMup_Tchannel_M500'  'HeavyNeutrinoToMupMup_Tchannel_M600'  'HeavyNeutrinoToMupMup_Tchannel_M700'  'HeavyNeutrinoToMupMup_Tchannel_M800'  
'HeavyNeutrinoToMupMup_Tchannel_M900'  'HeavyNeutrinoToMupMup_Tchannel_M1000'  'HeavyNeutrinoToMupMup_Tchannel_M1100'  'HeavyNeutrinoToMupMup_Tchannel_M1200'  
'HeavyNeutrinoToMupMup_Tchannel_M1300'  'HeavyNeutrinoToMupMup_Tchannel_M1400'  'HeavyNeutrinoToMupMup_Tchannel_M1500'  'HeavyNeutrinoToMupMup_Tchannel_M1700'  
'HeavyNeutrinoToMupMup_Tchannel_M2000'  'HeavyNeutrinoToMumMum_Tchannel_M300'  'HeavyNeutrinoToMumMum_Tchannel_M400'  'HeavyNeutrinoToMumMum_Tchannel_M500'  
'HeavyNeutrinoToMumMum_Tchannel_M600'  'HeavyNeutrinoToMumMum_Tchannel_M700'  'HeavyNeutrinoToMumMum_Tchannel_M800'  'HeavyNeutrinoToMumMum_Tchannel_M900'  
'HeavyNeutrinoToMumMum_Tchannel_M1000'  'HeavyNeutrinoToMumMum_Tchannel_M1100'  'HeavyNeutrinoToMumMum_Tchannel_M1200'  'HeavyNeutrinoToMumMum_Tchannel_M1300'  
'HeavyNeutrinoToMumMum_Tchannel_M1400'  'HeavyNeutrinoToMumMum_Tchannel_M1500'  'HeavyNeutrinoToMumMum_Tchannel_M1700'  'HeavyNeutrinoToMumMum_Tchannel_M2000'  
'HeavyNeutrinoToMupEp_Tchannel_M300'  'HeavyNeutrinoToMupEp_Tchannel_M400'  'HeavyNeutrinoToMupEp_Tchannel_M500'  'HeavyNeutrinoToMupEp_Tchannel_M600'  
'HeavyNeutrinoToMupEp_Tchannel_M700'  'HeavyNeutrinoToMupEp_Tchannel_M800'  'HeavyNeutrinoToMupEp_Tchannel_M900'  'HeavyNeutrinoToMupEp_Tchannel_M1000'  
'HeavyNeutrinoToMupEp_Tchannel_M1100'  'HeavyNeutrinoToMupEp_Tchannel_M1200'  'HeavyNeutrinoToMupEp_Tchannel_M1300'  'HeavyNeutrinoToMupEp_Tchannel_M1400'  
'HeavyNeutrinoToMupEp_Tchannel_M1500'  'HeavyNeutrinoToMupEp_Tchannel_M1700'  'HeavyNeutrinoToMupEp_Tchannel_M2000'  'HeavyNeutrinoToMumEm_Tchannel_M300'  
'HeavyNeutrinoToMumEm_Tchannel_M400'  'HeavyNeutrinoToMumEm_Tchannel_M500'  'HeavyNeutrinoToMumEm_Tchannel_M600'  'HeavyNeutrinoToMumEm_Tchannel_M700'  
'HeavyNeutrinoToMumEm_Tchannel_M800'  'HeavyNeutrinoToMumEm_Tchannel_M900'  'HeavyNeutrinoToMumEm_Tchannel_M1000'  'HeavyNeutrinoToMumEm_Tchannel_M1100'  
'HeavyNeutrinoToMumEm_Tchannel_M1200'  'HeavyNeutrinoToMumEm_Tchannel_M1300'  'HeavyNeutrinoToMumEm_Tchannel_M1400'  'HeavyNeutrinoToMumEm_Tchannel_M1500'  
'HeavyNeutrinoToMumEm_Tchannel_M1700'  'HeavyNeutrinoToMumEm_Tchannel_M2000'
)

declare -a privTch1=(
'HeavyNeutrinoToEpEp_Tchannel_M300' 'HeavyNeutrinoToEpEp_Tchannel_M400'  
'HeavyNeutrinoToEpEp_Tchannel_M500'  'HeavyNeutrinoToEpEp_Tchannel_M600'  'HeavyNeutrinoToEpEp_Tchannel_M700'  'HeavyNeutrinoToEpEp_Tchannel_M800'  
'HeavyNeutrinoToEpEp_Tchannel_M900'  'HeavyNeutrinoToEpEp_Tchannel_M1000'  'HeavyNeutrinoToEpEp_Tchannel_M1100'  'HeavyNeutrinoToEpEp_Tchannel_M1200'  
'HeavyNeutrinoToEpEp_Tchannel_M1300'  'HeavyNeutrinoToEpEp_Tchannel_M1400'  'HeavyNeutrinoToEpEp_Tchannel_M1500'  'HeavyNeutrinoToEpEp_Tchannel_M1700'  
'HeavyNeutrinoToEpEp_Tchannel_M2000'  'HeavyNeutrinoToEmEm_Tchannel_M300'  'HeavyNeutrinoToEmEm_Tchannel_M400'  'HeavyNeutrinoToEmEm_Tchannel_M500'  
'HeavyNeutrinoToEmEm_Tchannel_M600'  'HeavyNeutrinoToEmEm_Tchannel_M700'  'HeavyNeutrinoToEmEm_Tchannel_M800'  'HeavyNeutrinoToEmEm_Tchannel_M900'  
'HeavyNeutrinoToEmEm_Tchannel_M1000'  'HeavyNeutrinoToEmEm_Tchannel_M1100'  'HeavyNeutrinoToEmEm_Tchannel_M1200'  'HeavyNeutrinoToEmEm_Tchannel_M1300'  
'HeavyNeutrinoToEmEm_Tchannel_M1400'  'HeavyNeutrinoToEmEm_Tchannel_M1500'  'HeavyNeutrinoToEmEm_Tchannel_M1700'  'HeavyNeutrinoToEmEm_Tchannel_M2000'  
'HeavyNeutrinoToEpMup_Tchannel_M300'  'HeavyNeutrinoToEpMup_Tchannel_M400'  'HeavyNeutrinoToEpMup_Tchannel_M500'  'HeavyNeutrinoToEpMup_Tchannel_M600'  
'HeavyNeutrinoToEpMup_Tchannel_M700'  'HeavyNeutrinoToEpMup_Tchannel_M800'  'HeavyNeutrinoToEpMup_Tchannel_M900'  'HeavyNeutrinoToEpMup_Tchannel_M1000'  
'HeavyNeutrinoToEpMup_Tchannel_M1100'  'HeavyNeutrinoToEpMup_Tchannel_M1200'  'HeavyNeutrinoToEpMup_Tchannel_M1300'  'HeavyNeutrinoToEpMup_Tchannel_M1400'  
'HeavyNeutrinoToEpMup_Tchannel_M1500'  'HeavyNeutrinoToEpMup_Tchannel_M1700'  'HeavyNeutrinoToEpMup_Tchannel_M2000'  'HeavyNeutrinoToEmMum_Tchannel_M300'  
)
declare -a privTch2=(
'HeavyNeutrinoToEmMum_Tchannel_M400'  'HeavyNeutrinoToEmMum_Tchannel_M500'  'HeavyNeutrinoToEmMum_Tchannel_M600'  'HeavyNeutrinoToEmMum_Tchannel_M700'  
'HeavyNeutrinoToEmMum_Tchannel_M800'  'HeavyNeutrinoToEmMum_Tchannel_M900'  'HeavyNeutrinoToEmMum_Tchannel_M1000'  'HeavyNeutrinoToEmMum_Tchannel_M1100'  
'HeavyNeutrinoToEmMum_Tchannel_M1200'  'HeavyNeutrinoToEmMum_Tchannel_M1300'  'HeavyNeutrinoToEmMum_Tchannel_M1400'  'HeavyNeutrinoToEmMum_Tchannel_M1500'  
'HeavyNeutrinoToEmMum_Tchannel_M1700'  'HeavyNeutrinoToEmMum_Tchannel_M2000'  'HeavyNeutrinoToMupMup_Tchannel_M300'  'HeavyNeutrinoToMupMup_Tchannel_M400'  
'HeavyNeutrinoToMupMup_Tchannel_M500'  'HeavyNeutrinoToMupMup_Tchannel_M600'  'HeavyNeutrinoToMupMup_Tchannel_M700'  'HeavyNeutrinoToMupMup_Tchannel_M800'  
'HeavyNeutrinoToMupMup_Tchannel_M900'  'HeavyNeutrinoToMupMup_Tchannel_M1000'  'HeavyNeutrinoToMupMup_Tchannel_M1100'  'HeavyNeutrinoToMupMup_Tchannel_M1200'  
'HeavyNeutrinoToMupMup_Tchannel_M1300'  'HeavyNeutrinoToMupMup_Tchannel_M1400'  'HeavyNeutrinoToMupMup_Tchannel_M1500'  'HeavyNeutrinoToMupMup_Tchannel_M1700'  
'HeavyNeutrinoToMupMup_Tchannel_M2000'  'HeavyNeutrinoToMumMum_Tchannel_M300'  'HeavyNeutrinoToMumMum_Tchannel_M400'  'HeavyNeutrinoToMumMum_Tchannel_M500'  
'HeavyNeutrinoToMumMum_Tchannel_M600'  'HeavyNeutrinoToMumMum_Tchannel_M700'  'HeavyNeutrinoToMumMum_Tchannel_M800'  'HeavyNeutrinoToMumMum_Tchannel_M900'  
'HeavyNeutrinoToMumMum_Tchannel_M1000'  'HeavyNeutrinoToMumMum_Tchannel_M1100'  'HeavyNeutrinoToMumMum_Tchannel_M1200'  'HeavyNeutrinoToMumMum_Tchannel_M1300'  
'HeavyNeutrinoToMumMum_Tchannel_M1400'  'HeavyNeutrinoToMumMum_Tchannel_M1500'  'HeavyNeutrinoToMumMum_Tchannel_M1700'  'HeavyNeutrinoToMumMum_Tchannel_M2000'  
)
declare -a privTch3=(
'HeavyNeutrinoToMupEp_Tchannel_M300'  'HeavyNeutrinoToMupEp_Tchannel_M400'  'HeavyNeutrinoToMupEp_Tchannel_M500'  'HeavyNeutrinoToMupEp_Tchannel_M600'  
'HeavyNeutrinoToMupEp_Tchannel_M700'  'HeavyNeutrinoToMupEp_Tchannel_M800'  'HeavyNeutrinoToMupEp_Tchannel_M900'  'HeavyNeutrinoToMupEp_Tchannel_M1000'  
'HeavyNeutrinoToMupEp_Tchannel_M1100'  'HeavyNeutrinoToMupEp_Tchannel_M1200'  'HeavyNeutrinoToMupEp_Tchannel_M1300'  'HeavyNeutrinoToMupEp_Tchannel_M1400'  
'HeavyNeutrinoToMupEp_Tchannel_M1500'  'HeavyNeutrinoToMupEp_Tchannel_M1700'  'HeavyNeutrinoToMupEp_Tchannel_M2000'  'HeavyNeutrinoToMumEm_Tchannel_M300'  
'HeavyNeutrinoToMumEm_Tchannel_M400'  'HeavyNeutrinoToMumEm_Tchannel_M500'  'HeavyNeutrinoToMumEm_Tchannel_M600'  'HeavyNeutrinoToMumEm_Tchannel_M700'  
'HeavyNeutrinoToMumEm_Tchannel_M800'  'HeavyNeutrinoToMumEm_Tchannel_M900'  'HeavyNeutrinoToMumEm_Tchannel_M1000'  'HeavyNeutrinoToMumEm_Tchannel_M1100'  
'HeavyNeutrinoToMumEm_Tchannel_M1200'  'HeavyNeutrinoToMumEm_Tchannel_M1300'  'HeavyNeutrinoToMumEm_Tchannel_M1400'  'HeavyNeutrinoToMumEm_Tchannel_M1500'  
'HeavyNeutrinoToMumEm_Tchannel_M1700'  'HeavyNeutrinoToMumEm_Tchannel_M2000'
)

declare -a mass7585=(
'HeavyNeutrinoToEpEp_Schannel_M75' 'HeavyNeutrinoToEpEp_Schannel_M85' 'HeavyNeutrinoToEmEm_Schannel_M75' 'HeavyNeutrinoToEmEm_Schannel_M85' 'HeavyNeutrinoToEpMup_Schannel_M75' 'HeavyNeutrinoToEpMup_Schannel_M85' 'HeavyNeutrinoToEmMum_Schannel_M75' 'HeavyNeutrinoToEmMum_Schannel_M85' 'HeavyNeutrinoToMupMup_Schannel_M75' 'HeavyNeutrinoToMupMup_Schannel_M85' 'HeavyNeutrinoToMumMum_Schannel_M75' 'HeavyNeutrinoToMumMum_Schannel_M85' 'HeavyNeutrinoToMupEp_Schannel_M75' 'HeavyNeutrinoToMupEp_Schannel_M85' 'HeavyNeutrinoToMumEm_Schannel_M75' 'HeavyNeutrinoToMumEm_Schannel_M85'
)

declare -a dilepton_Moriond=(
'HNMoriondLLMupMup_50' 'HNMoriondLLMupMup_100' 'HNMoriondLLMupMup_200' 'HNMoriondLLMupMup_500' 'HNMoriondLLMupMup_1100' 
'HNMoriondLLMumMum_50' 'HNMoriondLLMumMum_100' 'HNMoriondLLMumMum_200' 'HNMoriondLLMumMum_500' 'HNMoriondLLMumMum_1100' 
'HNMoriondLLEpEp_50' 'HNMoriondLLEpEp_100' 'HNMoriondLLEpEp_200' 'HNMoriondLLEpEp_500' 'HNMoriondLLEpEp_1100' 
'HNMoriondLLEmEm_50' 'HNMoriondLLEmEm_100' 'HNMoriondLLEmEm_200' 'HNMoriondLLEmEm_500' 'HNMoriondLLEmEm_1100' 
'HNMoriondLLMupEp_50' 'HNMoriondLLMupEp_100' 'HNMoriondLLMupEp_200' 'HNMoriondLLMupEp_500' 'HNMoriondLLMupEp_1100' 
'HNMoriondLLMumEm_50' 'HNMoriondLLMumEm_100' 'HNMoriondLLMumEm_200' 'HNMoriondLLMumEm_500' 'HNMoriondLLMumEm_1100' 
'HNMoriondLLEpMup_50' 'HNMoriondLLEpMup_100' 'HNMoriondLLEpMup_200' 'HNMoriondLLEpMup_500' 'HNMoriondLLEpMup_1100' 
'HNMoriondLLEmMum_50' 'HNMoriondLLEmMum_100' 'HNMoriondLLEmMum_200' 'HNMoriondLLEmMum_500' 'HNMoriondLLEmMum_1100' 
)

declare -a shpriv=(
'HeavyNeutrinoToEE_Schannel_M20v2'
'HeavyNeutrinoToEMu_Schannel_M20v2'
'HeavyNeutrinoToMuE_Schannel_M20v2'
'HeavyNeutrinoToMuMu_Schannel_M20v2'
'HeavyNeutrinoToEE_Schannel_M30v2'
'HeavyNeutrinoToEMu_Schannel_M30v2'
'HeavyNeutrinoToMuE_Schannel_M30v2'
'HeavyNeutrinoToMuMu_Schannel_M30v2'
)

declare -a newTchannel=(
'HNDilepton_MuMu_Tchannel_M300'
'HNDilepton_MuMu_Tchannel_M600'
'HNDilepton_MuMu_Tchannel_M800'
'HNDilepton_MuMu_Tchannel_M1000'
'HNDilepton_MuMu_Tchannel_M1200'
'HNDilepton_EE_Tchannel_M300'
'HNDilepton_EE_Tchannel_M600'
'HNDilepton_EE_Tchannel_M800'
'HNDilepton_EE_Tchannel_M1000'
'HNDilepton_EE_Tchannel_M1200'
'HNDilepton_MuE_Tchannel_M300'
'HNDilepton_MuE_Tchannel_M600'
'HNDilepton_MuE_Tchannel_M800'
'HNDilepton_MuE_Tchannel_M1000'
'HNDilepton_MuE_Tchannel_M1200'
'HNDilepton_EMu_Tchannel_M300'
'HNDilepton_EMu_Tchannel_M600'
'HNDilepton_EMu_Tchannel_M800'
'HNDilepton_EMu_Tchannel_M1000'
'HNDilepton_EMu_Tchannel_M1200'
)

declare -a newTchannel1=(
'HNDilepton_MuMu_Tchannel_M300'
'HNDilepton_MuMu_Tchannel_M600'
'HNDilepton_MuMu_Tchannel_M800'
'HNDilepton_MuMu_Tchannel_M1000'
)
declare -a newTchannel2=(
'HNDilepton_MuMu_Tchannel_M1200'
'HNDilepton_EE_Tchannel_M300'
'HNDilepton_EE_Tchannel_M600'
'HNDilepton_EE_Tchannel_M800'
)
declare -a newTchannel3=(
'HNDilepton_EE_Tchannel_M1000'
'HNDilepton_EE_Tchannel_M1200'
'HNDilepton_MuE_Tchannel_M300'
'HNDilepton_MuE_Tchannel_M600'
)
declare -a newTchannel4=(
'HNDilepton_MuE_Tchannel_M800'
'HNDilepton_MuE_Tchannel_M1000'
'HNDilepton_MuE_Tchannel_M1200'
'HNDilepton_EMu_Tchannel_M300'
'HNDilepton_EMu_Tchannel_M600'
'HNDilepton_EMu_Tchannel_M800'
'HNDilepton_EMu_Tchannel_M1000'
'HNDilepton_EMu_Tchannel_M1200'
)
declare -a newTchannel5=(
'HNDilepton_EE_Tchannel_M1500'
'HNDilepton_MuMu_Tchannel_M1500'
'HNDilepton_MuE_Tchannel_M1500'
'HNDilepton_EMu_Tchannel_M1500'
)

declare -a nonprompt=(
'TT_powheg' 'WJets' 'DYJets'
)

########################
### SAMPLE LIST ########## 
#######################
declare -a tchannel_hn=('HNMoriondLL_Tchannel_EpEp_100' 'HNMoriondLL_Tchannel_EpEp_200' 'HNMoriondLL_Tchannel_EpEp_500' 'HNMoriondLL_Tchannel_EpEp_1100' 'HNMoriondLL_Tchannel_MupMup_100' 'HNMoriondLL_Tchannel_MupMup_200' 'HNMoriondLL_Tchannel_MupMup_500' 'HNMoriondLL_Tchannel_MupMup_1100' 'HNMumMum_40' 'HNMumMum_50' 'HNMumMum_200' 'HNMumMum_500' 'HNMumMum_1500'  'HNEmEm_40' 'HNEmEm_50' 'HNEmEm_200' 'HNEmEm_500' 'HNEmEm_1500'  'HNMupMup_40' 'HNMupMup_50' 'HNMupMup_200' 'HNMupMup_500' 'HNMupMup_1500' 'HNMupMup_100' 'HNMupMup_700' 'HNMupMup_1000' )

### FR ###
declare -a trilep_fake_bkg=(
'DYJets_10to50' 'DYJets'
'SingleTop_s' 'SingleTop_t' 'SingleTbar_t' 'SingleTop_tW' 'SingleTbar_tW'
'WJets'
'WZ' 'ZZ' 'WW'
'TT_powheg'
'ZGto2LG' 'WGtoLNuG'
)
declare -a trilep_fake_bkg_short=(
'SingleTbar_t' 'SingleTop_tW' 'SingleTbar_tW' 'SingleTop_s'
'DYJets_10to50'
'WZ' 'ZZ' 'WW'
)
declare -a trilep_fake_bkg_long=(
'SingleTop_t'
'DYJets'
'WJets'
'TT_powheg'
)
declare -a QCD_FR=(
'QCD_Pt-1000toInf_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-15to20_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-80to120_MuEnriched'
)

declare -a QCD_FR_EM=(
'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-30to50_EMEnriched' 'QCD_Pt-50to80_EMEnriched' 'QCD_Pt-80to120_EMEnriched'
'QCD_Pt-120to170_EMEnriched' 'QCD_Pt-170to300_EMEnriched' 'QCD_Pt-300toInf_EMEnriched'
'qcd_15to20_bctoe' 'qcd_170to250_bctoe' 'qcd_20to30_bctoe' 'qcd_250toinf_bctoe' 'qcd_30to80_bctoe' 'qcd_80to170_bctoe'
)
declare -a QCD_FR_EMEnriched=(
'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-30to50_EMEnriched' 'QCD_Pt-50to80_EMEnriched' 'QCD_Pt-80to120_EMEnriched'
'QCD_Pt-120to170_EMEnriched' 'QCD_Pt-170to300_EMEnriched' 'QCD_Pt-300toInf_EMEnriched'
)
declare -a QCD_FR_bctoe=(
'qcd_15to20_bctoe' 'qcd_170to250_bctoe' 'qcd_20to30_bctoe' 'qcd_250toinf_bctoe' 'qcd_30to80_bctoe' 'qcd_80to170_bctoe'
)


### FR MC Closure ###
declare -a FR_MC_Closure=(
'TT_powheg'
)
declare -a FR_MC_Closure_amc=(
'DYJets' 'WJets'
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
'TG' 'TTG'
'WWW' 'WWZ' 'WZZ' 'ZZZ'
'ttW' 'ttZ' 'ttH_nonbb'
'WWTo2L2Nu_DS' 'WpWpEWK' 'WpWpQCD'
)

declare -a test2=(
'HeavyNeutrinoToMupMup_Tchannel_M300'  'HeavyNeutrinoToMupMup_Tchannel_M400'  
'HeavyNeutrinoToMupMup_Tchannel_M500'  'HeavyNeutrinoToMupMup_Tchannel_M600'  'HeavyNeutrinoToMupMup_Tchannel_M700'  'HeavyNeutrinoToMupMup_Tchannel_M800'  
'HeavyNeutrinoToMupMup_Tchannel_M900'  'HeavyNeutrinoToMupMup_Tchannel_M1000'  'HeavyNeutrinoToMupMup_Tchannel_M1100'  'HeavyNeutrinoToMupMup_Tchannel_M1200'  
'HeavyNeutrinoToMupMup_Tchannel_M1300'  'HeavyNeutrinoToMupMup_Tchannel_M1400'  'HeavyNeutrinoToMupMup_Tchannel_M1500'  'HeavyNeutrinoToMupMup_Tchannel_M1700'  
'HeavyNeutrinoToMupMup_Tchannel_M2000'  'HeavyNeutrinoToMumMum_Tchannel_M300'  'HeavyNeutrinoToMumMum_Tchannel_M400'  'HeavyNeutrinoToMumMum_Tchannel_M500'  
'HeavyNeutrinoToMumMum_Tchannel_M600'  'HeavyNeutrinoToMumMum_Tchannel_M700'  'HeavyNeutrinoToMumMum_Tchannel_M800'  'HeavyNeutrinoToMumMum_Tchannel_M900'  
'HeavyNeutrinoToMumMum_Tchannel_M1000'  'HeavyNeutrinoToMumMum_Tchannel_M1100'  'HeavyNeutrinoToMumMum_Tchannel_M1200'  'HeavyNeutrinoToMumMum_Tchannel_M1300'  
'HeavyNeutrinoToMumMum_Tchannel_M1400'  'HeavyNeutrinoToMumMum_Tchannel_M1500'  'HeavyNeutrinoToMumMum_Tchannel_M1700'  'HeavyNeutrinoToMumMum_Tchannel_M2000'
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

declare dilebkg_long=(
'DYJets' 'TT_powheg' 'WJets'
)
declare dilepbkg_all=(
'DYJets_10to50'
'WWTo2L2Nu_DS' 'ww_ds'
'WWTo2L2Nu' 'ggWWto2L2Nu' 'ggHtoWW'
'WZTo3LNu_powheg' 'WZto2L2Q_amcatnlo' 'WZTo3LNu_amcatnlo'
'ZZTo4L_powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau' 'ggZZto2mu2nu' 'ggZZto2mu2tau' 'ggZZto4e' 'ggZZto4mu' 'ggZZto4tau'
'ZZTo2L2Nu_Powheg'
'ZZTo2L2Q_Powheg'
'vbhHtoZZ' 'ggHtoZZ'
'ZGto2LG'
'WGtoLNuG' 'WgstarToLNuEE' 'WgstarToLNuMuMu'
'TG'  'TTG'
'ttWToLNu'
'ttH_nonbb' 'ttZ' 'ttW' 'ttZToLL_M-1to10'
'WWW' 'WWZ' 'WZZ' 'ZZZ'
)

declare VVlepbkg=(
'WZTo3LNu_powheg'
'ZZTo4L_powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau' 'ggZZto2mu2nu' 'ggZZto2mu2tau' 'ggZZto4e' 'ggZZto4mu' 'ggZZto4tau' 'ggHtoZZ'
'TG' 'TTG' 'WGtoLNuG' 'ZGto2LG'
'WWW' 'WWZ' 'WZZ' 'ZZZ'
'ttW' 'ttZ' 'ttH_nonbb'
)

declare dilepbkg_SS=(
'WZTo3LNu_powheg'
'ZZTo4L_powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau' 'ggZZto2mu2nu' 'ggZZto2mu2tau' 'ggZZto4e' 'ggZZto4mu' 'ggZZto4tau' 'ggHtoZZ'
'TG' 'TTG' 'WGtoLNuG' 'ZGto2LG'
'WWW' 'WWZ' 'WZZ' 'ZZZ'
'ttW' 'ttZ' 'ttH_nonbb'
'WWTo2L2Nu_DS' 'WpWpEWK' 'WpWpQCD'
)
declare -a dilep_analysisCR_bkg=(
'WZTo3LNu_powheg'
'ZZTo4L_powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau' 'ggZZto2mu2nu' 'ggZZto2mu2tau' 'ggZZto4e' 'ggZZto4mu' 'ggZZto4tau' 'ggHtoZZ'
'TG' 'TTG' 'WGtoLNuG' 'ZGto2LG'
'WWW' 'WZZ' 'WWZ' 'ZZZ'
'ttW' 'ttZ' 'ttH_nonbb'
)


## HNEmEm_250
declare -a DiMuSig_Moriond=(
'HNMoriondLLMumMum_1100'
'HNMoriondLLMupMup_1100'
'HNMoriondLLMumMum_500'
'HNMoriondLLMupMup_500'
'HNMoriondLLMumMum_200'
'HNMoriondLLMupMup_200'
'HNMoriondLLMumMum_100'
'HNMoriondLLMupMup_100'
'HNMoriondLLMumMum_50'
'HNMoriondLLMupMup_50'
'HNMoriondLLEmEm_1100'
'HNMoriondLLEpEp_1100'
'HNMoriondLLEmEm_500'
'HNMoriondLLEpEp_500'
'HNMoriondLLEmEm_200'
'HNMoriondLLEpEp_200'
'HNMoriondLLEmEm_100'
'HNMoriondLLEpEp_100'
'HNMoriondLLEmEm_50'
'HNMoriondLLEpEp_50'
)

declare -a DiMuSig_all=(
'HNMupMup_40'  'HNMupMup_50'  'HNMupMup_60'  'HNMupMup_70'  'HNMupMup_80'  'HNMupMup_90'  'HNMupMup_100'  'HNMupMup_125'  'HNMupMup_150'  'HNMupMup_200'  'HNMupMup_250'  'HNMupMup_300'  'HNMupMup_400'  'HNMupMup_500'  'HNMupMup_600'  'HNMupMup_700'  'HNMupMup_800'  'HNMupMup_900'  'HNMupMup_1000'  'HNMupMup_1100'  'HNMupMup_1200'  'HNMupMup_1300'  'HNMupMup_1400'  'HNMupMup_1500'
'HNMumMum_40'  'HNMumMum_50'  'HNMumMum_60'  'HNMumMum_70'  'HNMumMum_80'  'HNMumMum_90'  'HNMumMum_100'  'HNMumMum_125'  'HNMumMum_150'  'HNMumMum_200'  'HNMumMum_250'  'HNMumMum_300'  'HNMumMum_400'  'HNMumMum_500'  'HNMumMum_600'  'HNMumMum_700'  'HNMumMum_800'  'HNMumMum_900'  'HNMumMum_1000'  'HNMumMum_1100'  'HNMumMum_1200'  'HNMumMum_1300'  'HNMumMum_1400'  'HNMumMum_1500'
)
#'HNMupMup_Tchannel_300'  'HNMupMup_Tchannel_400'  'HNMupMup_Tchannel_600'  'HNMupMup_Tchannel_700'  'HNMupMup_Tchannel_800'  'HNMupMup_Tchannel_900'  'HNMupMup_Tchannel_1000'  'HNMupMup_Tchannel_1200'  'HNMupMup_Tchannel_1300'  'HNMupMup_Tchannel_1400'  'HNMupMup_Tchannel_1500' 'HNMoriondLL_Tchannel_MupMup_100'  'HNMoriondLL_Tchannel_MupMup_200'  'HNMoriondLL_Tchannel_MupMup_500'  'HNMoriondLL_Tchannel_MupMup_1100'
#'HNMumMum_Tchannel_300'  'HNMumMum_Tchannel_400'  'HNMumMum_Tchannel_600'  'HNMumMum_Tchannel_700'  'HNMumMum_Tchannel_800'  'HNMumMum_Tchannel_900'  'HNMumMum_Tchannel_1000'  'HNMumMum_Tchannel_1200'  'HNMumMum_Tchannel_1300'  'HNMumMum_Tchannel_1400'  'HNMumMum_Tchannel_1500' 'HNMoriondLL_Tchannel_MumMum_100'  'HNMoriondLL_Tchannel_MumMum_200'  'HNMoriondLL_Tchannel_MumMum_500'  'HNMoriondLL_Tchannel_MumMum_1100'
#)

declare -a DiElSig_all=(
'HNEpEp_40'  'HNEpEp_50'  'HNEpEp_60'  'HNEpEp_70'  'HNEpEp_80'  'HNEpEp_90'  'HNEpEp_100'  'HNEpEp_125'  'HNEpEp_150'  'HNEpEp_200'  'HNEpEp_250'  'HNEpEp_300'  'HNEpEp_400'  'HNEpEp_500'  'HNEpEp_600'  'HNEpEp_700'  'HNEpEp_800'  'HNEpEp_900'  'HNEpEp_1000'  'HNEpEp_1100'  'HNEpEp_1200'  'HNEpEp_1300'  'HNEpEp_1400'  'HNEpEp_1500'
'HNEmEm_40'  'HNEmEm_50'  'HNEmEm_60'  'HNEmEm_70'  'HNEmEm_80'  'HNEmEm_90'  'HNEmEm_100'  'HNEmEm_125'  'HNEmEm_150'  'HNEmEm_200' 'HNEmEm_250' 'HNEmEm_300'  'HNEmEm_400'  'HNEmEm_500'  'HNEmEm_600'  'HNEmEm_700'  'HNEmEm_800'  'HNEmEm_900'  'HNEmEm_1000'  'HNEmEm_1100'  'HNEmEm_1200'  'HNEmEm_1300'  'HNEmEm_1400'  'HNEmEm_1500'
)
#'HNEpEp_Tchannel_300'  'HNEpEp_Tchannel_400'  'HNEpEp_Tchannel_600'  'HNEpEp_Tchannel_700'  'HNEpEp_Tchannel_800'  'HNEpEp_Tchannel_900'  'HNEpEp_Tchannel_1000'  'HNEpEp_Tchannel_1200'  'HNEpEp_Tchannel_1300'  'HNEpEp_Tchannel_1400'  'HNEpEp_Tchannel_1500' 'HNMoriondLL_Tchannel_EpEp_100'  'HNMoriondLL_Tchannel_EpEp_200'  'HNMoriondLL_Tchannel_EpEp_500'  'HNMoriondLL_Tchannel_EpEp_1100'
#'HNEmEm_Tchannel_300'  'HNEmEm_Tchannel_400'  'HNEmEm_Tchannel_600'  'HNEmEm_Tchannel_700'  'HNEmEm_Tchannel_800'  'HNEmEm_Tchannel_900'  'HNEmEm_Tchannel_1000'  'HNEmEm_Tchannel_1200'  'HNEmEm_Tchannel_1300'  'HNEmEm_Tchannel_1400'  'HNEmEm_Tchannel_1500' 'HNMoriondLL_Tchannel_EmEm_100'  'HNMoriondLL_Tchannel_EmEm_200'  'HNMoriondLL_Tchannel_EmEm_500'  'HNMoriondLL_Tchannel_EmEm_1100'
#)

declare -a EMuSig_all=(
'HNMupEp_40'  'HNMupEp_50'  'HNMupEp_60'  'HNMupEp_70'  'HNMupEp_80'  'HNMupEp_90'  'HNMupEp_100'  'HNMupEp_125'  'HNMupEp_150'  'HNMupEp_200'  'HNMupEp_250'  'HNMupEp_300'  'HNMupEp_400'  'HNMupEp_500'  'HNMupEp_600'  'HNMupEp_700'  'HNMupEp_800'  'HNMupEp_900'  'HNMupEp_1000'  'HNMupEp_1100'  'HNMupEp_1200'  'HNMupEp_1300'  'HNMupEp_1400'  'HNMupEp_1500'  
'HNMumEm_40'  'HNMumEm_50'  'HNMumEm_60'  'HNMumEm_70'  'HNMumEm_80'  'HNMumEm_90'  'HNMumEm_100'  'HNMumEm_125'  'HNMumEm_150'  'HNMumEm_200'  'HNMumEm_250'  'HNMumEm_300'  'HNMumEm_400'  'HNMumEm_500'  'HNMumEm_600'  'HNMumEm_700'  'HNMumEm_800'  'HNMumEm_900'  'HNMumEm_1000'  'HNMumEm_1100'  'HNMumEm_1200'  'HNMumEm_1300'  'HNMumEm_1400'  'HNMumEm_1500'  
'HNEpMup_40'  'HNEpMup_50'  'HNEpMup_60'  'HNEpMup_70'  'HNEpMup_80'  'HNEpMup_90'  'HNEpMup_100'  'HNEpMup_125'  'HNEpMup_150'  'HNEpMup_200'  'HNEpMup_250'  'HNEpMup_300'  'HNEpMup_400'  'HNEpMup_500'  'HNEpMup_600'  'HNEpMup_700'  'HNEpMup_800'  'HNEpMup_900'  'HNEpMup_1000'  'HNEpMup_1100'  'HNEpMup_1200'  'HNEpMup_1300'  'HNEpMup_1400'  'HNEpMup_1500'  
'HNEmMum_40'  'HNEmMum_50'  'HNEmMum_60'  'HNEmMum_70'  'HNEmMum_80'  'HNEmMum_90'  'HNEmMum_100'  'HNEmMum_125'  'HNEmMum_150'  'HNEmMum_200'  'HNEmMum_250'  'HNEmMum_300'  'HNEmMum_400'  'HNEmMum_500'  'HNEmMum_600'  'HNEmMum_700'  'HNEmMum_800'  'HNEmMum_900'  'HNEmMum_1000'  'HNEmMum_1100'  'HNEmMum_1200'  'HNEmMum_1300'  'HNEmMum_1400'  'HNEmMum_1500'
)
#'HNMupEp_Tchannel_300'  'HNMupEp_Tchannel_400'  'HNMupEp_Tchannel_600'  'HNMupEp_Tchannel_700'  'HNMupEp_Tchannel_800'  'HNMupEp_Tchannel_900'  'HNMupEp_Tchannel_1000'  'HNMupEp_Tchannel_1200'  'HNMupEp_Tchannel_1300'  'HNMupEp_Tchannel_1400'  'HNMupEp_Tchannel_1500'  
#'HNMoriondLL_Tchannel_MupEp_100'  'HNMoriondLL_Tchannel_MupEp_200'  'HNMoriondLL_Tchannel_MupEp_500'  'HNMoriondLL_Tchannel_MupEp_1100'  
#'HNMumEm_Tchannel_300'  'HNMumEm_Tchannel_400'  'HNMumEm_Tchannel_600'  'HNMumEm_Tchannel_700'  'HNMumEm_Tchannel_800'  'HNMumEm_Tchannel_900'  'HNMumEm_Tchannel_1000'  'HNMumEm_Tchannel_1200'  'HNMumEm_Tchannel_1300'  'HNMumEm_Tchannel_1400'  'HNMumEm_Tchannel_1500'  
#'HNMoriondLL_Tchannel_MumEm_100'  'HNMoriondLL_Tchannel_MumEm_200'  'HNMoriondLL_Tchannel_MumEm_500'  'HNMoriondLL_Tchannel_MumEm_1100'  
#'HNEpMup_Tchannel_300'  'HNEpMup_Tchannel_400'  'HNEpMup_Tchannel_600'  'HNEpMup_Tchannel_700'  'HNEpMup_Tchannel_800'  'HNEpMup_Tchannel_900'  'HNEpMup_Tchannel_1000'  'HNEpMup_Tchannel_1200'  'HNEpMup_Tchannel_1300'  'HNEpMup_Tchannel_1400'  'HNEpMup_Tchannel_1500'  
#'HNMoriondLL_Tchannel_EpMup_100'  'HNMoriondLL_Tchannel_EpMup_200'  'HNMoriondLL_Tchannel_EpMup_500'  'HNMoriondLL_Tchannel_EpMup_1100'  
#'HNEmMum_Tchannel_300'  'HNEmMum_Tchannel_400'  'HNEmMum_Tchannel_600'  'HNEmMum_Tchannel_700'  'HNEmMum_Tchannel_800'  'HNEmMum_Tchannel_900'  'HNEmMum_Tchannel_1000'  'HNEmMum_Tchannel_1200'  'HNEmMum_Tchannel_1300'  'HNEmMum_Tchannel_1400'  'HNEmMum_Tchannel_1500'  
#'HNMoriondLL_Tchannel_EmMum_100'  'HNMoriondLL_Tchannel_EmMum_200'  'HNMoriondLL_Tchannel_EmMum_500'  'HNMoriondLL_Tchannel_EmMum_1100'  
#)

declare -a EMuSig_mue=(
'HNMupEp_40'  'HNMupEp_50'  'HNMupEp_60'  'HNMupEp_70'  'HNMupEp_80'  'HNMupEp_90'  'HNMupEp_100'  'HNMupEp_125'  'HNMupEp_150'  'HNMupEp_200'  'HNMupEp_250'  'HNMupEp_300'  'HNMupEp_400'  'HNMupEp_500'  'HNMupEp_600'  'HNMupEp_700'  'HNMupEp_800'  'HNMupEp_900'  'HNMupEp_1000'  'HNMupEp_1100'  'HNMupEp_1200'  'HNMupEp_1300'  'HNMupEp_1400'  'HNMupEp_1500'  
'HNMumEm_40'  'HNMumEm_50'  'HNMumEm_60'  'HNMumEm_70'  'HNMumEm_80'  'HNMumEm_90'  'HNMumEm_100'  'HNMumEm_125'  'HNMumEm_150'  'HNMumEm_200'  'HNMumEm_250'  'HNMumEm_300'  'HNMumEm_400'  'HNMumEm_500'  'HNMumEm_600'  'HNMumEm_700'  'HNMumEm_800'  'HNMumEm_900'  'HNMumEm_1000'  'HNMumEm_1100'  'HNMumEm_1200'  'HNMumEm_1300'  'HNMumEm_1400'  'HNMumEm_1500'  
)
#'HNMupEp_Tchannel_300'  'HNMupEp_Tchannel_400'  'HNMupEp_Tchannel_600'  'HNMupEp_Tchannel_700'  'HNMupEp_Tchannel_800'  'HNMupEp_Tchannel_900'  'HNMupEp_Tchannel_1000'  'HNMupEp_Tchannel_1200'  'HNMupEp_Tchannel_1300'  'HNMupEp_Tchannel_1400'  'HNMupEp_Tchannel_1500'  
#'HNMoriondLL_Tchannel_MupEp_100'  'HNMoriondLL_Tchannel_MupEp_200'  'HNMoriondLL_Tchannel_MupEp_500'  'HNMoriondLL_Tchannel_MupEp_1100'  
#'HNMumEm_Tchannel_300'  'HNMumEm_Tchannel_400'  'HNMumEm_Tchannel_600'  'HNMumEm_Tchannel_700'  'HNMumEm_Tchannel_800'  'HNMumEm_Tchannel_900'  'HNMumEm_Tchannel_1000'  'HNMumEm_Tchannel_1200'  'HNMumEm_Tchannel_1300'  'HNMumEm_Tchannel_1400'  'HNMumEm_Tchannel_1500'  
#'HNMoriondLL_Tchannel_MumEm_100'  'HNMoriondLL_Tchannel_MumEm_200'  'HNMoriondLL_Tchannel_MumEm_500'  'HNMoriondLL_Tchannel_MumEm_1100'  
#)

declare -a EMuSig_emu=(
'HNEpMup_40'  'HNEpMup_50'  'HNEpMup_60'  'HNEpMup_70'  'HNEpMup_80'  'HNEpMup_90'  'HNEpMup_100'  'HNEpMup_125'  'HNEpMup_150'  'HNEpMup_200'  'HNEpMup_250'  'HNEpMup_300'  'HNEpMup_400'  'HNEpMup_500'  'HNEpMup_600'  'HNEpMup_700'  'HNEpMup_800'  'HNEpMup_900'  'HNEpMup_1000'  'HNEpMup_1100'  'HNEpMup_1200'  'HNEpMup_1300'  'HNEpMup_1400'  'HNEpMup_1500'  
'HNEmMum_40'  'HNEmMum_50'  'HNEmMum_60'  'HNEmMum_70'  'HNEmMum_80'  'HNEmMum_90'  'HNEmMum_100'  'HNEmMum_125'  'HNEmMum_150'  'HNEmMum_200'  'HNEmMum_250'  'HNEmMum_300'  'HNEmMum_400'  'HNEmMum_500'  'HNEmMum_600'  'HNEmMum_700'  'HNEmMum_800'  'HNEmMum_900'  'HNEmMum_1000'  'HNEmMum_1100'  'HNEmMum_1200'  'HNEmMum_1300'  'HNEmMum_1400'  'HNEmMum_1500'
)
#'HNEpMup_Tchannel_300'  'HNEpMup_Tchannel_400'  'HNEpMup_Tchannel_600'  'HNEpMup_Tchannel_700'  'HNEpMup_Tchannel_800'  'HNEpMup_Tchannel_900'  'HNEpMup_Tchannel_1000'  'HNEpMup_Tchannel_1200'  'HNEpMup_Tchannel_1300'  'HNEpMup_Tchannel_1400'  'HNEpMup_Tchannel_1500'  
#'HNMoriondLL_Tchannel_EpMup_100'  'HNMoriondLL_Tchannel_EpMup_200'  'HNMoriondLL_Tchannel_EpMup_500'  'HNMoriondLL_Tchannel_EpMup_1100'  
#'HNEmMum_Tchannel_300'  'HNEmMum_Tchannel_400'  'HNEmMum_Tchannel_600'  'HNEmMum_Tchannel_700'  'HNEmMum_Tchannel_800'  'HNEmMum_Tchannel_900'  'HNEmMum_Tchannel_1000'  'HNEmMum_Tchannel_1200'  'HNEmMum_Tchannel_1300'  'HNEmMum_Tchannel_1400'  'HNEmMum_Tchannel_1500'  
#'HNMoriondLL_Tchannel_EmMum_100'  'HNMoriondLL_Tchannel_EmMum_200'  'HNMoriondLL_Tchannel_EmMum_500'  'HNMoriondLL_Tchannel_EmMum_1100'  
#)

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



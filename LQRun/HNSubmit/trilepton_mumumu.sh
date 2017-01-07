#!/bin/bash

##### SR #####

sktree_bkg -a trilepton_mumumu -S DoubleMuon -s SKTree_DiLepSkim -n 100
sktree_bkg -a trilepton_mumumu -list trilep_bkg -s SKTree_DiLepSkim -n 100
sktree_bkg -a trilepton_mumumu -list trilep_signal -s SKTree_LeptonSkim -n 100
sktree_bkg -a trilepton_mumumu_FR_method -S DoubleMuon -s SKTree_DiLepSkim -n 100
sktree_bkg -a trilepton_mumumu -list trilep_nonprompt_bkg -s SKTree_DiLepSkim -n 100 -o /data2/CAT_SKTreeOutput/JobOutPut/jskim/LQanalyzer/data/output/CAT/trilepton_mumumu/periodBtoG/non_prompt_MC/
sktree_bkg -a trilepton_mumumu -list trilep_diboson_had -s SKTree_DiLepSkim -n 100 -o /data2/CAT_SKTreeOutput/JobOutPut/jskim/LQanalyzer/data/output/CAT/trilepton_mumumu/periodBtoG/non_prompt_MC/ -userflag diboson_had 

##### CR #####

sktree_bkg -a trilepton_mumumu_CR -S DoubleMuon -s SKTree_DiLepSkim -n 100
sktree_bkg -a trilepton_mumumu_CR -list trilep_CR_bkg -s SKTree_DiLepSkim -n 100
sktree_bkg -a trilepton_mumumu_CR -list trilep_signal -s SKTree_LeptonSkim -n 100
sktree_bkg -a trilepton_mumumu_CR_FR_method -S DoubleMuon -s SKTree_DiLepSkim -n 100

##### FR MC Closure #####
sktree_bkg -a trilepton_mumumu_CR -list FR_MC_Closure -s SKTree_DiLepSkim -n 100 -userflag MCClosure -o /data2/CAT_SKTreeOutput/JobOutPut/jskim/LQanalyzer/data/output/CAT/FR_MC_Closure/periodBtoG/
sktree_bkg -a trilepton_mumumu_CR_FR_method -list FR_MC_Closure -s SKTree_DiLepSkim -n 100 -userflag MCClosure -o /data2/CAT_SKTreeOutput/JobOutPut/jskim/LQanalyzer/data/output/CAT/FR_MC_Closure/periodBtoG/
sktree_bkg -a trilepton_mumumu_syst_FR -list FR_MC_Closure -s SKTree_DiLepSkim -n 100 -userflag MCClosure -o /data2/CAT_SKTreeOutput/JobOutPut/jskim/LQanalyzer/data/output/CAT/FR_MC_Closure/periodBtoG/syst_FR/

##### Fake yield systematic #####
sktree_bkg -a trilepton_mumumu_syst_FR -S DoubleMuon -s SKTree_DiLepSkim -n 100

##### Cut Op #####
sktree_bkg -a trilepton_mumumu -S DoubleMuon -s SKTree_DiLepSkim -n 100 -userflag cutop -o /data2/CAT_SKTreeOutput/JobOutPut/jskim/LQanalyzer/data/output/CAT/CutOp/periodBtoG/Data/
sktree_bkg -a trilepton_mumumu -list trilep_bkg -s SKTree_DiLepSkim -n 100 -userflag cutop -o /data2/CAT_SKTreeOutput/JobOutPut/jskim/LQanalyzer/data/output/CAT/CutOp/periodBtoG/
sktree_bkg -a trilepton_mumumu -list trilep_signal -s SKTree_LeptonSkim -n 100 -userflag cutop -o /data2/CAT_SKTreeOutput/JobOutPut/jskim/LQanalyzer/data/output/CAT/CutOp/periodBtoG/
sktree_bkg -a trilepton_mumumu_FR_method -S DoubleMuon -s SKTree_DiLepSkim -n 100 -userflag cutop -o /data2/CAT_SKTreeOutput/JobOutPut/jskim/LQanalyzer/data/output/CAT/CutOp/periodBtoG/Fake/

##### Syst (up down) #####
sktree_bkg -a trilepton_mumumu_ntp -S DoubleMuon -s SKTree_DiLepSkim -n 100
sktree_bkg -a trilepton_mumumu_ntp -list trilep_bkg -s SKTree_DiLepSkim -n 100
sktree_bkg -a trilepton_mumumu_ntp -list trilep_signal -s SKTree_LeptonSkim -n 100
sktree_bkg -a trilepton_mumumu_ntp_FR_method -S DoubleMuon -s SKTree_DiLepSkim -n 100

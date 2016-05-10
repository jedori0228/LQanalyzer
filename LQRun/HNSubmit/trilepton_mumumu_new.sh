#!/bin/bash
sktree -a trilepton_mumumu -S DoubleMuon -s SKTree_DiLepSkim -n 15
sktree -a trilepton_mumumu -list trilep_bkg -s SKTree_DiLepSkim -n 15
sktree -a trilepton_mumumu -list trilep_signal -s SKTree_LeptonSkim -n 15


sktree -a trilepton_mumumu -list trilep_bkg -s SKTree_TriLepSkim -n 15 -o /data2/CAT_SKTreeOutput/JobOutPut/jskim/LQanalyzer/data/output/CAT/trilepton_mumumu/periodCtoD/dXY_0p01_dZ_0p5/

sktree -a trilepton_mumumu -list trilep_signal -s SKTree_LeptonSkim -n 15 -o /data2/CAT_SKTreeOutput/JobOutPut/jskim/LQanalyzer/data/output/CAT/trilepton_mumumu/periodCtoD/dXY_0p01_dZ_0p5/

sktree -a trilepton_mumumu -list hn_alp_lll_mm -s SKTree_LeptonSkim -n 15 -o /data2/CAT_SKTreeOutput/JobOutPut/jskim/LQanalyzer/data/output/CAT/trilepton_mumumu/periodCtoD/dXY_0p01_dZ_0p5/

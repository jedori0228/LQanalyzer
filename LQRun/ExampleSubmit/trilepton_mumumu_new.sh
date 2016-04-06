#!/bin/bash
sktree -a trilepton_mumumu -S DoubleMuon -s SKTree_DiLepSkim -n 15
sktree -a trilepton_mumumu -list trilep_bkg -s SKTree_DiLepSkim -n 15
sktree -a trilepton_mumumu -list trilep_signal -s SKTree_LeptonSkim -n 15

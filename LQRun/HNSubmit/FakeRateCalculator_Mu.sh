#!/bin/bash
sktree_bkg -a FakeRateCalculator_Mu -S DoubleMuon -s SKTree_LeptonSkim -n 100
sktree_bkg -a FakeRateCalculator_Mu -list trilep_fake_bkg -s SKTree_LeptonSkim -n 100
sktree_bkg -a FakeRateCalculator_Mu -list QCD_FR -s SKTree_LeptonSkim -n 100
sktree_bkg -a FakeRateCalculator_Mu -list trilep_signal -s SKTree_LeptonSkim -n 100

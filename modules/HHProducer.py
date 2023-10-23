import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os
import numpy as np
from numpy import sign
from numpy import argsort

class HHProducer(Module):
  def __init__(self , year):
    self.year = year
  def beginJob(self):
    pass
  def endJob(self):
    pass
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch("HLT_passEle32WPTight", "I")
    self.out.branch("lhe_nlepton", "I")
    self.out.branch("nHad_tau", "I")
    self.out.branch("met_user","F")
    self.out.branch("met_phi_user","F")

    self.out.branch("TightFatJet_id","I",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_initid","I",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_bbvsQCD","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_ccvsQCD","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_qqvsQCD","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_drl1","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_drl2","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_pt","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_eta","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_phi","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_mass","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_dRinit","F",lenVar="nTightFatJet")
    self.out.branch("HZ_2F0R", "B")
    self.out.branch("HZ_2F0R_case", "I")
    self.out.branch("HZ_1F2R", "B")
    self.out.branch("HZ_1F2R_case", "I")
    self.out.branch("HZ_0F4R", "B")
    self.out.branch("HZ_0F4R_case", "I")
    self.out.branch("H_AK8Jet_id", "I")
    self.out.branch("H_AK8Jet_initid", "I")
    self.out.branch("H_AK8Jet_pt", "F")
    self.out.branch("H_AK8Jet_eta", "F")
    self.out.branch("H_AK8Jet_phi", "F")
    self.out.branch("H_AK8Jet_mass", "F")
    self.out.branch("H_AK8Jet_PNmass", "F")
    self.out.branch("H_AK8Jet_SDmass", "F")
    self.out.branch("H_AK8Jet_drl1", "F")
    self.out.branch("H_AK8Jet_drl2", "F")
    self.out.branch("Z_AK8Jet_id", "I")
    self.out.branch("Z_AK8Jet_initid", "I")
    self.out.branch("Z_AK8Jet_pt", "F")
    self.out.branch("Z_AK8Jet_eta", "F")
    self.out.branch("Z_AK8Jet_phi", "F")
    self.out.branch("Z_AK8Jet_mass", "F")
    self.out.branch("Z_AK8Jet_drl1", "F")
    self.out.branch("Z_AK8Jet_drl2", "F")

    self.out.branch("TightAK4JetCRall_id","I",lenVar="nTightAK4JetCRall")
    self.out.branch("TightAK4JetCRall_b_DeepJetL_id","I",lenVar="nTightAK4JetCRall")
    self.out.branch("TightAK4JetCRall_b_DeepJetM_id","I",lenVar="nTightAK4JetCRall")
    self.out.branch("nTightAK4JetCRall_b_DeepJetL","I")
    self.out.branch("nTightAK4JetCRall_b_DeepJetM","I")

    self.out.branch("TT_region_Lbb","B")
    self.out.branch("TT_Lb1_pt","F")
    self.out.branch("TT_Lb1_eta","F")
    self.out.branch("TT_Lb1_phi","F")
    self.out.branch("TT_Lb1_mass","F")
    self.out.branch("TT_Lb2_pt","F")
    self.out.branch("TT_Lb2_eta","F")
    self.out.branch("TT_Lb2_phi","F")
    self.out.branch("TT_Lb2_mass","F")
    self.out.branch("TT_Lmbb","F")
    self.out.branch("TT_Lb1_drl1","F")
    self.out.branch("TT_Lb1_drl2","F")
    self.out.branch("TT_Lb2_drl1","F")
    self.out.branch("TT_Lb2_drl2","F")

    self.out.branch("TT_region_Mbb","B")
    self.out.branch("TT_Mb1_pt","F")
    self.out.branch("TT_Mb1_eta","F")
    self.out.branch("TT_Mb1_phi","F")
    self.out.branch("TT_Mb1_mass","F")
    self.out.branch("TT_Mb2_pt","F")
    self.out.branch("TT_Mb2_eta","F")
    self.out.branch("TT_Mb2_phi","F")
    self.out.branch("TT_Mb2_mass","F")
    self.out.branch("TT_Mmbb","F")
    self.out.branch("TT_Mb1_drl1","F")
    self.out.branch("TT_Mb1_drl2","F")
    self.out.branch("TT_Mb2_drl1","F")
    self.out.branch("TT_Mb2_drl2","F")

    self.out.branch("DY2J_region","B")
    self.out.branch("DY2J_region_nob","B")

    self.out.branch("DY_j1_pt","F")
    self.out.branch("DY_j1_eta","F")
    self.out.branch("DY_j1_phi","F")
    self.out.branch("DY_j1_mass","F")
    self.out.branch("DY_j2_pt","F")
    self.out.branch("DY_j2_eta","F")
    self.out.branch("DY_j2_phi","F")
    self.out.branch("DY_j2_mass","F")
    self.out.branch("DY_mjj","F")
    self.out.branch("DY_j1_drl1","F")
    self.out.branch("DY_j1_drl2","F")
    self.out.branch("DY_j2_drl1","F")
    self.out.branch("DY_j2_drl2","F")

    self.out.branch("TightAK4Jet_pt","F",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_eta","F",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_phi","F",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_mass","F",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_id","I",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_initid","I",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_nob_id","I",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_b_DeepJetM_id","I",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_b_DeepJetL_id","I",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_drl1","F",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_drl2","F",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_dRinit","F",lenVar="nTightAK4Jet")
    self.out.branch("n_bjet_DeepB_M", "I")
    self.out.branch("n_bjet_DeepB_L", "I")
    self.out.branch("n_tight_nob", "I")
    self.out.branch("HT", "F")
    self.out.branch("Had_tau_id","I",lenVar="nTau")

    self.out.branch("h_j1_pt", "F")
    self.out.branch("h_j1_eta", "F")
    self.out.branch("h_j1_phi", "F")
    self.out.branch("h_j1_mass", "F")
    self.out.branch("h_j1_id", "I")
    self.out.branch("h_j1_drl1", "F")
    self.out.branch("h_j1_drl2", "F")
    self.out.branch("h_j2_pt", "F")
    self.out.branch("h_j2_eta", "F")
    self.out.branch("h_j2_phi", "F")
    self.out.branch("h_j2_mass", "F")
    self.out.branch("h_j2_id", "I")
    self.out.branch("h_j2_drl1", "F")
    self.out.branch("h_j2_drl2", "F")
    self.out.branch("h_mjj", "F")
    self.out.branch("h_detajj", "F")
    self.out.branch("h_dRjj", "F")
    self.out.branch("h_dphijj", "F")
    self.out.branch("zhad_j1_pt", "F",lenVar="nzhad")
    self.out.branch("zhad_j1_eta", "F",lenVar="nzhad")
    self.out.branch("zhad_j1_phi", "F",lenVar="nzhad")
    self.out.branch("zhad_j1_mass", "F",lenVar="nzhad")
    self.out.branch("zhad_j1_id", "I",lenVar="nzhad")
    self.out.branch("zhad_j2_pt", "F",lenVar="nzhad")
    self.out.branch("zhad_j2_eta", "F",lenVar="nzhad")
    self.out.branch("zhad_j2_phi", "F",lenVar="nzhad")
    self.out.branch("zhad_j2_mass", "F",lenVar="nzhad")
    self.out.branch("zhad_j2_id", "I",lenVar="nzhad")
    self.out.branch("zhad_mjj", "F",lenVar="nzhad")
    self.out.branch("zhad_detajj", "F",lenVar="nzhad")
    self.out.branch("zhad_dRjj", "F",lenVar="nzhad")
    self.out.branch("zhad_dphijj", "F",lenVar="nzhad")
    self.out.branch("zhad_zlep_mass", "F",lenVar="nzhad")

    # For H/Y 4 momentum, in 2F0R case and ZF2R(Z-fatjet and 2 bjets from Higgs), H/Y 4-p can be determined directly. In HF2R case and 0F4R case, H/Y 4-p depends on the choise of two jets from Z boson
    self.out.branch("Y_2F0R_pt", "F")
    self.out.branch("Y_2F0R_eta", "F")
    self.out.branch("Y_2F0R_phi", "F")
    self.out.branch("Y_2F0R_mass", "F")
    self.out.branch("X_2F0R_pt", "F")
    self.out.branch("X_2F0R_eta", "F")
    self.out.branch("X_2F0R_phi", "F")
    self.out.branch("X_2F0R_mass", "F")

    self.out.branch("Y_ZF2R_pt", "F")
    self.out.branch("Y_ZF2R_eta", "F")
    self.out.branch("Y_ZF2R_phi", "F")
    self.out.branch("Y_ZF2R_mass", "F")
    self.out.branch("X_ZF2R_pt", "F")
    self.out.branch("X_ZF2R_eta", "F")
    self.out.branch("X_ZF2R_phi", "F")
    self.out.branch("X_ZF2R_mass", "F")

    self.out.branch("YZQQ_pt", "F",lenVar="nYZQQ")
    self.out.branch("YZQQ_eta", "F",lenVar="nYZQQ")
    self.out.branch("YZQQ_phi", "F",lenVar="nYZQQ")
    self.out.branch("YZQQ_mass", "F",lenVar="nYZQQ")
    self.out.branch("XZQQ_pt", "F",lenVar="nXZQQ")
    self.out.branch("XZQQ_eta", "F",lenVar="nXZQQ")
    self.out.branch("XZQQ_phi", "F",lenVar="nXZQQ")
    self.out.branch("XZQQ_mass", "F",lenVar="nXZQQ")

    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))
    self.is_lhe = bool(inputTree.GetBranch("nLHEPart"))
    self.is_signal = bool(inputTree.GetBranch("GEN_zj1_pdgid"))

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    
    # PV selection
    if (event.PV_npvsGood<1): return False
    if self.is_signal and abs(event.GEN_zj1_pdgid)==5: return False

    # trigger selection
    # special action for 2017 single ele HLT, https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2017DataRecommendations#Single_Electron_Triggers
    # TrigObj_filterBits, 1024 = 1e (32_L1DoubleEG_AND_L1SingleEGOr)
    HLT_passEle32WPTight=0
    if self.year=="2017":
      trgobjs=Collection(event, 'TrigObj')
      if event.HLT_Ele32_WPTight_Gsf_L1DoubleEG==1:
	for iobj in range(0,event.nTrigObj):
	  if trgobjs[iobj].id==11 and (trgobjs[iobj].filterBits & (1<<10))== (1<<10):
	    HLT_passEle32WPTight=1

    self.out.fillBranch("HLT_passEle32WPTight",HLT_passEle32WPTight)

    lhe_nlepton=0
    if self.is_lhe:
      lheparticle = Collection(event, 'LHEPart')
      for ilhe in range(0, event.nLHEPart):
        if lheparticle[ilhe].status==1 and (abs(lheparticle[ilhe].pdgId)==11 or abs(lheparticle[ilhe].pdgId)==13 or abs(lheparticle[ilhe].pdgId)==15):
          lhe_nlepton=lhe_nlepton+1

    self.out.fillBranch("lhe_nlepton", lhe_nlepton)

    l1v4_tmp=TLorentzVector()
    l2v4_tmp=TLorentzVector()
    l1v4_tmp.SetPtEtaPhiM(event.l1_pt,event.l1_eta,event.l1_phi,event.l1_mass)
    l2v4_tmp.SetPtEtaPhiM(event.l2_pt,event.l2_eta,event.l2_phi,event.l2_mass)

    tau_v4_temp=TLorentzVector()
    taus = Collection(event, 'Tau')
    nHad_tau=0
    Had_tau_id=[]
    for itau in range(0, event.nTau):
      tau_v4_temp.SetPtEtaPhiM(taus[itau].pt, taus[itau].eta, taus[itau].phi, taus[itau].mass)
      pass_tau_lep_Dr=1
      if taus[itau].pt>20 and abs(taus[itau].eta)<2.3 and taus[itau].idDecayModeOldDMs and taus[itau].idDeepTau2017v2p1VSe>=4 and taus[itau].idDeepTau2017v2p1VSjet>=4 and taus[itau].idDeepTau2017v2p1VSmu>=1:
        if tau_v4_temp.DeltaR(l1v4_tmp)<0.4:pass_tau_lep_Dr=0
        if tau_v4_temp.DeltaR(l2v4_tmp)<0.4:pass_tau_lep_Dr=0
	if pass_tau_lep_Dr:
	  nHad_tau=nHad_tau+1
	  Had_tau_id.append(itau)
    self.out.fillBranch("nHad_tau", nHad_tau)

    met_user=-99
    met_phi_user=-99

    if self.is_mc:
      met_user=event.MET_T1Smear_pt
      met_phi_user=event.MET_T1Smear_phi
    else:
      met_user=event.MET_T1_pt
      met_phi_user=event.MET_T1_phi

    self.out.fillBranch("met_user",met_user)
    self.out.fillBranch("met_phi_user",met_phi_user)


    fatjet_v4_temp=TLorentzVector()
    # Fat jet selection
    # due to the potential overlap between jet and lepton, Tight ID is used. (*10>=2)
    fatjets = Collection(event, 'FatJetNoVlep')
    fatjets_init = Collection(event, 'FatJet')

    TightFatJet_id = []
    TightFatJet_initid = []
    TightFatJet_bbvsQCD = []
    TightFatJet_ccvsQCD = []
    TightFatJet_qqvsQCD = []
    TightFatJet_drl1 = []
    TightFatJet_drl2 = []
    TightFatJet_pt = []
    TightFatJet_eta = []
    TightFatJet_phi = []
    TightFatJet_mass = []
    TightFatJet_dRinit = []
    TightFatJet_v4_all = []
    # local jet id, with dr(lepton)>0.8 or <0.8
    TightFatJet_drLa08_id = []
    TightFatJet_drSm08_id = []

    nFatjet_MP=-1

    for ijet in range(0, event.nFatJetNoVlep):
      if abs(fatjets[ijet].eta)>2.4 :continue
      if not (fatjets[ijet].jetId>1 and fatjets[ijet].pt>300):continue
      # default jet mass is FatJet_mass
      fatjet_v4_temp.SetPtEtaPhiM(fatjets[ijet].pt_nom,fatjets[ijet].eta,fatjets[ijet].phi,fatjets[ijet].mass_nom)
      # order of jet in LepSubtraction FatJet collection
      TightFatJet_id.append(ijet)
      TightFatJet_dRinit.append(fatjets[ijet].dRinit)
      # order of jet in initial NanoAOD FatJet collection
      TightFatJet_initid.append(fatjets[ijet].id)
      TightFatJet_pt.append(fatjets[ijet].pt_nom)
      TightFatJet_eta.append(fatjets[ijet].eta)
      TightFatJet_phi.append(fatjets[ijet].phi)
      TightFatJet_mass.append(fatjets[ijet].mass_nom)
      if (fatjets_init[fatjets[ijet].id].particleNetMD_Xbb + fatjets_init[fatjets[ijet].id].particleNetMD_QCD)>0:
        TightFatJet_bbvsQCD.append(fatjets_init[fatjets[ijet].id].particleNetMD_Xbb/(fatjets_init[fatjets[ijet].id].particleNetMD_Xbb + fatjets_init[fatjets[ijet].id].particleNetMD_QCD))
      else:
        TightFatJet_bbvsQCD.append(-1.)
      if (fatjets_init[fatjets[ijet].id].particleNetMD_Xcc + fatjets_init[fatjets[ijet].id].particleNetMD_QCD)>0:
        TightFatJet_ccvsQCD.append(fatjets_init[fatjets[ijet].id].particleNetMD_Xcc/(fatjets_init[fatjets[ijet].id].particleNetMD_Xcc + fatjets_init[fatjets[ijet].id].particleNetMD_QCD))
      else:
        TightFatJet_ccvsQCD.append(-1.)
      if (fatjets_init[fatjets[ijet].id].particleNetMD_Xcc + fatjets_init[fatjets[ijet].id].particleNetMD_Xqq + fatjets_init[fatjets[ijet].id].particleNetMD_QCD)>0:
        TightFatJet_qqvsQCD.append((fatjets_init[fatjets[ijet].id].particleNetMD_Xcc + fatjets_init[fatjets[ijet].id].particleNetMD_Xqq)/(fatjets_init[fatjets[ijet].id].particleNetMD_Xcc + fatjets_init[fatjets[ijet].id].particleNetMD_Xqq + fatjets_init[fatjets[ijet].id].particleNetMD_QCD))
      else:
        TightFatJet_qqvsQCD.append(-1.)
      # https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/RecoBTag/ONNXRuntime/python/pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags_cfi.py#L74
#      fatJets_qqvsQCD.append((fatjets[ijet].deepTagMD_ZvsQCD)
      TightFatJet_drl1.append(fatjets[ijet].drl1)
      TightFatJet_drl2.append(fatjets[ijet].drl2)
      TightFatJet_v4_all.append(fatjet_v4_temp.Clone())

    nFatjet_MP=len(TightFatJet_id)
    if nFatjet_MP>0:
      for ij_tmp in range(nFatjet_MP):
        # e.g., if TightFatJet_id=[0,2,3], it mean the 1th,3rd,4th jet in FatJetNoVLep pass the MP
        if TightFatJet_drl1[ij_tmp]>0.8 and TightFatJet_drl2[ij_tmp]>0.8:
          TightFatJet_drLa08_id.append(ij_tmp)
        else:
          TightFatJet_drSm08_id.append(ij_tmp)

    # check the fat jet property
    # for ParticleNetMD, WP for bbVsQCD, in BTV-22-001, HP,MP,LP correspond to signal eff 40%, 60% and 80%
    # 2016apv pre-VFP HP:0.9883, MP:0.9737, LP:0.9088
    # 2016 post-VFP HP:0.9883, MP:0.9735, LP:0.9137
    # 2017 HP:0.9870, MP:0.9714, LP:0.9105
    # 2018 HP:0.9880, MP:0.9734, LP:0.9172
    Hbb_th=-1
    if self.year=="2016apv":Hbb_threshold=0.9088
    if self.year=="2016":Hbb_threshold=0.9137
    if self.year=="2017":Hbb_threshold=0.9105
    if self.year=="2018":Hbb_threshold=0.9172

    # for ParticleNetMD, WP for ccVsQCD, in BTV-22-001, HP,MP,LP correspond to signal eff 40%, 60% and 80%
    # 2016apv pre-VFP HP:0.9909, MP:0.9751, LP:0.9252
    # 2016 post-VFP HP:0.9905, MP:0.9743, LP:0.9252
    # 2017 HP:0.9909, MP:0.9765, LP:0.9347
    # 2018 HP:0.9917, MP:0.9777 LP:0.9368
    Hcc_th=-1
    if self.year=="2016apv":Hcc_threshold=0.9252
    if self.year=="2016":Hcc_threshold=0.9252
    if self.year=="2017":Hcc_threshold=0.9252
    if self.year=="2018":Hcc_threshold=0.9368
    
    # four categories (for H and Z decay), 2F0R: 2 fat jets, 1F2R: 1 fat and 2 resolved jets, 0F4R: 4 resolved jets
    HZ_2F0R=False
    HZ_2F0R_case=-1
    HZ_1F2R=False
    HZ_1F2R_case=-1
    HZ_0F4R=False
    HZ_0F4R_case=-1

    H_AK8Jet_id=-1
    H_AK8Jet_initid=-1
    H_AK8Jet_pt=-99
    H_AK8Jet_eta=-99
    H_AK8Jet_phi=-99
    #FatJet_mass
    H_AK8Jet_mass=-99
    #FatJet_particleNet_mass
    H_AK8Jet_PNmass=-99
    #FatJet_msoftdrop
    H_AK8Jet_SDmass=-99
    H_AK8Jet_drl1=-99
    H_AK8Jet_drl2=-99

    Z_AK8Jet_id=-1
    Z_AK8Jet_initid=-1
    Z_AK8Jet_pt=-99
    Z_AK8Jet_eta=-99
    Z_AK8Jet_phi=-99
    Z_AK8Jet_mass=-99
    Z_AK8Jet_drl1=-99
    Z_AK8Jet_drl2=-99

    # no fat jet passing tight ID
    if nFatjet_MP==0: 
      HZ_0F4R=True
    
    # only one fat jet passing tight ID
    # if it pass the MP then it's a Higgs->bb fat jet, otherwise make it to be the Z->qq fat jet
    elif nFatjet_MP==1:
      HZ_1F2R=True
      if TightFatJet_drl1[0]>0.8 and TightFatJet_drl2[0]>0.8:
        if TightFatJet_bbvsQCD[0]>Hbb_threshold and TightFatJet_ccvsQCD[0]<Hcc_threshold:
          HZ_1F2R_case=0
          H_AK8Jet_id=TightFatJet_id[0]
          H_AK8Jet_initid=TightFatJet_initid[0]
          H_AK8Jet_pt=TightFatJet_v4_all[0].Pt()
          H_AK8Jet_eta=TightFatJet_v4_all[0].Eta()
          H_AK8Jet_phi=TightFatJet_v4_all[0].Phi()
          H_AK8Jet_mass=TightFatJet_v4_all[0].M()
          H_AK8Jet_PNmass=event.FatJet_particleNet_mass[H_AK8Jet_initid]
          H_AK8Jet_SDmass=event.FatJet_msoftdrop[H_AK8Jet_initid]
          H_AK8Jet_drl1=TightFatJet_drl1[0]
          H_AK8Jet_drl2=TightFatJet_drl2[0]
        elif event.FatJetNoVlep_bbvsccqq[TightFatJet_id[0]]<0.5:
          HZ_1F2R_case=1
          Z_AK8Jet_id=TightFatJet_id[0]
          Z_AK8Jet_initid=TightFatJet_initid[0]
          Z_AK8Jet_pt=TightFatJet_v4_all[0].Pt()
          Z_AK8Jet_eta=TightFatJet_v4_all[0].Eta()
          Z_AK8Jet_phi=TightFatJet_v4_all[0].Phi()
          Z_AK8Jet_mass=TightFatJet_v4_all[0].M()
          Z_AK8Jet_drl1=TightFatJet_drl1[0]
          Z_AK8Jet_drl2=TightFatJet_drl2[0]
      else:
        HZ_1F2R_case=2
        Z_AK8Jet_id=TightFatJet_id[0]
        Z_AK8Jet_initid=TightFatJet_initid[0]
        Z_AK8Jet_pt=TightFatJet_v4_all[0].Pt()
        Z_AK8Jet_eta=TightFatJet_v4_all[0].Eta()
        Z_AK8Jet_phi=TightFatJet_v4_all[0].Phi()
        Z_AK8Jet_mass=TightFatJet_v4_all[0].M()
        Z_AK8Jet_drl1=TightFatJet_drl1[0]
        Z_AK8Jet_drl2=TightFatJet_drl2[0]


    # two fat jets passing tight ID
    elif nFatjet_MP==2:
      # all fat jet with dr(lep)<0.8
      if len(TightFatJet_drLa08_id)==0:
        HZ_1F2R=True
        HZ_1F2R_case=3
        Z_AK8Jet_id=TightFatJet_id[0]
        Z_AK8Jet_initid=TightFatJet_initid[0]
        Z_AK8Jet_pt=TightFatJet_pt[0]
        Z_AK8Jet_eta=TightFatJet_eta[0]
        Z_AK8Jet_phi=TightFatJet_phi[0]
        Z_AK8Jet_mass=TightFatJet_mass[0]
        Z_AK8Jet_drl1=TightFatJet_drl1[0]
        Z_AK8Jet_drl2=TightFatJet_drl2[0]
      elif len(TightFatJet_drLa08_id)==1:
        H_AK8_tempid=TightFatJet_drLa08_id[0]
        Z_AK8_tempid=TightFatJet_drSm08_id[0]
        if TightFatJet_bbvsQCD[H_AK8_tempid]>Hbb_threshold and TightFatJet_ccvsQCD[H_AK8_tempid]<Hcc_threshold:
          HZ_2F0R=True
          HZ_2F0R_case=1
          H_AK8Jet_id=TightFatJet_id[H_AK8_tempid]
          H_AK8Jet_initid=TightFatJet_initid[H_AK8_tempid]
          H_AK8Jet_pt=TightFatJet_pt[H_AK8_tempid]
          H_AK8Jet_eta=TightFatJet_eta[H_AK8_tempid]
          H_AK8Jet_phi=TightFatJet_phi[H_AK8_tempid]
          H_AK8Jet_mass=TightFatJet_mass[H_AK8_tempid]
          H_AK8Jet_PNmass=event.FatJet_particleNet_mass[H_AK8Jet_initid]
          H_AK8Jet_SDmass=event.FatJet_msoftdrop[H_AK8Jet_initid]
          H_AK8Jet_drl1=TightFatJet_drl1[H_AK8_tempid]
          H_AK8Jet_drl2=TightFatJet_drl2[H_AK8_tempid]

          Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
          Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
          Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
          Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
          Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
          Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
          Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
          Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]
        else:
          HZ_1F2R=True
          if H_AK8_tempid>Z_AK8_tempid:
          # Z_AK8_tempid with higher pT
            HZ_1F2R_case=4
            Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
            Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
            Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
            Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
            Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
            Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
            Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
            Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]
          else:
            if event.FatJetNoVlep_bbvsccqq[TightFatJet_id[H_AK8_tempid]]<0.5:
              HZ_1F2R_case=5
              Z_AK8Jet_id=TightFatJet_id[H_AK8_tempid]
              Z_AK8Jet_initid=TightFatJet_initid[H_AK8_tempid]
              Z_AK8Jet_pt=TightFatJet_pt[H_AK8_tempid]
              Z_AK8Jet_eta=TightFatJet_eta[H_AK8_tempid]
              Z_AK8Jet_phi=TightFatJet_phi[H_AK8_tempid]
              Z_AK8Jet_mass=TightFatJet_mass[H_AK8_tempid]
              Z_AK8Jet_drl1=TightFatJet_drl1[H_AK8_tempid]
              Z_AK8Jet_drl2=TightFatJet_drl2[H_AK8_tempid]
            else:
              HZ_1F2R_case=6
              Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
              Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
              Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
              Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
              Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
              Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
              Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
              Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]
      else:
        # both fat jet with dr(lep)>0.8
        H_AK8_tempid=argsort(TightFatJet_bbvsQCD)[-1]
        arr_tmp=[0,1]
        arr_tmp.remove(H_AK8_tempid)
        Z_AK8_tempid=arr_tmp[0]
        if TightFatJet_bbvsQCD[H_AK8_tempid]>Hbb_threshold and TightFatJet_ccvsQCD[H_AK8_tempid]<Hcc_threshold:
          HZ_2F0R=True
          HZ_2F0R_case=2
          H_AK8Jet_id=TightFatJet_id[H_AK8_tempid]
          H_AK8Jet_initid=TightFatJet_initid[H_AK8_tempid]
          H_AK8Jet_pt=TightFatJet_pt[H_AK8_tempid]
          H_AK8Jet_eta=TightFatJet_eta[H_AK8_tempid]
          H_AK8Jet_phi=TightFatJet_phi[H_AK8_tempid]
          H_AK8Jet_mass=TightFatJet_mass[H_AK8_tempid]
          H_AK8Jet_PNmass=event.FatJet_particleNet_mass[H_AK8Jet_id]
          H_AK8Jet_SDmass=event.FatJet_msoftdrop[H_AK8Jet_id]
          H_AK8Jet_drl1=TightFatJet_drl1[H_AK8_tempid]
          H_AK8Jet_drl2=TightFatJet_drl2[H_AK8_tempid]

          Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
          Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
          Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
          Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
          Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
          Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
          Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
          Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]

        elif TightFatJet_bbvsQCD[Z_AK8_tempid]>Hbb_threshold and TightFatJet_ccvsQCD[Z_AK8_tempid]<Hcc_threshold:
          HZ_2F0R=True
          HZ_2F0R_case=3
          H_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
          H_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
          H_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
          H_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
          H_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
          H_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
          H_AK8Jet_PNmass=event.FatJet_particleNet_mass[H_AK8Jet_id]
          H_AK8Jet_SDmass=event.FatJet_msoftdrop[H_AK8Jet_id]
          H_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
          H_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]

          Z_AK8Jet_id=TightFatJet_id[H_AK8_tempid]
          Z_AK8Jet_initid=TightFatJet_initid[H_AK8_tempid]
          Z_AK8Jet_pt=TightFatJet_pt[H_AK8_tempid]
          Z_AK8Jet_eta=TightFatJet_eta[H_AK8_tempid]
          Z_AK8Jet_phi=TightFatJet_phi[H_AK8_tempid]
          Z_AK8Jet_mass=TightFatJet_mass[H_AK8_tempid]
          Z_AK8Jet_drl1=TightFatJet_drl1[H_AK8_tempid]
          Z_AK8Jet_drl2=TightFatJet_drl2[H_AK8_tempid]

        else:
          HZ_1F2R=True
          if event.FatJetNoVlep_bbvsccqq[TightFatJet_id[0]]<0.5:
            HZ_1F2R_case=7
            Z_AK8Jet_id=TightFatJet_id[0]
            Z_AK8Jet_initid=TightFatJet_initid[0]
            Z_AK8Jet_pt=TightFatJet_pt[0]
            Z_AK8Jet_eta=TightFatJet_eta[0]
            Z_AK8Jet_phi=TightFatJet_phi[0]
            Z_AK8Jet_mass=TightFatJet_mass[0]
            Z_AK8Jet_drl1=TightFatJet_drl1[0]
            Z_AK8Jet_drl2=TightFatJet_drl2[0]
          elif event.FatJetNoVlep_bbvsccqq[TightFatJet_id[1]]<0.5:
            HZ_1F2R_case=8
            Z_AK8Jet_id=TightFatJet_id[1]
            Z_AK8Jet_initid=TightFatJet_initid[1]
            Z_AK8Jet_pt=TightFatJet_pt[1]
            Z_AK8Jet_eta=TightFatJet_eta[1]
            Z_AK8Jet_phi=TightFatJet_phi[1]
            Z_AK8Jet_mass=TightFatJet_mass[1]
            Z_AK8Jet_drl1=TightFatJet_drl1[1]
            Z_AK8Jet_drl2=TightFatJet_drl2[1]

    # at least three fat jets passing tight ID
    else:
      # all fat jet with dr(lep)<0.8
      if len(TightFatJet_drLa08_id)==0:
        HZ_1F2R=True
        HZ_1F2R_case=9
        Z_AK8Jet_id=TightFatJet_id[0]
        Z_AK8Jet_initid=TightFatJet_initid[0]
        Z_AK8Jet_pt=TightFatJet_pt[0]
        Z_AK8Jet_eta=TightFatJet_eta[0]
        Z_AK8Jet_phi=TightFatJet_phi[0]
        Z_AK8Jet_mass=TightFatJet_mass[0]
        Z_AK8Jet_drl1=TightFatJet_drl1[0]
        Z_AK8Jet_drl2=TightFatJet_drl2[0]

      elif len(TightFatJet_drLa08_id)==1:
        H_AK8_tempid=TightFatJet_drLa08_id[0]
        arr_tmp=[x for x in range(nFatjet_MP)]
        if TightFatJet_bbvsQCD[H_AK8_tempid]>Hbb_threshold and TightFatJet_ccvsQCD[H_AK8_tempid]<Hcc_threshold:
          HZ_2F0R=True
          HZ_2F0R_case=4
          H_AK8Jet_id=TightFatJet_id[H_AK8_tempid]
          H_AK8Jet_initid=TightFatJet_initid[H_AK8_tempid]
          H_AK8Jet_pt=TightFatJet_pt[H_AK8_tempid]
          H_AK8Jet_eta=TightFatJet_eta[H_AK8_tempid]
          H_AK8Jet_phi=TightFatJet_phi[H_AK8_tempid]
          H_AK8Jet_mass=TightFatJet_mass[H_AK8_tempid]
          H_AK8Jet_PNmass=event.FatJet_particleNet_mass[H_AK8Jet_id]
          H_AK8Jet_SDmass=event.FatJet_msoftdrop[H_AK8Jet_id]
          H_AK8Jet_drl1=TightFatJet_drl1[H_AK8_tempid]
          H_AK8Jet_drl2=TightFatJet_drl2[H_AK8_tempid]
        
          arr_tmp.remove(H_AK8_tempid)
          Z_AK8_tempid=arr_tmp[0]
          Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
          Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
          Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
          Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
          Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
          Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
          Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
          Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]

        else:
          HZ_1F2R=True
          arr_tmp.remove(H_AK8_tempid)
          Z_AK8_tempid=arr_tmp[0]
          if event.FatJetNoVlep_bbvsccqq[TightFatJet_id[H_AK8_tempid]]<0.5 and H_AK8_tempid<Z_AK8_tempid:
            HZ_1F2R_case=10
            Z_AK8Jet_id=TightFatJet_id[H_AK8_tempid]
            Z_AK8Jet_initid=TightFatJet_initid[H_AK8_tempid]
            Z_AK8Jet_pt=TightFatJet_pt[H_AK8_tempid]
            Z_AK8Jet_eta=TightFatJet_eta[H_AK8_tempid]
            Z_AK8Jet_phi=TightFatJet_phi[H_AK8_tempid]
            Z_AK8Jet_mass=TightFatJet_mass[H_AK8_tempid]
            Z_AK8Jet_drl1=TightFatJet_drl1[H_AK8_tempid]
            Z_AK8Jet_drl2=TightFatJet_drl2[H_AK8_tempid]
          else:
            HZ_1F2R_case=11
            Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
            Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
            Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
            Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
            Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
            Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
            Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
            Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]

      elif len(TightFatJet_drLa08_id)==2:
        H_AK8_tempid_1=TightFatJet_drLa08_id[0]
        H_AK8_tempid_2=TightFatJet_drLa08_id[1]
        if TightFatJet_bbvsQCD[H_AK8_tempid_1]>TightFatJet_bbvsQCD[H_AK8_tempid_2]:
          if TightFatJet_bbvsQCD[H_AK8_tempid_1]>Hbb_threshold and TightFatJet_ccvsQCD[H_AK8_tempid_1]<Hcc_threshold:
            arr_tmp=[x for x in range(nFatjet_MP)]
            HZ_2F0R=True
            HZ_2F0R_case=5
            H_AK8Jet_id=TightFatJet_id[H_AK8_tempid_1]
            H_AK8Jet_initid=TightFatJet_initid[H_AK8_tempid_1]
            H_AK8Jet_pt=TightFatJet_pt[H_AK8_tempid_1]
            H_AK8Jet_eta=TightFatJet_eta[H_AK8_tempid_1]
            H_AK8Jet_phi=TightFatJet_phi[H_AK8_tempid_1]
            H_AK8Jet_mass=TightFatJet_mass[H_AK8_tempid_1]
            H_AK8Jet_PNmass=event.FatJet_particleNet_mass[H_AK8Jet_initid]
            H_AK8Jet_SDmass=event.FatJet_msoftdrop[H_AK8Jet_initid]
            H_AK8Jet_drl1=TightFatJet_drl1[H_AK8_tempid_1]
            H_AK8Jet_drl2=TightFatJet_drl2[H_AK8_tempid_1]
        
            arr_tmp.remove(H_AK8_tempid_1)
            Z_AK8_tempid=arr_tmp[0]
            Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
            Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
            Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
            Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
            Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
            Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
            Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
            Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]
          elif TightFatJet_bbvsQCD[H_AK8_tempid_2]>Hbb_threshold and TightFatJet_ccvsQCD[H_AK8_tempid_2]<Hcc_threshold:
            arr_tmp=[x for x in range(nFatjet_MP)]
            HZ_2F0R=True
            HZ_2F0R_case=6
            H_AK8Jet_id=TightFatJet_id[H_AK8_tempid_2]
            H_AK8Jet_initid=TightFatJet_initid[H_AK8_tempid_2]
            H_AK8Jet_pt=TightFatJet_pt[H_AK8_tempid_2]
            H_AK8Jet_eta=TightFatJet_eta[H_AK8_tempid_2]
            H_AK8Jet_phi=TightFatJet_phi[H_AK8_tempid_2]
            H_AK8Jet_mass=TightFatJet_mass[H_AK8_tempid_2]
            H_AK8Jet_PNmass=event.FatJet_particleNet_mass[H_AK8Jet_initid]
            H_AK8Jet_SDmass=event.FatJet_msoftdrop[H_AK8Jet_initid]
            H_AK8Jet_drl1=TightFatJet_drl1[H_AK8_tempid_2]
            H_AK8Jet_drl2=TightFatJet_drl2[H_AK8_tempid_2]
        
            arr_tmp.remove(H_AK8_tempid_2)
            Z_AK8_tempid=arr_tmp[0]
            Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
            Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
            Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
            Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
            Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
            Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
            Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
            Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]
          else:
            HZ_1F2R=True
            if not (H_AK8_tempid_1==0 or H_AK8_tempid_2==0):
              HZ_1F2R_case=12
              Z_AK8Jet_id=TightFatJet_id[0]
              Z_AK8Jet_initid=TightFatJet_initid[0]
              Z_AK8Jet_pt=TightFatJet_pt[0]
              Z_AK8Jet_eta=TightFatJet_eta[0]
              Z_AK8Jet_phi=TightFatJet_phi[0]
              Z_AK8Jet_mass=TightFatJet_mass[0]
              Z_AK8Jet_drl1=TightFatJet_drl1[0]
              Z_AK8Jet_drl2=TightFatJet_drl2[0]
            elif H_AK8_tempid_1==0 and event.FatJetNoVlep_bbvsccqq[TightFatJet_id[H_AK8_tempid_1]]<0.5:
              HZ_1F2R_case=13
              Z_AK8Jet_id=TightFatJet_id[H_AK8_tempid_1]
              Z_AK8Jet_initid=TightFatJet_initid[H_AK8_tempid_1]
              Z_AK8Jet_pt=TightFatJet_pt[H_AK8_tempid_1]
              Z_AK8Jet_eta=TightFatJet_eta[H_AK8_tempid_1]
              Z_AK8Jet_phi=TightFatJet_phi[H_AK8_tempid_1]
              Z_AK8Jet_mass=TightFatJet_mass[H_AK8_tempid_1]
              Z_AK8Jet_drl1=TightFatJet_drl1[H_AK8_tempid_1]
              Z_AK8Jet_drl2=TightFatJet_drl2[H_AK8_tempid_1]
            elif H_AK8_tempid_2==0 and event.FatJetNoVlep_bbvsccqq[TightFatJet_id[H_AK8_tempid_2]]<0.5:
              HZ_1F2R_case=14
              Z_AK8Jet_id=TightFatJet_id[H_AK8_tempid_2]
              Z_AK8Jet_initid=TightFatJet_initid[H_AK8_tempid_2]
              Z_AK8Jet_pt=TightFatJet_pt[H_AK8_tempid_2]
              Z_AK8Jet_eta=TightFatJet_eta[H_AK8_tempid_2]
              Z_AK8Jet_phi=TightFatJet_phi[H_AK8_tempid_2]
              Z_AK8Jet_mass=TightFatJet_mass[H_AK8_tempid_2]
              Z_AK8Jet_drl1=TightFatJet_drl1[H_AK8_tempid_2]
              Z_AK8Jet_drl2=TightFatJet_drl2[H_AK8_tempid_2]
            else:
              HZ_1F2R_case=15
              arr_tmp=[x for x in range(nFatjet_MP)]
              arr_tmp.remove(H_AK8_tempid_1)
              arr_tmp.remove(H_AK8_tempid_2)
              Z_AK8_tempid=arr_tmp[0]
              Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
              Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
              Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
              Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
              Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
              Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
              Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
              Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]

        else:
          if TightFatJet_bbvsQCD[H_AK8_tempid_2]>Hbb_threshold and TightFatJet_ccvsQCD[H_AK8_tempid_2]<Hcc_threshold:
            arr_tmp=[x for x in range(nFatjet_MP)]
            HZ_2F0R=True
            HZ_2F0R_case=7
            H_AK8Jet_id=TightFatJet_id[H_AK8_tempid_2]
            H_AK8Jet_initid=TightFatJet_initid[H_AK8_tempid_2]
            H_AK8Jet_pt=TightFatJet_pt[H_AK8_tempid_2]
            H_AK8Jet_eta=TightFatJet_eta[H_AK8_tempid_2]
            H_AK8Jet_phi=TightFatJet_phi[H_AK8_tempid_2]
            H_AK8Jet_mass=TightFatJet_mass[H_AK8_tempid_2]
            H_AK8Jet_PNmass=event.FatJet_particleNet_mass[H_AK8Jet_initid]
            H_AK8Jet_SDmass=event.FatJet_msoftdrop[H_AK8Jet_initid]
            H_AK8Jet_drl1=TightFatJet_drl1[H_AK8_tempid_2]
            H_AK8Jet_drl2=TightFatJet_drl2[H_AK8_tempid_2]
        
            arr_tmp.remove(H_AK8_tempid_2)
            Z_AK8_tempid=arr_tmp[0]
            Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
            Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
            Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
            Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
            Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
            Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
            Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
            Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]
          elif TightFatJet_bbvsQCD[H_AK8_tempid_1]>Hbb_threshold and TightFatJet_ccvsQCD[H_AK8_tempid_1]<Hcc_threshold:
            arr_tmp=[x for x in range(nFatjet_MP)]
            HZ_2F0R=True
            HZ_2F0R_case=8
            H_AK8Jet_id=TightFatJet_id[H_AK8_tempid_1]
            H_AK8Jet_initid=TightFatJet_initid[H_AK8_tempid_1]
            H_AK8Jet_pt=TightFatJet_pt[H_AK8_tempid_1]
            H_AK8Jet_eta=TightFatJet_eta[H_AK8_tempid_1]
            H_AK8Jet_phi=TightFatJet_phi[H_AK8_tempid_1]
            H_AK8Jet_mass=TightFatJet_mass[H_AK8_tempid_1]
            H_AK8Jet_PNmass=event.FatJet_particleNet_mass[H_AK8Jet_initid]
            H_AK8Jet_SDmass=event.FatJet_msoftdrop[H_AK8Jet_initid]
            H_AK8Jet_drl1=TightFatJet_drl1[H_AK8_tempid_1]
            H_AK8Jet_drl2=TightFatJet_drl2[H_AK8_tempid_1]
        
            arr_tmp.remove(H_AK8_tempid_1)
            Z_AK8_tempid=arr_tmp[0]
            Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
            Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
            Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
            Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
            Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
            Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
            Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
            Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]
          else:
            HZ_1F2R=True
            if not (H_AK8_tempid_1==0 or H_AK8_tempid_2==0):
              HZ_1F2R_case=16
              Z_AK8Jet_id=TightFatJet_id[0]
              Z_AK8Jet_initid=TightFatJet_initid[0]
              Z_AK8Jet_pt=TightFatJet_pt[0]
              Z_AK8Jet_eta=TightFatJet_eta[0]
              Z_AK8Jet_phi=TightFatJet_phi[0]
              Z_AK8Jet_mass=TightFatJet_mass[0]
              Z_AK8Jet_drl1=TightFatJet_drl1[0]
              Z_AK8Jet_drl2=TightFatJet_drl2[0]
            elif H_AK8_tempid_1==0 and event.FatJetNoVlep_bbvsccqq[TightFatJet_id[H_AK8_tempid_1]]<0.5:
              HZ_1F2R_case=17
              Z_AK8Jet_id=TightFatJet_id[H_AK8_tempid_1]
              Z_AK8Jet_initid=TightFatJet_initid[H_AK8_tempid_1]
              Z_AK8Jet_pt=TightFatJet_pt[H_AK8_tempid_1]
              Z_AK8Jet_eta=TightFatJet_eta[H_AK8_tempid_1]
              Z_AK8Jet_phi=TightFatJet_phi[H_AK8_tempid_1]
              Z_AK8Jet_mass=TightFatJet_mass[H_AK8_tempid_1]
              Z_AK8Jet_drl1=TightFatJet_drl1[H_AK8_tempid_1]
              Z_AK8Jet_drl2=TightFatJet_drl2[H_AK8_tempid_1]
            elif H_AK8_tempid_2==0 and event.FatJetNoVlep_bbvsccqq[TightFatJet_id[H_AK8_tempid_2]]<0.5:
              HZ_1F2R_case=18
              Z_AK8Jet_id=TightFatJet_id[H_AK8_tempid_2]
              Z_AK8Jet_initid=TightFatJet_initid[H_AK8_tempid_2]
              Z_AK8Jet_pt=TightFatJet_pt[H_AK8_tempid_2]
              Z_AK8Jet_eta=TightFatJet_eta[H_AK8_tempid_2]
              Z_AK8Jet_phi=TightFatJet_phi[H_AK8_tempid_2]
              Z_AK8Jet_mass=TightFatJet_mass[H_AK8_tempid_2]
              Z_AK8Jet_drl1=TightFatJet_drl1[H_AK8_tempid_2]
              Z_AK8Jet_drl2=TightFatJet_drl2[H_AK8_tempid_2]
            else:
              HZ_1F2R_case=19
              arr_tmp=[x for x in range(nFatjet_MP)]
              arr_tmp.remove(H_AK8_tempid_1)
              arr_tmp.remove(H_AK8_tempid_2)
              Z_AK8_tempid=arr_tmp[0]
              Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
              Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
              Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
              Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
              Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
              Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
              Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
              Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]

      else:
        H_AK8_tempid_1=TightFatJet_drLa08_id[0]
        H_AK8_tempid_2=TightFatJet_drLa08_id[1]
        H_AK8_tempid_3=TightFatJet_drLa08_id[2]
        local_arr=[H_AK8_tempid_1, H_AK8_tempid_2, H_AK8_tempid_3]
        bbqcd_tmp=[TightFatJet_bbvsQCD[H_AK8_tempid_1],TightFatJet_bbvsQCD[H_AK8_tempid_2],TightFatJet_bbvsQCD[H_AK8_tempid_3]]
        # sort the three fatjet by bbvsqcd
        tmpid_bbqcd_local_1=local_arr[argsort(bbqcd_tmp)[-1]]
        tmpid_bbqcd_local_2=local_arr[argsort(bbqcd_tmp)[-2]]
        tmpid_bbqcd_local_3=local_arr[argsort(bbqcd_tmp)[-3]]
        if TightFatJet_bbvsQCD[tmpid_bbqcd_local_1]>Hbb_threshold and TightFatJet_ccvsQCD[tmpid_bbqcd_local_1]<Hcc_threshold:
          arr_tmp=[x for x in range(nFatjet_MP)]
          HZ_2F0R=True
          HZ_2F0R_case=9
          H_AK8Jet_id=TightFatJet_id[tmpid_bbqcd_local_1]
          H_AK8Jet_initid=TightFatJet_initid[tmpid_bbqcd_local_1]
          H_AK8Jet_pt=TightFatJet_pt[tmpid_bbqcd_local_1]
          H_AK8Jet_eta=TightFatJet_eta[tmpid_bbqcd_local_1]
          H_AK8Jet_phi=TightFatJet_phi[tmpid_bbqcd_local_1]
          H_AK8Jet_mass=TightFatJet_mass[tmpid_bbqcd_local_1]
          H_AK8Jet_PNmass=event.FatJet_particleNet_mass[H_AK8Jet_initid]
          H_AK8Jet_SDmass=event.FatJet_msoftdrop[H_AK8Jet_initid]
          H_AK8Jet_drl1=TightFatJet_drl1[tmpid_bbqcd_local_1]
          H_AK8Jet_drl2=TightFatJet_drl2[tmpid_bbqcd_local_1]
        
          arr_tmp.remove(tmpid_bbqcd_local_1)
          Z_AK8_tempid=arr_tmp[0]
          Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
          Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
          Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
          Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
          Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
          Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
          Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
          Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]
        elif TightFatJet_bbvsQCD[tmpid_bbqcd_local_2]>Hbb_threshold and TightFatJet_ccvsQCD[tmpid_bbqcd_local_2]<Hcc_threshold:
          arr_tmp=[x for x in range(nFatjet_MP)]
          HZ_2F0R=True
          HZ_2F0R_case=10
          H_AK8Jet_id=TightFatJet_id[tmpid_bbqcd_local_2]
          H_AK8Jet_initid=TightFatJet_initid[tmpid_bbqcd_local_2]
          H_AK8Jet_pt=TightFatJet_pt[tmpid_bbqcd_local_2]
          H_AK8Jet_eta=TightFatJet_eta[tmpid_bbqcd_local_2]
          H_AK8Jet_phi=TightFatJet_phi[tmpid_bbqcd_local_2]
          H_AK8Jet_mass=TightFatJet_mass[tmpid_bbqcd_local_2]
          H_AK8Jet_PNmass=event.FatJet_particleNet_mass[H_AK8Jet_initid]
          H_AK8Jet_SDmass=event.FatJet_msoftdrop[H_AK8Jet_initid]
          H_AK8Jet_drl1=TightFatJet_drl1[tmpid_bbqcd_local_2]
          H_AK8Jet_drl2=TightFatJet_drl2[tmpid_bbqcd_local_2]
        
          arr_tmp.remove(tmpid_bbqcd_local_2)
          Z_AK8_tempid=arr_tmp[0]
          Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
          Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
          Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
          Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
          Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
          Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
          Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
          Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]
        elif TightFatJet_bbvsQCD[tmpid_bbqcd_local_3]>Hbb_threshold and TightFatJet_ccvsQCD[tmpid_bbqcd_local_3]<Hcc_threshold:
          arr_tmp=[x for x in range(nFatjet_MP)]
          HZ_2F0R=True
          HZ_2F0R_case=11
          H_AK8Jet_id=TightFatJet_id[tmpid_bbqcd_local_3]
          H_AK8Jet_initid=TightFatJet_initid[tmpid_bbqcd_local_3]
          H_AK8Jet_pt=TightFatJet_pt[tmpid_bbqcd_local_3]
          H_AK8Jet_eta=TightFatJet_eta[tmpid_bbqcd_local_3]
          H_AK8Jet_phi=TightFatJet_phi[tmpid_bbqcd_local_3]
          H_AK8Jet_mass=TightFatJet_mass[tmpid_bbqcd_local_3]
          H_AK8Jet_PNmass=event.FatJet_particleNet_mass[H_AK8Jet_initid]
          H_AK8Jet_SDmass=event.FatJet_msoftdrop[H_AK8Jet_initid]
          H_AK8Jet_drl1=TightFatJet_drl1[tmpid_bbqcd_local_3]
          H_AK8Jet_drl2=TightFatJet_drl2[tmpid_bbqcd_local_3]
        
          arr_tmp.remove(tmpid_bbqcd_local_3)
          Z_AK8_tempid=arr_tmp[0]
          Z_AK8Jet_id=TightFatJet_id[Z_AK8_tempid]
          Z_AK8Jet_initid=TightFatJet_initid[Z_AK8_tempid]
          Z_AK8Jet_pt=TightFatJet_pt[Z_AK8_tempid]
          Z_AK8Jet_eta=TightFatJet_eta[Z_AK8_tempid]
          Z_AK8Jet_phi=TightFatJet_phi[Z_AK8_tempid]
          Z_AK8Jet_mass=TightFatJet_mass[Z_AK8_tempid]
          Z_AK8Jet_drl1=TightFatJet_drl1[Z_AK8_tempid]
          Z_AK8Jet_drl2=TightFatJet_drl2[Z_AK8_tempid]
        else:
          HZ_1F2R=True
          if event.FatJetNoVlep_bbvsccqq[TightFatJet_id[0]]<0.5:
            HZ_1F2R_case=20
            Z_AK8Jet_id=TightFatJet_id[0]
            Z_AK8Jet_initid=TightFatJet_initid[0]
            Z_AK8Jet_pt=TightFatJet_pt[0]
            Z_AK8Jet_eta=TightFatJet_eta[0]
            Z_AK8Jet_phi=TightFatJet_phi[0]
            Z_AK8Jet_mass=TightFatJet_mass[0]
            Z_AK8Jet_drl1=TightFatJet_drl1[0]
            Z_AK8Jet_drl2=TightFatJet_drl2[0]
          elif event.FatJetNoVlep_bbvsccqq[TightFatJet_id[1]]<0.5:
            HZ_1F2R_case=21
            Z_AK8Jet_id=TightFatJet_id[1]
            Z_AK8Jet_initid=TightFatJet_initid[1]
            Z_AK8Jet_pt=TightFatJet_pt[1]
            Z_AK8Jet_eta=TightFatJet_eta[1]
            Z_AK8Jet_phi=TightFatJet_phi[1]
            Z_AK8Jet_mass=TightFatJet_mass[1]
            Z_AK8Jet_drl1=TightFatJet_drl1[1]
            Z_AK8Jet_drl2=TightFatJet_drl2[1]
          elif event.FatJetNoVlep_bbvsccqq[TightFatJet_id[2]]<0.5:
            HZ_1F2R_case=22
            Z_AK8Jet_id=TightFatJet_id[2]
            Z_AK8Jet_initid=TightFatJet_initid[2]
            Z_AK8Jet_pt=TightFatJet_pt[2]
            Z_AK8Jet_eta=TightFatJet_eta[2]
            Z_AK8Jet_phi=TightFatJet_phi[2]
            Z_AK8Jet_mass=TightFatJet_mass[2]
            Z_AK8Jet_drl1=TightFatJet_drl1[2]
            Z_AK8Jet_drl2=TightFatJet_drl2[2]

    self.out.fillBranch("TightFatJet_id", TightFatJet_id)
    self.out.fillBranch("TightFatJet_bbvsQCD", TightFatJet_bbvsQCD)
    self.out.fillBranch("TightFatJet_ccvsQCD", TightFatJet_ccvsQCD)
    self.out.fillBranch("TightFatJet_qqvsQCD", TightFatJet_qqvsQCD)
    self.out.fillBranch("TightFatJet_drl1", TightFatJet_drl1)
    self.out.fillBranch("TightFatJet_drl2", TightFatJet_drl2)
    self.out.fillBranch("TightFatJet_pt", TightFatJet_pt)
    self.out.fillBranch("TightFatJet_eta", TightFatJet_eta)
    self.out.fillBranch("TightFatJet_phi", TightFatJet_phi)
    self.out.fillBranch("TightFatJet_mass", TightFatJet_mass)
    self.out.fillBranch("TightFatJet_dRinit", TightFatJet_dRinit)
    self.out.fillBranch("HZ_2F0R",HZ_2F0R)
    self.out.fillBranch("HZ_2F0R_case",HZ_2F0R_case)
    self.out.fillBranch("HZ_1F2R",HZ_1F2R)
    self.out.fillBranch("HZ_1F2R_case",HZ_1F2R_case)
    self.out.fillBranch("HZ_0F4R",HZ_0F4R)
    self.out.fillBranch("H_AK8Jet_id",H_AK8Jet_id)
    self.out.fillBranch("H_AK8Jet_initid",H_AK8Jet_initid)
    self.out.fillBranch("H_AK8Jet_pt",H_AK8Jet_pt)
    self.out.fillBranch("H_AK8Jet_eta",H_AK8Jet_eta)
    self.out.fillBranch("H_AK8Jet_phi",H_AK8Jet_phi)
    self.out.fillBranch("H_AK8Jet_mass",H_AK8Jet_mass)
    self.out.fillBranch("H_AK8Jet_PNmass",H_AK8Jet_PNmass)
    self.out.fillBranch("H_AK8Jet_SDmass",H_AK8Jet_SDmass)
    self.out.fillBranch("H_AK8Jet_drl1",H_AK8Jet_drl1)
    self.out.fillBranch("H_AK8Jet_drl2",H_AK8Jet_drl2)
    self.out.fillBranch("Z_AK8Jet_id",Z_AK8Jet_id)
    self.out.fillBranch("Z_AK8Jet_initid",Z_AK8Jet_initid)
    self.out.fillBranch("Z_AK8Jet_pt",Z_AK8Jet_pt)
    self.out.fillBranch("Z_AK8Jet_eta",Z_AK8Jet_eta)
    self.out.fillBranch("Z_AK8Jet_phi",Z_AK8Jet_phi)
    self.out.fillBranch("Z_AK8Jet_mass",Z_AK8Jet_mass)
    self.out.fillBranch("Z_AK8Jet_drl1",Z_AK8Jet_drl1)
    self.out.fillBranch("Z_AK8Jet_drl2",Z_AK8Jet_drl2)

    HYv4_tmp=TLorentzVector()
    Xv4_tmp=TLorentzVector()

    Y_2F0R_pt=-99
    Y_2F0R_eta=-99
    Y_2F0R_phi=-99
    Y_2F0R_mass=-99
    X_2F0R_pt=-99
    X_2F0R_eta=-99
    X_2F0R_phi=-99
    X_2F0R_mass=-99
    ZFatv4_tmp=TLorentzVector()
    HFatv4_tmp=TLorentzVector()
    if HZ_2F0R:
      ZFatv4_tmp.SetPtEtaPhiM(Z_AK8Jet_pt,Z_AK8Jet_eta,Z_AK8Jet_phi,Z_AK8Jet_mass)
      HFatv4_tmp.SetPtEtaPhiM(H_AK8Jet_pt,H_AK8Jet_eta,H_AK8Jet_phi,H_AK8Jet_mass)
      HYv4_tmp=l1v4_tmp+l2v4_tmp+ZFatv4_tmp
      Xv4_tmp=HYv4_tmp+HFatv4_tmp
      Y_2F0R_pt=HYv4_tmp.Pt()
      Y_2F0R_eta=HYv4_tmp.Eta()
      Y_2F0R_phi=HYv4_tmp.Phi()
      Y_2F0R_mass=HYv4_tmp.M()
      X_2F0R_pt=Xv4_tmp.Pt()
      X_2F0R_eta=Xv4_tmp.Eta()
      X_2F0R_phi=Xv4_tmp.Phi()
      X_2F0R_mass=Xv4_tmp.M()

    self.out.fillBranch("Y_2F0R_pt",Y_2F0R_pt)
    self.out.fillBranch("Y_2F0R_eta",Y_2F0R_eta)
    self.out.fillBranch("Y_2F0R_phi",Y_2F0R_phi)
    self.out.fillBranch("Y_2F0R_mass",Y_2F0R_mass)
    self.out.fillBranch("X_2F0R_pt",X_2F0R_pt)
    self.out.fillBranch("X_2F0R_eta",X_2F0R_eta)
    self.out.fillBranch("X_2F0R_phi",X_2F0R_phi)
    self.out.fillBranch("X_2F0R_mass",X_2F0R_mass)

    # https://btv-wiki.docs.cern.ch/ScaleFactors/#sf-campaigns
    # tight PF jets (ak4), *10>=2
    # medium B-tag WP
    # DeepFlavor=(nanoaod btagDeepFlavB) 
    # 2016apv loose: 0.0508, medium: 0.2598, tight: 0.6502
    # 2016 loose: 0.0480, medium: 0.2489, tight: 0.6377
    # 2017 loose: 0.0532, medium: 0.3040, tight: 0.7476
    # 2018 loose: 0.0490, medium: 0.2783, tight: 0.7100

    jets = Collection(event, 'JetNoVlep')
    jets_init = Collection(event, 'Jet')

    # don't remove overlap with selected Fatjet, used to define TTbar and DY CR
    TightAK4JetCRall_id = []
    TightAK4JetCRall_b_DeepJetL_id = []
    TightAK4JetCRall_b_DeepJetM_id = []

    # remove overlap with selected Fatjet, used for SR
    TightAK4Jet_id = []
    TightAK4Jet_initid = []
    TightAK4Jet_pt = []
    TightAK4Jet_eta = []
    TightAK4Jet_phi = []
    TightAK4Jet_mass = []
    TightAK4Jet_drl1 = []
    TightAK4Jet_drl2 = []
    TightAK4Jet_dRinit = []

    TightAK4Jet_nob_id = []
    TightAK4Jet_b_DeepJetM_id = []
    TightAK4Jet_b_DeepJetL_id = []
    nobjet_v4_passdrlep_id = []
    nobjet_v4_faildrlep_id = []

    jet_v4_all = []
    bjet_v4_all = []
    nobjet_v4_all = []
    nobjet_v4_passdrlep = []
    nobjet_v4_faildrlep = []

    # deepflavB, default is 2016apv WP
    medium_Bcut = 0.2598
    loose_Bcut = 0.0508
    if self.year=="2016":
      medium_Bcut = 0.2489
      loose_Bcut= 0.048
    if self.year=="2017":
      medium_Bcut = 0.3040
      loose_Bcut= 0.0532
    if self.year=="2018":
      medium_Bcut = 0.2783
      loose_Bcut= 0.0490

    jet_v4_temp=TLorentzVector()

    for ijet in range(0, event.nJetNoVlep):

      if abs(jets[ijet].eta)>2.4:continue

      jet_is_tau=0
      if nHad_tau>0:
        for ita in Had_tau_id:
          if ijet==event.Tau_jetIdx[ita]:jet_is_tau=1
      if jet_is_tau:continue

      # pass tight ID, i.e., not require dR(jet,lep) here
      if not (jets[ijet].jetId>1 and jets[ijet].pt_nom>30):continue 
      jet_v4_temp.SetPtEtaPhiM(jets[ijet].pt_nom,jets[ijet].eta,jets[ijet].phi,jets[ijet].mass_nom)

      # don't remove overlap with FatJet
      TightAK4JetCRall_id.append(ijet)
      if jets[ijet].drl1>0.4 and jets[ijet].drl2>0.4:
        if jets_init[jets[ijet].id].btagDeepFlavB > medium_Bcut:
          TightAK4JetCRall_b_DeepJetM_id.append(ijet)
          TightAK4JetCRall_b_DeepJetL_id.append(ijet)
        elif jets_init[jets[ijet].id].btagDeepFlavB > loose_Bcut:
          TightAK4JetCRall_b_DeepJetL_id.append(ijet)

      # remove overlap with FatJet
      if H_AK8Jet_id>-1:
        fatjetv4_tmp_=TLorentzVector()
        fatjetv4_tmp_.SetPtEtaPhiM(H_AK8Jet_pt,H_AK8Jet_eta,H_AK8Jet_phi,H_AK8Jet_mass)
        if jet_v4_temp.DeltaR(fatjetv4_tmp_)<0.8:continue
      if Z_AK8Jet_id>-1:
        fatjetv4_tmp_=TLorentzVector()
        fatjetv4_tmp_.SetPtEtaPhiM(Z_AK8Jet_pt,Z_AK8Jet_eta,Z_AK8Jet_phi,Z_AK8Jet_mass)
        if jet_v4_temp.DeltaR(fatjetv4_tmp_)<0.8:continue

      TightAK4Jet_id.append(ijet)
      TightAK4Jet_initid.append(jets[ijet].id)
      TightAK4Jet_pt.append(jets[ijet].pt_nom)
      TightAK4Jet_eta.append(jets[ijet].eta)
      TightAK4Jet_phi.append(jets[ijet].phi)
      TightAK4Jet_mass.append(jets[ijet].mass_nom)
      TightAK4Jet_drl1.append(jets[ijet].drl1)
      TightAK4Jet_drl2.append(jets[ijet].drl2)
      TightAK4Jet_dRinit.append(jets[ijet].dRinit)
      jet_v4_all.append(jet_v4_temp.Clone())
 
      # only use flavor tag infor for those jet with dR(lep)>0.4, i.e., the lepton subtraction will change the flavor tag infor and can't be used anymore after lepton subtraction
      if jets[ijet].drl1>0.4 and jets[ijet].drl2>0.4:
        if not jets_init[jets[ijet].id].btagDeepFlavB > loose_Bcut:
          nobjet_v4_passdrlep_id.append(ijet)
          nobjet_v4_passdrlep.append(jet_v4_temp.Clone())
        else:
          TightAK4Jet_b_DeepJetL_id.append(ijet)
          bjet_v4_all.append(jet_v4_temp.Clone())
          if jets_init[jets[ijet].id].btagDeepFlavB > medium_Bcut:
            TightAK4Jet_b_DeepJetM_id.append(ijet)
      else:
        nobjet_v4_faildrlep_id.append(ijet)
        nobjet_v4_faildrlep.append(jet_v4_temp.Clone())

    nTightAK4JetCRall_b_DeepJetL = len(TightAK4JetCRall_b_DeepJetL_id)
    nTightAK4JetCRall_b_DeepJetM = len(TightAK4JetCRall_b_DeepJetM_id)

    # TTbar region definition: mll>105
    TTCR_b1v4=TLorentzVector()
    TTCR_b2v4=TLorentzVector()
    TTCR_b1id_tmp=-1
    TTCR_b2id_tmp=-1

    TT_region_Lbb=False
    TT_Lb1_pt=-99
    TT_Lb1_eta=-99
    TT_Lb1_phi=-99
    TT_Lb1_mass=-99
    TT_Lb2_pt=-99
    TT_Lb2_eta=-99
    TT_Lb2_phi=-99
    TT_Lb2_mass=-99
    TT_Lmbb=-99
    TT_Lb1_drl1=-99
    TT_Lb1_drl2=-99
    TT_Lb2_drl1=-99
    TT_Lb2_drl2=-99

    TT_region_Mbb=False
    TT_Mb1_pt=-99
    TT_Mb1_eta=-99
    TT_Mb1_phi=-99
    TT_Mb1_mass=-99
    TT_Mb2_pt=-99
    TT_Mb2_eta=-99
    TT_Mb2_phi=-99
    TT_Mb2_mass=-99
    TT_Mmbb=-99
    TT_Mb1_drl1=-99
    TT_Mb1_drl2=-99
    TT_Mb2_drl1=-99
    TT_Mb2_drl2=-99

    if event.IsSS and (l1v4_tmp+l2v4_tmp).M()>105 and len(TightAK4JetCRall_b_DeepJetL_id)>1:
      TT_region_Lbb=True
    if TT_region_Lbb:
      TTCR_b1id_tmp=TightAK4JetCRall_b_DeepJetL_id[0]
      TTCR_b2id_tmp=TightAK4JetCRall_b_DeepJetL_id[1]
      TTCR_b1v4.SetPtEtaPhiM(jets[TTCR_b1id_tmp].pt_nom,jets[TTCR_b1id_tmp].eta,jets[TTCR_b1id_tmp].phi,jets[TTCR_b1id_tmp].mass_nom)
      TTCR_b2v4.SetPtEtaPhiM(jets[TTCR_b2id_tmp].pt_nom,jets[TTCR_b2id_tmp].eta,jets[TTCR_b2id_tmp].phi,jets[TTCR_b2id_tmp].mass_nom)
      TT_Lb1_pt=TTCR_b1v4.Pt()
      TT_Lb1_eta=TTCR_b1v4.Eta()
      TT_Lb1_phi=TTCR_b1v4.Phi()
      TT_Lb1_mass=TTCR_b1v4.M()
      TT_Lb2_pt=TTCR_b2v4.Pt()
      TT_Lb2_eta=TTCR_b2v4.Eta()
      TT_Lb2_phi=TTCR_b2v4.Phi()
      TT_Lb2_mass=TTCR_b2v4.M()
      TT_Lmbb=(TTCR_b1v4+TTCR_b2v4).M()
      TT_Lb1_drl1=TTCR_b1v4.DeltaR(l1v4_tmp)
      TT_Lb1_drl2=TTCR_b1v4.DeltaR(l2v4_tmp)
      TT_Lb2_drl1=TTCR_b2v4.DeltaR(l1v4_tmp)
      TT_Lb2_drl2=TTCR_b2v4.DeltaR(l2v4_tmp)
      
    if event.IsSS and (l1v4_tmp+l2v4_tmp).M()>105 and len(TightAK4JetCRall_b_DeepJetM_id)>1:
      TT_region_Mbb=True
    if TT_region_Mbb:
      TTCR_b1id_tmp=TightAK4JetCRall_b_DeepJetM_id[0]
      TTCR_b2id_tmp=TightAK4JetCRall_b_DeepJetM_id[1]
      TTCR_b1v4.SetPtEtaPhiM(jets[TTCR_b1id_tmp].pt_nom,jets[TTCR_b1id_tmp].eta,jets[TTCR_b1id_tmp].phi,jets[TTCR_b1id_tmp].mass_nom)
      TTCR_b2v4.SetPtEtaPhiM(jets[TTCR_b2id_tmp].pt_nom,jets[TTCR_b2id_tmp].eta,jets[TTCR_b2id_tmp].phi,jets[TTCR_b2id_tmp].mass_nom)
      TT_Mb1_pt=TTCR_b1v4.Pt()
      TT_Mb1_eta=TTCR_b1v4.Eta()
      TT_Mb1_phi=TTCR_b1v4.Phi()
      TT_Mb1_mass=TTCR_b1v4.M()
      TT_Mb2_pt=TTCR_b2v4.Pt()
      TT_Mb2_eta=TTCR_b2v4.Eta()
      TT_Mb2_phi=TTCR_b2v4.Phi()
      TT_Mb2_mass=TTCR_b2v4.M()
      TT_Mmbb=(TTCR_b1v4+TTCR_b2v4).M()
      TT_Mb1_drl1=TTCR_b1v4.DeltaR(l1v4_tmp)
      TT_Mb1_drl2=TTCR_b1v4.DeltaR(l2v4_tmp)
      TT_Mb2_drl1=TTCR_b2v4.DeltaR(l1v4_tmp)
      TT_Mb2_drl2=TTCR_b2v4.DeltaR(l2v4_tmp)
      
    # DY region definition: 75<mll<105, no b or (mbb>150 or mbb <100)
    DY2J_region=False
    DY2J_region_nob=False
    DYCR_j1v4=TLorentzVector()
    DYCR_j2v4=TLorentzVector()
    DYCR_j1id_tmp=-1
    DYCR_j2id_tmp=-1

    DY_j1_pt=-99
    DY_j1_eta=-99
    DY_j1_phi=-99
    DY_j1_mass=-99
    DY_j2_pt=-99
    DY_j2_eta=-99
    DY_j2_phi=-99
    DY_j2_mass=-99
    DY_mjj=-99
    DY_j1_drl1=-99
    DY_j1_drl2=-99
    DY_j2_drl1=-99
    DY_j2_drl2=-99

    if event.IsSS and (l1v4_tmp+l2v4_tmp).M()<105 and (l1v4_tmp+l2v4_tmp).M()>75 and len(TightAK4JetCRall_id)>1:
      DY2J_region=True
      DYCR_j1id_tmp=TightAK4JetCRall_id[0]
      DYCR_j2id_tmp=TightAK4JetCRall_id[1]
      DYCR_j1v4.SetPtEtaPhiM(jets[DYCR_j1id_tmp].pt_nom,jets[DYCR_j1id_tmp].eta,jets[DYCR_j1id_tmp].phi,jets[DYCR_j1id_tmp].mass_nom)
      DYCR_j2v4.SetPtEtaPhiM(jets[DYCR_j2id_tmp].pt_nom,jets[DYCR_j2id_tmp].eta,jets[DYCR_j2id_tmp].phi,jets[DYCR_j2id_tmp].mass_nom)
      DY_j1_pt=DYCR_j1v4.Pt()
      DY_j1_eta=DYCR_j1v4.Eta()
      DY_j1_phi=DYCR_j1v4.Phi()
      DY_j1_mass=DYCR_j1v4.M()
      DY_j2_pt=DYCR_j2v4.Pt()
      DY_j2_eta=DYCR_j2v4.Eta()
      DY_j2_phi=DYCR_j2v4.Phi()
      DY_j2_mass=DYCR_j2v4.M()
      DY_mjj=(DYCR_j1v4+DYCR_j2v4).M()
      DY_j1_drl1=DYCR_j1v4.DeltaR(l1v4_tmp)
      DY_j1_drl2=DYCR_j1v4.DeltaR(l2v4_tmp)
      DY_j2_drl1=DYCR_j2v4.DeltaR(l1v4_tmp)
      DY_j2_drl2=DYCR_j2v4.DeltaR(l2v4_tmp)

      if len(TightAK4JetCRall_b_DeepJetL_id)<1:
        DY2J_region_nob=True

    self.out.fillBranch("TT_region_Lbb",TT_region_Lbb)
    self.out.fillBranch("TT_Lb1_pt",TT_Lb1_pt)
    self.out.fillBranch("TT_Lb1_eta",TT_Lb1_eta)
    self.out.fillBranch("TT_Lb1_phi",TT_Lb1_phi)
    self.out.fillBranch("TT_Lb1_mass",TT_Lb1_mass)
    self.out.fillBranch("TT_Lb2_pt",TT_Lb2_pt)
    self.out.fillBranch("TT_Lb2_eta",TT_Lb2_eta)
    self.out.fillBranch("TT_Lb2_phi",TT_Lb2_phi)
    self.out.fillBranch("TT_Lb2_mass",TT_Lb2_mass)
    self.out.fillBranch("TT_Lmbb",TT_Lmbb)
    self.out.fillBranch("TT_Lb1_drl1",TT_Lb1_drl1)
    self.out.fillBranch("TT_Lb1_drl2",TT_Lb1_drl2)
    self.out.fillBranch("TT_Lb2_drl1",TT_Lb2_drl1)
    self.out.fillBranch("TT_Lb2_drl2",TT_Lb2_drl2)

    self.out.fillBranch("TT_region_Mbb",TT_region_Mbb)
    self.out.fillBranch("TT_Mb1_pt",TT_Mb1_pt)
    self.out.fillBranch("TT_Mb1_eta",TT_Mb1_eta)
    self.out.fillBranch("TT_Mb1_phi",TT_Mb1_phi)
    self.out.fillBranch("TT_Mb1_mass",TT_Mb1_mass)
    self.out.fillBranch("TT_Mb2_pt",TT_Mb2_pt)
    self.out.fillBranch("TT_Mb2_eta",TT_Mb2_eta)
    self.out.fillBranch("TT_Mb2_phi",TT_Mb2_phi)
    self.out.fillBranch("TT_Mb2_mass",TT_Mb2_mass)
    self.out.fillBranch("TT_Mmbb",TT_Mmbb)
    self.out.fillBranch("TT_Mb1_drl1",TT_Mb1_drl1)
    self.out.fillBranch("TT_Mb1_drl2",TT_Mb1_drl2)
    self.out.fillBranch("TT_Mb2_drl1",TT_Mb2_drl1)
    self.out.fillBranch("TT_Mb2_drl2",TT_Mb2_drl2)

    self.out.fillBranch("DY_region",DY_region)
    self.out.fillBranch("DY_region_nob",DY_region_nob)
    self.out.fillBranch("DY_j1_pt",DY_j1_pt)
    self.out.fillBranch("DY_j1_eta",DY_j1_eta)
    self.out.fillBranch("DY_j1_phi",DY_j1_phi)
    self.out.fillBranch("DY_j1_mass",DY_j1_mass)
    self.out.fillBranch("DY_j2_pt",DY_j2_pt)
    self.out.fillBranch("DY_j2_eta",DY_j2_eta)
    self.out.fillBranch("DY_j2_phi",DY_j2_phi)
    self.out.fillBranch("DY_j2_mass",DY_j2_mass)
    self.out.fillBranch("DY_mjj",DY_mjj)
    self.out.fillBranch("DY_j1_drl1",DY_j1_drl1)
    self.out.fillBranch("DY_j1_drl2",DY_j1_drl2)
    self.out.fillBranch("DY_j2_drl1",DY_j2_drl1)
    self.out.fillBranch("DY_j2_drl2",DY_j2_drl2)

    HT=0
    for ijet in TightAK4Jet_id:
      HT=HT + jets[ijet].pt_nom
    self.out.fillBranch("HT",HT)

    TightAK4Jet_nob_id = [x for x in TightAK4Jet_id if (x not in TightAK4Jet_b_DeepJetL_id)]
    nobjet_v4_all = [x for x in jet_v4_all if x not in bjet_v4_all]

    n_tight_jet=len(TightAK4Jet_id)

    n_bjet_DeepB_M = len(TightAK4Jet_b_DeepJetM_id)
    n_bjet_DeepB_L = len(TightAK4Jet_b_DeepJetL_id)
    n_tight_nob = len(TightAK4Jet_nob_id)
    self.out.fillBranch("n_bjet_DeepB_M",n_bjet_DeepB_M)
    self.out.fillBranch("n_bjet_DeepB_L",n_bjet_DeepB_L)
    self.out.fillBranch("n_tight_nob",n_tight_nob)

    Had_tau_id.extend(np.zeros(event.nTau-len(Had_tau_id),int)-1)
    self.out.fillBranch("Had_tau_id", Had_tau_id)
    
    # start the resolved decay of Z and H
    h_j1_pt=-99
    h_j1_eta=-99
    h_j1_phi=-99
    h_j1_mass=-99
    h_j1_id=-99
    h_j1_drl1=-99
    h_j1_drl2=-99
    h_j2_pt=-99
    h_j2_eta=-99
    h_j2_phi=-99
    h_j2_mass=-99
    h_j2_id=-99
    h_j2_drl1=-99
    h_j2_drl2=-99
    h_mjj=-99
    h_detajj=-99
    h_dRjj=-99
    h_dphijj=-99

    zhad_j1_pt=[]
    zhad_j1_eta=[]
    zhad_j1_phi=[]
    zhad_j1_mass=[]
    zhad_j1_id=[]
    zhad_j2_pt=[]
    zhad_j2_eta=[]
    zhad_j2_phi=[]
    zhad_j2_mass=[]
    zhad_j2_id=[]
    zhad_mjj=[]
    zhad_detajj=[]
    zhad_dRjj=[]
    zhad_dphijj=[]
    zhad_zlep_mass=[]

    z_j1_v4_temp = []
    z_j2_v4_temp = []

    # 4 resolved jets category
    if HZ_0F4R:
      # at least four good jets, and if only one b-tag jet, must be medium ID
      if n_tight_jet>3 and (n_bjet_DeepB_L>1 or n_bjet_DeepB_M>0):

        hbb_mass_threshold=99
        hbb_v4_temp=-1
  
        if n_bjet_DeepB_L==1:
          if len(nobjet_v4_passdrlep)>0:
            for ij in range(0,len(nobjet_v4_passdrlep)):
              if abs((bjet_v4_all[0]+nobjet_v4_passdrlep[ij]).M()-125.)<hbb_mass_threshold:
                hbb_mass_threshold=abs((bjet_v4_all[0]+nobjet_v4_passdrlep[ij]).M()-125.)
                hbb_v4_temp=ij
    
            if bjet_v4_all[0].Pt()>nobjet_v4_passdrlep[hbb_v4_temp].Pt():
              HZ_0F4R_case=0
              h_j1_pt=bjet_v4_all[0].Pt()
              h_j1_eta=bjet_v4_all[0].Eta()
              h_j1_phi=bjet_v4_all[0].Phi()
              h_j1_mass=bjet_v4_all[0].M()
              h_j1_id=TightAK4Jet_b_DeepJetL_id[0]
              h_j1_drl1=jets[h_j1_id].drl1
              h_j1_drl2=jets[h_j1_id].drl2
              h_j2_pt=nobjet_v4_passdrlep[hbb_v4_temp].Pt()
              h_j2_eta=nobjet_v4_passdrlep[hbb_v4_temp].Eta()
              h_j2_phi=nobjet_v4_all[hbb_v4_temp].Phi()
              h_j2_mass=nobjet_v4_all[hbb_v4_temp].M()
              h_j2_id=nobjet_v4_passdrlep_id[hbb_v4_temp]
              h_j2_drl1=jets[h_j2_id].drl1
              h_j2_drl2=jets[h_j2_id].drl2
            else:
              HZ_0F4R_case=1
              h_j1_pt=nobjet_v4_all[hbb_v4_temp].Pt()
              h_j1_eta=nobjet_v4_all[hbb_v4_temp].Eta()
              h_j1_phi=nobjet_v4_all[hbb_v4_temp].Phi()
              h_j1_mass=nobjet_v4_all[hbb_v4_temp].M()
              h_j1_id=nobjet_v4_passdrlep_id[hbb_v4_temp]
              h_j1_drl1=jets[h_j1_id].drl1
              h_j1_drl2=jets[h_j1_id].drl2
              h_j2_pt=bjet_v4_all[0].Pt()
              h_j2_eta=bjet_v4_all[0].Eta()
              h_j2_phi=bjet_v4_all[0].Phi()
              h_j2_mass=bjet_v4_all[0].M()
              h_j2_id=TightAK4Jet_b_DeepJetL_id[0]
              h_j2_drl1=jets[h_j2_id].drl1
              h_j2_drl2=jets[h_j2_id].drl2

            h_mjj=(bjet_v4_all[0]+nobjet_v4_passdrlep[hbb_v4_temp]).M()
            h_detajj=abs(h_j1_eta-h_j2_eta)
            h_dRjj=bjet_v4_all[0].DeltaR(nobjet_v4_passdrlep[hbb_v4_temp])
            h_dphijj=bjet_v4_all[0].DeltaPhi(nobjet_v4_passdrlep[hbb_v4_temp])
            
        if n_bjet_DeepB_L==2:
          HZ_0F4R_case=2
          h_j1_pt=bjet_v4_all[0].Pt()
          h_j1_eta=bjet_v4_all[0].Eta()
          h_j1_phi=bjet_v4_all[0].Phi()
          h_j1_mass=bjet_v4_all[0].M()
          h_j1_id=TightAK4Jet_b_DeepJetL_id[0]
          h_j1_drl1=jets[h_j1_id].drl1
          h_j1_drl2=jets[h_j1_id].drl2
          h_j2_pt=bjet_v4_all[1].Pt()
          h_j2_eta=bjet_v4_all[1].Eta()
          h_j2_phi=bjet_v4_all[1].Phi()
          h_j2_mass=bjet_v4_all[1].M()
          h_j2_id=TightAK4Jet_b_DeepJetL_id[1]
          h_j2_drl1=jets[h_j2_id].drl1
          h_j2_drl2=jets[h_j2_id].drl2
          h_mjj=(bjet_v4_all[0]+bjet_v4_all[1]).M()
          h_detajj=abs(h_j1_eta-h_j2_eta)
          h_dRjj=bjet_v4_all[0].DeltaR(bjet_v4_all[1])
          h_dphijj=bjet_v4_all[0].DeltaPhi(bjet_v4_all[1])
            
        hbb_id_item_temp=[]
        hbb_mass_temp=[]
        hbb_min_item_temp=-99
        if n_bjet_DeepB_L>2:
          for b1_id_temp,b1_v4_temp in enumerate(bjet_v4_all):
            for b2_id_temp,b2_v4_temp in enumerate(bjet_v4_all):
              if b1_id_temp<b2_id_temp:
                hbb_id_item_temp.append((b1_id_temp,b2_id_temp))
                hbb_mass_temp.append(abs((b1_v4_temp+b2_v4_temp).M()-125.))
          hbb_min_item_temp=hbb_mass_temp.index(min(hbb_mass_temp))
          h_j1_v4_temp=bjet_v4_all[hbb_id_item_temp[hbb_min_item_temp][0]]
          h_j2_v4_temp=bjet_v4_all[hbb_id_item_temp[hbb_min_item_temp][1]]
    
          HZ_0F4R_case=3
          h_j1_pt=h_j1_v4_temp.Pt()
          h_j1_eta=h_j1_v4_temp.Eta()
          h_j1_phi=h_j1_v4_temp.Phi()
          h_j1_mass=h_j1_v4_temp.M()
          h_j1_id=TightAK4Jet_b_DeepJetL_id[hbb_id_item_temp[hbb_min_item_temp][0]]
          h_j1_drl1=jets[h_j1_id].drl1
          h_j1_drl2=jets[h_j1_id].drl2
          h_j2_pt=h_j2_v4_temp.Pt()
          h_j2_eta=h_j2_v4_temp.Eta()
          h_j2_phi=h_j2_v4_temp.Phi()
          h_j2_mass=h_j2_v4_temp.M()
          h_j2_id=TightAK4Jet_b_DeepJetL_id[hbb_id_item_temp[hbb_min_item_temp][1]]
          h_j2_drl1=jets[h_j2_id].drl1
          h_j2_drl2=jets[h_j2_id].drl2
          h_mjj=(h_j1_v4_temp+h_j2_v4_temp).M()
          h_detajj=abs(h_j1_eta-h_j2_eta)
          h_dRjj=h_j1_v4_temp.DeltaR(h_j2_v4_temp)
          h_dphijj=h_j1_v4_temp.DeltaPhi(h_j2_v4_temp)

        zqq_id_item_temp=[]

        zqq_mass_temp_60=[]
        zqq_mass_temp_70=[]
        zqq_mass_temp_80=[]
        zqq_mass_temp_90=[]
        zqq_mass_temp_100=[]
        zqq_mass_temp_125=[]
        zqq_mass_temp_150=[]
        zqq_mass_temp_250=[]
        zqq_mass_temp_300=[]
        zqq_mass_temp_400=[]
        zqq_mass_temp_500=[]
        zqq_mass_temp_600=[]
        zqq_mass_temp_700=[]
        zqq_mass_temp_800=[]
        zqq_mass_temp_900=[]
        zqq_mass_temp_1000=[]
        zqq_mass_temp_1100=[]
        zqq_mass_temp_1200=[]
        zqq_mass_temp_1300=[]
        zqq_mass_temp_1400=[]
        zqq_mass_temp_1600=[]
        zqq_mass_temp_1800=[]
        zqq_mass_temp_2000=[]
        zqq_mass_temp_2200=[]
        zqq_mass_temp_2400=[]
        zqq_mass_temp_2500=[]
        zqq_mass_temp_2600=[]
        zqq_mass_temp_2800=[]

        zqq_min_item_temp_60=-99
        zqq_min_item_temp_70=-99
        zqq_min_item_temp_80=-99
        zqq_min_item_temp_90=-99
        zqq_min_item_temp_100=-99
        zqq_min_item_temp_125=-99
        zqq_min_item_temp_150=-99
        zqq_min_item_temp_250=-99
        zqq_min_item_temp_300=-99
        zqq_min_item_temp_400=-99
        zqq_min_item_temp_500=-99
        zqq_min_item_temp_600=-99
        zqq_min_item_temp_700=-99
        zqq_min_item_temp_800=-99
        zqq_min_item_temp_900=-99
        zqq_min_item_temp_1000=-99
        zqq_min_item_temp_1100=-99
        zqq_min_item_temp_1200=-99
        zqq_min_item_temp_1300=-99
        zqq_min_item_temp_1400=-99
        zqq_min_item_temp_1600=-99
        zqq_min_item_temp_1800=-99
        zqq_min_item_temp_2000=-99
        zqq_min_item_temp_2200=-99
        zqq_min_item_temp_2400=-99
        zqq_min_item_temp_2500=-99
        zqq_min_item_temp_2600=-99
        zqq_min_item_temp_2800=-99

        for q1_id_temp,q1_v4_temp in enumerate(jet_v4_all):
          if h_j1_id==TightAK4Jet_id[q1_id_temp] or h_j2_id==TightAK4Jet_id[q1_id_temp]:continue
          for q2_id_temp,q2_v4_temp in enumerate(jet_v4_all):
            if h_j1_id==TightAK4Jet_id[q1_id_temp] or h_j2_id==TightAK4Jet_id[q1_id_temp]:continue
            if q1_id_temp<q2_id_temp:
              zqq_id_item_temp.append((q1_id_temp,q2_id_temp))
              zqq_mass_temp_60.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-60.))
              zqq_mass_temp_70.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-70.))
              zqq_mass_temp_80.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-80.))
              zqq_mass_temp_90.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-90.))
              zqq_mass_temp_100.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-100.))
              zqq_mass_temp_125.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-125.))
              zqq_mass_temp_150.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-150.))
              zqq_mass_temp_250.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-250.))
              zqq_mass_temp_300.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-300.))
              zqq_mass_temp_400.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-400.))
              zqq_mass_temp_500.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-500.))
              zqq_mass_temp_600.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-600.))
              zqq_mass_temp_700.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-700.))
              zqq_mass_temp_800.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-800.))
              zqq_mass_temp_900.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-900.))
              zqq_mass_temp_1000.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-1000.))
              zqq_mass_temp_1100.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-1100.))
              zqq_mass_temp_1200.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-1200.))
              zqq_mass_temp_1300.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-1300.))
              zqq_mass_temp_1400.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-1400.))
              zqq_mass_temp_1600.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-1600.))
              zqq_mass_temp_1800.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-1800.))
              zqq_mass_temp_2000.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-2000.))
              zqq_mass_temp_2200.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-2200.))
              zqq_mass_temp_2400.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-2400.))
              zqq_mass_temp_2500.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-2500.))
              zqq_mass_temp_2600.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-2600.))
              zqq_mass_temp_2800.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-2800.))

        zqq_min_item_temp_60=zqq_mass_temp_60.index(min(zqq_mass_temp_60))
        zqq_min_item_temp_70=zqq_mass_temp_70.index(min(zqq_mass_temp_70))
        zqq_min_item_temp_80=zqq_mass_temp_80.index(min(zqq_mass_temp_80))
        zqq_min_item_temp_90=zqq_mass_temp_90.index(min(zqq_mass_temp_90))
        zqq_min_item_temp_100=zqq_mass_temp_100.index(min(zqq_mass_temp_100))
        zqq_min_item_temp_125=zqq_mass_temp_125.index(min(zqq_mass_temp_125))
        zqq_min_item_temp_150=zqq_mass_temp_150.index(min(zqq_mass_temp_150))
        zqq_min_item_temp_250=zqq_mass_temp_250.index(min(zqq_mass_temp_250))
        zqq_min_item_temp_300=zqq_mass_temp_300.index(min(zqq_mass_temp_300))
        zqq_min_item_temp_400=zqq_mass_temp_400.index(min(zqq_mass_temp_400))
        zqq_min_item_temp_500=zqq_mass_temp_500.index(min(zqq_mass_temp_500))
        zqq_min_item_temp_600=zqq_mass_temp_600.index(min(zqq_mass_temp_600))
        zqq_min_item_temp_700=zqq_mass_temp_700.index(min(zqq_mass_temp_700))
        zqq_min_item_temp_800=zqq_mass_temp_800.index(min(zqq_mass_temp_800))
        zqq_min_item_temp_900=zqq_mass_temp_900.index(min(zqq_mass_temp_900))
        zqq_min_item_temp_1000=zqq_mass_temp_1000.index(min(zqq_mass_temp_1000))
        zqq_min_item_temp_1100=zqq_mass_temp_1100.index(min(zqq_mass_temp_1100))
        zqq_min_item_temp_1200=zqq_mass_temp_1200.index(min(zqq_mass_temp_1200))
        zqq_min_item_temp_1300=zqq_mass_temp_1300.index(min(zqq_mass_temp_1300))
        zqq_min_item_temp_1400=zqq_mass_temp_1400.index(min(zqq_mass_temp_1400))
        zqq_min_item_temp_1600=zqq_mass_temp_1600.index(min(zqq_mass_temp_1600))
        zqq_min_item_temp_1800=zqq_mass_temp_1800.index(min(zqq_mass_temp_1800))
        zqq_min_item_temp_2000=zqq_mass_temp_2000.index(min(zqq_mass_temp_2000))
        zqq_min_item_temp_2200=zqq_mass_temp_2200.index(min(zqq_mass_temp_2200))
        zqq_min_item_temp_2400=zqq_mass_temp_2400.index(min(zqq_mass_temp_2400))
        zqq_min_item_temp_2500=zqq_mass_temp_2500.index(min(zqq_mass_temp_2500))
        zqq_min_item_temp_2600=zqq_mass_temp_2600.index(min(zqq_mass_temp_2600))
        zqq_min_item_temp_2800=zqq_mass_temp_2800.index(min(zqq_mass_temp_2800))


        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_60][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_70][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_80][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_90][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_100][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_125][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_150][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_250][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_300][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_400][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_500][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_600][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_700][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_800][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_900][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1000][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1100][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1200][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1300][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1400][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1600][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1800][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2000][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2200][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2400][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2500][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2600][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2800][0]].Clone())
                                                              
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_60][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_70][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_80][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_90][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_100][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_125][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_150][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_250][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_300][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_400][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_500][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_600][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_700][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_800][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_900][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1000][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1100][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1200][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1300][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1400][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1600][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1800][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2000][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2200][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2400][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2500][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2600][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2800][1]].Clone())

        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_60][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_70][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_80][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_90][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_100][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_125][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_150][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_250][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_300][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_400][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_500][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_600][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_700][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_800][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_900][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1000][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1100][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1200][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1300][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1400][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1600][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1800][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2000][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2200][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2400][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2500][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2600][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2800][0]])
        
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_60][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_70][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_80][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_90][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_100][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_125][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_150][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_250][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_300][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_400][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_500][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_600][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_700][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_800][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_900][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1000][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1100][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1200][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1300][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1400][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1600][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1800][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2000][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2200][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2400][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2500][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2600][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2800][1]])

        for iqjet in range(0,len(z_j2_v4_temp)):
          zhad_j1_pt.append(z_j1_v4_temp[iqjet].Pt())
          zhad_j1_eta.append(z_j1_v4_temp[iqjet].Eta())
          zhad_j1_phi.append(z_j1_v4_temp[iqjet].Phi())
          zhad_j1_mass.append(z_j1_v4_temp[iqjet].M())
          zhad_j2_pt.append(z_j2_v4_temp[iqjet].Pt())
          zhad_j2_eta.append(z_j2_v4_temp[iqjet].Eta())
          zhad_j2_phi.append(z_j2_v4_temp[iqjet].Phi())
          zhad_j2_mass.append(z_j2_v4_temp[iqjet].M())
          zhad_mjj.append((z_j1_v4_temp[iqjet]+z_j2_v4_temp[iqjet]).M())
          zhad_detajj.append(abs(zhad_j1_eta[iqjet]-zhad_j2_eta[iqjet]))
          zhad_dRjj.append(z_j1_v4_temp[iqjet].DeltaR(z_j2_v4_temp[iqjet]))
          zhad_dphijj.append(z_j1_v4_temp[iqjet].DeltaPhi(z_j2_v4_temp[iqjet]))
          zhad_zlep_mass.append((l1v4_tmp+l2v4_tmp+z_j1_v4_temp[iqjet]+z_j2_v4_temp[iqjet]).M())

    # 2 resolved jets category, at least two good jets
    if HZ_1F2R and n_tight_jet>1:

      # if only one b-tag jet, must be medium ID 
      if Z_AK8Jet_pt>0 and (n_bjet_DeepB_M>0 or n_bjet_DeepB_L>1):

        hbb_mass_threshold=99
        hbb_v4_temp=-1
  
        if n_bjet_DeepB_L==1:
          if len(nobjet_v4_passdrlep)>0:
            for ij in range(0,len(nobjet_v4_passdrlep)):
              if abs((bjet_v4_all[0]+nobjet_v4_passdrlep[ij]).M()-125.)<hbb_mass_threshold:
                hbb_mass_threshold=abs((bjet_v4_all[0]+nobjet_v4_passdrlep[ij]).M()-125.)
                hbb_v4_temp=ij
  
            if bjet_v4_all[0].Pt()>nobjet_v4_passdrlep[hbb_v4_temp].Pt():
              h_j1_pt=bjet_v4_all[0].Pt()
              h_j1_eta=bjet_v4_all[0].Eta()
              h_j1_phi=bjet_v4_all[0].Phi()
              h_j1_mass=bjet_v4_all[0].M()
              h_j1_id=TightAK4Jet_b_DeepJetL_id[0]
              h_j1_drl1=jets[h_j1_id].drl1
              h_j1_drl2=jets[h_j1_id].drl2
              h_j2_pt=nobjet_v4_passdrlep[hbb_v4_temp].Pt()
              h_j2_eta=nobjet_v4_passdrlep[hbb_v4_temp].Eta()
              h_j2_phi=nobjet_v4_passdrlep[hbb_v4_temp].Phi()
              h_j2_mass=nobjet_v4_passdrlep[hbb_v4_temp].M()
              h_j2_id=nobjet_v4_passdrlep_id[hbb_v4_temp]
              h_j2_drl1=jets[h_j2_id].drl1
              h_j2_drl2=jets[h_j2_id].drl2
            else:
              h_j1_pt=nobjet_v4_passdrlep[hbb_v4_temp].Pt()
              h_j1_eta=nobjet_v4_passdrlep[hbb_v4_temp].Eta()
              h_j1_phi=nobjet_v4_passdrlep[hbb_v4_temp].Phi()
              h_j1_mass=nobjet_v4_passdrlep[hbb_v4_temp].M()
              h_j1_id=nobjet_v4_passdrlep_id[hbb_v4_temp]
              h_j1_drl1=jets[h_j1_id].drl1
              h_j1_drl2=jets[h_j1_id].drl2
              h_j2_pt=bjet_v4_all[0].Pt()
              h_j2_eta=bjet_v4_all[0].Eta()
              h_j2_phi=bjet_v4_all[0].Phi()
              h_j2_mass=bjet_v4_all[0].M()
              h_j2_id=TightAK4Jet_b_DeepJetL_id[0]
              h_j2_drl1=jets[h_j2_id].drl1
              h_j2_drl2=jets[h_j2_id].drl2

            h_mjj=(bjet_v4_all[0]+nobjet_v4_passdrlep[hbb_v4_temp]).M()
            h_detajj=abs(h_j1_eta-h_j2_eta)
            h_dRjj=bjet_v4_all[0].DeltaR(nobjet_v4_passdrlep[hbb_v4_temp])
            h_dphijj=bjet_v4_all[0].DeltaPhi(nobjet_v4_passdrlep[hbb_v4_temp])
            
        if n_bjet_DeepB_L==2:
          h_j1_pt=bjet_v4_all[0].Pt()
          h_j1_eta=bjet_v4_all[0].Eta()
          h_j1_phi=bjet_v4_all[0].Phi()
          h_j1_mass=bjet_v4_all[0].M()
          h_j1_id=TightAK4Jet_b_DeepJetL_id[0]
          h_j1_drl1=jets[h_j1_id].drl1
          h_j1_drl2=jets[h_j1_id].drl2
          h_j2_pt=bjet_v4_all[1].Pt()
          h_j2_eta=bjet_v4_all[1].Eta()
          h_j2_phi=bjet_v4_all[1].Phi()
          h_j2_mass=bjet_v4_all[1].M()
          h_j2_id=TightAK4Jet_b_DeepJetL_id[1]
          h_j2_drl1=jets[h_j2_id].drl1
          h_j2_drl2=jets[h_j2_id].drl2
          h_mjj=(bjet_v4_all[0]+bjet_v4_all[1]).M()
          h_detajj=abs(h_j1_eta-h_j2_eta)
          h_dRjj=bjet_v4_all[0].DeltaR(bjet_v4_all[1])
          h_dphijj=bjet_v4_all[0].DeltaPhi(bjet_v4_all[1])
            
        hbb_id_item_temp=[]
        hbb_mass_temp=[]
        hbb_min_item_temp=-99
        if n_bjet_DeepB_L>2:
          for b1_id_temp,b1_v4_temp in enumerate(bjet_v4_all):
            for b2_id_temp,b2_v4_temp in enumerate(bjet_v4_all):
              if b1_id_temp<b2_id_temp:
                hbb_id_item_temp.append((b1_id_temp,b2_id_temp))
                hbb_mass_temp.append(abs((b1_v4_temp+b2_v4_temp).M()-125.))
          hbb_min_item_temp=hbb_mass_temp.index(min(hbb_mass_temp))
          h_j1_v4_temp=bjet_v4_all[hbb_id_item_temp[hbb_min_item_temp][0]]
          h_j2_v4_temp=bjet_v4_all[hbb_id_item_temp[hbb_min_item_temp][1]]
    
          h_j1_pt=h_j1_v4_temp.Pt()
          h_j1_eta=h_j1_v4_temp.Eta()
          h_j1_phi=h_j1_v4_temp.Phi()
          h_j1_mass=h_j1_v4_temp.M()
          h_j1_id=TightAK4Jet_b_DeepJetL_id[hbb_id_item_temp[hbb_min_item_temp][0]]
          h_j1_drl1=jets[h_j1_id].drl1
          h_j1_drl2=jets[h_j1_id].drl2
          h_j2_pt=h_j2_v4_temp.Pt()
          h_j2_eta=h_j2_v4_temp.Eta()
          h_j2_phi=h_j2_v4_temp.Phi()
          h_j2_mass=h_j2_v4_temp.M()
          h_j2_id=TightAK4Jet_b_DeepJetL_id[hbb_id_item_temp[hbb_min_item_temp][1]]
          h_j2_drl1=jets[h_j2_id].drl1
          h_j2_drl2=jets[h_j2_id].drl2
          h_mjj=(h_j1_v4_temp+h_j2_v4_temp).M()
          h_detajj=abs(h_j1_eta-h_j2_eta)
          h_dRjj=h_j1_v4_temp.DeltaR(h_j2_v4_temp)
          h_dphijj=h_j1_v4_temp.DeltaPhi(h_j2_v4_temp)

      if H_AK8Jet_pt>0:
        zqq_id_item_temp=[]

        zqq_mass_temp_60=[]
        zqq_mass_temp_70=[]
        zqq_mass_temp_80=[]
        zqq_mass_temp_90=[]
        zqq_mass_temp_100=[]
        zqq_mass_temp_125=[]
        zqq_mass_temp_150=[]
        zqq_mass_temp_250=[]
        zqq_mass_temp_300=[]
        zqq_mass_temp_400=[]
        zqq_mass_temp_500=[]
        zqq_mass_temp_600=[]
        zqq_mass_temp_700=[]
        zqq_mass_temp_800=[]
        zqq_mass_temp_900=[]
        zqq_mass_temp_1000=[]
        zqq_mass_temp_1100=[]
        zqq_mass_temp_1200=[]
        zqq_mass_temp_1300=[]
        zqq_mass_temp_1400=[]
        zqq_mass_temp_1600=[]
        zqq_mass_temp_1800=[]
        zqq_mass_temp_2000=[]
        zqq_mass_temp_2200=[]
        zqq_mass_temp_2400=[]
        zqq_mass_temp_2500=[]
        zqq_mass_temp_2600=[]
        zqq_mass_temp_2800=[]

        zqq_min_item_temp_60=-99
        zqq_min_item_temp_70=-99
        zqq_min_item_temp_80=-99
        zqq_min_item_temp_90=-99
        zqq_min_item_temp_100=-99
        zqq_min_item_temp_125=-99
        zqq_min_item_temp_150=-99
        zqq_min_item_temp_250=-99
        zqq_min_item_temp_300=-99
        zqq_min_item_temp_400=-99
        zqq_min_item_temp_500=-99
        zqq_min_item_temp_600=-99
        zqq_min_item_temp_700=-99
        zqq_min_item_temp_800=-99
        zqq_min_item_temp_900=-99
        zqq_min_item_temp_1000=-99
        zqq_min_item_temp_1100=-99
        zqq_min_item_temp_1200=-99
        zqq_min_item_temp_1300=-99
        zqq_min_item_temp_1400=-99
        zqq_min_item_temp_1600=-99
        zqq_min_item_temp_1800=-99
        zqq_min_item_temp_2000=-99
        zqq_min_item_temp_2200=-99
        zqq_min_item_temp_2400=-99
        zqq_min_item_temp_2500=-99
        zqq_min_item_temp_2600=-99
        zqq_min_item_temp_2800=-99

        for q1_id_temp,q1_v4_temp in enumerate(jet_v4_all):
          for q2_id_temp,q2_v4_temp in enumerate(jet_v4_all):
            if q1_id_temp<q2_id_temp:
              #print('reco zjpt:',q1_v4_temp.Pt(),q2_v4_temp.Pt())
              zqq_id_item_temp.append((q1_id_temp,q2_id_temp))
              #print('mass:',(l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M())
              zqq_mass_temp_60.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-60.))
              zqq_mass_temp_70.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-70.))
              zqq_mass_temp_80.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-80.))
              zqq_mass_temp_90.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-90.))
              zqq_mass_temp_100.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-100.))
              zqq_mass_temp_125.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-125.))
              zqq_mass_temp_150.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-150.))
              zqq_mass_temp_250.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-250.))
              zqq_mass_temp_300.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-300.))
              zqq_mass_temp_400.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-400.))
              zqq_mass_temp_500.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-500.))
              zqq_mass_temp_600.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-600.))
              zqq_mass_temp_700.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-700.))
              zqq_mass_temp_800.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-800.))
              zqq_mass_temp_900.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-900.))
              zqq_mass_temp_1000.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-1000.))
              zqq_mass_temp_1100.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-1100.))
              zqq_mass_temp_1200.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-1200.))
              zqq_mass_temp_1300.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-1300.))
              zqq_mass_temp_1400.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-1400.))
              zqq_mass_temp_1600.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-1600.))
              zqq_mass_temp_1800.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-1800.))
              zqq_mass_temp_2000.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-2000.))
              zqq_mass_temp_2200.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-2200.))
              zqq_mass_temp_2400.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-2400.))
              zqq_mass_temp_2500.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-2500.))
              zqq_mass_temp_2600.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-2600.))
              zqq_mass_temp_2800.append(abs((l1v4_tmp+l2v4_tmp+q1_v4_temp+q2_v4_temp).M()-2800.))

        zqq_min_item_temp_60=zqq_mass_temp_60.index(min(zqq_mass_temp_60))
        zqq_min_item_temp_70=zqq_mass_temp_70.index(min(zqq_mass_temp_70))
        zqq_min_item_temp_80=zqq_mass_temp_80.index(min(zqq_mass_temp_80))
        zqq_min_item_temp_90=zqq_mass_temp_90.index(min(zqq_mass_temp_90))
        zqq_min_item_temp_100=zqq_mass_temp_100.index(min(zqq_mass_temp_100))
        zqq_min_item_temp_125=zqq_mass_temp_125.index(min(zqq_mass_temp_125))
        zqq_min_item_temp_150=zqq_mass_temp_150.index(min(zqq_mass_temp_150))
        zqq_min_item_temp_250=zqq_mass_temp_250.index(min(zqq_mass_temp_250))
        zqq_min_item_temp_300=zqq_mass_temp_300.index(min(zqq_mass_temp_300))
        zqq_min_item_temp_400=zqq_mass_temp_400.index(min(zqq_mass_temp_400))
        zqq_min_item_temp_500=zqq_mass_temp_500.index(min(zqq_mass_temp_500))
        zqq_min_item_temp_600=zqq_mass_temp_600.index(min(zqq_mass_temp_600))
        zqq_min_item_temp_700=zqq_mass_temp_700.index(min(zqq_mass_temp_700))
        zqq_min_item_temp_800=zqq_mass_temp_800.index(min(zqq_mass_temp_800))
        zqq_min_item_temp_900=zqq_mass_temp_900.index(min(zqq_mass_temp_900))
        zqq_min_item_temp_1000=zqq_mass_temp_1000.index(min(zqq_mass_temp_1000))
        zqq_min_item_temp_1100=zqq_mass_temp_1100.index(min(zqq_mass_temp_1100))
        zqq_min_item_temp_1200=zqq_mass_temp_1200.index(min(zqq_mass_temp_1200))
        zqq_min_item_temp_1300=zqq_mass_temp_1300.index(min(zqq_mass_temp_1300))
        zqq_min_item_temp_1400=zqq_mass_temp_1400.index(min(zqq_mass_temp_1400))
        zqq_min_item_temp_1600=zqq_mass_temp_1600.index(min(zqq_mass_temp_1600))
        zqq_min_item_temp_1800=zqq_mass_temp_1800.index(min(zqq_mass_temp_1800))
        zqq_min_item_temp_2000=zqq_mass_temp_2000.index(min(zqq_mass_temp_2000))
        zqq_min_item_temp_2200=zqq_mass_temp_2200.index(min(zqq_mass_temp_2200))
        zqq_min_item_temp_2400=zqq_mass_temp_2400.index(min(zqq_mass_temp_2400))
        zqq_min_item_temp_2500=zqq_mass_temp_2500.index(min(zqq_mass_temp_2500))
        zqq_min_item_temp_2600=zqq_mass_temp_2600.index(min(zqq_mass_temp_2600))
        zqq_min_item_temp_2800=zqq_mass_temp_2800.index(min(zqq_mass_temp_2800))


        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_60][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_70][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_80][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_90][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_100][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_125][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_150][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_250][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_300][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_400][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_500][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_600][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_700][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_800][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_900][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1000][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1100][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1200][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1300][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1400][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1600][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1800][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2000][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2200][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2400][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2500][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2600][0]].Clone())
        z_j1_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2800][0]].Clone())

                                                              
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_60][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_70][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_80][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_90][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_100][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_125][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_150][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_250][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_300][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_400][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_500][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_600][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_700][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_800][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_900][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1000][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1100][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1200][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1300][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1400][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1600][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_1800][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2000][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2200][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2400][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2500][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2600][1]].Clone())
        z_j2_v4_temp.append(jet_v4_all[zqq_id_item_temp[zqq_min_item_temp_2800][1]].Clone())

        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_60][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_70][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_80][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_90][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_100][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_125][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_150][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_250][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_300][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_400][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_500][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_600][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_700][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_800][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_900][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1000][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1100][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1200][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1300][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1400][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1600][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1800][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2000][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2200][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2400][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2500][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2600][0]])
        zhad_j1_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2800][0]])
        
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_60][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_70][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_80][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_90][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_100][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_125][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_150][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_250][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_300][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_400][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_500][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_600][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_700][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_800][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_900][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1000][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1100][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1200][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1300][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1400][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1600][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_1800][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2000][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2200][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2400][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2500][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2600][1]])
        zhad_j2_id.append(TightAK4Jet_id[zqq_id_item_temp[zqq_min_item_temp_2800][1]])

        for iqjet in range(0,len(z_j2_v4_temp)):
          zhad_j1_pt.append(z_j1_v4_temp[iqjet].Pt())
          zhad_j1_eta.append(z_j1_v4_temp[iqjet].Eta())
          zhad_j1_phi.append(z_j1_v4_temp[iqjet].Phi())
          zhad_j1_mass.append(z_j1_v4_temp[iqjet].M())
          zhad_j2_pt.append(z_j2_v4_temp[iqjet].Pt())
          zhad_j2_eta.append(z_j2_v4_temp[iqjet].Eta())
          zhad_j2_phi.append(z_j2_v4_temp[iqjet].Phi())
          zhad_j2_mass.append(z_j2_v4_temp[iqjet].M())
          zhad_mjj.append((z_j1_v4_temp[iqjet]+z_j2_v4_temp[iqjet]).M())
          zhad_detajj.append(abs(zhad_j1_eta[iqjet]-zhad_j2_eta[iqjet]))
          zhad_dRjj.append(z_j1_v4_temp[iqjet].DeltaR(z_j2_v4_temp[iqjet]))
          zhad_dphijj.append(z_j1_v4_temp[iqjet].DeltaPhi(z_j2_v4_temp[iqjet]))
          zhad_zlep_mass.append((l1v4_tmp+l2v4_tmp+z_j1_v4_temp[iqjet]+z_j2_v4_temp[iqjet]).M())

    self.out.fillBranch("h_j1_pt",h_j1_pt)
    self.out.fillBranch("h_j1_eta",h_j1_eta)
    self.out.fillBranch("h_j1_phi",h_j1_phi)
    self.out.fillBranch("h_j1_mass",h_j1_mass)
    self.out.fillBranch("h_j1_id",h_j1_id)
    self.out.fillBranch("h_j1_drl1",h_j1_drl1)
    self.out.fillBranch("h_j1_drl2",h_j1_drl2)
    self.out.fillBranch("h_j2_pt",h_j2_pt)
    self.out.fillBranch("h_j2_eta",h_j2_eta)
    self.out.fillBranch("h_j2_phi",h_j2_phi)
    self.out.fillBranch("h_j2_mass",h_j2_mass)
    self.out.fillBranch("h_j2_id",h_j2_id)
    self.out.fillBranch("h_j2_drl1",h_j2_drl1)
    self.out.fillBranch("h_j2_drl2",h_j2_drl2)
    self.out.fillBranch("h_mjj",h_mjj)
    self.out.fillBranch("h_detajj",h_detajj)
    self.out.fillBranch("h_dRjj",h_dRjj)
    self.out.fillBranch("h_dphijj",h_dphijj)

    self.out.fillBranch("zhad_j1_pt",zhad_j1_pt)
    self.out.fillBranch("zhad_j1_eta",zhad_j1_eta)
    self.out.fillBranch("zhad_j1_phi",zhad_j1_phi)
    self.out.fillBranch("zhad_j1_mass",zhad_j1_mass)
    self.out.fillBranch("zhad_j1_id",zhad_j1_id)
    self.out.fillBranch("zhad_j2_pt",zhad_j2_pt)
    self.out.fillBranch("zhad_j2_eta",zhad_j2_eta)
    self.out.fillBranch("zhad_j2_phi",zhad_j2_phi)
    self.out.fillBranch("zhad_j2_mass",zhad_j2_mass)
    self.out.fillBranch("zhad_j2_id",zhad_j2_id)
    self.out.fillBranch("zhad_mjj",zhad_mjj)
    self.out.fillBranch("zhad_detajj",zhad_detajj)
    self.out.fillBranch("zhad_dRjj",zhad_dRjj)
    self.out.fillBranch("zhad_dphijj",zhad_dphijj)
    self.out.fillBranch("zhad_zlep_mass",zhad_zlep_mass)

    TightAK4JetCRall_b_DeepJetL_id.extend(np.zeros(len(TightAK4JetCRall_id)-len(TightAK4JetCRall_b_DeepJetL_id),int)-1)
    TightAK4JetCRall_b_DeepJetM_id.extend(np.zeros(len(TightAK4JetCRall_id)-len(TightAK4JetCRall_b_DeepJetM_id),int)-1)
    TightAK4Jet_nob_id.extend(np.zeros(n_tight_jet-n_tight_nob,int)-1)
    TightAK4Jet_b_DeepJetM_id.extend(np.zeros(n_tight_jet-n_bjet_DeepB_M,int)-1)
    TightAK4Jet_b_DeepJetL_id.extend(np.zeros(n_tight_jet-n_bjet_DeepB_L,int)-1)

    self.out.fillBranch("TightAK4JetCRall_id",TightAK4JetCRall_id)
    self.out.fillBranch("TightAK4JetCRall_b_DeepJetL_id",TightAK4JetCRall_b_DeepJetL_id)
    self.out.fillBranch("nTightAK4JetCRall_b_DeepJetL",nTightAK4JetCRall_b_DeepJetL)
    self.out.fillBranch("nTightAK4JetCRall_b_DeepJetM",nTightAK4JetCRall_b_DeepJetM)
    self.out.fillBranch("TightAK4JetCRall_b_DeepJetM_id",TightAK4JetCRall_b_DeepJetM_id)

    self.out.fillBranch("TightAK4Jet_id",TightAK4Jet_id)
    self.out.fillBranch("TightAK4Jet_initid",TightAK4Jet_initid)
    self.out.fillBranch("TightAK4Jet_nob_id",TightAK4Jet_nob_id)
    self.out.fillBranch("TightAK4Jet_b_DeepJetM_id",TightAK4Jet_b_DeepJetM_id)
    self.out.fillBranch("TightAK4Jet_b_DeepJetL_id",TightAK4Jet_b_DeepJetL_id)
    self.out.fillBranch("TightAK4Jet_drl1",TightAK4Jet_drl1)
    self.out.fillBranch("TightAK4Jet_drl2",TightAK4Jet_drl2)
    self.out.fillBranch("TightAK4Jet_pt",TightAK4Jet_pt)
    self.out.fillBranch("TightAK4Jet_eta",TightAK4Jet_eta)
    self.out.fillBranch("TightAK4Jet_phi",TightAK4Jet_phi)
    self.out.fillBranch("TightAK4Jet_mass",TightAK4Jet_mass)
    self.out.fillBranch("TightAK4Jet_dRinit",TightAK4Jet_dRinit)

    Y_ZF2R_pt=-99
    Y_ZF2R_eta=-99
    Y_ZF2R_phi=-99
    Y_ZF2R_mass=-99
    X_ZF2R_pt=-99
    X_ZF2R_eta=-99
    X_ZF2R_phi=-99
    X_ZF2R_mass=-99
    Hb1v4_tmp=TLorentzVector()
    Hb2v4_tmp=TLorentzVector()
    if HZ_1F2R and Z_AK8Jet_pt>0:
      ZFatv4_tmp.SetPtEtaPhiM(Z_AK8Jet_pt,Z_AK8Jet_eta,Z_AK8Jet_phi,Z_AK8Jet_mass)
      Hb1v4_tmp.SetPtEtaPhiM(h_j1_pt,h_j1_eta,h_j1_phi,h_j1_mass)
      Hb2v4_tmp.SetPtEtaPhiM(h_j2_pt,h_j2_eta,h_j2_phi,h_j2_mass)
      HYv4_tmp=ZFatv4_tmp+l1v4_tmp+l2v4_tmp
      Xv4_tmp=HYv4_tmp+Hb1v4_tmp+Hb2v4_tmp
      Y_ZF2R_pt=HYv4_tmp.Pt()
      Y_ZF2R_eta=HYv4_tmp.Eta()
      Y_ZF2R_phi=HYv4_tmp.Phi()
      Y_ZF2R_mass=HYv4_tmp.M()
      X_ZF2R_pt=Xv4_tmp.Pt()
      X_ZF2R_eta=Xv4_tmp.Eta()
      X_ZF2R_phi=Xv4_tmp.Phi()
      X_ZF2R_mass=Xv4_tmp.M()

    self.out.fillBranch("Y_ZF2R_pt",Y_ZF2R_pt)
    self.out.fillBranch("Y_ZF2R_eta",Y_ZF2R_eta)
    self.out.fillBranch("Y_ZF2R_phi",Y_ZF2R_phi)
    self.out.fillBranch("Y_ZF2R_mass",Y_ZF2R_mass)
    self.out.fillBranch("X_ZF2R_pt",X_ZF2R_pt)
    self.out.fillBranch("X_ZF2R_eta",X_ZF2R_eta)
    self.out.fillBranch("X_ZF2R_phi",X_ZF2R_phi)
    self.out.fillBranch("X_ZF2R_mass",X_ZF2R_mass)
      
    YZQQ_pt=[]
    YZQQ_eta=[]
    YZQQ_phi=[]
    YZQQ_mass=[]
    XZQQ_pt=[]
    XZQQ_eta=[]
    XZQQ_phi=[]
    XZQQ_mass=[]
    Zj1v4_tmp=TLorentzVector()
    Zj2v4_tmp=TLorentzVector()
    Hbbv4_tmp=TLorentzVector()
    if HZ_1F2R and H_AK8Jet_pt>0 and len(zhad_j1_pt)>0:
      Hbbv4_tmp.SetPtEtaPhiM(H_AK8Jet_pt,H_AK8Jet_eta,H_AK8Jet_phi,H_AK8Jet_mass)
      for ijj in range(0,len(zhad_j1_pt)):
        Zj1v4_tmp.SetPtEtaPhiM(zhad_j1_pt[ijj],zhad_j1_eta[ijj],zhad_j1_phi[ijj],zhad_j1_mass[ijj])
        Zj2v4_tmp.SetPtEtaPhiM(zhad_j2_pt[ijj],zhad_j2_eta[ijj],zhad_j2_phi[ijj],zhad_j2_mass[ijj])
        HYv4_tmp=Zj1v4_tmp+Zj2v4_tmp+l1v4_tmp+l2v4_tmp
        Xv4_tmp=HYv4_tmp+Hbbv4_tmp
        YZQQ_pt.append(HYv4_tmp.Pt())
        YZQQ_eta.append(HYv4_tmp.Eta())
        YZQQ_phi.append(HYv4_tmp.Phi())
        YZQQ_mass.append(HYv4_tmp.M())
        XZQQ_pt.append(Xv4_tmp.Pt())
        XZQQ_eta.append(Xv4_tmp.Eta())
        XZQQ_phi.append(Xv4_tmp.Phi())
        XZQQ_mass.append(Xv4_tmp.M())

    if HZ_0F4R and len(zhad_j1_pt)>0:
      Hb1v4_tmp.SetPtEtaPhiM(h_j1_pt,h_j1_eta,h_j1_phi,h_j1_mass)
      Hb2v4_tmp.SetPtEtaPhiM(h_j2_pt,h_j2_eta,h_j2_phi,h_j2_mass)
      for ijj in range(0,len(zhad_j1_pt)):
        Zj1v4_tmp.SetPtEtaPhiM(zhad_j1_pt[ijj],zhad_j1_eta[ijj],zhad_j1_phi[ijj],zhad_j1_mass[ijj])
        Zj2v4_tmp.SetPtEtaPhiM(zhad_j2_pt[ijj],zhad_j2_eta[ijj],zhad_j2_phi[ijj],zhad_j2_mass[ijj])
        HYv4_tmp=Zj1v4_tmp+Zj2v4_tmp+l1v4_tmp+l2v4_tmp
        Xv4_tmp=HYv4_tmp+Hb1v4_tmp+Hb2v4_tmp
        YZQQ_pt.append(HYv4_tmp.Pt())
        YZQQ_eta.append(HYv4_tmp.Eta())
        YZQQ_phi.append(HYv4_tmp.Phi())
        YZQQ_mass.append(HYv4_tmp.M())
        XZQQ_pt.append(Xv4_tmp.Pt())
        XZQQ_eta.append(Xv4_tmp.Eta())
        XZQQ_phi.append(Xv4_tmp.Phi())
        XZQQ_mass.append(Xv4_tmp.M())

    self.out.fillBranch("YZQQ_pt",YZQQ_pt)
    self.out.fillBranch("YZQQ_eta",YZQQ_eta)
    self.out.fillBranch("YZQQ_phi",YZQQ_phi)
    self.out.fillBranch("YZQQ_mass",YZQQ_mass)
    self.out.fillBranch("XZQQ_pt",XZQQ_pt)
    self.out.fillBranch("XZQQ_eta",XZQQ_eta)
    self.out.fillBranch("XZQQ_phi",XZQQ_phi)
    self.out.fillBranch("XZQQ_mass",XZQQ_mass)

    return True

HH2016apv = lambda: HHProducer("2016apv")
HH2016 = lambda: HHProducer("2016")
HH2017 = lambda: HHProducer("2017")
HH2018 = lambda: HHProducer("2018")

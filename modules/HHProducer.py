import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os
import numpy as np
from numpy import sign

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
    self.out.branch("n_tight_muon", "I")
    self.out.branch("n_loose_muon", "I")
    self.out.branch("n_tight_ele", "I")
    self.out.branch("n_loose_ele", "I")
    self.out.branch("n_tight_jet", "I")
    self.out.branch("n_bjet_DeepB_M", "I")
    self.out.branch("n_bjet_DeepB_L", "I")
    self.out.branch("n_tight_nob", "I")
    self.out.branch("HT", "F")
    self.out.branch("nHad_tau", "I")

    self.out.branch("l1_pt", "F")
    self.out.branch("l1_eta", "F")
    self.out.branch("l1_phi", "F")
    self.out.branch("l1_mass", "F")
    self.out.branch("l1_id", "I")
    self.out.branch("l1_pdgid", "I")
    self.out.branch("l2_pt", "F")
    self.out.branch("l2_eta", "F")
    self.out.branch("l2_phi", "F")
    self.out.branch("l2_mass", "F")
    self.out.branch("l2_id", "I")
    self.out.branch("l2_pdgid", "I")
    self.out.branch("mll", "F")
    self.out.branch("zlep_pt", "F")
    self.out.branch("zlep_eta", "F")
    self.out.branch("zlep_phi", "F")

    self.out.branch("tightJets_id","I",lenVar="nJet")
    self.out.branch("tightJets_nob_id","I",lenVar="nJet")
    self.out.branch("tightJets_b_DeepCSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightJets_b_DeepCSVloose_id","I",lenVar="nJet")
    self.out.branch("tightElectrons_id","I",lenVar="nElectron")
    self.out.branch("additional_vetoElectrons_id","I",lenVar="nElectron")
    self.out.branch("tightMuons_id","I",lenVar="nMuon")
    self.out.branch("additional_looseMuons_id","I",lenVar="nMuon")
    self.out.branch("Had_tau_id","I",lenVar="nTau")

    self.out.branch("met_user","F")
    self.out.branch("met_phi_user","F")

    self.out.branch("h_j1_pt", "F")
    self.out.branch("h_j1_eta", "F")
    self.out.branch("h_j1_phi", "F")
    self.out.branch("h_j1_mass", "F")
    self.out.branch("h_j1_id", "I")
    self.out.branch("h_j2_pt", "F")
    self.out.branch("h_j2_eta", "F")
    self.out.branch("h_j2_phi", "F")
    self.out.branch("h_j2_mass", "F")
    self.out.branch("h_j2_id", "I")
    self.out.branch("h_mjj", "F")
    self.out.branch("h_detajj", "F")
    self.out.branch("h_dRjj", "F")
    self.out.branch("h_dphijj", "F")
    self.out.branch("zhad_j1_pt", "F")
    self.out.branch("zhad_j1_eta", "F")
    self.out.branch("zhad_j1_phi", "F")
    self.out.branch("zhad_j1_mass", "F")
    self.out.branch("zhad_j1_id", "I")
    self.out.branch("zhad_j2_pt", "F")
    self.out.branch("zhad_j2_eta", "F")
    self.out.branch("zhad_j2_phi", "F")
    self.out.branch("zhad_j2_mass", "F")
    self.out.branch("zhad_j2_id", "I")
    self.out.branch("zhad_mjj", "F")
    self.out.branch("zhad_detajj", "F")
    self.out.branch("zhad_dRjj", "F")
    self.out.branch("zhad_dphijj", "F")
    self.out.branch("mass_zhad_zlep", "F")

    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))
    self.is_lhe = bool(inputTree.GetBranch("nLHEPart"))

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    
    # PV selection
    if (event.PV_npvsGood<1): return False

    # trigger selection
    # special action for 2017 single ele HLT, https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2017DataRecommendations#Single_Electron_Triggers
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

    # total number of ele+muon, currently require at least 1 leptons
    if ((event.nMuon + event.nElectron) < 2): return False
    if event.nJet<4: return False

    # Muon selection: tight cut-based ID + tight PF iso, or loose cut-based ID + loose PF iso, with pt > 20 GeV
    muons = Collection(event, 'Muon')
    muon_v4_temp=TLorentzVector()
    muon_v4_temp_raw=TLorentzVector()
    tightMuons = []
    tightMuons_raw = []
    tightMuons_pdgid = []
    tightMuons_id = []
    additional_looseMuons = []
    additional_looseMuons_pdgid = []
    additional_looseMuons_id = []

    jet_v4_temp=TLorentzVector()

    for imu in range(0, event.nMuon):
      # following cuts are preseletion for MVA muon ID
      if abs(muons[imu].eta)>2.4 or muons[imu].pfRelIso04_all>0.4: continue
      if not (muons[imu].looseId and event.Muon_corrected_pt[imu]>10):continue
      
      if muons[imu].mediumId:
        if muons[imu].pfRelIso04_all<0.25:
          muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
          muon_v4_temp_raw.SetPtEtaPhiM(muons[imu].pt, muons[imu].eta, muons[imu].phi, muons[imu].mass)
          tightMuons.append(muon_v4_temp.Clone())
          tightMuons_raw.append(muon_v4_temp_raw.Clone())
          tightMuons_pdgid.append(muons[imu].pdgId)
          tightMuons_id.append(imu)
        else:
          muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
          additional_looseMuons.append(muon_v4_temp.Clone())
          additional_looseMuons_pdgid.append(muons[imu].pdgId)
          additional_looseMuons_id.append(imu)
      else:
        muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
        additional_looseMuons.append(muon_v4_temp.Clone())
        additional_looseMuons_pdgid.append(muons[imu].pdgId)
        additional_looseMuons_id.append(imu)
        
    n_tight_muon = len(tightMuons)
    n_loose_muon = len(additional_looseMuons)

    self.out.fillBranch("n_tight_muon", n_tight_muon)
    self.out.fillBranch("n_loose_muon", n_loose_muon)
    tightMuons_id.extend(np.zeros(event.nMuon-len(tightMuons_id),int)-1)
    additional_looseMuons_id.extend(np.zeros(event.nMuon-len(additional_looseMuons_id),int)-1)
    self.out.fillBranch("tightMuons_id", tightMuons_id)
    self.out.fillBranch("additional_looseMuons_id", additional_looseMuons_id)

    # electron selection: tight (veto) cut-based ID + impact parameter cut, with pt > 15 GeV
    electrons = Collection(event, 'Electron')
    electron_v4_temp=TLorentzVector()
    electron_v4_temp_raw=TLorentzVector()
    tightElectrons = []
    tightElectrons_raw = []
    tightElectrons_pdgid = []
    tightElectrons_id = []
    additional_vetoElectrons = []
    additional_vetoElectrons_pdgid = []
    additional_vetoElectrons_id = []


    for iele in range(0, event.nElectron):
      if not (abs(electrons[iele].eta)<2.5 and electrons[iele].pt>10):continue
      if not electrons[iele].mvaFall17V2Iso_WPL: continue
      if electrons[iele].mvaFall17V2Iso_WP90 and electrons[iele].pt>15:
        electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
        electron_v4_temp_raw.SetPtEtaPhiM(electrons[iele].pt/electrons[iele].eCorr, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass/electrons[iele].eCorr)
        tightElectrons.append(electron_v4_temp.Clone())
        tightElectrons_raw.append(electron_v4_temp_raw.Clone())
        tightElectrons_pdgid.append(electrons[iele].pdgId)
        tightElectrons_id.append(iele)
      else:
        electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
        additional_vetoElectrons.append(electron_v4_temp.Clone())
        additional_vetoElectrons_pdgid.append(electrons[iele].pdgId)
        additional_vetoElectrons_id.append(iele)
        
    n_tight_ele = len(tightElectrons)
    n_loose_ele = len(additional_vetoElectrons)
    self.out.fillBranch("n_tight_ele", n_tight_ele)
    self.out.fillBranch("n_loose_ele", n_loose_ele)
    tightElectrons_id.extend(np.zeros(event.nElectron-len(tightElectrons_id),int)-1)
    additional_vetoElectrons_id.extend(np.zeros(event.nElectron-len(additional_vetoElectrons_id),int)-1)
    self.out.fillBranch("tightElectrons_id", tightElectrons_id)
    self.out.fillBranch("additional_vetoElectrons_id", additional_vetoElectrons_id)

    # tight leptons and additional loose leptons collection
    tightLeptons = tightMuons + tightElectrons
    tightLeptons.sort(key=lambda x: x.Pt(), reverse=True)
    tightLeptons_raw = tightMuons_raw + tightElectrons_raw
    tightLeptons_raw.sort(key=lambda x: x.Pt(), reverse=True)
    looseLeptons = additional_looseMuons + additional_vetoElectrons
    looseLeptons.sort(key=lambda x: x.Pt(), reverse=True)

    # only two tight leptons
    if not len(tightLeptons)==2:return False
    # two leptons should be same flavor
    if not (len(tightMuons)==2 or len(tightElectrons)==2):return False
    # two leptons should be opposite charge
    if len(tightMuons)==2:
      if not (tightMuons_pdgid[0]+tightMuons_pdgid[1])==0:return False
    if len(tightElectrons)==2:
      if not (tightElectrons_pdgid[0]+tightElectrons_pdgid[1])==0:return False
    # no additional loose lepton
    if len(looseLeptons)>0:return False

    l1_pt=-99
    l1_eta=-99
    l1_phi=-99
    l1_mass=-99
    l1_id=-1
    l1_pdgid=-99
    l2_pt=-99
    l2_eta=-99
    l2_phi=-99
    l2_mass=-99
    l2_id=-1
    l2_pdgid=-99
    mll=-99
    zlep_pt=-99
    zlep_eta=-99
    zlep_phi=-99
    
    if len(tightMuons)==2:
      l1_pt=tightMuons[0].Pt()
      l1_eta=tightMuons[0].Eta()
      l1_phi=tightMuons[0].Phi()
      l1_mass=tightMuons[0].M()
      l1_id=tightMuons_id[0]
      l1_pdgid=tightMuons_pdgid[0]
      l2_pt=tightMuons[1].Pt()
      l2_eta=tightMuons[1].Eta()
      l2_phi=tightMuons[1].Phi()
      l2_mass=tightMuons[1].M()
      l2_id=tightMuons_id[1]
      l2_pdgid=tightMuons_pdgid[1]
      mll=(tightMuons[0]+tightMuons[1]).M()
      zlep_pt=(tightMuons[0]+tightMuons[1]).Pt()
      zlep_eta=(tightMuons[0]+tightMuons[1]).Eta()
      zlep_phi=(tightMuons[0]+tightMuons[1]).Phi()
    if len(tightElectrons)==2:
      l1_pt=tightElectrons[0].Pt()
      l1_eta=tightElectrons[0].Eta()
      l1_phi=tightElectrons[0].Phi()
      l1_mass=tightElectrons[0].M()
      l1_id=tightElectrons_id[0]
      l1_pdgid=tightElectrons_pdgid[0]
      l2_pt=tightElectrons[1].Pt()
      l2_eta=tightElectrons[1].Eta()
      l2_phi=tightElectrons[1].Phi()
      l2_mass=tightElectrons[1].M()
      l2_id=tightElectrons_id[1]
      l2_pdgid=tightElectrons_pdgid[1]
      mll=(tightElectrons[0]+tightElectrons[1]).M()
      zlep_pt=(tightElectrons[0]+tightElectrons[1]).Pt()
      zlep_eta=(tightElectrons[0]+tightElectrons[1]).Eta()
      zlep_phi=(tightElectrons[0]+tightElectrons[1]).Phi()

    self.out.fillBranch("l1_pt",l1_pt)
    self.out.fillBranch("l1_eta",l1_eta)
    self.out.fillBranch("l1_phi",l1_phi)
    self.out.fillBranch("l1_mass",l1_mass)
    self.out.fillBranch("l1_id",l1_id)
    self.out.fillBranch("l1_pdgid",l1_pdgid)
    self.out.fillBranch("l2_pt",l2_pt)
    self.out.fillBranch("l2_eta",l2_eta)
    self.out.fillBranch("l2_phi",l2_phi)
    self.out.fillBranch("l2_mass",l2_mass)
    self.out.fillBranch("l2_id",l2_id)
    self.out.fillBranch("l2_pdgid",l2_pdgid)
    self.out.fillBranch("mll",mll)
    self.out.fillBranch("zlep_pt",zlep_pt)
    self.out.fillBranch("zlep_eta",zlep_eta)
    self.out.fillBranch("zlep_phi",zlep_phi)

    tau_v4_temp=TLorentzVector()
    taus = Collection(event, 'Tau')
    nHad_tau=0
    Had_tau_id=[]
    for itau in range(0, event.nTau):
      tau_v4_temp.SetPtEtaPhiM(taus[itau].pt, taus[itau].eta, taus[itau].phi, taus[itau].mass)
      pass_tau_lep_Dr=1
      if taus[itau].pt>20 and abs(taus[itau].eta)<2.3 and taus[itau].idDecayModeOldDMs and taus[itau].idDeepTau2017v2p1VSe>=4 and taus[itau].idDeepTau2017v2p1VSjet>=4 and taus[itau].idDeepTau2017v2p1VSmu>=1:
	for ilep in range(0,len(tightLeptons)):
          if tau_v4_temp.DeltaR(tightLeptons[ilep])<0.4:pass_tau_lep_Dr=0
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


    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
    # tightLepVeto PF jets (ak4), 2016 (111=7), 2017/2018 (110=6), medium B-tag WP
    # DeepCSV=(nanoaod btagDeepB) loose: 0.1355, medium: 0.4506, tight: 0.7738
    # DeepFlavor=(nanoaod btagDeepFlavB) loose: 0.0532, medium: 0.3040, tight: 0.7476

    # c-jet tag is based on two-D cuts, medium DeepJet WP:
    # CvsL=btagDeepFlavCvL: 0.085, CvsB=btagDeepFlavCvB: 0.34
    # c-tag not available in NANOAOD yet

    jets = Collection(event, 'Jet')

    tightJets_id = []
    tightJets_nob_id = []

    tightJets_b_DeepCSVmedium_id = []
    tightJets_b_DeepCSVloose_id = []

    # require DeltaR between Jets and tight leptons greater than 0.4
    jet_v4_all = []
    bjet_v4_all = []
    nobjet_v4_all = []

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

    for ijet in range(0, event.nJet):

      if abs(jets[ijet].eta)>2.4:continue

      jet_is_tau=0
      if nHad_tau>0:
        for ita in Had_tau_id:
          if ijet==event.Tau_jetIdx[ita]:jet_is_tau=1
      if jet_is_tau:continue

      pass_jet_lep_Dr=1
      jet_v4_temp.SetPtEtaPhiM(event.Jet_pt_nom[ijet],event.Jet_eta[ijet],event.Jet_phi[ijet],event.Jet_mass_nom[ijet])
      for ilep in range(0,len(tightLeptons)):
	if jet_v4_temp.DeltaR(tightLeptons[ilep])<0.4:pass_jet_lep_Dr=0

      if not (pass_jet_lep_Dr>0):continue
      if not (jets[ijet].jetId==6 and event.Jet_pt_nom[ijet]>30):continue 

      tightJets_id.append(ijet)
      jet_v4_all.append(jet_v4_temp.Clone())

      if jets[ijet].btagDeepFlavB > loose_Bcut:
        tightJets_b_DeepCSVloose_id.append(ijet)
        bjet_v4_all.append(jet_v4_temp.Clone())
      if jets[ijet].btagDeepFlavB > medium_Bcut:
        tightJets_b_DeepCSVmedium_id.append(ijet)


    HT=0
    for ijet in range(0,len(tightJets_id)):
      HT=HT+event.Jet_pt_nom[tightJets_id[ijet]]
    self.out.fillBranch("HT",HT)

    tightJets_nob_id = [x for x in tightJets_id if (x not in tightJets_b_DeepCSVloose_id)]
    nobjet_v4_all = [x for x in jet_v4_all if x not in bjet_v4_all]

    n_tight_jet = len(tightJets_id)
    n_bjet_DeepB_M = len(tightJets_b_DeepCSVmedium_id)
    n_bjet_DeepB_L = len(tightJets_b_DeepCSVloose_id)
    n_tight_nob = len(tightJets_nob_id)
    self.out.fillBranch("n_tight_jet",n_tight_jet)
    self.out.fillBranch("n_bjet_DeepB_M",n_bjet_DeepB_M)
    self.out.fillBranch("n_bjet_DeepB_L",n_bjet_DeepB_L)
    self.out.fillBranch("n_tight_nob",n_tight_nob)

    Had_tau_id.extend(np.zeros(event.nTau-len(Had_tau_id),int)-1)
    self.out.fillBranch("Had_tau_id", Had_tau_id)
    

    # at least four good jets
    if n_tight_jet<4:return False
    # if only one b-tag jet, must be medium ID 
    if n_bjet_DeepB_L<1:return False
    if n_bjet_DeepB_L==1 and n_bjet_DeepB_M<1:return False

    h_j1_pt=-99
    h_j1_eta=-99
    h_j1_phi=-99
    h_j1_mass=-99
    h_j1_id=-99
    h_j2_pt=-99
    h_j2_eta=-99
    h_j2_phi=-99
    h_j2_mass=-99
    h_j2_id=-99
    h_mjj=-99
    h_detajj=-99
    h_dRjj=-99
    h_dphijj=-99

    hbb_mass_threshold=99
    hbb_v4_temp=-1

    if n_bjet_DeepB_L==1:
      for ij in range(0,n_tight_nob):
        if abs((bjet_v4_all[0]+nobjet_v4_all[ij]).M()-125.)<hbb_mass_threshold:
          hbb_mass_threshold=abs((bjet_v4_all[0]+nobjet_v4_all[ij]).M()-125.)
          hbb_v4_temp=ij

      if bjet_v4_all[0].Pt()>nobjet_v4_all[hbb_v4_temp].Pt():
        h_j1_pt=bjet_v4_all[0].Pt()
        h_j1_eta=bjet_v4_all[0].Eta()
        h_j1_phi=bjet_v4_all[0].Phi()
        h_j1_mass=bjet_v4_all[0].M()
        h_j1_id=tightJets_b_DeepCSVloose_id[0]
        h_j2_pt=nobjet_v4_all[hbb_v4_temp].Pt()
        h_j2_eta=nobjet_v4_all[hbb_v4_temp].Eta()
        h_j2_phi=nobjet_v4_all[hbb_v4_temp].Phi()
        h_j2_mass=nobjet_v4_all[hbb_v4_temp].M()
        h_j2_id=tightJets_nob_id[hbb_v4_temp]
      else:
        h_j1_pt=nobjet_v4_all[hbb_v4_temp].Pt()
        h_j1_eta=nobjet_v4_all[hbb_v4_temp].Eta()
        h_j1_phi=nobjet_v4_all[hbb_v4_temp].Phi()
        h_j1_mass=nobjet_v4_all[hbb_v4_temp].M()
        h_j1_id=tightJets_nob_id[hbb_v4_temp]
        h_j2_pt=bjet_v4_all[0].Pt()
        h_j2_eta=bjet_v4_all[0].Eta()
        h_j2_phi=bjet_v4_all[0].Phi()
        h_j2_mass=bjet_v4_all[0].M()
        h_j2_id=tightJets_b_DeepCSVloose_id[0]
      h_mjj=(bjet_v4_all[0]+nobjet_v4_all[hbb_v4_temp]).M()
      h_detajj=abs(h_j1_eta-h_j2_eta)
      h_dRjj=bjet_v4_all[0].DeltaR(nobjet_v4_all[hbb_v4_temp])
      h_dphijj=bjet_v4_all[0].DeltaPhi(nobjet_v4_all[hbb_v4_temp])
        
    if n_bjet_DeepB_L==2:
      h_j1_pt=bjet_v4_all[0].Pt()
      h_j1_eta=bjet_v4_all[0].Eta()
      h_j1_phi=bjet_v4_all[0].Phi()
      h_j1_mass=bjet_v4_all[0].M()
      h_j1_id=tightJets_b_DeepCSVloose_id[0]
      h_j2_pt=bjet_v4_all[1].Pt()
      h_j2_eta=bjet_v4_all[1].Eta()
      h_j2_phi=bjet_v4_all[1].Phi()
      h_j2_mass=bjet_v4_all[1].M()
      h_j2_id=tightJets_b_DeepCSVloose_id[1]
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
      h_j1_id=tightJets_b_DeepCSVloose_id[hbb_id_item_temp[hbb_min_item_temp][0]]
      h_j2_pt=h_j2_v4_temp.Pt()
      h_j2_eta=h_j2_v4_temp.Eta()
      h_j2_phi=h_j2_v4_temp.Phi()
      h_j2_mass=h_j2_v4_temp.M()
      h_j2_id=tightJets_b_DeepCSVloose_id[hbb_id_item_temp[hbb_min_item_temp][1]]
      h_mjj=(h_j1_v4_temp+h_j2_v4_temp).M()
      h_detajj=abs(h_j1_eta-h_j2_eta)
      h_dRjj=h_j1_v4_temp.DeltaR(h_j2_v4_temp)
      h_dphijj=h_j1_v4_temp.DeltaPhi(h_j2_v4_temp)

    self.out.fillBranch("h_j1_pt",h_j1_pt)
    self.out.fillBranch("h_j1_eta",h_j1_eta)
    self.out.fillBranch("h_j1_phi",h_j1_phi)
    self.out.fillBranch("h_j1_mass",h_j1_mass)
    self.out.fillBranch("h_j1_id",h_j1_id)
    self.out.fillBranch("h_j2_pt",h_j2_pt)
    self.out.fillBranch("h_j2_eta",h_j2_eta)
    self.out.fillBranch("h_j2_phi",h_j2_phi)
    self.out.fillBranch("h_j2_mass",h_j2_mass)
    self.out.fillBranch("h_j2_id",h_j2_id)
    self.out.fillBranch("h_mjj",h_mjj)
    self.out.fillBranch("h_detajj",h_detajj)
    self.out.fillBranch("h_dRjj",h_dRjj)
    self.out.fillBranch("h_dphijj",h_dphijj)

    zhad_j1_pt=-99
    zhad_j1_eta=-99
    zhad_j1_phi=-99
    zhad_j1_mass=-99
    zhad_j1_id=-99
    zhad_j2_pt=-99
    zhad_j2_eta=-99
    zhad_j2_phi=-99
    zhad_j2_mass=-99
    zhad_j2_id=-99
    zhad_mjj=-99
    zhad_detajj=-99
    zhad_dRjj=-99
    zhad_dphijj=-99
    mass_zhad_zlep=-99

    zqq_id_item_temp=[]
    zqq_mass_temp=[]
    zqq_min_item_temp=-99
    for q1_id_temp,q1_v4_temp in enumerate(jet_v4_all):
      if h_j1_id==tightJets_id[q1_id_temp] or h_j2_id==tightJets_id[q1_id_temp]:continue
      for q2_id_temp,q2_v4_temp in enumerate(jet_v4_all):
        if h_j1_id==tightJets_id[q1_id_temp] or h_j2_id==tightJets_id[q1_id_temp]:continue
        if q1_id_temp<q2_id_temp:
          zqq_id_item_temp.append((q1_id_temp,q2_id_temp))
          zqq_mass_temp.append(abs((tightLeptons[0]+tightLeptons[1]+q1_v4_temp+q2_v4_temp).M()-125.))
    zqq_min_item_temp=zqq_mass_temp.index(min(zqq_mass_temp))
    z_j1_v4_temp=jet_v4_all[zqq_id_item_temp[zqq_min_item_temp][0]]
    z_j2_v4_temp=jet_v4_all[zqq_id_item_temp[zqq_min_item_temp][1]]

    zhad_j1_pt=z_j1_v4_temp.Pt()
    zhad_j1_eta=z_j1_v4_temp.Eta()
    zhad_j1_phi=z_j1_v4_temp.Phi()
    zhad_j1_mass=z_j1_v4_temp.M()
    zhad_j1_id=tightJets_id[zqq_id_item_temp[zqq_min_item_temp][0]]
    zhad_j2_pt=z_j2_v4_temp.Pt()
    zhad_j2_eta=z_j2_v4_temp.Eta()
    zhad_j2_phi=z_j2_v4_temp.Phi()
    zhad_j2_mass=z_j2_v4_temp.M()
    zhad_j2_id=tightJets_id[zqq_id_item_temp[zqq_min_item_temp][1]]
    zhad_mjj=(z_j1_v4_temp+z_j2_v4_temp).M()
    zhad_detajj=abs(zhad_j1_eta-zhad_j2_eta)
    zhad_dRjj=z_j1_v4_temp.DeltaR(z_j2_v4_temp)
    zhad_dphijj=z_j1_v4_temp.DeltaPhi(z_j2_v4_temp)
    mass_zhad_zlep=(tightLeptons[0]+tightLeptons[1]+z_j1_v4_temp+z_j2_v4_temp).M()

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
    self.out.fillBranch("mass_zhad_zlep",mass_zhad_zlep)

    tightJets_id.extend(np.zeros(event.nJet-len(tightJets_id),int)-1)
    tightJets_nob_id.extend(np.zeros(event.nJet-len(tightJets_nob_id),int)-1)
    tightJets_b_DeepCSVmedium_id.extend(np.zeros(event.nJet-len(tightJets_b_DeepCSVmedium_id),int)-1)
    tightJets_b_DeepCSVloose_id.extend(np.zeros(event.nJet-len(tightJets_b_DeepCSVloose_id),int)-1)

    self.out.fillBranch("tightJets_id",tightJets_id)
    self.out.fillBranch("tightJets_nob_id",tightJets_nob_id)
    self.out.fillBranch("tightJets_b_DeepCSVmedium_id",tightJets_b_DeepCSVmedium_id)
    self.out.fillBranch("tightJets_b_DeepCSVloose_id",tightJets_b_DeepCSVloose_id)

    if mll<20:return False

    return True

HH2016apv = lambda: HHProducer("2016apv")
HH2016 = lambda: HHProducer("2016")
HH2017 = lambda: HHProducer("2017")
HH2018 = lambda: HHProducer("2018")

import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class GetGenPart(Module):

  def __init__(self):
    pass

  def beginJob(self):
    pass
  def endJob(self):
    pass

  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))
    if self.is_mc:
      self.out.branch("GEN_y_pt", "F")
      self.out.branch("GEN_y_eta", "F")
      self.out.branch("GEN_y_phi", "F")
      self.out.branch("GEN_y_mass", "F")
      self.out.branch("GEN_hbb_pt", "F")
      self.out.branch("GEN_hbb_eta", "F")
      self.out.branch("GEN_hbb_phi", "F")
      self.out.branch("GEN_hbb_mass", "F")
      self.out.branch("GEN_hzz_pt", "F")
      self.out.branch("GEN_hzz_eta", "F")
      self.out.branch("GEN_hzz_phi", "F")
      self.out.branch("GEN_hzz_mass", "F")
      self.out.branch("GEN_b1_pt", "F")
      self.out.branch("GEN_b1_eta", "F")
      self.out.branch("GEN_b1_phi", "F")
      self.out.branch("GEN_b1_mass", "F")
      self.out.branch("GEN_b2_pt", "F")
      self.out.branch("GEN_b2_eta", "F")
      self.out.branch("GEN_b2_phi", "F")
      self.out.branch("GEN_b2_mass", "F")
      self.out.branch("GEN_dr_b1b2", "F")
      self.out.branch("GEN_zll_pt", "F")
      self.out.branch("GEN_zll_eta", "F")
      self.out.branch("GEN_zll_phi", "F")
      self.out.branch("GEN_zll_mass", "F")
      self.out.branch("GEN_zqq_pt", "F")
      self.out.branch("GEN_zqq_eta", "F")
      self.out.branch("GEN_zqq_phi", "F")
      self.out.branch("GEN_zqq_mass", "F")
      self.out.branch("GEN_dr_z1z2", "F")
      self.out.branch("GEN_zl1_pt", "F")
      self.out.branch("GEN_zl1_eta", "F")
      self.out.branch("GEN_zl1_phi", "F")
      self.out.branch("GEN_zl1_mass", "F")
      self.out.branch("GEN_zl2_pt", "F")
      self.out.branch("GEN_zl2_eta", "F")
      self.out.branch("GEN_zl2_phi", "F")
      self.out.branch("GEN_zl2_mass", "F")
      self.out.branch("GEN_dr_l1l2", "F")
      self.out.branch("GEN_zj1_pt", "F")
      self.out.branch("GEN_zj1_pdgid", "I")
      self.out.branch("GEN_zj1_eta", "F")
      self.out.branch("GEN_zj1_phi", "F")
      self.out.branch("GEN_zj1_mass", "F")
      self.out.branch("GEN_zj2_pt", "F")
      self.out.branch("GEN_zj2_pdgid", "I")
      self.out.branch("GEN_zj2_eta", "F")
      self.out.branch("GEN_zj2_phi", "F")
      self.out.branch("GEN_zj2_mass", "F")
      self.out.branch("GEN_dr_j1j2", "F")

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    if self.is_mc:
      GEN_y_pt = -99
      GEN_y_eta = -99
      GEN_y_phi = -99
      GEN_y_mass = -99
      GEN_hbb_pt = -99
      GEN_hbb_eta = -99
      GEN_hbb_phi = -99
      GEN_hbb_mass = -99
      GEN_hzz_pt = -99
      GEN_hzz_eta = -99
      GEN_hzz_phi = -99
      GEN_hzz_mass = -99
      GEN_b1_pt = -99
      GEN_b1_eta = -99
      GEN_b1_phi = -99
      GEN_b1_mass = -99
      GEN_b2_pt = -99
      GEN_b2_eta = -99
      GEN_b2_phi = -99
      GEN_b2_mass = -99
      GEN_dr_b1b2 = -99
      GEN_zll_pt = -99
      GEN_zll_eta = -99
      GEN_zll_phi = -99
      GEN_zll_mass = -99
      GEN_zqq_pt = -99
      GEN_zqq_eta = -99
      GEN_zqq_phi = -99
      GEN_zqq_mass = -99
      GEN_dr_z1z2 = -99
      GEN_zl1_pt = -99
      GEN_zl1_eta = -99
      GEN_zl1_phi = -99
      GEN_zl1_mass = -99
      GEN_zl2_pt = -99
      GEN_zl2_eta = -99
      GEN_zl2_phi = -99
      GEN_zl2_mass = -99
      GEN_dr_l1l2 = -99
      GEN_zj1_pt = -99
      GEN_zj1_pdgid = -99
      GEN_zj1_eta = -99
      GEN_zj1_phi = -99
      GEN_zj1_mass = -99
      GEN_zj2_pt = -99
      GEN_zj2_pdgid = -99
      GEN_zj2_eta = -99
      GEN_zj2_phi = -99
      GEN_zj2_mass = -99
      GEN_dr_j1j2 = -99

      tmp1_v4=ROOT.TLorentzVector()
      tmp2_v4=ROOT.TLorentzVector()

      GenPart = Collection(event, 'GenPart')
      GenPart_pdgId=[]
      GenPart_genPartIdxMother=[]
      for ig in range(event.nGenPart):
        GenPart_pdgId.append(event.GenPart_pdgId[ig])
        GenPart_genPartIdxMother.append(event.GenPart_genPartIdxMother[ig])

      # for signal X->HY, there is only one Y, for signal Y->HH, there are two Ys
      # obtain array id for Y
      temp_gen=np.array(list(range(0,event.nGenPart,1)))
      temp_y1=np.array(GenPart_pdgId)-35
      temp_y2=np.nonzero(temp_y1)[0]
      index_y=sorted(list(set(temp_gen)-set(temp_y2)))
      if len(index_y)==1:
        GEN_y_pt=GenPart[index_y[0]].pt
        GEN_y_eta=GenPart[index_y[0]].eta
        GEN_y_phi=GenPart[index_y[0]].phi
        GEN_y_mass=GenPart[index_y[0]].mass
      if len(index_y)==2:
        GEN_y_pt=GenPart[index_y[-1]].pt
        GEN_y_eta=GenPart[index_y[-1]].eta
        GEN_y_phi=GenPart[index_y[-1]].phi
        GEN_y_mass=GenPart[index_y[-1]].mass

      # obtain array id for Higgs bosons
      temp_h1=np.array(GenPart_pdgId)-25
      temp_h2=np.nonzero(temp_h1)[0]
      index_h1=-1
      index_h2=-1
      index_h=sorted(list(set(temp_gen)-set(temp_h2)))

      # obtain array id for b-jet
      temp_bs=np.abs(np.array(GenPart_pdgId))-5
      temp_bs_=np.nonzero(temp_bs)[0]
      bs_index_=sorted(list(set(temp_gen)-set(temp_bs_)))
      bs_index=[]
      # remove ISR
      for ib in bs_index_:
        if GenPart[ib].status==23: bs_index.append(ib)
      b1_index=bs_index[0]
      b2_index=bs_index[1]

      # obtain array id for two Z bosons
      temp_z1=np.array(GenPart_pdgId)-23
      temp_z2=np.nonzero(temp_z1)[0]
      index_z=sorted(list(set(temp_gen)-set(temp_z2)))

      # obtain array id for decay products for z1/z2
      temp_z1_decay1=np.array(GenPart_genPartIdxMother)-index_z[0]
      temp_z1_decay2=np.nonzero(temp_z1_decay1)[0]
      index_z1_decay=sorted(list(set(temp_gen)-set(temp_z1_decay2)))
      temp_z2_decay1=np.array(GenPart_genPartIdxMother)-index_z[1]
      temp_z2_decay2=np.nonzero(temp_z2_decay1)[0]
      index_z2_decay=sorted(list(set(temp_gen)-set(temp_z2_decay2)))

      if abs(GenPart_pdgId[index_z1_decay[0]])<10:
        GEN_zll_pt=GenPart[index_z[1]].pt
        GEN_zll_eta=GenPart[index_z[1]].eta
        GEN_zll_phi=GenPart[index_z[1]].phi
        GEN_zll_mass=GenPart[index_z[1]].mass
        GEN_zqq_pt=GenPart[index_z[0]].pt
        GEN_zqq_eta=GenPart[index_z[0]].eta
        GEN_zqq_phi=GenPart[index_z[0]].phi
        GEN_zqq_mass=GenPart[index_z[0]].mass
        tmp1_v4.SetPtEtaPhiM(GEN_zll_pt,GEN_zll_eta,GEN_zll_phi,GEN_zll_mass)
        tmp2_v4.SetPtEtaPhiM(GEN_zqq_pt,GEN_zqq_eta,GEN_zqq_phi,GEN_zqq_mass)
        GEN_dr_z1z2=tmp1_v4.DeltaR(tmp2_v4)
        GEN_zl1_pt=GenPart[index_z2_decay[0]].pt
        GEN_zl1_eta=GenPart[index_z2_decay[0]].eta
        GEN_zl1_phi=GenPart[index_z2_decay[0]].phi
        GEN_zl1_mass=GenPart[index_z2_decay[0]].mass
        GEN_zl2_pt=GenPart[index_z2_decay[1]].pt
        GEN_zl2_eta=GenPart[index_z2_decay[1]].eta
        GEN_zl2_phi=GenPart[index_z2_decay[1]].phi
        GEN_zl2_mass=GenPart[index_z2_decay[1]].mass
        tmp1_v4.SetPtEtaPhiM(GEN_zl1_pt,GEN_zl1_eta,GEN_zl1_phi,GEN_zl1_mass)
        tmp2_v4.SetPtEtaPhiM(GEN_zl2_pt,GEN_zl2_eta,GEN_zl2_phi,GEN_zl2_mass)
        GEN_dr_l1l2=tmp1_v4.DeltaR(tmp2_v4)
        GEN_zj1_pt=GenPart[index_z1_decay[0]].pt
        GEN_zj1_pdgid=GenPart[index_z1_decay[0]].pdgId
        GEN_zj1_eta=GenPart[index_z1_decay[0]].eta
        GEN_zj1_phi=GenPart[index_z1_decay[0]].phi
        GEN_zj1_mass=GenPart[index_z1_decay[0]].mass
        GEN_zj2_pt=GenPart[index_z1_decay[1]].pt
        GEN_zj2_pdgid=GenPart[index_z1_decay[1]].pdgId
        GEN_zj2_eta=GenPart[index_z1_decay[1]].eta
        GEN_zj2_phi=GenPart[index_z1_decay[1]].phi
        GEN_zj2_mass=GenPart[index_z1_decay[1]].mass
        tmp1_v4.SetPtEtaPhiM(GEN_zj1_pt,GEN_zj1_eta,GEN_zj1_phi,GEN_zj1_mass)
        tmp2_v4.SetPtEtaPhiM(GEN_zj2_pt,GEN_zj2_eta,GEN_zj2_phi,GEN_zj2_mass)
        GEN_dr_j1j2=tmp1_v4.DeltaR(tmp2_v4)
      else:
        GEN_zqq_pt=GenPart[index_z[1]].pt
        GEN_zqq_eta=GenPart[index_z[1]].eta
        GEN_zqq_phi=GenPart[index_z[1]].phi
        GEN_zqq_mass=GenPart[index_z[1]].mass
        GEN_zll_pt=GenPart[index_z[0]].pt
        GEN_zll_eta=GenPart[index_z[0]].eta
        GEN_zll_phi=GenPart[index_z[0]].phi
        GEN_zll_mass=GenPart[index_z[0]].mass
        tmp1_v4.SetPtEtaPhiM(GEN_zll_pt,GEN_zll_eta,GEN_zll_phi,GEN_zll_mass)
        tmp2_v4.SetPtEtaPhiM(GEN_zqq_pt,GEN_zqq_eta,GEN_zqq_phi,GEN_zqq_mass)
        GEN_dr_z1z2=tmp1_v4.DeltaR(tmp2_v4)
        GEN_zl1_pt=GenPart[index_z1_decay[0]].pt
        GEN_zl1_eta=GenPart[index_z1_decay[0]].eta
        GEN_zl1_phi=GenPart[index_z1_decay[0]].phi
        GEN_zl1_mass=GenPart[index_z1_decay[0]].mass
        GEN_zl2_pt=GenPart[index_z1_decay[1]].pt
        GEN_zl2_eta=GenPart[index_z1_decay[1]].eta
        GEN_zl2_phi=GenPart[index_z1_decay[1]].phi
        GEN_zl2_mass=GenPart[index_z1_decay[1]].mass
        tmp1_v4.SetPtEtaPhiM(GEN_zl1_pt,GEN_zl1_eta,GEN_zl1_phi,GEN_zl1_mass)
        tmp2_v4.SetPtEtaPhiM(GEN_zl2_pt,GEN_zl2_eta,GEN_zl2_phi,GEN_zl2_mass)
        GEN_dr_l1l2=tmp1_v4.DeltaR(tmp2_v4)
        GEN_zj1_pt=GenPart[index_z2_decay[0]].pt
        GEN_zj1_pdgid=GenPart[index_z2_decay[0]].pdgId
        GEN_zj1_eta=GenPart[index_z2_decay[0]].eta
        GEN_zj1_phi=GenPart[index_z2_decay[0]].phi
        GEN_zj1_mass=GenPart[index_z2_decay[0]].mass
        GEN_zj2_pt=GenPart[index_z2_decay[1]].pt
        GEN_zj2_pdgid=GenPart[index_z2_decay[1]].pdgId
        GEN_zj2_eta=GenPart[index_z2_decay[1]].eta
        GEN_zj2_phi=GenPart[index_z2_decay[1]].phi
        GEN_zj2_mass=GenPart[index_z2_decay[1]].mass
        tmp1_v4.SetPtEtaPhiM(GEN_zj1_pt,GEN_zj1_eta,GEN_zj1_phi,GEN_zj1_mass)
        tmp2_v4.SetPtEtaPhiM(GEN_zj2_pt,GEN_zj2_eta,GEN_zj2_phi,GEN_zj2_mass)
        GEN_dr_j1j2=tmp1_v4.DeltaR(tmp2_v4)

      # in X->YH case, only one Higgs boson
      if len(index_y)==1:
        temp_h_decay1=np.array(GenPart_genPartIdxMother)-index_h[0]
        temp_h_decay2=np.nonzero(temp_h_decay1)[0]
        index_h_decay=sorted(list(set(temp_gen)-set(temp_h_decay2)))

        # Higgs 4-m
        GEN_hbb_pt=GenPart[index_h[0]].pt
        GEN_hbb_eta=GenPart[index_h[0]].eta
        GEN_hbb_phi=GenPart[index_h[0]].phi
        GEN_hbb_mass=GenPart[index_h[0]].mass
        GEN_b1_pt=GenPart[index_h_decay[0]].pt
        GEN_b1_eta=GenPart[index_h_decay[0]].eta
        GEN_b1_phi=GenPart[index_h_decay[0]].phi
        GEN_b1_mass=GenPart[index_h_decay[0]].mass
        GEN_b2_pt=GenPart[index_h_decay[1]].pt
        GEN_b2_eta=GenPart[index_h_decay[1]].eta
        GEN_b2_phi=GenPart[index_h_decay[1]].phi
        GEN_b2_mass=GenPart[index_h_decay[1]].mass
        tmp1_v4.SetPtEtaPhiM(GEN_b1_pt,GEN_b1_eta,GEN_b1_phi,GEN_b1_mass)
        tmp2_v4.SetPtEtaPhiM(GEN_b2_pt,GEN_b2_eta,GEN_b2_phi,GEN_b2_mass)
        GEN_dr_b1b2=tmp1_v4.DeltaR(tmp2_v4)

      # in X->HH case, two Higgs bosons
      if len(index_y)==2:
        # Higgs 4-m
        GEN_hbb_pt=GenPart[GenPart_genPartIdxMother[b1_index]].pt
        GEN_hbb_eta=GenPart[GenPart_genPartIdxMother[b1_index]].eta
        GEN_hbb_phi=GenPart[GenPart_genPartIdxMother[b1_index]].phi
        GEN_hbb_mass=GenPart[GenPart_genPartIdxMother[b1_index]].mass
        GEN_hzz_pt=GenPart[GenPart_genPartIdxMother[index_z[0]]].pt
        GEN_hzz_eta=GenPart[GenPart_genPartIdxMother[index_z[0]]].eta
        GEN_hzz_phi=GenPart[GenPart_genPartIdxMother[index_z[0]]].phi
        GEN_hzz_mass=GenPart[GenPart_genPartIdxMother[index_z[0]]].mass
        GEN_b1_pt=GenPart[b1_index].pt
        GEN_b1_eta=GenPart[b1_index].eta
        GEN_b1_phi=GenPart[b1_index].phi
        GEN_b1_mass=GenPart[b1_index].mass
        GEN_b2_pt=GenPart[b2_index].pt
        GEN_b2_eta=GenPart[b2_index].eta
        GEN_b2_phi=GenPart[b2_index].phi
        GEN_b2_mass=GenPart[b2_index].mass
        tmp1_v4.SetPtEtaPhiM(GEN_b1_pt,GEN_b1_eta,GEN_b1_phi,GEN_b1_mass)
        tmp2_v4.SetPtEtaPhiM(GEN_b2_pt,GEN_b2_eta,GEN_b2_phi,GEN_b2_mass)
        GEN_dr_b1b2=tmp1_v4.DeltaR(tmp2_v4)

      self.out.fillBranch("GEN_y_pt", GEN_y_pt)
      self.out.fillBranch("GEN_y_eta", GEN_y_eta)
      self.out.fillBranch("GEN_y_phi", GEN_y_phi)
      self.out.fillBranch("GEN_y_mass", GEN_y_mass)
      self.out.fillBranch("GEN_hbb_pt", GEN_hbb_pt)
      self.out.fillBranch("GEN_hbb_eta", GEN_hbb_eta)
      self.out.fillBranch("GEN_hbb_phi", GEN_hbb_phi)
      self.out.fillBranch("GEN_hbb_mass", GEN_hbb_mass)
      self.out.fillBranch("GEN_hzz_pt", GEN_hzz_pt)
      self.out.fillBranch("GEN_hzz_eta", GEN_hzz_eta)
      self.out.fillBranch("GEN_hzz_phi", GEN_hzz_phi)
      self.out.fillBranch("GEN_hzz_mass", GEN_hzz_mass)
      self.out.fillBranch("GEN_b1_pt", GEN_b1_pt)
      self.out.fillBranch("GEN_b1_eta", GEN_b1_eta)
      self.out.fillBranch("GEN_b1_phi", GEN_b1_phi)
      self.out.fillBranch("GEN_b1_mass", GEN_b1_mass)
      self.out.fillBranch("GEN_b2_pt", GEN_b2_pt)
      self.out.fillBranch("GEN_b2_eta", GEN_b2_eta)
      self.out.fillBranch("GEN_b2_phi", GEN_b2_phi)
      self.out.fillBranch("GEN_b2_mass", GEN_b2_mass)
      self.out.fillBranch("GEN_dr_b1b2", GEN_dr_b1b2)
      self.out.fillBranch("GEN_zll_pt", GEN_zll_pt)
      self.out.fillBranch("GEN_zll_eta", GEN_zll_eta)
      self.out.fillBranch("GEN_zll_phi", GEN_zll_phi)
      self.out.fillBranch("GEN_zll_mass", GEN_zll_mass)
      self.out.fillBranch("GEN_zqq_pt", GEN_zqq_pt)
      self.out.fillBranch("GEN_zqq_eta", GEN_zqq_eta)
      self.out.fillBranch("GEN_zqq_phi", GEN_zqq_phi)
      self.out.fillBranch("GEN_zqq_mass", GEN_zqq_mass)
      self.out.fillBranch("GEN_dr_z1z2", GEN_dr_z1z2)
      self.out.fillBranch("GEN_zl1_pt", GEN_zl1_pt)
      self.out.fillBranch("GEN_zl1_eta", GEN_zl1_eta)
      self.out.fillBranch("GEN_zl1_phi", GEN_zl1_phi)
      self.out.fillBranch("GEN_zl1_mass", GEN_zl1_mass)
      self.out.fillBranch("GEN_zl2_pt", GEN_zl2_pt)
      self.out.fillBranch("GEN_zl2_eta", GEN_zl2_eta)
      self.out.fillBranch("GEN_zl2_phi", GEN_zl2_phi)
      self.out.fillBranch("GEN_zl2_mass", GEN_zl2_mass)
      self.out.fillBranch("GEN_dr_l1l2", GEN_dr_l1l2)
      self.out.fillBranch("GEN_zj1_pt", GEN_zj1_pt)
      self.out.fillBranch("GEN_zj1_pdgid", GEN_zj1_pdgid)
      self.out.fillBranch("GEN_zj1_eta", GEN_zj1_eta)
      self.out.fillBranch("GEN_zj1_phi", GEN_zj1_phi)
      self.out.fillBranch("GEN_zj1_mass", GEN_zj1_mass)
      self.out.fillBranch("GEN_zj2_pt", GEN_zj2_pt)
      self.out.fillBranch("GEN_zj2_pdgid", GEN_zj2_pdgid)
      self.out.fillBranch("GEN_zj2_eta", GEN_zj2_eta)
      self.out.fillBranch("GEN_zj2_phi", GEN_zj2_phi)
      self.out.fillBranch("GEN_zj2_mass", GEN_zj2_mass)
      self.out.fillBranch("GEN_dr_j1j2", GEN_dr_j1j2)

      return True

GenPartProducer = lambda: GetGenPart()

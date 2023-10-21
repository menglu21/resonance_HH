import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class ElectronSelector(Module):

  def __init__(self):
    pass

  def beginJob(self):
    pass
  def endJob(self):
    pass
    
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch("GoodElectron_pt", "F", lenVar="nGoodElectron")
    self.out.branch("GoodElectron_rawpt", "F", lenVar="nGoodElectron")
    self.out.branch("GoodElectron_eta", "F", lenVar="nGoodElectron")
    self.out.branch("GoodElectron_phi", "F", lenVar="nGoodElectron")
    self.out.branch("GoodElectron_mass", "F", lenVar="nGoodElectron")
    self.out.branch("GoodElectron_rawmass", "F", lenVar="nGoodElectron")
    self.out.branch("GoodElectron_id", "I", lenVar="nGoodElectron")
    self.out.branch("GoodElectron_pdgid", "I", lenVar="nGoodElectron")
    self.out.branch("FakeElectron_pt", "F", lenVar="nFakeElectron")
    self.out.branch("FakeElectron_rawpt", "F", lenVar="nFakeElectron")
    self.out.branch("FakeElectron_eta", "F", lenVar="nFakeElectron")
    self.out.branch("FakeElectron_phi", "F", lenVar="nFakeElectron")
    self.out.branch("FakeElectron_mass", "F", lenVar="nFakeElectron")
    self.out.branch("FakeElectron_rawmass", "F", lenVar="nFakeElectron")
    self.out.branch("FakeElectron_id", "I", lenVar="nFakeElectron")
    self.out.branch("FakeElectron_pdgid", "I", lenVar="nFakeElectron")
    self.out.branch("LooseElectron_pt", "F", lenVar="nLooseElectron")
    self.out.branch("LooseElectron_eta", "F", lenVar="nLooseElectron")
    self.out.branch("LooseElectron_phi", "F", lenVar="nLooseElectron")
    self.out.branch("LooseElectron_mass", "F", lenVar="nLooseElectron")
    self.out.branch("LooseElectron_id", "I", lenVar="nLooseElectron")
    self.out.branch("LooseElectron_pdgid", "I", lenVar="nLooseElectron")

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):

    electrons = Collection(event, 'Electron')

    goodElectrons_pt = []
    goodElectrons_rawpt = []
    goodElectrons_eta = []
    goodElectrons_phi = []
    goodElectrons_mass = []
    goodElectrons_rawmass = []
    goodElectrons_pdgid = []
    goodElectrons_id = []
    fakeElectrons_pt = []
    fakeElectrons_rawpt = []
    fakeElectrons_eta = []
    fakeElectrons_phi = []
    fakeElectrons_mass = []
    fakeElectrons_rawmass = []
    fakeElectrons_pdgid = []
    fakeElectrons_id = []
    looseElectrons_pt = []
    looseElectrons_eta = []
    looseElectrons_phi = []
    looseElectrons_mass = []
    looseElectrons_pdgid = []
    looseElectrons_id = []

    for iele in range(0, event.nElectron):
      if abs(electrons[iele].eta+electrons[iele].deltaEtaSC)>2.5:continue 
      if abs(electrons[iele].eta)>1.44 and abs(electrons[iele].eta<1.57):continue
      if electrons[iele].pt<10:continue
      if not electrons[iele].mvaFall17V2noIso_WPL: continue

      if electrons[iele].mvaFall17V2noIso_WP90 and electrons[iele].pt>15 and electrons[iele].sip3d<8 and abs(electrons[iele].dxy)<0.05 and abs(electrons[iele].dz)<0.1:
        # here need to get the pt before and after Egamma correction, the former one is used for jet-subtraction, the latter one is used for signal electron
        if electrons[iele].miniPFRelIso_all<0.2:
          goodElectrons_pt.append(electrons[iele].pt)
          goodElectrons_rawpt.append(electrons[iele].pt/electrons[iele].eCorr)
          goodElectrons_eta.append(electrons[iele].eta)
          goodElectrons_phi.append(electrons[iele].phi)
          goodElectrons_mass.append(electrons[iele].mass)
          goodElectrons_rawmass.append(electrons[iele].mass/electrons[iele].eCorr)
          goodElectrons_pdgid.append(electrons[iele].pdgId)
          goodElectrons_id.append(iele)
        elif electrons[iele].miniPFRelIso_all<0.4:
          fakeElectrons_pt.append(electrons[iele].pt)
          fakeElectrons_rawpt.append(electrons[iele].pt/electrons[iele].eCorr)
          fakeElectrons_eta.append(electrons[iele].eta)
          fakeElectrons_phi.append(electrons[iele].phi)
          fakeElectrons_mass.append(electrons[iele].mass)
          fakeElectrons_rawmass.append(electrons[iele].mass/electrons[iele].eCorr)
          fakeElectrons_pdgid.append(electrons[iele].pdgId)
          fakeElectrons_id.append(iele)
      else:
        looseElectrons_pt.append(electrons[iele].pt)
        looseElectrons_eta.append(electrons[iele].eta)
        looseElectrons_phi.append(electrons[iele].phi)
        looseElectrons_mass.append(electrons[iele].mass)
        looseElectrons_pdgid.append(electrons[iele].pdgId)
        looseElectrons_id.append(iele)

    self.out.fillBranch("GoodElectron_pt", goodElectrons_pt)
    self.out.fillBranch("GoodElectron_rawpt", goodElectrons_rawpt)
    self.out.fillBranch("GoodElectron_eta", goodElectrons_eta)
    self.out.fillBranch("GoodElectron_phi", goodElectrons_phi)
    self.out.fillBranch("GoodElectron_mass", goodElectrons_mass)
    self.out.fillBranch("GoodElectron_rawmass", goodElectrons_rawmass)
    self.out.fillBranch("GoodElectron_id", goodElectrons_id)
    self.out.fillBranch("GoodElectron_pdgid", goodElectrons_pdgid)
    self.out.fillBranch("FakeElectron_pt", fakeElectrons_pt)
    self.out.fillBranch("FakeElectron_rawpt", fakeElectrons_rawpt)
    self.out.fillBranch("FakeElectron_eta", fakeElectrons_eta)
    self.out.fillBranch("FakeElectron_phi", fakeElectrons_phi)
    self.out.fillBranch("FakeElectron_mass", fakeElectrons_mass)
    self.out.fillBranch("FakeElectron_rawmass", fakeElectrons_rawmass)
    self.out.fillBranch("FakeElectron_id", fakeElectrons_id)
    self.out.fillBranch("FakeElectron_pdgid", fakeElectrons_pdgid)
    self.out.fillBranch("LooseElectron_pt", looseElectrons_pt)
    self.out.fillBranch("LooseElectron_eta", looseElectrons_eta)
    self.out.fillBranch("LooseElectron_phi", looseElectrons_phi)
    self.out.fillBranch("LooseElectron_mass", looseElectrons_mass)
    self.out.fillBranch("LooseElectron_id", looseElectrons_id)
    self.out.fillBranch("LooseElectron_pdgid", looseElectrons_pdgid)

    return True

ElectronProducer = lambda: ElectronSelector()

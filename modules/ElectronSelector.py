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
    self.out.branch("GoodElectron_eta", "F", lenVar="nGoodElectron")
    self.out.branch("GoodElectron_phi", "F", lenVar="nGoodElectron")
    self.out.branch("GoodElectron_mass", "F", lenVar="nGoodElectron")
    self.out.branch("GoodElectron_id", "I", lenVar="nGoodElectron")
    self.out.branch("GoodElectron_pdgid", "I", lenVar="nGoodElectron")
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
    goodElectrons_eta = []
    goodElectrons_phi = []
    goodElectrons_mass = []
    goodElectrons_pdgid = []
    goodElectrons_id = []
    looseElectrons_pt = []
    looseElectrons_eta = []
    looseElectrons_phi = []
    looseElectrons_mass = []
    looseElectrons_pdgid = []
    looseElectrons_id = []

    for iele in range(0, event.nElectron):
      if not (abs(electrons[iele].eta)<2.5 and electrons[iele].pt>10):continue
      if not electrons[iele].mvaFall17V2noIso_WPL: continue

      if electrons[iele].mvaFall17V2noIso_WP90 and electrons[iele].pt>15:
        goodElectrons_pt.append(electrons[iele].pt)
        goodElectrons_eta.append(electrons[iele].eta)
        goodElectrons_phi.append(electrons[iele].phi)
        goodElectrons_mass.append(electrons[iele].mass)
        goodElectrons_pdgid.append(electrons[iele].pdgId)
        goodElectrons_id.append(iele)
      else:
        looseElectrons_pt.append(electrons[iele].pt)
        looseElectrons_eta.append(electrons[iele].eta)
        looseElectrons_phi.append(electrons[iele].phi)
        looseElectrons_mass.append(electrons[iele].mass)
        looseElectrons_pdgid.append(electrons[iele].pdgId)
        looseElectrons_id.append(iele)

    self.out.fillBranch("GoodElectron_pt", goodElectrons_pt)
    self.out.fillBranch("GoodElectron_eta", goodElectrons_eta)
    self.out.fillBranch("GoodElectron_phi", goodElectrons_phi)
    self.out.fillBranch("GoodElectron_mass", goodElectrons_mass)
    self.out.fillBranch("GoodElectron_id", goodElectrons_id)
    self.out.fillBranch("GoodElectron_pdgid", goodElectrons_pdgid)
    self.out.fillBranch("LooseElectron_pt", looseElectrons_pt)
    self.out.fillBranch("LooseElectron_eta", looseElectrons_eta)
    self.out.fillBranch("LooseElectron_phi", looseElectrons_phi)
    self.out.fillBranch("LooseElectron_mass", looseElectrons_mass)
    self.out.fillBranch("LooseElectron_id", looseElectrons_id)
    self.out.fillBranch("LooseElectron_pdgid", looseElectrons_pdgid)

    return True

ElectronProducer = lambda: ElectronSelector()

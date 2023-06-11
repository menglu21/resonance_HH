import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class MuonSelector(Module):

  def __init__(self):
    pass

  def beginJob(self):
    pass
  def endJob(self):
    pass
    
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch("GoodMuon_pt", "F", lenVar="nGoodMuon")
    self.out.branch("GoodMuon_rawpt", "F", lenVar="nGoodMuon")
    self.out.branch("GoodMuon_eta", "F", lenVar="nGoodMuon")
    self.out.branch("GoodMuon_phi", "F", lenVar="nGoodMuon")
    self.out.branch("GoodMuon_mass", "F", lenVar="nGoodMuon")
    self.out.branch("GoodMuon_id", "I", lenVar="nGoodMuon")
    self.out.branch("GoodMuon_pdgid", "I", lenVar="nGoodMuon")
    self.out.branch("LooseMuon_pt", "F", lenVar="nLooseMuon")
    self.out.branch("LooseMuon_eta", "F", lenVar="nLooseMuon")
    self.out.branch("LooseMuon_phi", "F", lenVar="nLooseMuon")
    self.out.branch("LooseMuon_mass", "F", lenVar="nLooseMuon")
    self.out.branch("LooseMuon_id", "I", lenVar="nLooseMuon")
    self.out.branch("LooseMuon_pdgid", "I", lenVar="nLooseMuon")

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):

    muons = Collection(event, 'Muon')

    goodMuons_pt = []
    goodMuons_rawpt = []
    goodMuons_eta = []
    goodMuons_phi = []
    goodMuons_mass = []
    goodMuons_pdgid = []
    goodMuons_id = []
    looseMuons_pt = []
    looseMuons_eta = []
    looseMuons_phi = []
    looseMuons_mass = []
    looseMuons_pdgid = []
    looseMuons_id = []

    for imu in range(0, event.nMuon):
      if abs(muons[imu].eta)>2.4: continue
      if not (muons[imu].looseId and event.Muon_corrected_pt[imu]>10):continue

      if muons[imu].mediumId:
        goodMuons_pt.append(event.Muon_corrected_pt[imu])
        goodMuons_rawpt.append(muons[imu].pt)
        goodMuons_eta.append(muons[imu].eta)
        goodMuons_phi.append(muons[imu].phi)
        goodMuons_mass.append(muons[imu].mass)
        goodMuons_pdgid.append(muons[imu].pdgId)
        goodMuons_id.append(imu)
      else:
        looseMuons_pt.append(event.Muon_corrected_pt[imu])
        looseMuons_eta.append(muons[imu].eta)
        looseMuons_phi.append(muons[imu].phi)
        looseMuons_mass.append(muons[imu].mass)
        looseMuons_pdgid.append(muons[imu].pdgId)
        looseMuons_id.append(imu)

    self.out.fillBranch("GoodMuon_pt", goodMuons_pt)
    self.out.fillBranch("GoodMuon_rawpt", goodMuons_rawpt)
    self.out.fillBranch("GoodMuon_eta", goodMuons_eta)
    self.out.fillBranch("GoodMuon_phi", goodMuons_phi)
    self.out.fillBranch("GoodMuon_mass", goodMuons_mass)
    self.out.fillBranch("GoodMuon_id", goodMuons_id)
    self.out.fillBranch("GoodMuon_pdgid", goodMuons_pdgid)
    self.out.fillBranch("LooseMuon_pt", looseMuons_pt)
    self.out.fillBranch("LooseMuon_eta", looseMuons_eta)
    self.out.fillBranch("LooseMuon_phi", looseMuons_phi)
    self.out.fillBranch("LooseMuon_mass", looseMuons_mass)
    self.out.fillBranch("LooseMuon_id", looseMuons_id)
    self.out.fillBranch("LooseMuon_pdgid", looseMuons_pdgid)

    return True

MuonProducer = lambda: MuonSelector()

import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True
from numpy import sqrt

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class ShuffleSubGenJetMerger(Module):

  def __init__(self):
    pass

  def beginJob(self):
    pass
  def endJob(self):
    pass
    
  def dr_etaphi(eta1, phi1, eta2, phi2):
    if (phi1 - phi2)>3.14159: 
      dphi = phi1 - phi2 - 2*3.14159
    elif (phi1 - phi2)<-3.14159:
      dphi = 2*3.14159 - (phi1 - phi2)
    else:
      dphi = phi1 - phi2
    deta = eta1 - eta2
    return sqrt(deta*deta + dphi*dphi)

  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    if self.is_mc:
      SubGenJet=Collection(event,"SubGenJetAK8NoVlep")
      GenJetAK8=Collection(event,"GenJetAK8NoVlep")
      _pt=[]
      _eta=[]
      _phi=[]
      _mass=[]
      _dRinit=[]

  
      if len(GenJetAK8)>0 and len(SubGenJet)>0:
        for ij8 in range(0,len(GenJetAK8)):
          _subid1=GenJetAK8[ij8].subj1ID
          _subid2=GenJetAK8[ij8].subj2ID
          if _subid1>-1 and _subid2>-1:
            if SubGenJet[_subid1].pt > SubGenJet[_subid2].pt:
              _pt.append(SubGenJet[_subid1].pt)
              _pt.append(SubGenJet[_subid2].pt)
              _eta.append(SubGenJet[_subid1].eta)
              _eta.append(SubGenJet[_subid2].eta)
              _phi.append(SubGenJet[_subid1].phi)
              _phi.append(SubGenJet[_subid2].phi)
              _mass.append(SubGenJet[_subid1].mass)
              _mass.append(SubGenJet[_subid2].mass)
              _dRinit.append(SubGenJet[_subid1].dRinit)
              _dRinit.append(SubGenJet[_subid2].dRinit)
            else:
              _pt.append(SubGenJet[_subid2].pt)
              _pt.append(SubGenJet[_subid1].pt)
              _eta.append(SubGenJet[_subid2].eta)
              _eta.append(SubGenJet[_subid1].eta)
              _phi.append(SubGenJet[_subid2].phi)
              _phi.append(SubGenJet[_subid1].phi)
              _mass.append(SubGenJet[_subid2].mass)
              _mass.append(SubGenJet[_subid1].mass)
              _dRinit.append(SubGenJet[_subid2].dRinit)
              _dRinit.append(SubGenJet[_subid1].dRinit)
          elif _subid1<0 and _subid2>-1:
            _pt.append(SubGenJet[_subid2].pt)
            _eta.append(SubGenJet[_subid2].eta)
            _phi.append(SubGenJet[_subid2].phi)
            _mass.append(SubGenJet[_subid2].mass)
            _dRinit.append(SubGenJet[_subid2].dRinit)
          elif _subid1>-1 and _subid2<0:
            _pt.append(SubGenJet[_subid1].pt)
            _eta.append(SubGenJet[_subid1].eta)
            _phi.append(SubGenJet[_subid1].phi)
            _mass.append(SubGenJet[_subid1].mass)
            _dRinit.append(SubGenJet[_subid1].dRinit)
        
        self.out.fillBranch("SubGenJetAK8NoVlep_pt", _pt)
        self.out.fillBranch("SubGenJetAK8NoVlep_eta", _eta)
        self.out.fillBranch("SubGenJetAK8NoVlep_phi", _phi)
        self.out.fillBranch("SubGenJetAK8NoVlep_mass", _mass)
        self.out.fillBranch("SubGenJetAK8NoVlep_dRinit", _dRinit)

    return True

ShuffleSubGenJet = lambda: ShuffleSubGenJetMerger()

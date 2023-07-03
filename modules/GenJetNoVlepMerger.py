import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class GenJetNoVlepMerger(Module):

  def __init__(self,inputs,output):
    self.inputs = inputs
    self.output = output
    pass

  def beginJob(self):
    pass
  def endJob(self):
    pass
    
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))
    if self.is_mc:
      self.out.branch("%s_pt"%(self.output), "F", lenVar="n%s"%self.output)
      self.out.branch("%s_eta"%(self.output), "F", lenVar="n%s"%self.output)
      self.out.branch("%s_phi"%(self.output), "F", lenVar="n%s"%self.output)
      self.out.branch("%s_mass"%(self.output), "F", lenVar="n%s"%self.output)
      self.out.branch("%s_dRinit"%(self.output), "F", lenVar="n%s"%self.output)
      if "GenJet"==self.inputs or "GenJetAK8"==self.inputs:
        self.out.branch("%s_initID"%(self.output), "I", lenVar="n%s"%self.output)
      if "GenJetAK8"==self.inputs:
        self.out.branch("%s_subj1ID"%(self.output), "I", lenVar="n%s"%self.output)
        self.out.branch("%s_subj2ID"%(self.output), "I", lenVar="n%s"%self.output)

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    if self.is_mc:
      jetcoll=Collection(event,"%s"%self.inputs)
      Gendressedlep=Collection(event,"GenDressedLepton")
      _pt=[]
      _eta=[]
      _phi=[]
      _mass=[]
      _initID=[]
      _dRinit=[]
      if "GenJetAK8"==self.inputs:
        SubGenJet=Collection(event,"SubGenJetAK8")
  
      dr=100.
      if "GenJet"==self.inputs or "SubGenJetAK8"==self.inputs:
        dr=0.4

      # firstly get the subjet of the corresponding AK8Jet
      if "GenJetAK8"==self.inputs:
        dr=0.8
        _subj1ID=[]
        _subj2ID=[]
        for ij in range(len(jetcoll)):
          _id_tmp1=-1
          _id_tmp2=-1
          tmp_ak8_v4=TLorentzVector()
          tmp_ak8_v4.SetPtEtaPhiM(jetcoll[ij].pt,jetcoll[ij].eta,jetcoll[ij].phi,jetcoll[ij].mass)
          for ijj in range(len(SubGenJet)):
            tmp_subak8_v4=TLorentzVector()
            tmp_subak8_v4.SetPtEtaPhiM(SubGenJet[ijj].pt,SubGenJet[ijj].eta,SubGenJet[ijj].phi,SubGenJet[ijj].mass)
            if tmp_ak8_v4.DeltaR(tmp_subak8_v4)<0.8:
              if _id_tmp1<0:_id_tmp1=ijj
              elif _id_tmp2<0:_id_tmp2=ijj
          _subj1ID.append(_id_tmp1)
          _subj2ID.append(_id_tmp2)
   
          
  
      lepv4_arr=[]
      lepv4_=TLorentzVector()
  
      if len(Gendressedlep)>0:
        for il in range(0,len(Gendressedlep)):
          lepv4_.SetPtEtaPhiM(Gendressedlep[il].pt,Gendressedlep[il].eta,Gendressedlep[il].phi,Gendressedlep[il].mass)
          lepv4_arr.append(lepv4_.Clone())
  
   
      jetv4_=TLorentzVector()
      jetv4_arr=[]
      jetv4_arr_init=[]

      if len(jetcoll)>0:
        for ij in range(0,len(jetcoll)):
          jetv4_.SetPtEtaPhiM(jetcoll[ij].pt,jetcoll[ij].eta,jetcoll[ij].phi,jetcoll[ij].mass) 
          jetv4_arr.append(jetv4_.Clone())
          jetv4_arr_init.append(jetv4_.Clone())


      if len(Gendressedlep)>0:
        if len(jetcoll)>1:
          for il in range(0,len(Gendressedlep)):
            lep_dr=[]
            lep_dr_copy=[]
            for ij in range(0,len(jetcoll)):
              lep_dr.append(jetv4_arr[ij].DeltaR(lepv4_arr[il]))
              lep_dr_copy.append(jetv4_arr[ij].DeltaR(lepv4_arr[il]))
     
            # get the minimum deltaR between current lepton and all jets
            min_index=lep_dr.index(min(lep_dr))
            lep_dr_copy.remove(min(lep_dr))
            submin_index=lep_dr.index(min(lep_dr_copy))
            
            # subtract lepton from the jet if dR<dr_set. N.B., one lepton can be subtracted from both jet and fatjet, because jet can be the one reconstructed as fatjet
            if min(lep_dr)<dr:
              if (jetv4_arr[min_index] - lepv4_arr[il]).M()>0:
                jetv4_arr[min_index]=jetv4_arr[min_index] - lepv4_arr[il]
              elif min(lep_dr_copy)<dr and (jetv4_arr[submin_index] - lepv4_arr[il]).M()>0:
                jetv4_arr[submin_index]=jetv4_arr[submin_index] - lepv4_arr[il]

        elif len(jetcoll)>0:
          for il in range(0,len(Gendressedlep)):
            if lepv4_arr[il].DeltaR(jetv4_arr[0])<dr:
              if (jetv4_arr[0]-lepv4_arr[il]).M()>0:
                jetv4_arr[0] = jetv4_arr[0]-lepv4_arr[il]

      for ij in range(0,len(jetcoll)):
        _pt.append(jetv4_arr[ij].Pt())
        _eta.append(jetv4_arr[ij].Eta())
        _phi.append(jetv4_arr[ij].Phi())
        _mass.append(jetv4_arr[ij].M())
        _initID.append(ij)
        _dRinit.append(jetv4_arr[ij].DeltaR(jetv4_arr_init[ij]))

      # sort those lepton-subtracted GenJets according to their pT
      if "GenJet"==self.inputs:
        jet_sortedArr = sorted(zip(_pt, _eta, _phi, _mass, _initID, _dRinit), key=lambda x: x[0], reverse=True)

        _pt=[x[0] for x in jet_sortedArr]
        _eta=[x[1] for x in jet_sortedArr]
        _phi=[x[2] for x in jet_sortedArr]
        _mass=[x[3] for x in jet_sortedArr]
        _initID=[x[4] for x in jet_sortedArr]
        _dRinit=[x[5] for x in jet_sortedArr]

      if "GenJetAK8"==self.inputs:
        jet_sortedArr = sorted(zip(_pt, _eta, _phi, _mass, _initID, _dRinit, _subj1ID, _subj2ID), key=lambda x: x[0], reverse=True)

        _pt=[x[0] for x in jet_sortedArr]
        _eta=[x[1] for x in jet_sortedArr]
        _phi=[x[2] for x in jet_sortedArr]
        _mass=[x[3] for x in jet_sortedArr]
        _initID=[x[4] for x in jet_sortedArr]
        _dRinit=[x[5] for x in jet_sortedArr]
        _subj1ID=[x[6] for x in jet_sortedArr]
        _subj2ID=[x[7] for x in jet_sortedArr]
  
  
      self.out.fillBranch("%s_pt"%(self.output), _pt)
      self.out.fillBranch("%s_eta"%(self.output), _eta)
      self.out.fillBranch("%s_phi"%(self.output), _phi)
      self.out.fillBranch("%s_mass"%(self.output), _mass)
      self.out.fillBranch("%s_dRinit"%(self.output), _dRinit)
      if "GenJet"==self.inputs or "GenJetAK8"==self.inputs:
        self.out.fillBranch("%s_initID"%(self.output), _initID)
      if "GenJetAK8"==self.inputs:
        self.out.fillBranch("%s_subj1ID"%(self.output), _subj1ID)
        self.out.fillBranch("%s_subj2ID"%(self.output), _subj2ID)

    return True

GenJetAK4NoVlepMerger = lambda: GenJetNoVlepMerger("GenJet","GenJetNoVlep")
GenJetAK8NoVlepMerger = lambda: GenJetNoVlepMerger("GenJetAK8","GenJetAK8NoVlep")
SubGenJetAK8NoVlepMerger = lambda: GenJetNoVlepMerger("SubGenJetAK8","SubGenJetAK8NoVlep")

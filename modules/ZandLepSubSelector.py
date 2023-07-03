import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class ZandLepSubSelector(Module):

  def __init__(self):
    pass

  def beginJob(self):
    pass
  def endJob(self):
    pass
    
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree

    self.out.branch("JetNoVlep_pt", "F", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_eta", "F", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_phi", "F", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_mass", "F", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_rawFactor", "F", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_area", "F", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_id", "I", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_jetId", "I", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_muonIdx1", "I", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_muonIdx2", "I", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_muonSubtrFactor", "F", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_neEmEF", "F", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_chEmEF", "F", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_drl1", "F", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_drl2", "F", lenVar="nJetNoVlep")
    self.out.branch("JetNoVlep_dRinit", "F", lenVar="nJetNoVlep")
    self.out.branch("SubJetNoVlep_pt", "F", lenVar="nSubJetNoVlep")
    self.out.branch("SubJetNoVlep_eta", "F", lenVar="nSubJetNoVlep")
    self.out.branch("SubJetNoVlep_phi", "F", lenVar="nSubJetNoVlep")
    self.out.branch("SubJetNoVlep_mass", "F", lenVar="nSubJetNoVlep")
    self.out.branch("SubJetNoVlep_rawFactor", "F", lenVar="nSubJetNoVlep")
    self.out.branch("SubJetNoVlep_dRinit", "F", lenVar="nSubJetNoVlep")
    self.out.branch("SubJetNoVlep_id", "I", lenVar="nSubJetNoVlep")
    self.out.branch("FatJetNoVlep_pt", "F", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_eta", "F", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_phi", "F", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_mass", "F", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_PNmass", "F", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_SDmass", "F", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_rawFactor", "F", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_area", "F", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_id", "I", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_jetId", "I", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_overlapl1", "I", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_overlapl2", "I", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_subJetIdx1", "I", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_subJetIdx2", "I", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_drl1", "F", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_drl2", "F", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_dRinit", "F", lenVar="nFatJetNoVlep")
    self.out.branch("FatJetNoVlep_bbvsccqq", "F", lenVar="nFatJetNoVlep")
    self.out.branch("mu_channel", "B")
    self.out.branch("l1_pt", "F")
    self.out.branch("l1_rawpt", "F")
    self.out.branch("l1_eta", "F")
    self.out.branch("l1_phi", "F")
    self.out.branch("l1_mass", "F")
    self.out.branch("l1_rawmass", "F")
    self.out.branch("l1_relPtl2", "F")
    self.out.branch("l2_pt", "F")
    self.out.branch("l2_rawpt", "F")
    self.out.branch("l2_eta", "F")
    self.out.branch("l2_phi", "F")
    self.out.branch("l2_mass", "F")
    self.out.branch("l2_rawmass", "F")
    self.out.branch("l2_relPtl1", "F")
    self.out.branch("drll", "F")
    self.out.branch("detall", "F")
    self.out.branch("dphill", "F")
    self.out.branch("zlep_pt", "F")
    self.out.branch("zlep_eta", "F")
    self.out.branch("zlep_phi", "F")
    self.out.branch("zlep_mass", "F")

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):

    muons = Collection(event, 'Muon')
    goodmuons = Collection(event, 'GoodMuon')
    goodeles = Collection(event, 'GoodElectron')
  
    recojets = Collection(event, 'Jet')
    recosubjets = Collection(event, 'SubJet')
    recofatjets = Collection(event, 'FatJet')

    mu_channel=False
    JetNoVlep_pt=[]
    JetNoVlep_eta=[]
    JetNoVlep_phi=[]
    JetNoVlep_mass=[]
    JetNoVlep_rawFactor=[]
    JetNoVlep_area=[]
    JetNoVlep_id=[]
    JetNoVlep_jetId=[]
    JetNoVlep_muonIdx1=[]
    JetNoVlep_muonIdx2=[]
    JetNoVlep_muonSubtrFactor=[]
    JetNoVlep_neEmEF=[]
    JetNoVlep_chEmEF=[]
    JetNoVlep_drl1=[]
    JetNoVlep_drl2=[]
    JetNoVlep_dRinit=[]
    SubJetNoVlep_pt=[]
    SubJetNoVlep_eta=[]
    SubJetNoVlep_phi=[]
    SubJetNoVlep_mass=[]
    SubJetNoVlep_dRinit=[]
    SubJetNoVlep_rawFactor=[]
    SubJetNoVlep_id=[]
    FatJetNoVlep_pt=[]
    FatJetNoVlep_eta=[]
    FatJetNoVlep_phi=[]
    FatJetNoVlep_mass=[]
    FatJetNoVlep_PNmass=[]
    FatJetNoVlep_SDmass=[]
    FatJetNoVlep_rawFactor=[]
    FatJetNoVlep_area=[]
    FatJetNoVlep_id=[]
    FatJetNoVlep_jetId=[]
    FatJetNoVlep_overlapl1=[]
    FatJetNoVlep_overlapl2=[]
    FatJetNoVlep_subJetIdx1=[]
    FatJetNoVlep_subJetIdx2=[]
    FatJetNoVlep_drl1=[]
    FatJetNoVlep_drl2=[]
    FatJetNoVlep_dRinit=[]
    FatJetNoVlep_bbvsccqq=[]
    l1_pt=-99
    l1_rawpt=-99
    l1_eta=-99
    l1_phi=-99
    l1_mass=-99
    l1_rawmass=-99
    l1_relPtl2=-99
    l2_pt=-99
    l2_rawpt=-99
    l2_eta=-99
    l2_phi=-99
    l2_mass=-99
    l2_rawmass=-99
    l2_relPtl1=-99
    drll=-99
    detall=-99
    dphill=-99
    zlep_pt=-99
    zlep_eta=-99
    zlep_phi=-99
    zlep_mass=-99

    # only two tight leptons
    if not (event.nGoodMuon+event.nGoodElectron)==2:return False
    if not (event.nGoodMuon==2 or event.nGoodElectron==2):return False

    l1v4_tmp=TLorentzVector()
    l2v4_tmp=TLorentzVector()
    l1v4raw_tmp=TLorentzVector()
    l2v4raw_tmp=TLorentzVector()

    if event.nGoodMuon==2:
      mu_channel=True
      l1_pt=goodmuons[0].pt
      l1_rawpt=goodmuons[0].rawpt
      l1_eta=goodmuons[0].eta
      l1_phi=goodmuons[0].phi
      l1_mass=goodmuons[0].mass
      l1_rawmass=goodmuons[0].rawmass
      l2_pt=goodmuons[1].pt
      l2_rawpt=goodmuons[1].rawpt
      l2_eta=goodmuons[1].eta
      l2_phi=goodmuons[1].phi
      l2_mass=goodmuons[1].mass
      l2_rawmass=goodmuons[1].rawmass
    else:
      mu_channel=False
      l1_pt=goodeles[0].pt
      l1_rawpt=goodeles[0].rawpt
      l1_eta=goodeles[0].eta
      l1_phi=goodeles[0].phi
      l1_mass=goodeles[0].mass
      l1_rawmass=goodeles[0].rawmass
      l2_pt=goodeles[1].pt
      l2_rawpt=goodeles[1].rawpt
      l2_eta=goodeles[1].eta
      l2_phi=goodeles[1].phi
      l2_mass=goodeles[1].mass
      l2_rawmass=goodeles[1].rawmass

    l1v4_tmp.SetPtEtaPhiM(l1_pt,l1_eta,l1_phi,l1_mass)
    l2v4_tmp.SetPtEtaPhiM(l2_pt,l2_eta,l2_phi,l2_mass)
    l1v4raw_tmp.SetPtEtaPhiM(l1_rawpt,l1_eta,l1_phi,l1_rawmass)
    l2v4raw_tmp.SetPtEtaPhiM(l2_rawpt,l2_eta,l2_phi,l2_rawmass)
    l1_relPtl2=l1v4_tmp.Vect().Cross(l2v4_tmp.Vect().Unit()).Mag()
    l2_relPtl1=l2v4_tmp.Vect().Cross(l1v4_tmp.Vect().Unit()).Mag()
    drll=l1v4_tmp.DeltaR(l2v4_tmp)
    detall=abs(l1_eta - l2_eta)
    dphill=l1v4_tmp.DeltaPhi(l2v4_tmp)
    zlep_pt=(l1v4_tmp + l2v4_tmp).Pt()
    zlep_eta=(l1v4_tmp + l2v4_tmp).Eta()
    zlep_phi=(l1v4_tmp + l2v4_tmp).Phi()
    zlep_mass=(l1v4_tmp + l2v4_tmp).M()

    # re-calculate jet information after leptons subtraction
    if event.nJet>0:
      ak4jet_drl1_=[]
      ak4jet_drl2_=[]
      ak4jet_arr_=[]
      ak4jet_arr2_=[]
      chEmEF_total=[]
      neEmEF_total=[]
      rawenergy_total=[]
      for ij in range(0,event.nJet):
        JetNoVlep_area.append(recojets[ij].area)
        JetNoVlep_jetId.append(recojets[ij].jetId)
        jetv4_tmp=TLorentzVector()
        jetv4_tmp.SetPtEtaPhiM((1-recojets[ij].rawFactor)*recojets[ij].pt,recojets[ij].eta,recojets[ij].phi,(1-recojets[ij].rawFactor)*recojets[ij].mass)
        ak4jet_arr_.append(jetv4_tmp.Clone())
        ak4jet_arr2_.append(jetv4_tmp.Clone())
        rawenergy_total.append(jetv4_tmp.E())
        
        # recalculate chEmEF and neEmEF
        # https://github.com/cms-sw/cmssw/blob/b1e36dcb824e36bbbb8cb30963c988890b0b28c1/PhysicsTools/Heppy/python/physicsobjects/Jet.py#L78
        # https://github.com/cms-sw/cmssw/blob/b1e36dcb824e36bbbb8cb30963c988890b0b28c1/RecoJets/JetProducers/src/JetSpecific.cc#L254-L300
        chEmEF_total.append(jetv4_tmp.E()*recojets[ij].chEmEF)
        neEmEF_total.append(jetv4_tmp.E()*recojets[ij].neEmEF)

        jetmatched_mu1=-1
        jetmatched_mu2=-1

        # re-select the matched muons except the signal muon overlap with Jets
        if event.nMuon>0:
          muv4_tmp=TLorentzVector()
          for imu in range(0,event.nMuon):
            muv4_tmp.SetPtEtaPhiM(event.Muon_pt[imu],event.Muon_eta[imu],event.Muon_phi[imu],event.Muon_mass[imu])
            if jetv4_tmp.DeltaR(muv4_tmp)<0.4:
              if mu_channel:
                if not (imu==event.GoodMuon_id[0] or imu==event.GoodMuon_id[1]):
                  if jetmatched_mu1<0:jetmatched_mu1=imu
                  elif jetmatched_mu2<0:jetmatched_mu2=imu
              else:
                if jetmatched_mu1<0:jetmatched_mu1=imu
                elif jetmatched_mu2<0:jetmatched_mu2=imu
        JetNoVlep_muonIdx1.append(jetmatched_mu1)
        JetNoVlep_muonIdx2.append(jetmatched_mu2)

        ak4jet_drl1_.append(l1v4_tmp.DeltaR(jetv4_tmp))
        ak4jet_drl2_.append(l2v4_tmp.DeltaR(jetv4_tmp))

      min_index_ak4jet_drl1=ak4jet_drl1_.index(min(ak4jet_drl1_))
      min_index_ak4jet_drl2=ak4jet_drl2_.index(min(ak4jet_drl2_))

      # be careful of the case: lepton overlap with jet but not clustered in the jet, which mean even they have dR smaller than 0.4, the subtraction is not needed otherwise it will cause negative Jet mass
      
      if min(ak4jet_drl1_)<0.4:
        if min(ak4jet_drl2_)<0.4:
          # the closest jet of l1 and l2 are the same one
          if min_index_ak4jet_drl1==min_index_ak4jet_drl2:
            if (ak4jet_arr_[min_index_ak4jet_drl1]-l1v4raw_tmp-l2v4raw_tmp).M()>0:
              ak4jet_arr_[min_index_ak4jet_drl1]=ak4jet_arr_[min_index_ak4jet_drl1]-l1v4raw_tmp-l2v4raw_tmp
            elif (ak4jet_arr_[min_index_ak4jet_drl1]-l1v4raw_tmp).M()>0:
              ak4jet_arr_[min_index_ak4jet_drl1]=ak4jet_arr_[min_index_ak4jet_drl1]-l1v4raw_tmp
            elif (ak4jet_arr_[min_index_ak4jet_drl1]-l2v4raw_tmp).M()>0:
              ak4jet_arr_[min_index_ak4jet_drl1]=ak4jet_arr_[min_index_ak4jet_drl1]-l2v4raw_tmp
          # the closest jet of l1 and l2 are not the same one
          else:
            if (ak4jet_arr_[min_index_ak4jet_drl1]-l1v4raw_tmp).M()>0:
              ak4jet_arr_[min_index_ak4jet_drl1]=ak4jet_arr_[min_index_ak4jet_drl1]-l1v4raw_tmp
            if (ak4jet_arr_[min_index_ak4jet_drl2]-l2v4raw_tmp).M()>0:
              ak4jet_arr_[min_index_ak4jet_drl2]=ak4jet_arr_[min_index_ak4jet_drl2]-l2v4raw_tmp
        else:
          if (ak4jet_arr_[min_index_ak4jet_drl1]-l1v4raw_tmp).M()>0:
            ak4jet_arr_[min_index_ak4jet_drl1]=ak4jet_arr_[min_index_ak4jet_drl1]-l1v4raw_tmp
      else:
        if min(ak4jet_drl2_)<0.4:
          if (ak4jet_arr_[min_index_ak4jet_drl2]-l2v4raw_tmp).M()>0:
            ak4jet_arr_[min_index_ak4jet_drl2]=ak4jet_arr_[min_index_ak4jet_drl2]-l2v4raw_tmp


      for ij in range(0,event.nJet):

        JetNoVlep_eta.append(ak4jet_arr_[ij].Eta())
        JetNoVlep_phi.append(ak4jet_arr_[ij].Phi())
        JetNoVlep_id.append(ij)

        # if jet undergo the lepton subtraction, store the raw pt and store 0 as the rawFactor
        if not ( ak4jet_arr_[ij].Pt()== ak4jet_arr2_[ij].Pt() ):
          JetNoVlep_rawFactor.append(0.0)
          JetNoVlep_pt.append(ak4jet_arr_[ij].Pt())
          JetNoVlep_mass.append(ak4jet_arr_[ij].M())
          rawenergy_total[ij] = ak4jet_arr_[ij].E()
        else:
          JetNoVlep_rawFactor.append(recojets[ij].rawFactor)
          JetNoVlep_pt.append(recojets[ij].pt)
          JetNoVlep_mass.append(recojets[ij].mass)

        JetNoVlep_neEmEF.append(1.*neEmEF_total[ij]/rawenergy_total[ij])
        JetNoVlep_chEmEF.append(1.*chEmEF_total[ij]/rawenergy_total[ij])
        JetNoVlep_drl1.append(ak4jet_arr_[ij].DeltaR(l1v4_tmp))
        JetNoVlep_drl2.append(ak4jet_arr_[ij].DeltaR(l2v4_tmp))
        JetNoVlep_dRinit.append(ak4jet_arr_[ij].DeltaR(ak4jet_arr2_[ij]))

        ak4_tmp_v4_=ak4jet_arr_[ij].Clone()
        #construct the muonSubtrFactor, refer to https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetmetUncertainties.py#L500-L519
        if JetNoVlep_muonIdx1[ij]>-1 and muons[JetNoVlep_muonIdx1[ij]].isGlobal:
          ak4_tmp_v4_ = ak4jet_arr_[ij] - muons[JetNoVlep_muonIdx1[ij]].p4()
        if JetNoVlep_muonIdx2[ij]>-1 and muons[JetNoVlep_muonIdx2[ij]].isGlobal:
          ak4_tmp_v4_ = ak4jet_arr_[ij] - muons[JetNoVlep_muonIdx2[ij]].p4()

        JetNoVlep_muonSubtrFactor.append(1. - ak4_tmp_v4_.Pt()/ak4jet_arr_[ij].Pt())


    JetNoVlep_array = sorted(zip(JetNoVlep_pt,JetNoVlep_eta,JetNoVlep_phi,JetNoVlep_mass,JetNoVlep_rawFactor,JetNoVlep_id,JetNoVlep_area,JetNoVlep_jetId,JetNoVlep_muonIdx1,JetNoVlep_muonIdx2,JetNoVlep_muonSubtrFactor,JetNoVlep_neEmEF,JetNoVlep_chEmEF,JetNoVlep_drl1,JetNoVlep_drl2,JetNoVlep_dRinit), key=lambda x: x[0], reverse=True)
    JetNoVlep_pt=[x[0] for x in JetNoVlep_array]
    JetNoVlep_eta=[x[1] for x in JetNoVlep_array]
    JetNoVlep_phi=[x[2] for x in JetNoVlep_array]
    JetNoVlep_mass=[x[3] for x in JetNoVlep_array]
    JetNoVlep_rawFactor=[x[4] for x in JetNoVlep_array]
    JetNoVlep_id=[x[5] for x in JetNoVlep_array]
    JetNoVlep_area=[x[6] for x in JetNoVlep_array]
    JetNoVlep_jetId=[x[7] for x in JetNoVlep_array]
    JetNoVlep_muonIdx1=[x[8] for x in JetNoVlep_array]
    JetNoVlep_muonIdx2=[x[9] for x in JetNoVlep_array]
    JetNoVlep_muonSubtrFactor=[x[10] for x in JetNoVlep_array]
    JetNoVlep_neEmEF=[x[11] for x in JetNoVlep_array]
    JetNoVlep_chEmEF=[x[12] for x in JetNoVlep_array]
    JetNoVlep_drl1=[x[13] for x in JetNoVlep_array]
    JetNoVlep_drl2=[x[14] for x in JetNoVlep_array]
    JetNoVlep_dRinit=[x[15] for x in JetNoVlep_array]


    # re-calculate fatjet information after leptons subtraction
    if event.nFatJet>0:
      for ij in range(0,event.nFatJet):
        if (recofatjets[ij].particleNetMD_Xbb+recofatjets[ij].particleNetMD_Xcc+recofatjets[ij].particleNetMD_Xqq)==0:
          FatJetNoVlep_bbvsccqq.append(0)
        else:
          FatJetNoVlep_bbvsccqq.append(recofatjets[ij].particleNetMD_Xbb/(recofatjets[ij].particleNetMD_Xbb+recofatjets[ij].particleNetMD_Xcc+recofatjets[ij].particleNetMD_Xqq))
        FatJetNoVlep_area.append(recofatjets[ij].area)
        FatJetNoVlep_jetId.append(recofatjets[ij].jetId)
        FatJetNoVlep_subJetIdx1.append(recofatjets[ij].subJetIdx1)
        FatJetNoVlep_subJetIdx2.append(recofatjets[ij].subJetIdx2)
        fatjetv4_tmp=TLorentzVector()
        fatjetv4_tmp2=TLorentzVector()
        fatjetv4_tmp.SetPtEtaPhiM((1-recofatjets[ij].rawFactor)*recofatjets[ij].pt,recofatjets[ij].eta,recofatjets[ij].phi,(1-recofatjets[ij].rawFactor)*recofatjets[ij].mass)
        fatjetv4_tmp2.SetPtEtaPhiM((1-recofatjets[ij].rawFactor)*recofatjets[ij].pt,recofatjets[ij].eta,recofatjets[ij].phi,(1-recofatjets[ij].rawFactor)*recofatjets[ij].mass)
        fatj_overlap_l1=0
        fatj_overlap_l2=0
        if l1v4_tmp.DeltaR(fatjetv4_tmp)<0.8:
          fatj_overlap_l1=1
          FatJetNoVlep_overlapl1.append(1)
        else:
          FatJetNoVlep_overlapl1.append(0)
        if l2v4_tmp.DeltaR(fatjetv4_tmp)<0.8:
          fatj_overlap_l2=1
          FatJetNoVlep_overlapl2.append(1)
        else:
          FatJetNoVlep_overlapl2.append(0)

        if fatj_overlap_l1>0:
          if fatj_overlap_l2>0:
            if (fatjetv4_tmp - l1v4raw_tmp - l2v4raw_tmp).M()>0:
              fatjetv4_tmp=fatjetv4_tmp - l1v4raw_tmp - l2v4raw_tmp
            elif (fatjetv4_tmp - l1v4raw_tmp).M()>0:
              fatjetv4_tmp=fatjetv4_tmp - l1v4raw_tmp
            elif (fatjetv4_tmp - l2v4raw_tmp).M()>0:
              fatjetv4_tmp=fatjetv4_tmp - l2v4raw_tmp

        FatJetNoVlep_eta.append(fatjetv4_tmp.Eta())
        FatJetNoVlep_phi.append(fatjetv4_tmp.Phi())
        FatJetNoVlep_id.append(ij)
        FatJetNoVlep_PNmass.append(event.FatJet_particleNet_mass[ij])
        FatJetNoVlep_SDmass.append(event.FatJet_msoftdrop[ij])

        if not (fatjetv4_tmp.Pt() == fatjetv4_tmp2.Pt()):
          FatJetNoVlep_rawFactor.append(0.0)
          FatJetNoVlep_pt.append(fatjetv4_tmp.Pt())
          FatJetNoVlep_mass.append(fatjetv4_tmp.M())
        else:
          FatJetNoVlep_rawFactor.append(recofatjets[ij].rawFactor)
          FatJetNoVlep_pt.append(recofatjets[ij].pt)
          FatJetNoVlep_mass.append(recofatjets[ij].mass)

        FatJetNoVlep_drl1.append(fatjetv4_tmp.DeltaR(l1v4_tmp))
        FatJetNoVlep_drl2.append(fatjetv4_tmp.DeltaR(l2v4_tmp))
        FatJetNoVlep_dRinit.append(fatjetv4_tmp.DeltaR(fatjetv4_tmp2))

    FatJetNoVlep_array = sorted(zip(FatJetNoVlep_pt,FatJetNoVlep_eta,FatJetNoVlep_phi,FatJetNoVlep_mass,FatJetNoVlep_rawFactor,FatJetNoVlep_id,FatJetNoVlep_overlapl1,FatJetNoVlep_overlapl2,FatJetNoVlep_area,FatJetNoVlep_jetId,FatJetNoVlep_subJetIdx1,FatJetNoVlep_subJetIdx2,FatJetNoVlep_drl1,FatJetNoVlep_drl2,FatJetNoVlep_dRinit,FatJetNoVlep_PNmass,FatJetNoVlep_SDmass,FatJetNoVlep_bbvsccqq), key=lambda x: x[0], reverse=True)
    FatJetNoVlep_pt=[x[0] for x in FatJetNoVlep_array]
    FatJetNoVlep_eta=[x[1] for x in FatJetNoVlep_array]
    FatJetNoVlep_phi=[x[2] for x in FatJetNoVlep_array]
    FatJetNoVlep_mass=[x[3] for x in FatJetNoVlep_array]
    FatJetNoVlep_rawFactor=[x[4] for x in FatJetNoVlep_array]
    FatJetNoVlep_id=[x[5] for x in FatJetNoVlep_array]
    FatJetNoVlep_overlapl1=[x[6] for x in FatJetNoVlep_array]
    FatJetNoVlep_overlapl2=[x[7] for x in FatJetNoVlep_array]
    FatJetNoVlep_area=[x[8] for x in FatJetNoVlep_array]
    FatJetNoVlep_jetId=[x[9] for x in FatJetNoVlep_array]
    FatJetNoVlep_subJetIdx1=[x[10] for x in FatJetNoVlep_array]
    FatJetNoVlep_subJetIdx2=[x[11] for x in FatJetNoVlep_array]
    FatJetNoVlep_drl1=[x[12] for x in FatJetNoVlep_array]
    FatJetNoVlep_drl2=[x[13] for x in FatJetNoVlep_array]
    FatJetNoVlep_dRinit=[x[14] for x in FatJetNoVlep_array]
    FatJetNoVlep_PNmass=[x[15] for x in FatJetNoVlep_array]
    FatJetNoVlep_SDmass=[x[16] for x in FatJetNoVlep_array]
    FatJetNoVlep_bbvsccqq=[x[17] for x in FatJetNoVlep_array]


    # re-calculate subjet information after leptons subtraction
    if event.nSubJet>0:
      subjet_drl1_=[]
      subjet_drl2_=[]
      subjetv4_arr_=[]
      subjetv4_arr2_=[]
      for ij in range(0,event.nSubJet):
        subjetv4_tmp=TLorentzVector()
        subjetv4_tmp.SetPtEtaPhiM((1-recosubjets[ij].rawFactor)*recosubjets[ij].pt,recosubjets[ij].eta,recosubjets[ij].phi,(1-recosubjets[ij].rawFactor)*recosubjets[ij].mass)
        subjet_drl1_.append(subjetv4_tmp.DeltaR(l1v4_tmp))
        subjet_drl2_.append(subjetv4_tmp.DeltaR(l2v4_tmp))
        subjetv4_arr_.append(subjetv4_tmp.Clone())
        subjetv4_arr2_.append(subjetv4_tmp.Clone())

      # get the second smallest dr
      if event.nSubJet>1:
        subjet_drl1_copy=subjet_drl1_[:]
        subjet_drl1_copy.remove(min(subjet_drl1_copy))
        subjet_drl2_copy=subjet_drl2_[:]
        subjet_drl2_copy.remove(min(subjet_drl2_copy))

      min_drl1=min(subjet_drl1_)
      min_drl2=min(subjet_drl2_)
      if event.nSubJet>1:
        submin_drl1=min(subjet_drl1_copy)
        submin_drl2=min(subjet_drl2_copy)

      min_index_subjet_drl1=subjet_drl1_.index(min_drl1)
      min_index_subjet_drl2=subjet_drl2_.index(min_drl2)
      if event.nSubJet>1:
        submin_index_subjet_drl1=subjet_drl1_.index(submin_drl1)
        submin_index_subjet_drl2=subjet_drl2_.index(submin_drl2)

      lep1subtracted_jetid_=-1
      lep2subtracted_jetid_=-1
      
      if event.nSubJet>1:
        if min_drl1<0.4:
          if min_drl2<0.4:
            # l1 and l2 overlap with the subjet
            if min_index_subjet_drl1==min_index_subjet_drl2:
              # l1 and l2 are closest to the same jet
              sub_tmptmp_v4=TLorentzVector()
              sub_tmptmp_v4=subjetv4_arr_[min_index_subjet_drl1]-l1v4raw_tmp-l2v4raw_tmp
              if sub_tmptmp_v4.M()<0:
                # the jet subtracting two leptons has negative mass, try other combination between subjets and leptons
                if (subjetv4_arr_[min_index_subjet_drl1]-l1v4raw_tmp).M()>0:
                  subjetv4_arr_[min_index_subjet_drl1]=subjetv4_arr_[min_index_subjet_drl1]-l1v4raw_tmp
                  lep1subtracted_jetid_=min_index_subjet_drl1
                  if submin_drl2<0.4 and (subjetv4_arr_[submin_index_subjet_drl2]-l2v4raw_tmp).M()>0:
                    subjetv4_arr_[submin_index_subjet_drl2]=subjetv4_arr_[submin_index_subjet_drl2]-l2v4raw_tmp
                    lep2subtracted_jetid_=submin_index_subjet_drl2

                elif (subjetv4_arr_[min_index_subjet_drl1]-l2v4raw_tmp).M()>0:
                  subjetv4_arr_[min_index_subjet_drl1]=subjetv4_arr_[min_index_subjet_drl1]-l2v4raw_tmp
                  lep2subtracted_jetid_=min_index_subjet_drl1
                  if submin_drl1<0.4 and (subjetv4_arr_[submin_index_subjet_drl1]-l1v4raw_tmp).M()>0:
                    subjetv4_arr_[submin_index_subjet_drl1]=subjetv4_arr_[submin_index_subjet_drl1]-l1v4raw_tmp
                    lep1subtracted_jetid_=submin_index_subjet_drl1

                else:
                  if submin_drl1<0.4 and (subjetv4_arr_[submin_index_subjet_drl1]-l1v4raw_tmp).M()>0:
                    subjetv4_arr_[submin_index_subjet_drl1]=subjetv4_arr_[submin_index_subjet_drl1]-l1v4raw_tmp
                    lep1subtracted_jetid_=submin_index_subjet_drl1
                  if submin_drl2<0.4 and (subjetv4_arr_[submin_index_subjet_drl2]-l2v4raw_tmp).M()>0:
                    subjetv4_arr_[submin_index_subjet_drl2]=subjetv4_arr_[submin_index_subjet_drl2]-l2v4raw_tmp
                    lep2subtracted_jetid_=submin_index_subjet_drl2

              else:
                subjetv4_arr_[min_index_subjet_drl1]=subjetv4_arr_[min_index_subjet_drl1]-l1v4raw_tmp-l2v4raw_tmp
                lep1subtracted_jetid_=min_index_subjet_drl1
                lep2subtracted_jetid_=min_index_subjet_drl1

            # l1 and l2 overlap with different jets
            else:
              if (subjetv4_arr_[min_index_subjet_drl1]-l1v4raw_tmp).M()>0:
                subjetv4_arr_[min_index_subjet_drl1]=subjetv4_arr_[min_index_subjet_drl1]-l1v4raw_tmp
                lep1subtracted_jetid_=min_index_subjet_drl1
              elif submin_drl1<0.4 and (subjetv4_arr_[submin_index_subjet_drl1]-l1v4raw_tmp).M()>0:
                subjetv4_arr_[submin_index_subjet_drl1]=subjetv4_arr_[submin_index_subjet_drl1]-l1v4raw_tmp
                lep1subtracted_jetid_=submin_index_subjet_drl1

              if (subjetv4_arr_[min_index_subjet_drl2]-l2v4raw_tmp).M()>0:
                subjetv4_arr_[min_index_subjet_drl2]=subjetv4_arr_[min_index_subjet_drl2]-l2v4raw_tmp
                lep2subtracted_jetid_=min_index_subjet_drl2
              elif submin_drl2<0.4 and (subjetv4_arr_[submin_index_subjet_drl2]-l2v4raw_tmp).M()>0:
                subjetv4_arr_[submin_index_subjet_drl2]=subjetv4_arr_[submin_index_subjet_drl2]-l2v4raw_tmp
                lep2subtracted_jetid_=submin_index_subjet_drl2

          else:
            # l2 doesn't overlap with jet
            if (subjetv4_arr_[min_index_subjet_drl1]-l1v4raw_tmp).M()<0:
              if submin_drl1<0.4 and (subjetv4_arr_[submin_index_subjet_drl1]-l1v4raw_tmp).M()>0:
                subjetv4_arr_[submin_index_subjet_drl1]=subjetv4_arr_[submin_index_subjet_drl1]-l1v4raw_tmp
                lep1subtracted_jetid_=submin_index_subjet_drl1
            else:
              subjetv4_arr_[min_index_subjet_drl1]=subjetv4_arr_[min_index_subjet_drl1]-l1v4raw_tmp
              lep1subtracted_jetid_=min_index_subjet_drl1
        else:
          if min_drl2<0.4:
            if (subjetv4_arr_[min_index_subjet_drl2]-l2v4raw_tmp).M()<0:
              if submin_drl2<0.4 and (subjetv4_arr_[submin_index_subjet_drl2]-l2v4raw_tmp).M()>0:
                subjetv4_arr_[submin_index_subjet_drl2]=subjetv4_arr_[submin_index_subjet_drl2]-l2v4raw_tmp
                lep2subtracted_jetid_=submin_index_subjet_drl2
            else:
              subjetv4_arr_[min_index_subjet_drl2]=subjetv4_arr_[min_index_subjet_drl2]-l2v4raw_tmp
              lep2subtracted_jetid_=min_index_subjet_drl2

      # only one subjet
      else:
        if min_drl1<0.4:
          if min_drl2<0.4:
            sub_tmptmp_v4=TLorentzVector()
            sub_tmptmp_v4=subjetv4_arr_[0]-l1v4raw_tmp-l2v4raw_tmp
            if sub_tmptmp_v4.M()<0:
              if (subjetv4_arr_[0]-l1v4raw_tmp).M()>0:
                subjetv4_arr_[0]=subjetv4_arr_[0]-l1v4raw_tmp
                lep1subtracted_jetid_=0
              elif (subjetv4_arr_[0]-l2v4raw_tmp).M()>0:
                subjetv4_arr_[0]=subjetv4_arr_[0]-l2v4raw_tmp
                lep2subtracted_jetid_=0
            else:
              subjetv4_arr_[0]=subjetv4_arr_[0]-l1v4raw_tmp-l2v4raw_tmp
              lep1subtracted_jetid_=0
              lep2subtracted_jetid_=0

          else:
            # l2 doesn't overlap with jet
            if (subjetv4_arr_[0]-l1v4raw_tmp).M()>0:
              subjetv4_arr_[0]=subjetv4_arr_[0]-l1v4raw_tmp
              lep1subtracted_jetid_=0
        else:
          if min_drl2<0.4:
            # l1 doesn't overlap with jet
            if (subjetv4_arr_[0]-l2v4raw_tmp).M()>0:
              subjetv4_arr_[0]=subjetv4_arr_[0]-l2v4raw_tmp


      for ij in range(0,event.nSubJet):
        SubJetNoVlep_eta.append(subjetv4_arr_[ij].Eta())
        SubJetNoVlep_phi.append(subjetv4_arr_[ij].Phi())
        SubJetNoVlep_dRinit.append(subjetv4_arr_[ij].DeltaR(subjetv4_arr2_[ij]))
        SubJetNoVlep_id.append(ij)
        if not ( subjetv4_arr_[ij].Pt()== subjetv4_arr2_[ij].Pt()):
          SubJetNoVlep_rawFactor.append(0.0)
          SubJetNoVlep_pt.append(subjetv4_arr_[ij].Pt())
          SubJetNoVlep_mass.append(subjetv4_arr_[ij].M())
        else:
          SubJetNoVlep_rawFactor.append(recosubjets[ij].rawFactor)
          SubJetNoVlep_pt.append(recosubjets[ij].pt)
          SubJetNoVlep_mass.append(recosubjets[ij].mass)

    # subjet no need to sort by pt, the default ID of subjet is needed to match the subjet id of Fatjet

    self.out.fillBranch("JetNoVlep_pt", JetNoVlep_pt)
    self.out.fillBranch("JetNoVlep_eta", JetNoVlep_eta)
    self.out.fillBranch("JetNoVlep_phi", JetNoVlep_phi)
    self.out.fillBranch("JetNoVlep_mass", JetNoVlep_mass)
    self.out.fillBranch("JetNoVlep_rawFactor", JetNoVlep_rawFactor)
    self.out.fillBranch("JetNoVlep_area", JetNoVlep_area)
    self.out.fillBranch("JetNoVlep_id", JetNoVlep_id)
    self.out.fillBranch("JetNoVlep_jetId", JetNoVlep_jetId)
    self.out.fillBranch("JetNoVlep_muonIdx1", JetNoVlep_muonIdx1)
    self.out.fillBranch("JetNoVlep_muonIdx2", JetNoVlep_muonIdx2)
    self.out.fillBranch("JetNoVlep_muonSubtrFactor", JetNoVlep_muonSubtrFactor)
    self.out.fillBranch("JetNoVlep_neEmEF", JetNoVlep_neEmEF)
    self.out.fillBranch("JetNoVlep_chEmEF", JetNoVlep_chEmEF)
    self.out.fillBranch("JetNoVlep_drl1", JetNoVlep_drl1)
    self.out.fillBranch("JetNoVlep_drl2", JetNoVlep_drl2)
    self.out.fillBranch("JetNoVlep_dRinit", JetNoVlep_dRinit)
    self.out.fillBranch("SubJetNoVlep_pt", SubJetNoVlep_pt)
    self.out.fillBranch("SubJetNoVlep_eta", SubJetNoVlep_eta)
    self.out.fillBranch("SubJetNoVlep_phi", SubJetNoVlep_phi)
    self.out.fillBranch("SubJetNoVlep_mass", SubJetNoVlep_mass)
    self.out.fillBranch("SubJetNoVlep_rawFactor", SubJetNoVlep_rawFactor)
    self.out.fillBranch("SubJetNoVlep_dRinit", SubJetNoVlep_dRinit)
    self.out.fillBranch("SubJetNoVlep_id", SubJetNoVlep_id)
    self.out.fillBranch("FatJetNoVlep_pt", FatJetNoVlep_pt)
    self.out.fillBranch("FatJetNoVlep_eta", FatJetNoVlep_eta)
    self.out.fillBranch("FatJetNoVlep_phi", FatJetNoVlep_phi)
    self.out.fillBranch("FatJetNoVlep_mass", FatJetNoVlep_mass)
    self.out.fillBranch("FatJetNoVlep_PNmass", FatJetNoVlep_PNmass)
    self.out.fillBranch("FatJetNoVlep_SDmass", FatJetNoVlep_SDmass)
    self.out.fillBranch("FatJetNoVlep_rawFactor", FatJetNoVlep_rawFactor)
    self.out.fillBranch("FatJetNoVlep_area", FatJetNoVlep_area)
    self.out.fillBranch("FatJetNoVlep_id", FatJetNoVlep_id)
    self.out.fillBranch("FatJetNoVlep_jetId", FatJetNoVlep_jetId)
    self.out.fillBranch("FatJetNoVlep_overlapl1", FatJetNoVlep_overlapl1)
    self.out.fillBranch("FatJetNoVlep_overlapl2", FatJetNoVlep_overlapl2)
    self.out.fillBranch("FatJetNoVlep_subJetIdx1", FatJetNoVlep_subJetIdx1)
    self.out.fillBranch("FatJetNoVlep_subJetIdx2", FatJetNoVlep_subJetIdx2)
    self.out.fillBranch("FatJetNoVlep_drl1", FatJetNoVlep_drl1)
    self.out.fillBranch("FatJetNoVlep_drl2", FatJetNoVlep_drl2)
    self.out.fillBranch("FatJetNoVlep_dRinit", FatJetNoVlep_dRinit)
    self.out.fillBranch("FatJetNoVlep_bbvsccqq", FatJetNoVlep_bbvsccqq)
    self.out.fillBranch("mu_channel", mu_channel)
    self.out.fillBranch("l1_pt", l1_pt)
    self.out.fillBranch("l1_rawpt", l1_rawpt)
    self.out.fillBranch("l1_eta", l1_eta)
    self.out.fillBranch("l1_phi", l1_phi)
    self.out.fillBranch("l1_mass", l1_mass)
    self.out.fillBranch("l1_rawmass", l1_rawmass)
    self.out.fillBranch("l1_relPtl2", l1_relPtl2)
    self.out.fillBranch("l2_pt", l2_pt)
    self.out.fillBranch("l2_rawpt", l2_rawpt)
    self.out.fillBranch("l2_eta", l2_eta)
    self.out.fillBranch("l2_phi", l2_phi)
    self.out.fillBranch("l2_mass", l2_mass)
    self.out.fillBranch("l2_rawmass", l2_rawmass)
    self.out.fillBranch("l2_relPtl1", l2_relPtl1)
    self.out.fillBranch("drll", drll)
    self.out.fillBranch("detall", detall)
    self.out.fillBranch("dphill", dphill)
    self.out.fillBranch("zlep_pt", zlep_pt)
    self.out.fillBranch("zlep_eta", zlep_eta)
    self.out.fillBranch("zlep_phi", zlep_phi)
    self.out.fillBranch("zlep_mass", zlep_mass)

    return True

ZandLepSubProducer = lambda: ZandLepSubSelector()

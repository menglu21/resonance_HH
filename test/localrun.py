import os
import sys
import optparse
import ROOT
import re

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.common.countHistogramsModule import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.eleRECOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.eleIDSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.muonScaleResProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.muonIDISOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.HHProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.ZandLepSubSelector import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.ShuffleSubGenJet import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.ElectronSelector import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.MuonSelector import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.GenJetNoVlepMerger import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.PrefireCorr import *
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles, runsAndLumis
### main python file to run ###

def main():

  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('--year', dest='year', help='which year sample', default='2018', type='string')
  parser.add_option('-m', dest='ismc', help='to apply sf correction or not', default=False, action='store_true')
  parser.add_option('-n','--nEve', dest='nEvent', help='number of event', type='int', action='store')
  parser.add_option('-i', '--in', dest='inputs', help='input directory with files', default=None, type='string')
  parser.add_option('-d', dest='ismc', help='to apply sf correction or not', action='store_false')
  parser.add_option('-o', '--out', dest='output', help='output directory with files', default=None, type='string')
  (opt, args) = parser.parse_args()

  if opt.ismc:
    if opt.year == "2016a":
      p = PostProcessor(opt.output, [opt.inputs], modules=[countHistogramsModule(),puAutoWeight_2016_preAPV(),GenJetAK4NoVlepMerger(),GenJetAK8NoVlepMerger(),PrefCorr2016(),muonIDISOSF2016(),muonScaleRes2016a(),eleRECOSF2016(),eleIDSF2016(), jmeCorrections_UL2016APVMC(), btagSF2016ULapv(), HH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2016b":
      p = PostProcessor(opt.output, [opt.inputs], modules=[countHistogramsModule(),puAutoWeight_2016_postAPV(),GenJetAK4NoVlepMerger(),GenJetAK8NoVlepMerger(),PrefCorr2016(),muonIDISOSF2016(),muonScaleRes2016b(),eleRECOSF2016(),eleIDSF2016(), jmeCorrections_UL2016MC(), btagSF2016UL(), HH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2017":
      p = PostProcessor(opt.output, [opt.inputs], modules=[countHistogramsModule(),puWeight_2017(),GenJetAK4NoVlepMerger(),GenJetAK8NoVlepMerger(),SubGenJetAK8NoVlepMerger(),ShuffleSubGenJet(),muonIDISOSF2017(),muonScaleRes2017(),eleRECOSF2017(),eleIDSF2017(),MuonProducer(),ElectronProducer(),ZandLepSubProducer(),jmeCorrections_UL2017MC(),fatjmeCorrections_UL2017MC(),btagSF2017UL(),HH2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
      #p = PostProcessor(opt.output, [opt.inputs], modules=[countHistogramsModule(),puWeight_2017(),GenJetAK4NoVlepMerger(),GenJetAK8NoVlepMerger(),PrefCorr(),muonIDISOSF2017(),muonScaleRes2017(),eleRECOSF2017(),eleIDSF2017(), jmeCorrections_UL2017MC(),btagSF2017UL(), HH2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2018":
      p = PostProcessor(opt.output, [opt.inputs], modules=[countHistogramsModule(),puWeight_2018(),muonIDISOSF2018(),muonScaleRes2018(),eleRECOSF2018(),eleIDSF2018(),jmeCorrections_UL2018MC(), btagSF2018UL(),HH2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)


# Sequence for data
  if not (opt.ismc):
    if opt.year == "2016b":
      print "run on 2016B data"
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016a(),jmeCorrections_UL2016B(),HH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt", maxEntries=opt.nEvent)
    if opt.year == "2016c":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016a(),jmeCorrections_UL2016C(),HH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt", maxEntries=opt.nEvent)  
    if opt.year == "2016d":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016a(),jmeCorrections_UL2016D(),HH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt", maxEntries=opt.nEvent)
    if opt.year == "2016e":
      print "run on 2016E data"
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016a(),jmeCorrections_UL2016E(),HH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt", maxEntries=opt.nEvent)
    if opt.year == "2016f_apv":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016a(),jmeCorrections_UL2016APVF(),HH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt", maxEntries=opt.nEvent)
    if opt.year == "2016f":
      print "run on 2016F data"
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016b(),jmeCorrections_UL2016F(),HH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt", maxEntries=opt.nEvent)
    if opt.year == "2016g":
      print "run on 2016G data"
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016b(),jmeCorrections_UL2016G(),HH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt", maxEntries=opt.nEvent)
    if opt.year == "2016h":
      print "run on 2016H data"
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016b(),jmeCorrections_UL2016H(),HH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt", maxEntries=opt.nEvent)
    if opt.year == "2017b":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2017(),jmeCorrections_UL2017B(),HH2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2017c":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2017(),jmeCorrections_UL2017C(),HH2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2017d":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2017(),jmeCorrections_UL2017D(),HH2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2017e":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2017(),jmeCorrections_UL2017E(),HH2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2017f":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2017(),jmeCorrections_UL2017F(),HH2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2018a":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2018(),jmeCorrections_UL2018A(),HH2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2018b":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2018(),jmeCorrections_UL2018B(),HH2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2018c":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2018(),jmeCorrections_UL2018C(),HH2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2018d":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2018(),jmeCorrections_UL2018D(),HH2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
  p.run()

if __name__ == "__main__":
    sys.exit(main())

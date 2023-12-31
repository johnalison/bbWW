import sys
import optparse
#from pyxrootd import client
import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
#sys.path.insert(0, 'ZZ4b/nTupleAnalysis/scripts/')
#from cfgHelper import *
sys.path.insert(0, 'nTupleAnalysis/python/') #https://github.com/patrickbryant/nTupleAnalysis
from commandLineHelpers import *



def xstr(s):
    if s is None:
        return ''
    return str(s)

print( "Input command")
print( " ".join(sys.argv))

parser = optparse.OptionParser()
parser.add_option('-d', '--debug',                dest="debug",         action="store_true", default=False, help="debug")
parser.add_option('-m', '--isMC',                 dest="isMC",          action="store_true", default=False, help="isMC")
parser.add_option('-y', '--year',                 dest="year",          default="2016", help="Year specifies trigger (and lumiMask for data)")
parser.add_option(      '--firstEvent',           default=0, help="First event in the data set to proccess")
parser.add_option('-l', '--lumi', type="float",   dest="lumi",          default=1.0,    help="Luminosity for MC normalization: units [pb]")
#parser.add_option(      '--bTagger',              dest="bTagger",       default="CSVv2", help="bTagging algorithm")
#parser.add_option('-b', '--bTag',                 dest="bTag",          default="0.8484", help="bTag cut value: default is medium WP for CSVv2 https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco")
parser.add_option(      '--bTagger',              dest="bTagger",       default="deepFlavB", help="bTagging algorithm")
parser.add_option('-b', '--bTag', type="float",   dest="bTag",          default="0.2770", help="bTag cut value: default is medium WP for deepFlavB (DeepJet?) https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X")
parser.add_option(      '--bTagSF',               dest="bTagSF",        action="store_true", default=False, help="Load btagging SFs")
parser.add_option(      '--bTagSyst',             dest="bTagSyst",      action="store_true", default=False, help="run btagging systematics")
parser.add_option(      '--JECSyst',              dest="JECSyst",       default="", help="Name of JEC Systematic uncertainty, examples: _jerDown, _jesTotalUp")
parser.add_option('-i', '--input',                dest="input",         default="ZZ4b/fileLists/data2016H.txt", help="Input file(s). If it ends in .txt, will treat it as a list of input files.")
parser.add_option(      '--friends',              dest="friends",       default=None, help="Extra friend files. comma separated list where each item replaces picoAOD in the input file, ie FvT,SvB for FvT.root stored in same location as picoAOD.root")
parser.add_option('-o', '--outputBase',           dest="outputBase",    default="/uscms/home/bryantp/nobackup/ZZ4b/", help="Base path for storing output histograms and picoAOD")
parser.add_option('-p', '--createPicoAOD',        dest="createPicoAOD", type="string", help="Create picoAOD with given name. Use NONE (case does not matter) to explicitly avoid creating any picoAOD")
parser.add_option('-f', '--fastSkim',             dest="fastSkim",      action="store_true", default=False, help="Do minimal computation to maximize event loop rate for picoAOD production")
parser.add_option('-n', '--nevents',              dest="nevents",       default="-1", help="Number of events to process. Default -1 for no limit.")
parser.add_option(      '--histDetailLevel',        dest="histDetailLevel", default="allEvents.HLTStudy.preSel.bTags", help="Hist Detail level.")
parser.add_option(   '--unBlind',    default=False, action="store_true",help="Flag to not blind (Needed for Mixed samples and eventually in real Data!)")
parser.add_option(      '--histFile',             dest="histFile",      default="hists.root", help="name of ouptut histogram file")
parser.add_option(   '--writeOutEventNumbers',      action="store_true", default=False, help="boolean  to toggle writing out event/run numbers")
#parser.add_option('-r', '--reweight',             dest="reweight",      default="", help="Reweight file containing TSpline3 of nTagClassifier ratio")
parser.add_option(      '--otherWeights',     type="string", default=None, help="other event weights to apply")
parser.add_option(   '--condor',   action="store_true", default=False,           help="currenty does nothing. Try to keep it that way")
o, a = parser.parse_args()


bjetSF = 'deepjet'+o.year
year = o.year.replace('_preVFP','').replace('_postVFP','')
if not o.isMC or not o.bTagSF:
    bjetSF = ''
btagVariations = 'central'
if 'jes' in o.JECSyst:
    if 'Down' in o.JECSyst:
        btagVariations = 'down'+o.JECSyst.replace('Down','')
    if 'Up' in o.JECSyst:
        btagVariations = 'up'+o.JECSyst.replace('Up','')
if o.bTagSyst:
    btagVariations += ' down_hfstats1 up_hfstats1'
    btagVariations += ' down_hfstats2 up_hfstats2'
    btagVariations += ' down_lfstats1 up_lfstats1'
    btagVariations += ' down_lfstats2 up_lfstats2'
    btagVariations += ' down_hf up_hf'
    btagVariations += ' down_lf up_lf'
    btagVariations += ' down_cferr1 up_cferr1'
    btagVariations += ' down_cferr2 up_cferr2'

#
# Basic Configuration
#
outputBase = o.outputBase + ('/' if o.outputBase[-1] != '/' else '') # make sure it ends with a slash
isData     = not o.isMC
blind      = True and isData and not o.unBlind
#https://cms-service-dqmdc.web.cern.ch/CAF/certification/
JSONfiles  = {'2015':'',
              '2016':        'bbWW/lumiMasks/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt', #Ultra Legacy
              '2016_preVFP': 'bbWW/lumiMasks/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt', #Ultra Legacy
              '2016_postVFP':'bbWW/lumiMasks/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt', #Ultra Legacy
              '2017':        'bbWW/lumiMasks/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt', #Ultra Legacy
              '2018':        'bbWW/lumiMasks/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'} #Ultra Legacy
# Calculated lumi per lumiBlock from brilcalc. See README
lumiData   = {'2015':'',
              '2016':        'bbWW/lumiMasks/brilcalc_2016_HLT_QuadJet45_TripleBTagCSV_p087.csv', 
              '2016_preVFP': 'bbWW/lumiMasks/brilcalc_2016_HLT_QuadJet45_TripleBTagCSV_p087.csv', 
              '2016_postVFP':'bbWW/lumiMasks/brilcalc_2016_HLT_QuadJet45_TripleBTagCSV_p087.csv', 
              '2017':        'bbWW/lumiMasks/brilcalc_2017_HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0.csv',
              '2018':        'bbWW/lumiMasks/brilcalc_2018_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5.csv'} 

# for MC we need to normalize the sample to the recommended cross section * BR times the target luminosity
## Higgs BRs https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR BR(h125->bb) = 0.5824 BR(h125->\tau\tau) = 0.06272 BR(Z->bb) = 0.1512, BR(Z->\tau\tau) = 0.03696
## ZH cross sections https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV#ZH_Process
## ZZ cross section 15.0 +0.7 -0.6 +/-0.2 (MCFM at NLO in QCD with additional contributions from LO gg -> ZZ diagrams) or 16.2 +0.6 -0.4 (calculated at NNLO in QCD via MATRIX) https://arxiv.org/pdf/1607.08834.pdf pg 10
## ZH->bb\tau\tau xs = (0.7612+0.1227)*(0.58*0.036+0.15*0.067) = 27/fb ~ 10x HH cross section
## HH->bb\tau\tau xs = 34*0.58*0.067*2 = 2.6/fb
## Higgs BR(mH=125.0) = 0.5824, BR(mH=125.09) = 0.5809: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
## Z BR = 0.1512+/-0.0005 from PDG
## store all process cross sections in pb. Can compute xs of sample with GenXsecAnalyzer. Example: 
## https://twiki.cern.ch/twiki/bin/viewauth/CMS/HowToGenXSecAnalyzer
## cd genproductions/Utilities/calculateXSectionAndFilterEfficiency; ./calculateXSectionAndFilterEfficiency.sh -f ../../../ZZ_dataset.txt -c RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1 -d MINIAODSIM -n -1 
## tt xs NNLO and measurement in dilep and semilep tt+jets, tt+bb: https://cds.cern.ch/record/2684606/files/TOP-18-002-paper-v19.pdf
xsDictionary = {'ggZH4b':  0.1227*0.5824*0.1512, #0.0432 from GenXsecAnalyzer, does not include BR for H, does include BR(Z->hadrons) = 0.69911. 0.0432/0.69911 = 0.0618, almost exactly half the LHCXSWG value... NNLO = 2x NLO??
                  'ZH4b':  0.7612*0.5824*0.1512, #0.5540 from GenXsecAnalyzer, does not include BR for H, does include BR(Z->hadrons) = 0.69911. 0.5540/0.69911 = 0.7924, 4% larger than the LHCXSWG value.
              'bothZH4b': (0.1227+0.7612)*0.5824*0.1512,
                  'ZZ4b': 15.5   *0.1512*0.1512,#0.3688 from GenXsecAnalyzer gives 16.13 dividing by BR^2. mcEventSumw/mcEventCount * FxFx Jet Matching eff. = 542638/951791 * 0.647 = 0.3688696216. Jet matching not included in genWeight!
                  'HH4b': 0.03105*0.5824**2, # (0.0457 2018, doesn't include BR, 0.009788 2016, does include BR...) https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWGHH recommends 31.05fb*BR^2=10.53fb
                'TTJets': 831.76, #749.5 get xs from GenXsecAnalyzer, McM is just wrong... TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8. Apply 4b scale k-factor 5.5/3.6=1.53 https://cds.cern.ch/record/2687373/files/TOP-18-011-paper-v15.pdf
                'TTToHadronic': 377.9607353256, #313.9 from McM. NNLO tt xs = 831.76, W hadronic BR = 0.6741 => NNLO = 831.76*0.6741^2 = 377.9607353256
                'TTToSemiLeptonic': 365.7826460496, #300.9 from McM. NNLO = 831.76*2*(1-0.6741)*0.6747 = 365.7826460496
                'TTTo2L2Nu': 88.3419033256, #72.1 from McM. NNLO = 831.76*(1-0.6741)^2 = 88.3419033256
                'WHHTo4B_CV_0_5_C2V_1_0_C3_1_0':2.870e-04*0.5824*0.5824,  # 2.870e-04from GenXsecAnalyzer, does not include BR for H 
                'WHHTo4B_CV_1_0_C2V_0_0_C3_1_0':1.491e-04*0.5824*0.5824,  # 1.491e-04from GenXsecAnalyzer, does not include BR for H 
                'WHHTo4B_CV_1_0_C2V_1_0_C3_0_0':2.371e-04*0.5824*0.5824,  # 2.371e-04from GenXsecAnalyzer, does not include BR for H 
                'WHHTo4B_CV_1_0_C2V_1_0_C3_1_0':4.152e-04*0.5824*0.5824,  # 4.152e-04from GenXsecAnalyzer, does not include BR for H 
                'WHHTo4B_CV_1_0_C2V_1_0_C3_2_0':6.880e-04*0.5824*0.5824,  # 6.880e-04from GenXsecAnalyzer, does not include BR for H 
                'WHHTo4B_CV_1_0_C2V_2_0_C3_1_0':1.115e-03*0.5824*0.5824,  # 1.115e-03from GenXsecAnalyzer, does not include BR for H 
                'WHHTo4B_CV_1_5_C2V_1_0_C3_1_0':8.902e-04*0.5824*0.5824,  # 8.902e-04from GenXsecAnalyzer, does not include BR for H 
                'WHHTo4B_CV_1_0_C2V_1_0_C3_20_0':2.158e-02*0.5824*0.5824, # 2.158e-02from GenXsecAnalyzer, does not include BR for H 
                'ZHHTo4B_CV_0_5_C2V_1_0_C3_1_0':1.663e-04*0.5824*0.5824,  # 1.663e-04from GenXsecAnalyzer, does not include BR for H 
                'ZHHTo4B_CV_1_0_C2V_0_0_C3_1_0':9.037e-05*0.5824*0.5824,  # 9.037e-05from GenXsecAnalyzer, does not include BR for H 
                'ZHHTo4B_CV_1_0_C2V_1_0_C3_0_0':1.544e-04*0.5824*0.5824,  # 1.544e-04from GenXsecAnalyzer, does not include BR for H 
                'ZHHTo4B_CV_1_0_C2V_1_0_C3_1_0':2.642e-04*0.5824*0.5824,  # 2.642e-04from GenXsecAnalyzer, does not include BR for H 
                'ZHHTo4B_CV_1_0_C2V_1_0_C3_2_0':4.255e-04*0.5824*0.5824,  # 4.255e-04from GenXsecAnalyzer, does not include BR for H 
                'ZHHTo4B_CV_1_0_C2V_2_0_C3_1_0':6.770e-04*0.5824*0.5824,  # 6.770e-04from GenXsecAnalyzer, does not include BR for H 
                'ZHHTo4B_CV_1_5_C2V_1_0_C3_1_0':5.738e-04*0.5824*0.5824,  # 5.738e-04from GenXsecAnalyzer, does not include BR for H 
                'ZHHTo4B_CV_1_0_C2V_1_0_C3_20_0':1.229e-02*0.5824*0.5824, # 1.229e-02from GenXsecAnalyzer, does not include BR for H 
                } 

## figure out what sample is being run from the name of the input
sample = ''
for key in xsDictionary.keys():
    if key in o.input:
        sample = key
        break

if sample == '':
    if 'TTJets' in o.input: sample = 'TTJets'
    elif 'TTToHadronic' in o.input: sample = 'TTToHadronic'
    elif 'TTToSemiLeptonic' in o.input: sample = 'TTToSemiLeptonic'
    elif 'TTTo2L2Nu' in o.input: sample = 'TTTo2L2Nu'
    elif 'ggZH' in o.input: sample = 'ggZH4b'
    elif 'bothZH' in o.input: sample = 'bothZH4b'
    elif 'ZH' in o.input: sample =   'ZH4b'
    elif 'HH' in o.input: sample =   'HH4b'
    elif 'ZZ' in o.input: sample =   'ZZ4b' #make sure this is last, ZZ in path name...
xs = 1
if o.isMC: 
    xs = xsDictionary[sample] if sample in xsDictionary else 1.0
    print( 'Simulated sample:',sample,'| xs =',xs)


fileNames = []
inputList=False
if '.txt' in o.input:
    inputList = True
    for line in open(o.input, 'r').readlines():
        line = line.replace('\n','').strip()
        if line    == '' : continue
        if line[0] == '#': continue
        fileNames.append(line.replace('\n',''))
else:
    fileNames.append(o.input)




useOtherPicoAOD = True if 'picoAOD' in o.input else False

pathOut = outputBase
if 'root://cmsxrootd-site.fnal.gov//store/' in pathOut: 
    pathOut = pathOut + fileNames[0].replace('root://cmsxrootd-site.fnal.gov//store/', '') #make it a local path
if useOtherPicoAOD:
    pathOut = o.input
pathOut = '/'.join(pathOut.split('/')[:-1])+'/' #remove <fileName>.root
    
if inputList: #use simplified directory structure based on grouping of filelists
    sampleDirectory = o.input.split('/')[-1].replace('.txt','/')
    pathOut = outputBase+sampleDirectory

mkpath(pathOut)

histOut = pathOut+o.histFile


#
# Automatic picoAOD logic
#
defaultPicoAOD = 'picoAOD.root'

#check if the defaultPicoAOD already exists, if it doesn't we probably want to make one.
defaultPicoAODExists = exists(pathOut + defaultPicoAOD)

#check if we are explicitly being asked by the user to make the default picoAOD
createDefaultPicoAOD = o.createPicoAOD == defaultPicoAOD

# Use the default picoAOD if it exists and we aren't explicitly being asked to make it or use a different picoAOD
useDefaultPicoAOD = defaultPicoAODExists and not createDefaultPicoAOD and not useOtherPicoAOD and xstr(o.createPicoAOD).lower() != 'none'
if useDefaultPicoAOD: fileNames = [pathOut+defaultPicoAOD]

# construct the picoAOD file path
picoAOD = pathOut+(o.createPicoAOD if o.createPicoAOD else 'picoAOD.root')

# create this picoAOD if the user specified one or if the default picoAOD.root does not exist and not explicitly overriding the auto-create
create = bool(o.createPicoAOD or not defaultPicoAODExists) and xstr(o.createPicoAOD).lower() != 'none'

# make sure the input and output files are not the same
print('create picoAOD:',create,picoAOD)
if fileNames[0] == picoAOD and create:
    print( 'ERROR: Trying to overwrite input picoAOD:',picoAOD)
    sys.exit()


friendFiles = []
if o.friends:
    if '.txt' in o.friends:
        for line in open(o.friends, 'r').readlines():
            line = line.replace('\n','').strip()
            if line    == '' : continue
            if line[0] == '#': continue
            friendFiles.append(line.replace('\n',''))
    else:
        friends = o.friends.split(',')
        for friend in friends:
            friendFileName = pathOut+friend+'.root'
            friendFiles.append(friendFileName)
            print('Friend:',friendFileName)
    



#
#  Logic to prepare the other weights
#
otherWeights = []
if o.otherWeights:
    otherWeights = o.otherWeights.split(',')


#
# ParameterSets for use in bin/<script>.cc 
#
process = cms.PSet()

#Setup framework lite input file object
process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(fileNames),
    maxEvents   = cms.int32(int(float(o.nevents))),
    )

# LumiMask
process.inputs = cms.PSet(
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    )
if isData:
    # get JSON file correctly parced
    myList = LumiList.LumiList(filename = JSONfiles[year]).getCMSSWString().split(',')
    process.inputs.lumisToProcess.extend(myList)

# Setup picoAOD
process.picoAOD = cms.PSet(
    fileName = cms.string(picoAOD),
    create   = cms.bool(create),
    fastSkim = cms.bool(o.fastSkim),
    )




# Setup framwork lite output file object
process.fwliteOutput = cms.PSet(
    fileName  = cms.string(histOut),
    )


#Setup event loop object
process.jpsiAnalysis = cms.PSet(
    debug   = cms.bool(o.debug),
    isMC    = cms.bool(o.isMC),
    blind   = cms.bool(blind),
    year    = cms.string(year),
    lumi    = cms.double(o.lumi),
    firstEvent  = cms.int32(int(o.firstEvent)),
    xs      = cms.double(xs),
    bTag    = cms.double(o.bTag),
    bTagger = cms.string(o.bTagger),
    bjetSF  = cms.string(bjetSF),
    btagVariations = cms.string(btagVariations),
    JECSyst = cms.string(o.JECSyst),
    friendFile = cms.string(fileNames[0].replace('.root','_Friend.root')),
    lumiData= cms.string(lumiData[year]),
    histDetailLevel = cms.string(o.histDetailLevel),
    writeOutEventNumbers  = cms.bool(o.writeOutEventNumbers),
    #friends          = cms.vstring(friendFiles),
    )

print('jpsiAnalysis_cfg.py done')

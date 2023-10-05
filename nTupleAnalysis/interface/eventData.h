// -*- C++ -*-

#if !defined(bbWW_eventData_H)
#define bbWW_eventData_H

#include <iostream>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include "nTupleAnalysis/baseClasses/interface/initBranch.h"
#include "nTupleAnalysis/baseClasses/interface/truthData.h"
#include "nTupleAnalysis/baseClasses/interface/jetData.h"
#include "nTupleAnalysis/baseClasses/interface/trackJetData.h"
#include "nTupleAnalysis/baseClasses/interface/muonData.h"
#include "nTupleAnalysis/baseClasses/interface/elecData.h"
#include "nTupleAnalysis/baseClasses/interface/eventData.h"
#include "nTupleAnalysis/baseClasses/interface/dimuon.h"
#include "nTupleAnalysis/baseClasses/interface/dijet.h"
#include "nTupleAnalysis/baseClasses/interface/trigData.h"

namespace bbWW {

  class eventData {

  public:

    // Member variables
    TChain* tree;
    bool isMC;
    float year;
    bool debug;
    bool printCurrentFile = true;
    bool fastSkim = false;
    
    nTupleAnalysis::eventData* eventDataTree = NULL;

    Float_t   bTagSF = 1;
    int       nTrueBJets = 0;
    float weight = 1.0;

    nTupleAnalysis::truthData* truth = NULL;
    nTupleAnalysis::truthParticle* truthJets = NULL;
    std::vector<nTupleAnalysis::particlePtr> genJets;//jets passing pt/eta and bTagging requirements

    float mbbWW;

    //Predefine btag sorting functions
    float       bTag    = 0.4168; //0.6;
    std::string bTagger = "deepFlavB";

    //triggers
    std::map<std::string, bool> L1_triggers;
    std::map<std::string, bool> HLT_triggers;
    std::map<std::string, std::map<std::string, bool*>> HLT_L1_seeds;
    bool passL1              = false;
    bool passHLT             = false;
    bool wStarlnu            = false;
    bool lepton25            = true;
    bool qjet25              = true;
    bool qjetEta             = true;
    bool qjetLead30          = true;

    std::map<std::string, bool*> L1_triggers_mon;

  public:

    float jetPtMin = 0;
    const float jetEtaMax= 2.4;
    const int puIdMin = 0b110;//7=tight, 6=medium, 4=loose working point
    const bool doJetCleaning=true;
     
    nTupleAnalysis::jetData* treeJets;
    std::vector<nTupleAnalysis::jetPtr> allJets;//all jets in nTuple
    std::vector<nTupleAnalysis::jetPtr> selJets;//jets passing pt/eta requirements
    std::vector<nTupleAnalysis::jetPtr> selJets30;//jets passing pt > 30/eta requirements
    std::vector<nTupleAnalysis::jetPtr> btagJets;//jets passing pt/eta and bTagging requirements

    //  
    // For the TTBar TandP calibration
    //
    std::vector< nTupleAnalysis::jetPtr > WqTagJets;
    std::vector< nTupleAnalysis::jetPtr > WqProbeJets;

    //std::vector<nTupleAnalysis::jetPtr> ctagJets;//jets passing pt/eta and cTagging requirements

    nTupleAnalysis::trackJetData* treeTrackJets;
    std::vector<nTupleAnalysis::trackJetPtr> allTrackJets;//all track jets in nTuple

    nTupleAnalysis::muonData* treeMuons;
    std::vector<nTupleAnalysis::muonPtr> allMuons;
    std::vector<nTupleAnalysis::muonPtr> preSelMuons;

    std::vector<nTupleAnalysis::dimuonPtr> preSelDiMuons;
    bool diMuonIsolated = true;

    std::shared_ptr<nTupleAnalysis::dijet> diBJets;
    std::shared_ptr<nTupleAnalysis::dijet> diCJets;



    nTupleAnalysis::elecData* treeElecs;
    std::vector<nTupleAnalysis::elecPtr> allElecs;


    nTupleAnalysis::trigData* treeTrig = NULL;

    // Constructors and member functions
    eventData(TChain* t, bool mc, std::string y, bool d, bool _fastSkim = false, std::string bjetSF = "", std::string btagVariations = "central",
	      std::string JECSyst = "");
        
    void update(long int);
    void buildEvent();
    void resetEvent();
    void doTTbarTandPSelection();

    void dump();
    ~eventData(); 

    std::string currentFile = "";


  };

}
#endif // bbWW_eventData_H

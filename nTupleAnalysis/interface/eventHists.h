// -*- C++ -*-
#if !defined(eventHists_H)
#define eventHists_H

#include <iostream>
#include <TH1F.h>
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "bbWW/nTupleAnalysis/interface/eventData.h"
#include "nTupleAnalysis/baseClasses/interface/fourVectorHists.h"
#include "nTupleAnalysis/baseClasses/interface/truthParticleHists.h"
#include "nTupleAnalysis/baseClasses/interface/jetHists.h"
#include "nTupleAnalysis/baseClasses/interface/muonHists.h"
#include "nTupleAnalysis/baseClasses/interface/dimuonHists.h"
#include "nTupleAnalysis/baseClasses/interface/dijetHists.h"
#include "nTupleAnalysis/baseClasses/interface/elecHists.h"


namespace bbWW {

  class eventHists {
  public:

    TFileDirectory dir;
    bool debug;

    // Truth
    nTupleAnalysis::fourVectorHists* vHH;

    nTupleAnalysis::truthParticleHists* vHbb;
    nTupleAnalysis::truthParticleHists* vbLead;
    nTupleAnalysis::truthParticleHists* vbSubl;

    nTupleAnalysis::truthParticleHists* vHWW;
    nTupleAnalysis::truthParticleHists* vW;
    nTupleAnalysis::truthParticleHists* vWstar;
    nTupleAnalysis::truthParticleHists* vqLead;
    nTupleAnalysis::truthParticleHists* vqLead_s;
    nTupleAnalysis::truthParticleHists* vqLead_c;

    nTupleAnalysis::truthParticleHists* vqSubl;
    nTupleAnalysis::truthParticleHists* vqSubl_c;
    nTupleAnalysis::truthParticleHists* vqSubl_s;


    nTupleAnalysis::truthParticleHists* vqLead_mJet;
    nTupleAnalysis::truthParticleHists* vqSubl_mJet;
    nTupleAnalysis::truthParticleHists* vqLead_mJet25;
    nTupleAnalysis::truthParticleHists* vqSubl_mJet25;

    nTupleAnalysis::truthParticleHists* vqLead_mTrkJet;
    nTupleAnalysis::truthParticleHists* vqSubl_mTrkJet;


    TH2F* qLeadPt_vs_qSublPt;
    TH2F* qLeadEta_vs_qSublPt;
    TH2F* mW_vs_qSublPt;
    TProfile* ave_mW_vs_qSublPt     ;
    TProfile* ave_qLeadPt_vs_qSublPt;
    TProfile* ave_qLeadEta_vs_qSublPt;

    nTupleAnalysis::truthParticleHists* vl;
    nTupleAnalysis::truthParticleHists* vnu;

    TH1F* dRHH;
    TH1F* dRbb;
    TH1F* dRWW;
    TH1F* dRqq;
    TH1F* WstarHadronic;


    // Object Level
    nTupleAnalysis::jetHists*  allJets;
    nTupleAnalysis::jetHists*  selJets;
    nTupleAnalysis::fourVectorHists* allTrackJets;
    //nTupleAnalysis::jetHists*  btagJets;
    //nTupleAnalysis::dijetHists* diBJets;

    //nTupleAnalysis::muonHists* allMuons;
    //nTupleAnalysis::muonHists* preSelMuons;
    //nTupleAnalysis::dimuonHists* preSelDiMuons;

    //nTupleAnalysis::elecHists* allElecs;


    eventHists(std::string, fwlite::TFileService&, bool isMC = false, bool blind = true, std::string histDetailLevel = "", bool _debug = false, eventData* event=NULL);
    void Fill(eventData*);
    //void Fill(eventData* event, std::vector<std::shared_ptr<eventView>> &views);
    ~eventHists(); 

  };

}
#endif // eventHists_H

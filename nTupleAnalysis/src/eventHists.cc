#include "bbWW/nTupleAnalysis/interface/eventHists.h"
#include "nTupleAnalysis/baseClasses/interface/helpers.h"

using namespace bbWW;

using std::cout; using std::endl;

eventHists::eventHists(std::string name, fwlite::TFileService& fs, bool isMC, bool blind, std::string histDetailLevel, bool _debug, eventData* event) {
  std::cout << "Initialize >> eventHists: " << name << " with detail level: " << histDetailLevel << std::endl;
  dir = fs.mkdir(name);
  debug = _debug;

  if(isMC){
    vHH = new nTupleAnalysis::fourVectorHists(name+"/HH", fs, "HH System");
    dRHH = dir.make<TH1F>("dRHH", (name+"/dRHH;  #DeltaR(H,H); Entries").c_str(), 50,0,5);

    vHbb = new nTupleAnalysis::truthParticleHists(name+"/Hbb", fs, "H->bb System");
    dRbb = dir.make<TH1F>("dRbb", (name+"/dRbb;  #DeltaR(b,b); Entries").c_str(), 50,0,5);
    vbLead = new nTupleAnalysis::truthParticleHists(name+"/blead", fs, "leading b");
    vbSubl = new nTupleAnalysis::truthParticleHists(name+"/bsubl", fs, "sub-leading b");
    
    vHWW   = new nTupleAnalysis::truthParticleHists(name+"/HWW", fs, "H->WW System");
    dRWW = dir.make<TH1F>("dRWW", (name+"/dRWW;  #DeltaR(W,W); Entries").c_str(), 50,0,5);

    WstarHadronic = dir.make<TH1F>("WstarHadronic", (name+"/WstarHadronic;  Is W^{*} Hadronic?; Entries").c_str(), 2,-0.5,1.5);

    vW     = new nTupleAnalysis::truthParticleHists(name+"/W", fs, "W System");
    dRqq = dir.make<TH1F>("dRqq", (name+"/dRqq;  #DeltaR(q,q); Entries").c_str(), 50,0,5);
    vWstar = new nTupleAnalysis::truthParticleHists(name+"/Wstar", fs, "W^* System");

    vqLead = new nTupleAnalysis::truthParticleHists(name+"/qlead", fs, "leading q (from W)");
    vqLead_s = new nTupleAnalysis::truthParticleHists(name+"/qlead_s", fs, "leading s-q (from W)");
    vqLead_c = new nTupleAnalysis::truthParticleHists(name+"/qlead_c", fs, "leading c-q (from W)");
    vqLead_mJet = new nTupleAnalysis::truthParticleHists(name+"/qlead_mJet", fs, "leading q (from W)");
    vqLead_mJet25 = new nTupleAnalysis::truthParticleHists(name+"/qlead_mJet25", fs, "leading q (from W)");
    vqLead_mTrkJet = new nTupleAnalysis::truthParticleHists(name+"/qlead_mTrkJet", fs, "leading q (from W)");
    
    qLeadPt_vs_qSublPt   = dir.make<TH2F>("qLeadPt_vs_qSublPt",  "qLeadPt_vs_qSublPt;Subl-Leading q P_{T}-gen [GeV]; Leading q P_{T}-gen [GeV]",50,0,100,50,0,100);
    qLeadEta_vs_qSublPt   = dir.make<TH2F>("qLeadEta_vs_qSublPt",  "qLeadEta_vs_qSublPt;Subl-Leading q P_{T}-gen [GeV]; Leading q |#eta|-gen",50,0,100,50,0,3);
    mW_vs_qSublPt        = dir.make<TH2F>("mW_vs_qSublPt",  "mW_vs_qSublPt;Sub-leading q P_{T}-gen [GeV]; mW(qq)-gen [GeV]",50,0,100,50,30,130);
//      2) if it looks feasible, you could then plot <mW> vs pT,jet2 for selection with pT,jet1>30 GeV (or another suitable cut). 
//        Probably good to also simultaneously plot <pT,jet1> and <|eta,jet1> vs pT,jet2 for pT,jet1>30 GeV to make sure the reference jet position is not drifting, so we can trust that the reference JEC is stable
    
    ave_mW_vs_qSublPt        = dir.make<TProfile>("ave_mW_vs_qSublPt",        "ave_mW_vs_qSublPt;Sub-leading q P_{T}-gen [GeV]; mW(qq)-gen [GeV]", 50, 0, 100, 0,  100,"s");  
    ave_qLeadPt_vs_qSublPt   = dir.make<TProfile>("ave_qLeadPt_vs_qSublPt"  , "ave_qLeadPt_vs_qSublPt;Sub-leading q P_{T}-gen [GeV];<Leading q P_{T}-gen> [GeV]", 50, 0, 100, 0,  200,"s");  
    ave_qLeadEta_vs_qSublPt   = dir.make<TProfile>("ave_qLeadEta_vs_qSublPt"  ,"ave_qLeadEta_vs_qSublPt;Sub-leading q P_{T}-gen [GeV];<Leading |#eta|-gen>", 50, 0, 100, 0,  3, "s");  

    vqSubl = new nTupleAnalysis::truthParticleHists(name+"/qsubl", fs, "sub-leading q (from W)");
    vqSubl_s = new nTupleAnalysis::truthParticleHists(name+"/qsubl_s", fs, "sub-leading s-q (from W)");
    vqSubl_c = new nTupleAnalysis::truthParticleHists(name+"/qsubl_c", fs, "sub-leading c-q (from W)");

    vqSubl_mJet = new nTupleAnalysis::truthParticleHists(name+"/qsubl_mJet", fs, "sub-leading q (from W)");
    vqSubl_mJet25 = new nTupleAnalysis::truthParticleHists(name+"/qsubl_mJet25", fs, "sub-leading q (from W)");
    vqSubl_mTrkJet = new nTupleAnalysis::truthParticleHists(name+"/qsubl_mTrkJet", fs, "sub-leading q (from W)");

    vl     = new nTupleAnalysis::truthParticleHists(name+"/l", fs, "l (from W)");
    vnu    = new nTupleAnalysis::truthParticleHists(name+"/nu", fs, "nu (from W)");
  }

  allJets = new nTupleAnalysis::jetHists(name+"/allJets", fs, "All Jets");
  selJets = new nTupleAnalysis::jetHists(name+"/selJets", fs, "Selected Jets");

  allTrackJets = new nTupleAnalysis::fourVectorHists(name+"/allTrackJets", fs, "All Track Jets");
//  btagJets = new nTupleAnalysis::jetHists(name+"/btagJets", fs, "b-Tagged Jets");
//  diBJets  = new nTupleAnalysis::dijetHists(name+"/diBJets", fs, " Di-BJets");
//
//  allMuons        = new nTupleAnalysis::muonHists(name+"/allMuons", fs, "All Muons");
//  preSelMuons     = new nTupleAnalysis::muonHists(name+"/preSelMuons", fs, "PreSel Muons");
//  preSelDiMuons   = new nTupleAnalysis::dimuonHists(name+"/preSelDiMuons", fs, "PreDiSel Muons");
//  selDiMuons      = new nTupleAnalysis::dimuonHists(name+"/selDiMuons", fs, "SelDiSel Muons");
//
//  allElecs        = new nTupleAnalysis::elecHists(name+"/allElecs", fs, "All Elecs");

} 

void eventHists::Fill(eventData* event){
  if(debug) cout << "eventHists::Fill " << endl;


  if(event->truth){

    if(event->truth->Hbbs.size() == 1 && event->truth->HWWs.size() == 1){ 
      
      const nTupleAnalysis::particlePtr& Hbb = event->truth->Hbbs[0];
      const nTupleAnalysis::particlePtr& HWW = event->truth->HWWs[0];

      TLorentzVector pHH = (Hbb->p + HWW->p);
      vHH->Fill(pHH, event->weight);
      dRHH->Fill(Hbb->p.DeltaR(HWW->p), event->weight);

      // To add:
      //  - delta R(H,H)
      //  - delta R(bb)
      //  - delta R(WW)
      //  - delta R(qq)

      //
      //  H-> bb system
      //
      vHbb->Fill(Hbb, event->weight);
      dRbb->Fill(Hbb->daughters[0]->p.DeltaR(Hbb->daughters[1]->p), event->weight);
      if(Hbb->daughters[0]->p.Pt() > Hbb->daughters[1]->p.Pt()){
      	vbLead->Fill(Hbb->daughters[0], event->weight);
      	vbSubl->Fill(Hbb->daughters[1], event->weight);
      }else{
      	vbLead->Fill(Hbb->daughters[1], event->weight);
      	vbSubl->Fill(Hbb->daughters[0], event->weight);
      }

      //
      //   H->WW System
      // 
      vHWW->Fill(HWW, event->weight);
      dRWW->Fill(HWW->daughters[0]->p.DeltaR(HWW->daughters[1]->p), event->weight);
      if(HWW->daughters[0]->p.M() > HWW->daughters[1]->p.M()){
      	vW    ->Fill(HWW->daughters[0], event->weight);
      	vWstar->Fill(HWW->daughters[1], event->weight);
      }else{
      	vW    ->Fill(HWW->daughters[1], event->weight);
      	vWstar->Fill(HWW->daughters[0], event->weight);
      }

      bool isWstarHadronic = (event->truth->Wqqs[0]->p.M() < event->truth->Wlnus[0]->p.M());
      WstarHadronic->Fill(isWstarHadronic, event->weight);
    
    }// BBWW


    //
    //  W->qq
    //
    if(event->truth->Wqqs.size() == 1){   
      const nTupleAnalysis::particlePtr& Wqq = event->truth->Wqqs[0];
      dRqq->Fill(Wqq->daughters[0]->p.DeltaR(Wqq->daughters[1]->p), event->weight);

      unsigned int leadq_idx = 0;
      unsigned int sublq_Idx = 1;
      if(Wqq->daughters[0]->p.Pt() < Wqq->daughters[1]->p.Pt()){
	leadq_idx = 1;
	sublq_Idx = 0;
      }
      const nTupleAnalysis::particlePtr& leadq = Wqq->daughters[leadq_idx];
      const nTupleAnalysis::particlePtr& sublq = Wqq->daughters[sublq_Idx];
      
      vqLead  ->Fill(leadq, event->weight);

      if(abs(leadq->pdgId) == 3)
	vqLead_s  ->Fill(leadq, event->weight);

      if(abs(leadq->pdgId) == 4)
	vqLead_c  ->Fill(leadq, event->weight);

      vqSubl  ->Fill(sublq, event->weight);
      if(abs(sublq->pdgId) == 3)
	vqSubl_s  ->Fill(sublq, event->weight);
      if(abs(sublq->pdgId) == 4)
	vqSubl_c  ->Fill(sublq, event->weight);

      float qSublPt = sublq->p.Pt();
      float qLeadPt = leadq->p.Pt();
      float qLeadEta = fabs(leadq->p.Eta());
      float mWqq = Wqq->p.M();
      qLeadPt_vs_qSublPt->Fill(qSublPt,qLeadPt, event->weight);
      qLeadEta_vs_qSublPt->Fill(qSublPt,qLeadEta, event->weight);
      mW_vs_qSublPt->Fill(qSublPt, mWqq,       event->weight);

      ave_mW_vs_qSublPt->Fill(qSublPt,  mWqq, event->weight);
      ave_qLeadPt_vs_qSublPt->Fill(qSublPt,  qLeadPt, event->weight);
      ave_qLeadEta_vs_qSublPt->Fill(qSublPt,  qLeadEta, event->weight);


      const nTupleAnalysis::jetPtr leadq_matchedJet = leadq->matchedJet.lock();
      if(leadq_matchedJet){
	vqLead_mJet  ->Fill(leadq, event->weight);
	if(leadq_matchedJet->p.Pt() > 25) vqLead_mJet25  ->Fill(leadq, event->weight);
      }

      if(leadq->matchedTrackJet.lock()){
	vqLead_mTrkJet  ->Fill(leadq, event->weight);
      }



      const nTupleAnalysis::jetPtr sublq_matchedJet = sublq->matchedJet.lock();
      if(sublq_matchedJet){
	vqSubl_mJet  ->Fill(sublq, event->weight);
	if(sublq_matchedJet->p.Pt() > 25) vqSubl_mJet25  ->Fill(sublq, event->weight);
      }


      if(sublq->matchedTrackJet.lock()){
	vqSubl_mTrkJet  ->Fill(sublq, event->weight);
      }


    }// W->qq

    //
    //  W->lnu
    //
    if(event->truth->Wlnus.size() == 1){
      const nTupleAnalysis::particlePtr& Wlnu = event->truth->Wlnus[0];
      if(abs(Wlnu->daughters[0]->pdgId) % 2 == 1){
      	vl   ->Fill(Wlnu->daughters[0], event->weight);
      	vnu  ->Fill(Wlnu->daughters[1], event->weight);
      }else{
      	vl   ->Fill(Wlnu->daughters[1], event->weight);
      	vnu  ->Fill(Wlnu->daughters[0], event->weight);
      }
    }// Wlni


    //
    //  TTbars
    //
    if(event->truth->TWbs.size() == 2){

      const nTupleAnalysis::particlePtr& T0 = event->truth->TWbs[0];
      nTupleAnalysis::particlePtr& bFromT0 = (abs(T0->daughters[0]->pdgId) == 5) ? T0->daughters[0] : T0->daughters[1];

      const nTupleAnalysis::particlePtr& T1 = event->truth->TWbs[1];
      nTupleAnalysis::particlePtr& bFromT1 = (abs(T1->daughters[0]->pdgId) == 5) ? T1->daughters[0] : T1->daughters[1];

      nTupleAnalysis::particlePtr& leadb = bFromT0;
      nTupleAnalysis::particlePtr& sublb = bFromT1;
      if(leadb->p.Pt() < sublb->p.Pt()){
	leadb = bFromT1;
	sublb = bFromT0;
      }

      vbLead->Fill(leadb, event->weight);
      vbSubl->Fill(sublb, event->weight);

    }// Ts



  
  }//event->truth
  


  //
  // Object Level
  //

  //
  // Jets 
  //
  if(debug) cout << "eventHists::fillJets " << endl;
  allJets->nJets->Fill(event->allJets.size(), event->weight);
  for(auto &jet: event->allJets) allJets->Fill(jet, event->weight);

  selJets->nJets->Fill(event->selJets.size(), event->weight);
  for(auto &jet: event->selJets) selJets->Fill(jet, event->weight);

  //
  //  Track Jets
  //
  //allJets->nJets->Fill(event->allJets.size(), event->weight);
  for(auto &trkJet: event->allTrackJets) allTrackJets->Fill(trkJet->p, event->weight);


//  unsigned int nBTagJets = event->btagJets.size()   ;
//  btagJets->nJets->Fill(nBTagJets, event->weight);
//  for(auto &jet: event->btagJets) btagJets->Fill(jet, event->weight);


//  //
//  //  Muons
//  //
//  if(debug) cout << "eventHists::fillMuons " << endl;
//  allMuons->nMuons->Fill(event->allMuons.size(), event->weight);
//  for(auto &muon: event->allMuons) allMuons->Fill(muon, event->weight);
//
//  preSelMuons->nMuons->Fill(event->preSelMuons.size(), event->weight);
//  for(auto &muon: event->preSelMuons) preSelMuons->Fill(muon, event->weight);
//
//  unsigned int nPreSelDiMuons = event->preSelDiMuons.size();
//  preSelDiMuons->nDiMuons->Fill(nPreSelDiMuons, event->weight);
//  for(auto &dimuon: event->preSelDiMuons) preSelDiMuons->Fill(dimuon, event->weight);
//
//  if(nPreSelDiMuons)
//    selDiMuons->Fill(event->preSelDiMuons[0], event->weight);
//


//  //
//  //  Electrons
//  //
//  if(debug) cout << "eventHists::fillElectrons " << endl;
//  allElecs->nElecs->Fill(event->allElecs.size(), event->weight);
//  for(auto &elec: event->allElecs)             allElecs->Fill(elec, event->weight);
//
//  //
//  //  di-BJets
//  // 
//  if(nPreSelDiMuons && nBTagJets > 1){
//    if(debug) cout << "eventHists::fillEvents " << endl;
//
//    diBJets->Fill(event->diBJets, event->weight);
//
//    //const nTupleAnalysis::dimuonPtr& jPsi = event->preSelDiMuons[0]; 
//
//
//  }


  if(debug) cout << "eventHists::Fill left " << endl;
  return;
}

eventHists::~eventHists(){} 


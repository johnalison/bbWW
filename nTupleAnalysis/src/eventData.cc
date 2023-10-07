#include "bbWW//nTupleAnalysis/interface/eventData.h"
#include "nTupleAnalysis/baseClasses/interface/helpers.h"

using namespace bbWW;

using std::cout; using std::endl; 
using std::vector; using std::string;
#include <cmath>

// Sorting functions
bool sortPt(       std::shared_ptr<nTupleAnalysis::jet>       &lhs, std::shared_ptr<nTupleAnalysis::jet>       &rhs){ return (lhs->pt        > rhs->pt   );     } // put largest  pt    first in list
bool sortDeepFlavB(std::shared_ptr<nTupleAnalysis::jet>       &lhs, std::shared_ptr<nTupleAnalysis::jet>       &rhs){ return (lhs->deepFlavB > rhs->deepFlavB); } // put largest  deepB first in list
bool sortDimuonMJPsi(std::shared_ptr<nTupleAnalysis::dimuon>  &lhs, std::shared_ptr<nTupleAnalysis::dimuon>    &rhs){ return (fabs(lhs->m - nTupleAnalysis::mJpsi) < fabs(rhs->m - nTupleAnalysis::mJpsi)); } 
bool sortDeepFlavCvLandCvB(std::shared_ptr<nTupleAnalysis::jet>       &lhs, std::shared_ptr<nTupleAnalysis::jet>       &rhs){ return ( lhs->XCvLCvB < rhs->XCvLCvB); } // put largest  deepB first in list


eventData::eventData(TChain* t, bool mc, std::string y, bool d, bool _fastSkim, std::string bjetSF, std::string btagVariations, std::string JECSyst){
  std::cout << "eventData::eventData()" << std::endl;
  tree  = t;
  isMC  = mc;
  year  = ::atof(y.c_str());
  debug = d;
  fastSkim = _fastSkim;


  eventDataTree = new nTupleAnalysis::eventData("", tree, true, isMC);

  if(isMC){
    //inputBranch(tree, "genWeight", genWeight);

    if(tree->FindBranch("nGenPart")){
      truth = new nTupleAnalysis::truthData(tree, debug, "GenPart", true, true);
    }else{
      cout << "No GenPart (missing branch 'nGenPart'). Will ignore ..." << endl;
    }

    if(tree->FindBranch("nGenJet")){
      truthJets = new nTupleAnalysis::truthParticle("GenJet", tree, true);
    }else{
      cout << "No GenPart (missing branch 'nGenPart'). Will ignore ..." << endl;
    }

  }


  //triggers https://twiki.cern.ch/twiki/bin/viewauth/CMS/HLTPathsRunIIList
  if(year==2016){

  }

  if(year==2017){

  }

  if(year==2018){

    HLT_triggers["HLT_Dimuon0_Jpsi3p5_Muon2"] = false;
    //HLT_triggers["HLT_Dimuon25_Jpsi_noCorrL1_v"] = false;
    HLT_triggers["HLT_Dimuon25_Jpsi"] = false;
  }

  for(auto &trigger:  L1_triggers)     inputBranch(tree, trigger.first, trigger.second);
  for(auto &trigger: HLT_triggers)     inputBranch(tree, trigger.first, trigger.second);
  //for(auto &trigger:  L1_triggers_mon){
  //  if(L1_triggers.find(trigger.first)!=L1_triggers.end()) continue; // don't initialize branch again!
  //  inputBranch(tree, trigger.first, trigger.second);
  //}


  std::cout << "eventData::eventData() Initialize jets" << std::endl;
  treeJets  = new  nTupleAnalysis::jetData(    "Jet", tree, true, isMC, "", "", bjetSF, btagVariations, JECSyst);
  std::cout << "eventData::eventData() Initialize track jets" << std::endl;
  treeTrackJets  = new  nTupleAnalysis::trackJetData(    "SoftActivityJet", tree, true);
  std::cout << "eventData::eventData() Initialize muons" << std::endl;
  treeMuons = new nTupleAnalysis::muonData(   "Muon", tree, true, isMC);
  std::cout << "eventData::eventData() Initialize elecs" << std::endl;
  treeElecs = new nTupleAnalysis::elecData(   "Electron", tree, true, isMC);
  std::cout << "eventData::eventData() Initialize TrigObj" << std::endl;
  //treeTrig  = new trigData("TrigObj", tree);
} 




void eventData::resetEvent(){
  if(debug) std::cout<<"Reset eventData"<<std::endl;
  preSelDiMuons.clear();


  WqqTandPPairs.clear();
  WqqTandPPairs_dR.clear();
  WqqTandPPairs_pT.clear();
  WqqTandPPairs_xWbW.clear();


  passL1  = false;
  passHLT = false;
  bTagSF = 1;
  treeJets->resetSFs();
  nTrueBJets = 0;
  qjet25 = false;
  qjetEta = false;
  lepton25 = false;

  WqTagJets  .clear();
}



void eventData::update(long int e){
  if(debug){
    std::cout<<"Get Entry "<<e<<std::endl;
    std::cout<<tree->GetCurrentFile()->GetName()<<std::endl;
    tree->Show(e);
  }

  // if(printCurrentFile && tree->GetCurrentFile()->GetName() != currentFile){
  //   currentFile = tree->GetCurrentFile()->GetName();
  //   std::cout<< std::endl << "Loading: " << currentFile << std::endl;
  // }

  Long64_t loadStatus = tree->LoadTree(e);
  if(loadStatus<0){
   std::cout << "Error "<<loadStatus<<" getting event "<<e<<std::endl; 
   return;
  }

  tree->GetEntry(e);
  if(debug) std::cout<<"Got Entry "<<e<<std::endl;


  //
  // Reset the derived data
  //
  resetEvent();

  if(truth) truth->update();



  if(truth){

    genJets   = truthJets->getParticles();

    if(truth->Wqqs.size() == 1 && truth->Wlnus.size() == 1  ){ 
      wStarlnu = (truth->Wqqs[0]->p.M() > truth->Wlnus[0]->p.M());

      if(abs(truth->Wlnus[0]->daughters[0]->pdgId) % 2 == 1){
      	lepton25 = (truth->Wlnus[0]->daughters[0]->p.Pt() > 25);
      }else{
	lepton25 = (truth->Wlnus[0]->daughters[1]->p.Pt() > 25);
      }

      qjet25 = ((truth->Wqqs[0]->daughters[0]->p.Pt() > 25) && (truth->Wqqs[0]->daughters[1]->p.Pt() > 25));

      qjetEta = ( (fabs(truth->Wqqs[0]->daughters[0]->p.Eta()) < 2.5) && (fabs(truth->Wqqs[0]->daughters[1]->p.Eta()) < 2.5));
      qjetLead30 = ((truth->Wqqs[0]->daughters[0]->p.Pt() > 30) || (truth->Wqqs[0]->daughters[1]->p.Pt() > 30));

    }

    //
    //  Calculate mbbWW
    //
    if(truth->Hbbs.size() == 1 && truth->HWWs.size() == 1){ 
      mbbWW = (truth->Hbbs[0]->p + truth->HWWs[0]->p).M();
    }else{
      mbbWW = -1;
    }
    

    //cout << " ------------------- " << endl;
    //truth->dump();
  }

  //Objects from ntuple
  if(debug) std::cout << "Get Jets\n";
  //getJets(float ptMin = -1e6, float ptMax = 1e6, float etaMax = 1e6, bool clean = false, float tagMin = -1e6, std::string tagger = "CSVv2", bool antiTag = false, int puIdMin = 0);

  if(debug) std::cout << "Get Muons\n";

  buildEvent();


  //
  //  Do q to jet matching
  //
  if(truth){

    for(const nTupleAnalysis::particlePtr& W :truth->Wqqs){
      for(const nTupleAnalysis::particlePtr& q : W->daughters){
	//cout << " " << q->p.Pt() << " " << q->p.Eta() <<  "  " << q->p.Phi() << endl;

	//
	//  Match Jets
	//
	float mindR = 99;
	nTupleAnalysis::jetPtr matchedJet  = nullptr;
	for(const nTupleAnalysis::jetPtr &jet: allJets){
	  float thisDr = jet->p.DeltaR(q->p);
	  
	  if(thisDr < mindR){
	    mindR = thisDr;
	    matchedJet = jet;
	  }

	}// jets

	if(mindR < 0.4){
	  q->matchedJet = matchedJet;
	}

	//
	//  Try to Match Track Jets
	//
	float mindR_trkJet = 99;
	nTupleAnalysis::trackJetPtr matchedTrackJet  = nullptr;
	for(const nTupleAnalysis::trackJetPtr &trackJet: allTrackJets){
	  float thisDr = trackJet->p.DeltaR(q->p);
	  
	  if(thisDr < mindR_trkJet){
	    mindR_trkJet = thisDr;
	    matchedTrackJet = trackJet;
	  }

	}// track jets

	if(mindR_trkJet < 0.4){
	  q->matchedTrackJet = matchedTrackJet;
	}


      }// qs
    }// Wqqs   
  }


  for(auto &trigger: HLT_triggers){
    ///bool pass_seed = boost::accumulate(HLT_L1_seeds[trigger.first] | boost::adaptors::map_values, false, [](bool pass, bool *seed){return pass||*seed;});//std::logical_or<bool>());
    //passL1  = passL1  || pass_seed;
    //passHLT = passHLT || (trigger.second && pass_seed);
    passHLT = passHLT || (trigger.second);
  }


  //
  // HACK for now
  //
  passHLT = true;
  
  if(debug) std::cout<<"eventData updated\n";
  return;
}

void eventData::buildEvent(){

  //
  // Select Jets
  //
  allJets   = treeJets->getJets(0.0, 1e6, 1e6, false, -1e6, bTagger, false, puIdMin);
  selJets   = treeJets->getJets(       allJets, jetPtMin,   1e6, jetEtaMax, doJetCleaning);
  selJets30 = treeJets->getJets(       selJets,   30,   1e6, jetEtaMax, doJetCleaning);
  btagJets  = treeJets->getJets(       selJets30, 30,   1e6, jetEtaMax, doJetCleaning, bTag,   bTagger);


  //
  //  Select Track Jets
  //
  allTrackJets  = treeTrackJets->getJets(0.0, 1e6, 1e6);
  
  // Tight cuts for now https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/#ak4-c-tagging
  //std::vector<nTupleAnalysis::jetPtr> notLJets = treeJets->getJets(       selJets, jetPtMin,   1e6, jetEtaMax, doJetCleaning, 0.282,   "deepFlavCvL");
  //ctagJets        = treeJets->getJets(       notLJets, jetPtMin,   1e6, jetEtaMax, doJetCleaning, 0.267,   "deepFlavCvB");


  //btag SF
  if(isMC){
    //for(auto &jet: selJets) bTagSF *= treeJets->getSF(jet->eta, jet->pt, jet->deepFlavB, jet->hadronFlavour);
    treeJets->updateSFs(selJets, debug);
    bTagSF = treeJets->m_btagSFs["central"];

    if(debug) std::cout << "eventData buildEvent bTagSF = " << bTagSF << std::endl;
    //weight *= bTagSF;
    for(auto &jet: allJets) nTrueBJets += jet->hadronFlavour == 5 ? 1 : 0;
  }
  
  //
  // Select Muons
  //
  allMuons         = treeMuons->getMuons();
  preSelMuons      = treeMuons->getMuons(3.0, 2.4, 5, false, -1);


  //
  //  Select Electrons
  //
  allElecs         = treeElecs->getElecs();

  
  

//  //
//  // If Di-Jets 
//  //
//  if(selJets.size() > 1){
//    std::sort(selJets.begin(),       selJets.end(),       sortDeepFlavB);
//    diBJets = std::make_shared<nTupleAnalysis::dijet>(nTupleAnalysis::dijet(selJets[0], selJets[1]));
//
//    std::sort(selJets.begin(),       selJets.end(),       sortDeepFlavCvLandCvB);
//    diCJets = std::make_shared<nTupleAnalysis::dijet>(nTupleAnalysis::dijet(selJets[0], selJets[1]));
//  }
//
//  //
//  // If Di-Jets and di-Muons
//  //
//  if(selJets.size() > 1 && preSelDiMuons.size()){
//    const nTupleAnalysis::dimuonPtr& jPsi = preSelDiMuons[0]; 
// 
//
//
//  }


  


  
  if(debug) std::cout<<"eventData buildEvent done\n";
  return;
}




void eventData::doTTbarTandPSelection(){

  //
  //  WqTag jets are 30 GeV non-bTagged jets
  //
  for(nTupleAnalysis::jetPtr& sJet30 : selJets30){
    if(find(btagJets.begin(), btagJets.end(), sJet30) != btagJets.end()){
      if(debug) cout << "fails WqTag-bTag overlap " <<endl;
      continue;
    }
    WqTagJets.push_back(sJet30);

    
    for(nTupleAnalysis::jetPtr& sJet : selJets){
      if(sJet == sJet30) continue;

      if(sJet->p.Pt() > sJet30->p.Pt()) continue;

      if(find(btagJets.begin(), btagJets.end(), sJet) != btagJets.end()){
	if(debug) cout << "fails WqProbe-bTag overlap " <<endl;
	continue;
      }

      std::shared_ptr<nTupleAnalysis::trijet> trijet_0 = std::make_shared<nTupleAnalysis::trijet>(nTupleAnalysis::trijet(btagJets.at(0),sJet30,sJet, truth));
      std::shared_ptr<nTupleAnalysis::trijet> trijet_1 = std::make_shared<nTupleAnalysis::trijet>(nTupleAnalysis::trijet(btagJets.at(1),sJet30,sJet, truth));
      
      if(trijet_0->xWt < trijet_1->xWt)
	WqqTandPPairs.push_back(trijet_0);
      else
	WqqTandPPairs.push_back(trijet_1);

    }
  }

  //
  // Cuts on the TandP probes
  //
  for(std::shared_ptr<nTupleAnalysis::trijet>& pair : WqqTandPPairs){
    if(pair->W->dR > 3.5) continue;      
    WqqTandPPairs_dR.push_back(pair);

    if((pair->W->lead->p.Pt() + 1.33*pair->W->subl->p.Pt()) < 70) continue;
    WqqTandPPairs_pT.push_back(pair);

    if(pair->xWbW  > 4) continue;
    WqqTandPPairs_xWbW.push_back(pair);

    
  }
// 
  


}




void eventData::dump(){

  std::cout << "   Run: " << eventDataTree->run    << std::endl;
  std::cout << " Event: " << eventDataTree->event  << std::endl;  
  //std::cout << "Weight: " << eventDataTree->weight << std::endl;
  //std::cout << "Trigger Weight : " << trigWeight << std::endl;
  //std::cout << "WeightNoTrig: " << weightNoTrigger << std::endl;
  std::cout << " allJets: " << allJets .size() << " |  selJets: " << selJets .size() << " | btagJets: " << btagJets.size() << std::endl;
  std::cout << "allMuons: " << allMuons.size() << " | preSelMuons: " << preSelMuons.size() << std::endl;

  cout << "All Jets" << endl;
  for(auto& jet : allJets){
    std::cout << "\t " << jet->pt << " (" << jet->pt_wo_bRegCorr << ") " <<  jet->eta << " " << jet->phi << " " << jet->deepB  << " " << jet->deepFlavB << " " << (jet->pt - 40) << std::endl;
  }

  cout << "Sel Jets" << endl;
  for(auto& jet : selJets){
    std::cout << "\t " << jet->pt << " " << jet->eta << " " << jet->phi << " " << jet->deepB  << " " << jet->deepFlavB << std::endl;
  }

  cout << "bTag Jets" << endl;
  for(auto& jet : btagJets){
    std::cout << "\t " << jet->pt << " " << jet->eta << " " << jet->phi << " " << jet->deepB  << " " << jet->deepFlavB << std::endl;
  }


  return;
}

eventData::~eventData(){} 


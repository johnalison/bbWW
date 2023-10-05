//#include "TChain.h"
#include "bbWW/nTupleAnalysis/interface/cutflowHists.h"

using namespace bbWW;

cutflowHists::cutflowHists(std::string name, fwlite::TFileService& fs, bool isMC, bool _debug) {
  debug = _debug;
  dir = fs.mkdir(name);
  unitWeight = dir.make<TH1I>("unitWeight", (name+"/unitWeight; ;Entries").c_str(),  1,1,2);
  unitWeight->SetCanExtend(1);
  
  weighted = dir.make<TH1D>("weighted", (name+"/weighted; ;Entries").c_str(),  1,1,2);
  weighted->SetCanExtend(1);

  AddCut("all");

} 

void cutflowHists::AddCut(std::string cut){
  unitWeight->GetXaxis()->FindBin(cut.c_str());  
  weighted->GetXaxis()->FindBin(cut.c_str());
}

void cutflowHists::BasicFill(const std::string& cut, eventData* event, float weight){
  if(debug) std::cout << "cutflowHists::BasicFill(const std::string& cut, eventData* event, float weight) " << cut << std::endl; 
  unitWeight->Fill(cut.c_str(), 1);
  weighted  ->Fill(cut.c_str(), weight);
  return;
}

void cutflowHists::BasicFill(const std::string& cut, eventData* event){
  BasicFill(cut, event, event->weight);
  return;
}

void cutflowHists::Fill(const std::string& cut, eventData* event){
  if(debug) std::cout << "cutflowHists::Fill(const std::string& cut, eventData* event) " << cut << std::endl;

  BasicFill(cut, event);

  return;
}

void cutflowHists::labelsDeflate(){
  unitWeight->LabelsDeflate("X");
  //unitWeight->LabelsOption("a");
  weighted  ->LabelsDeflate("X");
  //weighted  ->LabelsOption("a");
  return;  
}

cutflowHists::~cutflowHists(){}

// -*- C++ -*-
#if !defined(cutflowHists_H)
#define cutflowHists_H

#include <iostream>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptor/map.hpp>
#include <TH1F.h>
#include <TH2F.h>
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "bbWW/nTupleAnalysis/interface/eventData.h"

namespace bbWW {

  class cutflowHists {
  public:
    bool debug = false;
    TFileDirectory dir;
    
    TH1I* unitWeight;
    TH1D* weighted;

    cutflowHists(std::string, fwlite::TFileService&, bool, bool);
    void BasicFill(const std::string& cut, eventData* event);
    void BasicFill(const std::string& cut, eventData* event, float weight);
    void Fill(const std::string&, eventData*);

    void labelsDeflate();

    void AddCut(std::string cut);
    
    ~cutflowHists(); 

  };

}
#endif // cutflowHists_H

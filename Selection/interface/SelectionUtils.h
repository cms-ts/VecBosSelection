//#include "DataFormats/PatCandidates/interface/PATObject.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/Event.h"

using namespace std;
using namespace edm;

class SelectionUtils {

  public:

   bool DoWP80(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,bool removePU_);
   bool DoWP80Pf(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,bool removePU_);
   bool DoWP80Pf_NewHE(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,bool removePU_);
   bool DoHLTMatch(pat::ElectronCollection::const_iterator,edm::Event&);
   std::vector<bool> MakeEleIDAnalysis(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,bool removePU_);
   std::vector<bool> MakePfEleIDAnalysis(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,bool removePU_); 

//private:
		//bool removePU_;

};


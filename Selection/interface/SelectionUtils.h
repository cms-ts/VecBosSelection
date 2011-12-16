//#include "DataFormats/PatCandidates/interface/PATObject.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/Event.h"

using namespace std;
using namespace edm;

class SelectionUtils {

	public:

		bool DoWP80(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent);
		bool DoHLTMatch(pat::ElectronCollection::const_iterator,edm::Event&);


private:
      bool removePU_;

};

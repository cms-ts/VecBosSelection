#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include <TMath.h>

using namespace std;
using namespace edm;


typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > > IsoDepositMaps;
typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;

class SelectionUtils {

  public:

   bool DoWP80(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,bool removePU_);
   bool DoWP80Pf(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent);
   bool DoWP80PfGSF(reco::GsfElectronCollection::const_iterator recoElectron,edm::Event& iEvent);
   bool DoWP90Pf(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent);
   bool DoWP90PfGSF(reco::GsfElectronCollection::const_iterator recoElectron,edm::Event& iEvent);
   bool DoWP80Pf_NewHE(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,bool removePU_);
   bool DoHLTMatch(pat::ElectronCollection::const_iterator,edm::Event&);
   std::vector<bool> MakeEleIDAnalysis(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,bool removePU_);
   std::vector<bool> MakePfEleIDAnalysis(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,bool removePU_);

   bool DoMedSel2011(pat::ElectronCollection::const_iterator recoElectron,
		     edm::Event& iEvent,
		     const edm::Handle<reco::ConversionCollection> &conversions,
		     const reco::BeamSpot &beamspot,
		     const edm::Handle<reco::VertexCollection> &vtxs);
   //bool DoIso2011(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent, IsoDepositMaps &IsoMap, IsoDepositVals &IsoVal); 
   bool DoIso2011(pat::ElectronCollection::const_iterator recoElectron, edm::Event& iEvent, IsoDepositVals &IsoVals); 
   bool DoIso2011GSF(reco::GsfElectronRef myElectronRef, edm::Event& iEvent, IsoDepositVals &IsoVals); 
   std::vector<bool> MakePfEleNewIDAnalysis(pat::ElectronCollection::const_iterator recoElectron,edm::Event& iEvent,
					    bool useNewID_,bool doWP90_,const edm::Handle<reco::ConversionCollection> &conversions,
					    const reco::BeamSpot &beamspot,const edm::Handle<reco::VertexCollection> &vtxs, 
					    IsoDepositVals &IsoVals);

   static bool isGoodConversion(const reco::Conversion &conv, const math::XYZPoint &beamspot, float lxyMin=2.0, float probMin=1e-6, unsigned int nHitsBeforeVtxMax=1);   
   static bool matchesConversion(const pat::Electron &ele, const reco::Conversion &conv, bool allowCkfMatch=true);
   static bool  hasMatchedConversion(const pat::Electron &ele, const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot, bool allowCkfMatch=true, float lxyMin=2.0, float probMin=1e-6, unsigned int nHitsBeforeVtxMax=0);

};


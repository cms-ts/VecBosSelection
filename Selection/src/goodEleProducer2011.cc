#include "VecBosSelection/Selection/interface/goodEleProducer2011.h"

using namespace std;
using namespace edm;
using namespace reco;


//
// member functions
//

// ------------ method called for each event  ------------
void
goodEleProducer2011::produce(edm::Event & iEvent, edm::EventSetup const & iSetup)
{
   
   // get gfs Electron  
   auto_ptr< reco::PFCandidateCollection > 
      pOutput( new reco::PFCandidateCollection ); 

  if (Debug2) cout<<"------- NEW Event -----"<<endl;
  double delta = 0.001;
  int cont = 0;

  Handle < pat::ElectronCollection > electronCollection;  
  iEvent.getByLabel (theElectronCollectionLabel, electronCollection);

  Handle<reco::PFCandidateCollection> pfElecCollection;
  iEvent.getByLabel (pflowEleCollection_, pfElecCollection);
  
  // conversions
  edm::Handle<reco::ConversionCollection> conversions_h;
  iEvent.getByLabel(conversionsInputTag_, conversions_h);
  
  // iso deposits
  IsoDepositVals isoVals(isoValInputTags_.size());
  for (size_t j = 0; j < isoValInputTags_.size(); ++j) {
     iEvent.getByLabel(isoValInputTags_[j], isoVals[j]);
  }
  
  // beam spot
  edm::Handle<reco::BeamSpot> beamspot_h;
  iEvent.getByLabel(beamSpotInputTag_, beamspot_h);
  const reco::BeamSpot &beamSpot = *(beamspot_h.product());
  
  // vertices
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByLabel(primaryVertexInputTag_, vtx_h);

  if (electronCollection.isValid () && electronCollection->size()>1){

     int i=0;
     bool protection=false;
     int jj=0;
     int sizePat=electronCollection->size();
     protection=false;  
     passSelection=true;


     /// NEW DS
     // Cutting on WP80
     for (pat::ElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
			
	protection=true;
	jj++;
	//Perform checks on each ele ID criteria
	// Here you get a plot full of information. Each electron contributes with one entry (so total numer of entries = 3* #electrons)
	// To have the "%", each bin value shold be divided by total numer of entries/3
	std::vector<bool> result=SelectionUtils::MakePfEleNewIDAnalysis(recoElectron,iEvent,useNewID_,conversions_h,beamSpot,vtx_h, isoVals);
	passIDEleCriteria->SetBinContent(1,passIDEleCriteria->GetBinContent(1)+1);
	if (result[0]) passIDEleCriteria->SetBinContent(2,passIDEleCriteria->GetBinContent(2)+1);
	if (result[1]) passIDEleCriteria->SetBinContent(3,passIDEleCriteria->GetBinContent(3)+1);
	if (result[2]) passIDEleCriteria->SetBinContent(4,passIDEleCriteria->GetBinContent(4)+1);

	if ( (!doID_ || ((!useNewID_ && SelectionUtils::DoWP80Pf(recoElectron,iEvent))
			|| (useNewID_ && SelectionUtils::DoMedSel2011(recoElectron,iEvent,conversions_h,beamSpot,vtx_h))))
	     && (!doIsolation_ ||  SelectionUtils::DoIso2011(recoElectron, iEvent, isoVals))
	     //&& SelectionUtils::DoHLTMatch(recoElectron,iEvent) 
	     //&& recoElectron->pt()>secondEleEnThrhold
	   ){
	   for (reco::PFCandidateCollection::const_iterator pfElectron = pfElecCollection->begin (); pfElectron != pfElecCollection->end (); pfElectron++) {
	      if (fabs(pfElectron->pdgId())==11){
		 if ( fabs(pfElectron->pt() - recoElectron->pt())<delta && 
		      fabs(pfElectron->eta() - recoElectron->eta())<delta &&
		      fabs(pfElectron->phi() - recoElectron->phi())<delta ) {
		    pOutput->push_back(*pfElectron);
		    h_electronEn->Fill(pfElectron->energy());
		    h_electronEta->Fill(pfElectron->eta());
		    h_electronPt->Fill(pfElectron->pt());
		    cont++;
		    continue;
		 }
	      }
	   }
	}	
     }
     
     if (!protection) {
	cout<<"size pat is "<<sizePat<<" while jj is "<<jj<<" and protection "<<protection<<endl;
	cout<<"problems with PAT collection, in Efficiency.cc-->... Please check..."<<endl;    
	passSelection = false;
     }
      
  } // if electronCollection is valid
  
  eventAccept->Fill(cont); 
  iEvent.put( pOutput );
  
}


// ------------ method called once each job just before starting event loop  ------------
void
goodEleProducer2011::beginJob (){

//beginJob

}


// ------------ method called once each job just after ending the event loop  ------------
void
goodEleProducer2011::endJob ()
{

//endJob

}

DEFINE_FWK_MODULE (goodEleProducer2011);

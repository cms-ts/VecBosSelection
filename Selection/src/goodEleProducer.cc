#include "VecBosSelection/Selection/interface/goodEleProducer.h"

using namespace std;
using namespace edm;
using namespace reco;


//
// member functions
//

// ------------ method called for each event  ------------
void
goodEleProducer::produce(edm::Event & iEvent, edm::EventSetup const & iSetup)
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
	std::vector<bool> result=SelectionUtils::MakeEleIDAnalysis(recoElectron,iEvent,removePU_); 
	if (result[0]) passIDEleCriteria->SetBinContent(1,passIDEleCriteria->GetBinContent(1)+1);
	if (result[1]) passIDEleCriteria->SetBinContent(2,passIDEleCriteria->GetBinContent(2)+1);
	if (result[2]) passIDEleCriteria->SetBinContent(3,passIDEleCriteria->GetBinContent(3)+1);
	
	if ( SelectionUtils::DoWP80Pf(recoElectron,iEvent,removePU_) 
	     && SelectionUtils::DoHLTMatch(recoElectron,iEvent) 
	     //&& recoElectron->pt()>20.0
	   ){
	   for (reco::PFCandidateCollection::const_iterator pfElectron = pfElecCollection->begin (); 
		pfElectron != pfElecCollection->end (); pfElectron++) {
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
     
     //if (cont>2) h_zPt_3e->Fill(e_pair.Pt());
	   	   
	
  } // if electronCollection is valid
  
  eventAccept->Fill(cont); 
  iEvent.put( pOutput );
  
}


// ------------ method called once each job just before starting event loop  ------------
void
goodEleProducer::beginJob (){

//beginJob

}


// ------------ method called once each job just after ending the event loop  ------------
void
goodEleProducer::endJob ()
{

//endJob

}

DEFINE_FWK_MODULE (goodEleProducer);

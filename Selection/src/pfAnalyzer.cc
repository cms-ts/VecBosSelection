#include "VecBosSelection/Selection/interface/pfAnalyzer.h"

using namespace std;
using namespace edm;
using namespace reco;


//
// member functions
//

// ------------ method called for each event  ------------
void
pfAnalyzer::produce(edm::Event & iEvent, edm::EventSetup const & iSetup)
{

   // get gfs Electron  
   auto_ptr< reco::PFCandidateCollection > 
      pOutput( new reco::PFCandidateCollection ); 

  if (Debug2) cout<<"------- NEW Event -----"<<endl;
  double delta = 0.001;
  int cont = 0;

  Handle < pat::ElectronCollection > electronCollection;  
  iEvent.getByLabel (theElectronCollectionLabel, electronCollection);
  if (electronCollection.isValid () && electronCollection->size()>1){
 
     pat::ElectronCollection::const_iterator highestptele;
     pat::ElectronCollection::const_iterator secondptele;

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
	
	if ( SelectionUtils::DoWP80Pf(recoElectron,iEvent,removePU_) && SelectionUtils::DoHLTMatch(recoElectron,iEvent) && recoElectron->pt()>10.0){
	//if ( SelectionUtils::DoWP80(recoElectron,iEvent) && SelectionUtils::DoHLTMatch(recoElectron,iEvent)){
	   if (Debug2) cout<<"Tag is a WP80 electron..."<<endl;
	   
	   //Sort in Pt
	   if (Debug2) cout<<"Electron pt value ->"<<recoElectron->pt()<<endl;
	   if (Debug2) cout<<" MMMM ele trigger size "<<recoElectron->triggerObjectMatches().size()<<endl;
	   if (i==0) highestptele=recoElectron;
	   if (i==1){
	      if (highestptele->pt()<recoElectron->pt()){
		 secondptele=highestptele;
		 highestptele=recoElectron;
	      }
	      else{
		 secondptele=recoElectron;
	      }
	   }
	   if (i>1){
	      if (highestptele->pt()<recoElectron->pt()){
		 secondptele=highestptele;
		 highestptele=recoElectron;
	      }
	      else{
		 if (secondptele->pt()<recoElectron->pt()) secondptele=recoElectron;
	      }
	   }
	   i++;
	} 

	else{
	   if (Debug2) cout<<"Tag IS a NOT WP80 electron...Exit"<<endl;
	}
		
     }
     
     if (!protection) {
	cout<<"size pat is "<<sizePat<<" while jj is "<<jj<<" and protection "<<protection<<endl;
	cout<<"problems with PAT collection, in Efficiency.cc-->... Please check..."<<endl;    
	passSelection = false;
     }
     
     
     if(i>1 && highestptele->pt()>=20.0) { //you NEED at least two electrons :)
	
	//--------------
	// Match the HLT
	pat::ElectronCollection::const_iterator tag;
	pat::ElectronCollection::const_iterator probe;
	
	tag=highestptele; //Ã¨ solo un rinominare le cose... (non e' necessario solo retaggio di codice copiato)
	probe=secondptele;
	
	//Calculating Invariant Mass
	TLorentzVector tagv;
	tagv.SetPtEtaPhiM(tag->pt(),tag->eta(),tag->phi(), 0.0);
	TLorentzVector probev;
	probev.SetPtEtaPhiM(probe->pt(),probe->eta(),probe->phi(), 0.0);
	
	TLorentzVector e_pair = tagv + probev;
	double e_ee_invMass = e_pair.M ();
	
	
	//Cut on the tag and probe mass...
	if (e_ee_invMass>120 || e_ee_invMass<60) passSelection=false;
	
	
// searching for the corresponding GSF electron 
	if (passSelection){
	   Handle<reco::PFCandidateCollection> pfElecCollection;
	   iEvent.getByLabel (pflowEleCollection_, pfElecCollection);
	   
	   for (reco::PFCandidateCollection::const_iterator pfElectron = pfElecCollection->begin (); pfElectron != pfElecCollection->end (); pfElectron++) {
	      if (fabs(pfElectron->pdgId())==11){
		 if ( fabs(pfElectron->pt() - tag->pt())<delta && 
		      fabs(pfElectron->eta() - tag->eta())<delta &&
		      fabs(pfElectron->phi() - tag->phi())<delta ) {
		    pOutput->push_back(*pfElectron);
		    h_electronEn->Fill(pfElectron->energy());
		    h_electronEta->Fill(pfElectron->eta());
		    h_electronPt->Fill(pfElectron->pt());
		    cont++;
		    continue;
		 }
		 if ( fabs(pfElectron->pt() - probe->pt())<delta && 
		      fabs(pfElectron->eta() - probe->eta())<delta &&
		      fabs(pfElectron->phi() - probe->phi())<delta ) {
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
	
	h_electronInvMass->Fill(e_ee_invMass);
	if (cont==2) h_electronInvMassPass->Fill(e_ee_invMass); 
	
     }  // if at least two electron
  } // if electronCollection is valid
  
  eventAccept->Fill(cont); 
  iEvent.put( pOutput );
  
}


// ------------ method called once each job just before starting event loop  ------------
void
pfAnalyzer::beginJob (){

//beginJob

}


// ------------ method called once each job just after ending the event loop  ------------
void
pfAnalyzer::endJob ()
{

//endJob

}

DEFINE_FWK_MODULE (pfAnalyzer);

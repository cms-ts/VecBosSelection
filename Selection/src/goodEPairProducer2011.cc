#include "VecBosSelection/Selection/interface/goodEPairProducer2011.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

using namespace std;
using namespace edm;
using namespace reco;


//
// member functions
//

// ------------ method called for each event  ------------
void
goodEPairProducer2011::produce(edm::Event & iEvent, edm::EventSetup const & iSetup)
{
  Debug2=false;   
   eleSelStepByStep->SetBinContent(2,eleSelStepByStep->GetBinContent(2)+1); //Number of events in which hlt has fired (1)
   // get gfs Electron  
   auto_ptr< reco::PFCandidateCollection > 
     pOutput( new reco::PFCandidateCollection ); 
   
   //working std::auto_ptr<std::vector<TLorentzVector*> > pOutputmu( new std::vector<TLorentzVector*>);
   
  if (Debug2) cout<<"------- NEW Event -----"<<endl;
  double delta = 0.001;
  int cont = 0;

  if (isElectron_){
    Handle < pat::ElectronCollection > electronCollection;  
    iEvent.getByLabel (theElectronCollectionLabel, electronCollection);
    
    
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
 
      pat::ElectronCollection::const_iterator highestptele;
      pat::ElectronCollection::const_iterator secondptele;
      
      double secondEleEnThrhold=secondEleEnThrhold_;
      double firstEleEnThrhold=firstEleEnThrhold_;  
      double lowZmassLimit=lowZmassLimit_;
      double highZmassLimit=highZmassLimit_;
      double maxEtaForElectron=maxEtaForElectron_;
      
      //Set of counters to follow the Z selection history  and form the plot labels..
      int twoEleGoodEtaCount=0;
      int twoEleHLTCount=0;
      int WP80Count=0;
      int lowThrholdCount=0;
      int isID=0; 
      int isIso=0;
      int isConv=0;
      
      int i=0;
      bool protection=false;
      int jj=0;
      int sizePat=electronCollection->size();
      protection=false;  
      passSelection=true;
      
      eleSelStepByStep->SetBinContent(3,eleSelStepByStep->GetBinContent(3)+1); // (1) + at least 2 ele in the event (2)
      
      /// NEW DS
      // Cutting on WP80
      for (pat::ElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
	
	protection=true;
	jj++;
	//Check whether the electron is within the acceptance
	if (fabs(recoElectron ->superCluster()->eta()) > maxEtaForElectron) {
	  continue;
	}
	twoEleGoodEtaCount++;
	if (twoEleGoodEtaCount==2) {
	  eleSelStepByStep->SetBinContent(4,eleSelStepByStep->GetBinContent(4)+1); //(2) + 2 ele whithin eta acceptance (3)
	}
	if (SelectionUtils::DoHLTMatch(recoElectron,iEvent)) twoEleHLTCount++;
	if (twoEleHLTCount==2) {
	  eleSelStepByStep->SetBinContent(5,eleSelStepByStep->GetBinContent(5)+1); //(3) + 2 ele HLT matched (4)
	  twoEleHLTCount++; //To avoid multiple insertion...
	}
	//Perform checks on each ele ID criteria
	// Here you get a plot full of information. Each electron contributes with one entry (so total numer of entries = 3* #electrons)
	// To have the "%", each bin value shold be divided by total numer of entries/3
	std::vector<bool> result=SelectionUtils::MakePfEleNewIDAnalysis(recoElectron,iEvent,useNewID_,doWP90_,conversions_h,beamSpot,vtx_h, isoVals);
	passIDEleCriteria->SetBinContent(1,passIDEleCriteria->GetBinContent(1)+1);
	if (result[0]) {
	  passIDEleCriteria->SetBinContent(2,passIDEleCriteria->GetBinContent(2)+1);
	  if (isIso<2) isIso++;
	}
	if (result[1]) {
	  passIDEleCriteria->SetBinContent(3,passIDEleCriteria->GetBinContent(3)+1);
	  if (isID<2) isID++;
	}
	if (result[2]) {
	  passIDEleCriteria->SetBinContent(4,passIDEleCriteria->GetBinContent(4)+1);
	  if (isConv<2) isConv++;
	}
	if (isID==2) eleSelStepByStep->SetBinContent(6,eleSelStepByStep->GetBinContent(6)+1); // (5) + 2 ele pt > lowTh (6)
	if (isID==2 && isIso==2) eleSelStepByStep->SetBinContent(7,eleSelStepByStep->GetBinContent(7)+1); // (5) + 2 ele pt > lowTh (6)
	if (isID==2 && isIso==2 && isConv==2) eleSelStepByStep->SetBinContent(8,eleSelStepByStep->GetBinContent(8)+1); // (5) + 2 ele pt > lowTh (6)
	
	if (result[0] && result[1] && result[2]) WP80Count++;
	if (WP80Count==2) eleSelStepByStep->SetBinContent(9,eleSelStepByStep->GetBinContent(9)+1); //(4) + 2 ele WP80 (5)
	
	if ( (!doID_ || ( (!useNewID_ && ((doWP90_ || SelectionUtils::DoWP80Pf(recoElectron,iEvent)) && 
					  (!doWP90_ || SelectionUtils::DoWP90Pf(recoElectron,iEvent))) )
			  || (useNewID_ && SelectionUtils::DoMedSel2011(recoElectron,iEvent,conversions_h,beamSpot,vtx_h))) )
	     && ( !doIsolation_ || SelectionUtils::DoIso2011(recoElectron, iEvent, isoVals))
	     && SelectionUtils::DoHLTMatch(recoElectron,iEvent) && recoElectron->pt()>secondEleEnThrhold){
	  lowThrholdCount++;
	  if (lowThrholdCount==2)
	    eleSelStepByStep->SetBinContent(10,eleSelStepByStep->GetBinContent(10)+1); // (5) + 2 ele pt > lowTh (6)
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
      
      if(i>1 && highestptele->pt()>=firstEleEnThrhold) { //you NEED at least two electrons :)
	eleSelStepByStep->SetBinContent(11,eleSelStepByStep->GetBinContent(11)+1); // (8) + 1 ele pt > LowPt and + 1 ele pt > highPt (7)	
	//Check if the charge are opposite..
	if(highestptele->charge() != secondptele->charge()) {
	  eleSelStepByStep->SetBinContent(12,eleSelStepByStep->GetBinContent(12)+1); // (9) + Opposite Charge
	  
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
	  if (e_ee_invMass>highZmassLimit || e_ee_invMass<lowZmassLimit) passSelection=false;
	  
	  
	  // searching for the corresponding GSF electron 
	  if (passSelection){	      
	    eleSelStepByStep->SetBinContent(13,eleSelStepByStep->GetBinContent(13)+1);
	    
	    if (i>2) h_zPt_3e->Fill(e_pair.Pt());
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
	  
	}
      } // if at least two electron
    } // if electronCollection is valid
  }
  else{
     //Handle < reco::CompositeCandidateCollection > ZmumuCandidates; 
     //iEvent.getByLabel (ZmumuCandidates_, ZmumuCandidates);
     //const reco::CompositeCandidateCollection & Zmm = *(ZmumuCandidates.product());
     // if (Zmm.size()==0) return;

    // //const reco::CompositeCandidateCollection::const_iterator zmmIter = Zmm.begin();
//     //const reco::CompositeCandidate zmm = *zmmIter;

//     const reco::Candidate* lep0 = Zmm[0].daughter(0);   
//     const reco::Candidate* lep1 = Zmm[0].daughter(1); 
//     const pat::Muon* muon0 = dynamic_cast<const pat::Muon*>(&(*lep0->masterClone()));     
//     const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(&(*lep1->masterClone()));

//     reco::PFCandidate m0(lep0->charge(), lep0->p4() , reco::PFCandidate::mu);
//     reco::PFCandidate m1(lep1->charge(), lep1->p4() , reco::PFCandidate::mu);

    Handle < pat::MuonCollection > muonCollection;
    iEvent.getByLabel (theMuCollectionLabel, muonCollection);
    if (muonCollection.isValid () && muonCollection->size()>1) {
   
       pat::MuonCollection::const_iterator highestptmu;
       pat::MuonCollection::const_iterator secondptmu;
  
       double lowZmassLimit=lowZmassLimit_;
       double highZmassLimit=highZmassLimit_;

       int i=0;

       bool protection=false;
       int jj=0;
       int sizePat=muonCollection->size();
       passSelection=true;
       
       /// NEW DS
       for (pat::MuonCollection::const_iterator recoElectron = muonCollection->begin (); recoElectron != muonCollection->end (); recoElectron++) {
	  protection=true;
	  jj++;
	  
	  if (i==0) highestptmu=recoElectron;
	  if (i==1){
	     if (highestptmu->pt()<recoElectron->pt()){
		secondptmu=highestptmu;
		highestptmu=recoElectron;
	     }
	     else{
		secondptmu=recoElectron;
	     }
	  }
	  if (i>1){
	     if (highestptmu->pt()<recoElectron->pt()){
		secondptmu=highestptmu;
		highestptmu=recoElectron;
	     }
	     else{
		if (secondptmu->pt()<recoElectron->pt()) secondptmu=recoElectron;
	     }
	  }
	  i++;
       }
       
       if (!protection) {
	  cout<<"size pat is "<<sizePat<<" while jj is "<<jj<<" and protection "<<protection<<endl;
	  cout<<"problems with PAT collection, in ZpatFilter.cc-->... Please check..."<<endl;    
	  passSelection=false;
       }

       if(i>1                                                  //you NEED at least two electrons :)
	  && highestptmu->charge() != secondptmu->charge()) {  //Check if the charge are opposite..

	  //Calculating Invariant Mass
	  TLorentzVector tagv;
	  tagv.SetPtEtaPhiM(highestptmu->pt(),highestptmu->eta(),highestptmu->phi(), highestptmu->mass());
	  TLorentzVector probev;
	  probev.SetPtEtaPhiM(secondptmu->pt(),secondptmu->eta(),secondptmu->phi(), secondptmu->mass());
	  
	  TLorentzVector mu_pair = tagv + probev;
	  double mu_invMass = mu_pair.M ();
	  
	  
	  //Cut on the tag and probe mass...
	  if (mu_invMass>highZmassLimit || mu_invMass<lowZmassLimit) passSelection=false;
	  
	  
	  // searching for the corresponding GSF electron 
	  if (passSelection){	      
	    
	    if (i>2) h_zPt_3e->Fill(mu_pair.Pt());

	    Handle<reco::PFCandidateCollection> pfMuCollection;
	    iEvent.getByLabel (pflowMuCollection_, pfMuCollection);
	    
	    for (reco::PFCandidateCollection::const_iterator pfMu = pfMuCollection->begin (); pfMu != pfMuCollection->end (); pfMu++) {
	       if (fabs(pfMu->pdgId())==13){
		  if ( fabs(pfMu->pt() - highestptmu->pt())<delta && 
		       fabs(pfMu->eta() - highestptmu->eta())<delta &&
		       fabs(pfMu->phi() - highestptmu->phi())<delta ) {
		     pOutput->push_back(*pfMu);
		     h_electronEn->Fill(pfMu->energy());
		     h_electronEta->Fill(pfMu->eta());
		     h_electronPt->Fill(pfMu->pt());
		     cont++;
		     continue;
		  }
		  if ( fabs(pfMu->pt() - secondptmu->pt())<delta && 
		       fabs(pfMu->eta() - secondptmu->eta())<delta &&
		       fabs(pfMu->phi() - secondptmu->phi())<delta ) {
		     pOutput->push_back(*pfMu);
		     h_electronEn->Fill(pfMu->energy());
		     h_electronEta->Fill(pfMu->eta());
		     h_electronPt->Fill(pfMu->pt());
		     cont++;
		     continue;
		  }
	       }
	    }
	  
	    //pOutput->push_back(m0);
	    //pOutput->push_back(m1);
	    /* it's working, in case it does not work
	       TLorentzVector *l0,*l1;
	       l0->SetPtEtaPhiM(muon0->pt(),muon0->eta(),muon0->phi(),muon0->mass());
	       l1->SetPtEtaPhiM(muon1->pt(),muon1->eta(),muon1->phi(),muon1->mass());
	       pOutputmu->push_back(l0);
	       pOutputmu->push_back(l1);
	    */
	  }  // passSelection

	  h_electronInvMass->Fill(mu_invMass);
	  if (cont==2) h_electronInvMassPass->Fill(mu_invMass);

       }  // opposite charge
    }    // valid && at least 2 
  }  // isElectron
  eventAccept->Fill(cont); 
  iEvent.put( pOutput );
  
}


// ------------ method called once each job just before starting event loop  ------------
void
goodEPairProducer2011::beginJob (){
  cout<<endl;
  cout<<"##############################"<<endl;
  cout<<"#   Good Epair Producer Parameters   #"<<endl;
  cout<<"##############################"<<endl;
  cout<<endl; 
  cout<<"Transverse Energy cut on first Ele="<<firstEleEnThrhold_<<"GeV, and "<<secondEleEnThrhold_<<"GeV on the second"<<endl;
  cout<<"Z invariant mass limit: low="<<lowZmassLimit_<<"GeV, high="<<highZmassLimit_<<"GeV"<<endl;
  cout<<"Electron max acceptance="<<maxEtaForElectron_<<endl;

  cout<<endl;
  eleSelStepByStep->GetXaxis()->SetBinLabel(1,"Total # of Events");
  eleSelStepByStep->GetXaxis()->SetBinLabel(2,"event HLT Fired");
  eleSelStepByStep->GetXaxis()->SetBinLabel(3,">= 2 ele");
  eleSelStepByStep->GetXaxis()->SetBinLabel(4,"2 ele <= eta Acceptance");
  eleSelStepByStep->GetXaxis()->SetBinLabel(5,"2 ele HLT matched");
  eleSelStepByStep->GetXaxis()->SetBinLabel(6,"2 ele WP80");
  eleSelStepByStep->GetXaxis()->SetBinLabel(7,"2 ele pt > lowPt");
  eleSelStepByStep->GetXaxis()->SetBinLabel(8,"ele pt > lowPt + pt > HighPt");
  passIDEleCriteria->GetXaxis()->SetBinLabel(1,"TotEle");
  passIDEleCriteria->GetXaxis()->SetBinLabel(2,"Isolated");
  passIDEleCriteria->GetXaxis()->SetBinLabel(3,"ID");
  passIDEleCriteria->GetXaxis()->SetBinLabel(4,"NotConverted");

//beginJob

}


// ------------ method called once each job just after ending the event loop  ------------
void
goodEPairProducer2011::endJob ()
{

//endJob

}

DEFINE_FWK_MODULE (goodEPairProducer2011);

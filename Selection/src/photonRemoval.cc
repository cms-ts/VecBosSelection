#include "VecBosSelection/Selection/interface/photonRemoval.h"

using namespace std;
using namespace edm;
using namespace reco;


typedef pair<const reco::GenParticle*, const reco::GenParticle*> eleGamma;
typedef vector<eleGamma> vecEleGamma;

//
// member functions
//

class GreaterPte{
public:
  bool operator()( const reco::GenParticle* & a, const reco::GenParticle* & b) {
     return (a->pt() > b->pt());
  }
};

class GreaterPt{
public:
  bool operator()( eleGamma & a, eleGamma & b) {
     return ( (a.first)->pt() > (b.first)->pt() );
  }
};

// ------------ method called for each event  ------------
void
photonRemoval::produce(edm::Event & iEvent, edm::EventSetup const & iSetup)
{
   
   // get gfs Electron  
   auto_ptr< reco::GenParticleCollection > 
      pOutput( new reco::GenParticleCollection ); 

   Handle < reco::GenParticleCollection > particles;  
   iEvent.getByLabel (particleCollectionLabel, particles);
   Handle <reco::GenParticleRefVector> bottom;  
   iEvent.getByLabel (bottomCollectionLabel, bottom); 
   vecEleGamma gammaCone;
   gammaCone.clear();
   vector<const reco::GenParticle*> eVector;
   eVector.clear();

   // if (particles.isValid ()){
//       for (reco::GenParticleCollection::const_iterator it = particles->begin (); it != particles->end (); it++) {
// 	 if (fabs(it->pdgId()==11) && it->status()==1 && searchHistory(&(*it)) ) eVector.push_back(e1_test);
//       }
//    }
   
   if (particles.isValid ()){
      for (reco::GenParticleCollection::const_iterator it = particles->begin (); it != particles->end (); it++) {
	 if (fabs(it->pdgId())==11 && it->status()==1){
	    // save the electrons
	    eVector.push_back(&(*it));
	 }
      }
   }

   //stable_sort(eVector.begin(), eVector.end(), GreaterPte());

   if (particles.isValid ()){
      for (reco::GenParticleCollection::const_iterator it = particles->begin (); it != particles->end (); it++) {
	 if (it->pdgId()==22 && it->status()==1){
	    // save the photons
	    for (vector<const reco::GenParticle*>::const_iterator itEl = eVector.begin(); itEl != eVector.end(); itEl++){
	       double deltaPhi = fabs(it->phi()- (*itEl)->phi());
	       if (deltaPhi > acos(-1)) deltaPhi= 2*acos(-1) - deltaPhi;
	       double deltaR = sqrt( deltaPhi*deltaPhi  + pow(it->eta()-(*itEl)->eta(),2) );
	       if ((*itEl)->eta()< 1.479 && deltaR < barrelRCone) {
		  pair<const reco::GenParticle*, const reco::GenParticle*> coppia;
		  coppia = make_pair (*itEl,&(*it));
		  gammaCone.push_back(coppia);
		  break;
	       } else if ((*itEl)->eta()>= 1.479 && (*itEl)->eta() < 3.0 && deltaR < endcapRCone){
		  pair<const reco::GenParticle*, const reco::GenParticle*> coppia;
		  coppia = make_pair (*itEl,&(*it));
		  gammaCone.push_back(coppia);
		  break;
	       } 
	    }
	 }
      }
   }   
   
   //stable_sort(gammaCone.begin(), gammaCone.end(), GreaterPt());
   //reco::GenParticle* ePlusTmp, eMinusTmp;
   double ePlusPt=-1, eMinusPt=-1;
   double ePx=0, ePy=0;
   int ePlusI=-1, eMinusI=-1, eCont=0;
   TLorentzVector ePlusP4, eMinusP4, zP4;

   for (vector<const reco::GenParticle*>::const_iterator itEl = eVector.begin(); itEl != eVector.end(); itEl++){
      ePx=(*itEl)->px(); ePy=(*itEl)->py();
      for (vecEleGamma::const_iterator it = gammaCone.begin(); it!= gammaCone.end(); it++){
	 if (it->first == *itEl ){
	    ePx += (it->second)->px();
	    ePy += (it->second)->py();
	 }
      }
      if ( (*itEl)->pdgId()==11 && sqrt( ePx*ePx  + ePy*ePy)> ePlusPt ){
	 ePlusPt = sqrt( ePx*ePx  + ePy*ePy);
	 //ePlusTmp = *itEl;
	 ePlusP4.SetPtEtaPhiM((*itEl)->pt(),(*itEl)->eta(),(*itEl)->phi(),(*itEl)->mass());
	 ePlusI = eCont;
      }
      if ( (*itEl)->pdgId()==-11 && sqrt( ePx*ePx  + ePy*ePy)> eMinusPt ){
	 eMinusPt = sqrt( ePx*ePx  + ePy*ePy);
	 //eMinusTmp = *itEl;
	 eMinusP4.SetPtEtaPhiM((*itEl)->pt(),(*itEl)->eta(),(*itEl)->phi(),(*itEl)->mass());
	 eMinusI = eCont;
      }
      eCont++;
   }

   bool someElectrons = false;
   if (ePlusI!=-1 && eMinusI!=-1){
      zP4 = ePlusP4 + eMinusP4;
      if (zP4.M() > minZMass && zP4.M() < maxZMass){
	 someElectrons = true;
      } else {
	 eNotRemovedMass->Fill(zP4.M());
	 eNotRemovedPt->Fill(ePlusP4.Pt());
	 eNotRemovedPt->Fill(eMinusP4.Pt());
	 eNotRemovedEta->Fill(ePlusP4.Eta());
	 eNotRemovedEta->Fill(eMinusP4.Eta());
      }
   }

	 
   bool isRemoval = false;
   for  (reco::GenParticleRefVector::const_iterator itB = bottom->begin (); itB != bottom->end (); itB++) {
      isRemoval = false;
      const reco::GenParticle* partB = itB->get();
      if (someElectrons){
	 const reco::GenParticle* eleP = eVector[ePlusI];
	 const reco::GenParticle* eleM = eVector[eMinusI];
	 // there are electron and gamma to remove
	 if (partB->pdgId()==22 && partB->status()==1){
	    // check on the photons
	    for (vecEleGamma::const_iterator itG = gammaCone.begin (); itG != gammaCone.end (); itG++){
	       if ((itG->first == eleP || itG->first == eleM) && partB == itG->second) {
		  isRemoval = true;
		  gammaRemovedPt->Fill(partB->pt());
		  gammaRemovedEta->Fill(partB->eta());
		  break;
	       }
	    }
	 }
	 if (fabs(partB->pdgId())==11 && partB->status()==1){
	       if (partB == eleP || partB == eleM) {
	    // check on the electronsif ((itG->first == eleP || itG->first == eleM) && partB == itG->second) {
		  isRemoval = true;
		  eRemovedPt->Fill(partB->pt());
		  eRemovedEta->Fill(partB->eta());
		  //break;
	       }
	 }
      } // end someElectrons
      if (isRemoval == false) pOutput->push_back(*partB);
   }

   iEvent.put( pOutput );
   
}


// ------------ method called once each job just before starting event loop  ------------
void
photonRemoval::beginJob (){

//beginJob

}


// ------------ method called once each job just after ending the event loop  ------------
void
photonRemoval::endJob ()
{

//endJob

}

DEFINE_FWK_MODULE (photonRemoval);

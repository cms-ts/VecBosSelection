// -*- C++ -*-
//
// Package:    Efficiency
// Class:      Efficiency
// 
/**\class Efficiency Efficiency.cc Zmonitoring/Efficiency/src/Efficiency.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Vieri Candelise & Matteo Marone
//         Created:  Wed May 11 14:53:26 CESDo2011



// system include files
#include <memory>

// user include files

#include "VecBosSelection/Selection/interface/EfficiencyFilter.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include <string>
#include <cmath>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <iostream>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "HLTrigger/HLTcore/interface/TriggerSummaryAnalyzerAOD.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


using namespace std;
using namespace edm;
using namespace reco;

bool Debug=false; //Activate with true if you wonna have verbosity for Debug

//
// member functions
//

// ------------ method called for each event  ------------
bool
EfficiencyFilter::filter (edm::Event & iEvent, edm::EventSetup const & iSetup)
{
  if (Debug) cout<<"------- NEW Event -----"<<endl;
  using namespace edm;


  // Pick up SUPERCLUSTERS:

  Handle<reco::SuperClusterCollection> superClusters_EB_h;
  iEvent.getByLabel(superClusterCollection_EB_,superClusters_EB_h );
  if ( ! superClusters_EB_h.isValid() ) return false;

  Handle<reco::SuperClusterCollection> superClusters_EE_h;
  iEvent.getByLabel(superClusterCollection_EE_,superClusters_EE_h );
  if ( ! superClusters_EE_h.isValid() ) return false;
  
  reco::SuperClusterCollection::const_iterator tag_SC;
  reco::SuperClusterCollection::const_iterator probe_SC;
  bool SC_isPassingProbe=false;
  int l=0;
  
  //EB superclusters:
  for (reco::SuperClusterCollection::const_iterator superCluster = superClusters_EB_h->begin(); superCluster != superClusters_EB_h->end(); superCluster++) {
    if ((superCluster->energy()/cosh(superCluster->eta())) > 20.0 && fabs(superCluster->eta())<=1.4442) {
      if (l==0) tag_SC=superCluster;
      if (l==1){
	if (tag_SC->energy()/cosh(tag_SC->eta()) < superCluster->energy()/cosh(superCluster->eta())){
          probe_SC=tag_SC;
          tag_SC=superCluster;
        }
        else{
          probe_SC=superCluster;
	}
      }
      if (l>1){
	if (tag_SC->energy()/cosh(tag_SC->eta()) < superCluster->energy()/cosh(superCluster->eta())){
          probe_SC=tag_SC;
          tag_SC=superCluster;
        }
        else{
          if (probe_SC->energy()/cosh(probe_SC->eta()) < superCluster->energy()/cosh(superCluster->eta())) probe_SC=superCluster;
	}
      }
      l++;
    }
  }
  //EE superclusters:  
  for (reco::SuperClusterCollection::const_iterator superCluster = superClusters_EE_h->begin(); superCluster != superClusters_EE_h->end(); superCluster++) {
    if ((superCluster->energy()/cosh(superCluster->eta())) > 20.0 && (fabs(superCluster->eta())>=1.5660 && fabs(superCluster->eta())<=2.4000)) {
      if (l==0) tag_SC=superCluster;
      if (l==1){
	if (tag_SC->energy()/cosh(tag_SC->eta()) < superCluster->energy()/cosh(superCluster->eta())){
          probe_SC=tag_SC;
          tag_SC=superCluster;
        }
        else{
          probe_SC=superCluster;
	}
      }
      if (l>1){
	if (tag_SC->energy()/cosh(tag_SC->eta()) < superCluster->energy()/cosh(superCluster->eta())){
          probe_SC=tag_SC;
          tag_SC=superCluster;
        }
        else{
          if (probe_SC->energy()/cosh(probe_SC->eta()) < superCluster->energy()/cosh(superCluster->eta())) probe_SC=superCluster;
	}
      }
      l++;
    }
  }

  if (l<2) return false;


  // Mixing TAG and PROBE SuperClusters

  if (RECO_efficiency_) {
    int SCchoice = rand()%2;
    if (SCchoice == 1) {
      reco::SuperClusterCollection::const_iterator temp_SC;
      temp_SC = tag_SC;
      tag_SC = probe_SC;
      probe_SC = temp_SC;
    }
  }


  // Pick up ELECTRONS:

  Handle < pat::ElectronCollection > electronCollection;
  iEvent.getByLabel (theElectronCollectionLabel, electronCollection);
  if (!electronCollection.isValid ()) return false;

  pat::ElectronCollection::const_iterator highestptele;
  pat::ElectronCollection::const_iterator secondptele;

  // iso deposits
  IsoDepositVals isoVals(isoValInputTags_.size());
  for (size_t j = 0; j < isoValInputTags_.size(); ++j) {
     iEvent.getByLabel(isoValInputTags_[j], isoVals[j]);
  }

  int i=0;
  if (electronCollection->size()<=1) return false;
  bool protection=false;
  double deltaPhi = 0.0;
  
  for (pat::ElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
    protection=true;
    if (recoElectron->pt()>20.0 && (fabs(recoElectron->superCluster()->eta())<1.4442 || (fabs(recoElectron->superCluster()->eta())>1.566 && fabs(recoElectron->superCluster()->eta())<2.4)) && (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_)) {
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
    if (recoElectron->pt()>20.0 && (fabs(recoElectron->superCluster()->eta())<1.4442 || (fabs(recoElectron->superCluster()->eta())>1.566 && fabs(recoElectron->superCluster()->eta())<2.4)) && RECO_efficiency_) {
      if (fabs(probe_SC->phi()-recoElectron->superCluster()->phi()) < 3.1416) {
	deltaPhi = probe_SC->phi()-recoElectron->superCluster()->phi();
      } else {
	deltaPhi = fabs(probe_SC->phi()-recoElectron->superCluster()->phi()) - 6.2832;
      }
      if (sqrt((probe_SC->eta()-recoElectron->superCluster()->eta())*(probe_SC->eta()-recoElectron->superCluster()->eta())+(deltaPhi*deltaPhi)) < 0.2) {
	SC_isPassingProbe=true;
	continue;
      }
      if (i==0) highestptele=recoElectron;
      if (i>0){
	if (highestptele->pt()<recoElectron->pt()){
	  highestptele=recoElectron;
	}
      }
      i++;
    }
  }
  if (!protection) {
    cout<<"problems with PAT collection, in Efficiency.cc-->... Please check..."<<endl;    
    return false;
  }
  
  if (Debug) cout<<"Out of the SC and RECOele loops!"<<endl;  

  if(i<1) return false;
  if((WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) && (i<2 || highestptele->charge() == secondptele->charge())) return false;

// MC match

  if (matchMC_) {

    Handle<reco::GenParticleCollection> genPart;
    iEvent.getByLabel(genParticleCollection_, genPart);

    if (genPart.isValid()) {

      int nMatch = 0;

      for  (reco::GenParticleCollection::const_iterator itgen=genPart->begin(); itgen!=genPart->end(); itgen++) {

        if (itgen->status() == 3 && fabs(itgen->pdgId()) == 11 && itgen->mother()->pdgId() == 23) {

        double deltaPhi1, deltaR1;
          if (fabs(itgen->phi()-highestptele->phi()) < 3.1416) {
            deltaPhi1 = itgen->phi()-highestptele->phi();
          } else {
            deltaPhi1 = fabs(itgen->phi()-highestptele->phi()) - 6.2832;
          }
          deltaR1 = sqrt( pow(itgen->eta()-highestptele->eta(), 2) + deltaPhi1*deltaPhi1 );

          double deltaPhi2, deltaR2;
          if ( RECO_efficiency_) {
            if (fabs(itgen->phi()-probe_SC->phi()) < 3.1416) {
              deltaPhi2 = itgen->phi()-probe_SC->phi();
            } else {
              deltaPhi2 = fabs(itgen->phi()-probe_SC->phi()) - 6.2832;
            }
            deltaR2 = sqrt( pow(itgen->eta()-probe_SC->eta(), 2) + deltaPhi2*deltaPhi2 );
          } else {
            if (fabs(itgen->phi()-secondptele->phi()) < 3.1416) {
              deltaPhi2 = itgen->phi()-secondptele->phi();
            } else {
              deltaPhi2 = fabs(itgen->phi()-secondptele->phi()) - 6.2832;
            }
            deltaR2 = sqrt( pow(itgen->eta()-secondptele->eta(), 2) + deltaPhi2*deltaPhi2 );
          }

          if ( deltaR1 < 0.1 && highestptele->charge() == itgen->charge() ) {
            nMatch++;
            if (Debug) std::cout << "match 1" << std::endl;
          } else if ( deltaR2 < 0.1 && RECO_efficiency_ ) {
            if (Debug) std::cout << "match 2" << std::endl;
            nMatch++;
          } else if ( deltaR2 < 0.1 && !RECO_efficiency_ && secondptele->charge() == itgen->charge() ) {
            if (Debug) std::cout << "match 2" << std::endl;
            nMatch++;
          } else {
            if (Debug) std::cout << "no match" << std::endl;
          }

        }

      }

      if (Debug) std::cout << "# match = " << nMatch << std::endl;

      if (nMatch > 2) std::cout << "TOO MANY MATCHED ELECTRONS" << std::endl;

      if (nMatch != 2) return false;

    }

  }

// MC match


  // TRIGGER MATCHING (set flags to TRUE if matched):

  bool ProbeHLTmatch_highestptele = false;
  bool ProbeHLTmatch_secondptele = false;
  bool TagHLTmatch_highestptele = false;
  bool TagHLTmatch_secondptele = false;

  double deltaReles = 999.0;

  Handle < edm::RefToBaseVector<reco::GsfElectron> > TagHLTelectronCollection;
  iEvent.getByLabel (theTagHLTElectronCollectionLabel, TagHLTelectronCollection);
  if (!TagHLTelectronCollection.isValid ()) return false;

  for (edm::RefToBaseVector<reco::GsfElectron>::const_iterator TagHLTElectron = TagHLTelectronCollection->begin (); TagHLTElectron != TagHLTelectronCollection->end (); TagHLTElectron++) {

    if ((*TagHLTElectron)->pt()>20.0) {

      if (fabs(highestptele->phi()-(*TagHLTElectron)->phi()) < 3.1416) {
	deltaPhi = highestptele->phi()-(*TagHLTElectron)->phi();
      } else {
	deltaPhi = fabs(highestptele->phi()-(*TagHLTElectron)->phi()) - 6.2832;
      }
      deltaReles = sqrt((highestptele->eta()-(*TagHLTElectron)->eta())*
			(highestptele->eta()-(*TagHLTElectron)->eta())+
			(deltaPhi*deltaPhi));
      if (deltaReles < 0.2) TagHLTmatch_highestptele = true;
      
      if (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) {
	if (fabs(secondptele->phi()-(*TagHLTElectron)->phi()) < 3.1416) {
	  deltaPhi = secondptele->phi()-(*TagHLTElectron)->phi();
	} else {
	  deltaPhi = fabs(secondptele->phi()-(*TagHLTElectron)->phi()) - 6.2832;
	}
	deltaReles = sqrt((secondptele->eta()-(*TagHLTElectron)->eta())*
			  (secondptele->eta()-(*TagHLTElectron)->eta())+
			  (deltaPhi*deltaPhi));
	if (deltaReles < 0.2) TagHLTmatch_secondptele = true;
      }
    }
  }


  Handle < edm::RefToBaseVector<reco::GsfElectron> > ProbeHLTelectronCollection;
  iEvent.getByLabel (theProbeHLTElectronCollectionLabel, ProbeHLTelectronCollection);
  if (!ProbeHLTelectronCollection.isValid ()) return false;
  
  if (HLTele17_efficiency_ || HLTele8_efficiency_) {
    for (edm::RefToBaseVector<reco::GsfElectron>::const_iterator ProbeHLTElectron = ProbeHLTelectronCollection->begin (); ProbeHLTElectron != ProbeHLTelectronCollection->end (); ProbeHLTElectron++) {

      if ((*ProbeHLTElectron)->pt()>20.0) {
	
	if (fabs(highestptele->phi()-(*ProbeHLTElectron)->phi()) < 3.1416) {
	  deltaPhi = highestptele->phi()-(*ProbeHLTElectron)->phi();
	} else {
	  deltaPhi = fabs(highestptele->phi()-(*ProbeHLTElectron)->phi()) - 6.2832;
	}
	deltaReles = sqrt((highestptele->eta()-(*ProbeHLTElectron)->eta())*
			  (highestptele->eta()-(*ProbeHLTElectron)->eta())+
			  (deltaPhi*deltaPhi));
	if (deltaReles < 0.2) ProbeHLTmatch_highestptele = true;
	
	if (fabs(secondptele->phi()-(*ProbeHLTElectron)->phi()) < 3.1416) {
	  deltaPhi = secondptele->phi()-(*ProbeHLTElectron)->phi();
	} else {
	  deltaPhi = fabs(secondptele->phi()-(*ProbeHLTElectron)->phi()) - 6.2832;
	}
	deltaReles = sqrt((secondptele->eta()-(*ProbeHLTElectron)->eta())*
			  (secondptele->eta()-(*ProbeHLTElectron)->eta())+
			  (deltaPhi*deltaPhi));
	if (deltaReles < 0.2) ProbeHLTmatch_secondptele = true;
      }
    }
  }


  //Calculating tag & probe stuff
  TLorentzVector tagv;
  tagv.SetPtEtaPhiM(highestptele->pt(),highestptele->eta(),highestptele->phi(), 0.0);
  TLorentzVector probev;
  if (RECO_efficiency_)  {
    probev.SetXYZT(probe_SC->energy()*cos(probe_SC->phi())/cosh(probe_SC->eta()),
		   probe_SC->energy()*sin(probe_SC->phi())/cosh(probe_SC->eta()),
		   probe_SC->energy()*tanh(probe_SC->eta()),
		   probe_SC->energy());
  } else {
    probev.SetPtEtaPhiM(secondptele->pt(),secondptele->eta(),secondptele->phi(), 0.0);
  }


  TLorentzVector e_pair = tagv + probev;
  double e_ee_invMass = e_pair.M ();

  //cut on the tag and probe mass...
  if (e_ee_invMass>120.0 || e_ee_invMass<60.0) return false;


  //Estraggo il contenuto dei jets
  int nJet = 0;
  double leadingJet_pT=-1.0;
  double subleadingJet_pT=-1.0;
  double subsubleadingJet_pT=-1.0;
  double subsubsubleadingJet_pT=-1.0;
  double deltaR1jet = 999.0;
  double deltaR2jet = 999.0;

  Handle<PFJetCollection> pfJets;
  iEvent.getByLabel(theJetCollectionLabel_, pfJets);
  if (pfJets.isValid()) {
    for(reco::PFJetCollection::const_iterator jet = pfJets->begin();jet != pfJets->end();jet++) { 

      if (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) {
	if (fabs(jet->phi()-secondptele->superCluster()->phi()) < 3.1416) {
	  deltaPhi = jet->phi()-secondptele->superCluster()->phi();
	} else {
	  deltaPhi = fabs(jet->phi()-secondptele->superCluster()->phi()) - 6.2832;
	}
	deltaR1jet = sqrt( pow(jet->eta()-secondptele->superCluster()->eta(),2)+pow(deltaPhi,2) );
      } else {
	if (fabs(jet->phi()-probe_SC->phi()) < 3.1416) {
	  deltaPhi = jet->phi()-probe_SC->phi();
	} else {
	  deltaPhi = fabs(jet->phi()-probe_SC->phi()) - 6.2832;
	}
	deltaR1jet = sqrt( pow(jet->eta()-probe_SC->eta(),2)+pow(deltaPhi,2) );
      }
      if (fabs(jet->phi()-highestptele->superCluster()->phi()) < 3.1416) {
	deltaPhi = jet->phi()-highestptele->superCluster()->phi();
      } else {
	deltaPhi = fabs(jet->phi()-highestptele->superCluster()->phi()) - 6.2832;
      }
      deltaR2jet= sqrt( pow(jet->eta()-highestptele->superCluster()->eta(),2)+pow(deltaPhi,2) );

      if (deltaR1jet > 0.3 && deltaR2jet > 0.3 && fabs(jet->eta())<2.4 
	  && jet->pt()>30
	  && jet->chargedEmEnergyFraction() < 0.99
	  && jet->neutralHadronEnergyFraction() < 0.99
	  && jet->neutralEmEnergyFraction() < 0.99
	  && jet->chargedHadronEnergyFraction() > 0.0
	  && jet->chargedMultiplicity() > 0
	  ) {
	nJet++;
	if (jet->pt() > leadingJet_pT && jet->pt() > subleadingJet_pT && jet->pt() > subsubleadingJet_pT && jet->pt() > subsubsubleadingJet_pT){
	  subsubsubleadingJet_pT = subsubleadingJet_pT;
	  subsubleadingJet_pT = subleadingJet_pT;
	  subleadingJet_pT = leadingJet_pT;
	  leadingJet_pT = jet->pt();
	} else if (jet->pt() < leadingJet_pT && jet->pt() > subleadingJet_pT && jet->pt() > subsubleadingJet_pT && jet->pt() > subsubsubleadingJet_pT){
	  subsubsubleadingJet_pT = subsubleadingJet_pT;
	  subsubleadingJet_pT = subleadingJet_pT;
	  subleadingJet_pT = jet->pt();
	} else if (jet->pt() < leadingJet_pT && jet->pt() < subleadingJet_pT && jet->pt() > subsubleadingJet_pT && jet->pt() > subsubsubleadingJet_pT){
	  subsubsubleadingJet_pT = subsubleadingJet_pT;
	  subsubleadingJet_pT = jet->pt();
	} else if (jet->pt() < leadingJet_pT && jet->pt() < subleadingJet_pT && jet->pt() < subsubleadingJet_pT && jet->pt() > subsubsubleadingJet_pT) {
          subsubsubleadingJet_pT = jet->pt();
	}
	if (Debug) cout<<"Jet eta "<<jet->eta()<<" pt "<<jet->pt()<<endl;
      }
    }
  }
  else{
    cout<<"No valid Jets Collection"<<endl;
  }
  
  if (Debug) cout<<"This event has jets #->"<<nJet<<endl;
  
  
  // Filling TAP distributions:

  // WP80 & HLT:

  if (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) {
    // 1st leg WP80 efficiency
    if (SelectionUtils::DoWP80Pf(highestptele,iEvent) && SelectionUtils::DoIso2011(highestptele, iEvent, isoVals) && TagHLTmatch_highestptele){
      probeall_pt->Fill(secondptele->pt());
      probeall_eta->Fill(secondptele->superCluster()->eta());
      probeall_mee->Fill(e_ee_invMass);
      probeall_leadjetpt->Fill(leadingJet_pT);
      if ((WP80_efficiency_ && SelectionUtils::DoWP90Pf(secondptele,iEvent) && SelectionUtils::DoIso2011(secondptele, iEvent, isoVals)) || 
	  ((HLTele17_efficiency_ || HLTele8_efficiency_) && ProbeHLTmatch_secondptele)){
	probepass_pt->Fill(secondptele->pt());
	probepass_eta->Fill(secondptele->superCluster()->eta());
	probepass_mee->Fill(e_ee_invMass);
	if (nJet==0) probepass0jet->Fill(e_ee_invMass);
	if (nJet==1) probepass1jet->Fill(e_ee_invMass);
	if (nJet==2) probepass2jet->Fill(e_ee_invMass);
	if (nJet==3) probepass3jet->Fill(e_ee_invMass);
	if (nJet==4) probepass4jet->Fill(e_ee_invMass);
	if (nJet>4)  probepass5jet->Fill(e_ee_invMass);
	if (nJet > 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 1) {
	  if (subleadingJet_pT >=  30.0  && subleadingJet_pT <  40.0) probepass0subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  40.0  && subleadingJet_pT <  50.0) probepass1subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  50.0  && subleadingJet_pT <  70.0) probepass2subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  70.0  && subleadingJet_pT <  90.0) probepass3subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  90.0  && subleadingJet_pT < 150.0) probepass4subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >= 150.0) probepass5subleadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 2) {
	  if (subsubleadingJet_pT >=  30.0  && subsubleadingJet_pT <  50.0) probepass0subsubleadjetpt->Fill(e_ee_invMass);
	  if (subsubleadingJet_pT >=  50.0  && subsubleadingJet_pT <  150.0) probepass1subsubleadjetpt->Fill(e_ee_invMass);
	  if (subsubleadingJet_pT >= 150.0) probepass2subsubleadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 3) {
	  probepass0subsubsubleadjetpt->Fill(e_ee_invMass);
	}
      }
      else{
	probefail_pt->Fill(secondptele->pt());
	probefail_eta->Fill(secondptele->superCluster()->eta());
	probefail_mee->Fill(e_ee_invMass);
	if (nJet==0) probefail0jet->Fill(e_ee_invMass);
	if (nJet==1) probefail1jet->Fill(e_ee_invMass);
	if (nJet==2) probefail2jet->Fill(e_ee_invMass);
	if (nJet==3) probefail3jet->Fill(e_ee_invMass);
	if (nJet==4) probefail4jet->Fill(e_ee_invMass);
	if (nJet>4)  probefail5jet->Fill(e_ee_invMass);
	if (nJet > 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 1) {
	  if (subleadingJet_pT >=  30.0  && subleadingJet_pT <  40.0) probefail0subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  40.0  && subleadingJet_pT <  50.0) probefail1subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  50.0  && subleadingJet_pT <  70.0) probefail2subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  70.0  && subleadingJet_pT <  90.0) probefail3subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  90.0  && subleadingJet_pT < 150.0) probefail4subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >= 150.0) probefail5subleadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 2) {
	  if (subsubleadingJet_pT >=  30.0  && subsubleadingJet_pT <  50.0) probefail0subsubleadjetpt->Fill(e_ee_invMass);
	  if (subsubleadingJet_pT >=  50.0  && subsubleadingJet_pT <  150.0) probefail1subsubleadjetpt->Fill(e_ee_invMass);
	  if (subsubleadingJet_pT >= 150.0) probefail2subsubleadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 3) {
	  probefail0subsubsubleadjetpt->Fill(e_ee_invMass);
	}
      }
    }  
    if (SelectionUtils::DoWP80Pf(secondptele,iEvent) && SelectionUtils::DoIso2011(secondptele, iEvent, isoVals) && TagHLTmatch_secondptele){
      tagall_pt->Fill(highestptele->pt());
      tagall_eta->Fill(highestptele->superCluster()->eta());
      tagall_mee->Fill(e_ee_invMass);
      tagall_leadjetpt->Fill(leadingJet_pT);
      if ( (WP80_efficiency_ && SelectionUtils::DoWP90Pf(highestptele,iEvent) && SelectionUtils::DoIso2011(highestptele, iEvent, isoVals)) || 
	   ((HLTele17_efficiency_ || HLTele8_efficiency_) && ProbeHLTmatch_highestptele)){
	tagpass_pt->Fill(highestptele->pt());
	tagpass_eta->Fill(highestptele->superCluster()->eta());
	tagpass_mee->Fill(e_ee_invMass);
	if (nJet==0) tagpass0jet->Fill(e_ee_invMass);
	if (nJet==1) tagpass1jet->Fill(e_ee_invMass);
	if (nJet==2) tagpass2jet->Fill(e_ee_invMass);
	if (nJet==3) tagpass3jet->Fill(e_ee_invMass);
	if (nJet==4) tagpass4jet->Fill(e_ee_invMass);
	if (nJet>4)  tagpass5jet->Fill(e_ee_invMass);
	if (nJet > 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagpass0leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagpass1leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagpass2leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagpass3leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagpass4leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagpass5leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagpass6leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagpass7leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagpass8leadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 1) {
	  if (subleadingJet_pT >=  30.0  && subleadingJet_pT <  40.0) tagpass0subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  40.0  && subleadingJet_pT <  50.0) tagpass1subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  50.0  && subleadingJet_pT <  70.0) tagpass2subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  70.0  && subleadingJet_pT <  90.0) tagpass3subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  90.0  && subleadingJet_pT < 150.0) tagpass4subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >= 150.0) tagpass5subleadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 2) {
	  if (subsubleadingJet_pT >=  30.0  && subsubleadingJet_pT <  50.0) tagpass0subsubleadjetpt->Fill(e_ee_invMass);
	  if (subsubleadingJet_pT >=  50.0  && subsubleadingJet_pT <  150.0) tagpass1subsubleadjetpt->Fill(e_ee_invMass);
	  if (subsubleadingJet_pT >= 150.0) tagpass2subsubleadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 3) {
	  tagpass0subsubsubleadjetpt->Fill(e_ee_invMass);
	}
      }
      else{
	tagfail_pt->Fill(highestptele->pt());
	tagfail_eta->Fill(highestptele->superCluster()->eta());
	tagfail_mee->Fill(e_ee_invMass);
	if (nJet==0) tagfail0jet->Fill(e_ee_invMass);
	if (nJet==1) tagfail1jet->Fill(e_ee_invMass);
	if (nJet==2) tagfail2jet->Fill(e_ee_invMass);
	if (nJet==3) tagfail3jet->Fill(e_ee_invMass);
	if (nJet==4) tagfail4jet->Fill(e_ee_invMass);
	if (nJet>4)  tagfail5jet->Fill(e_ee_invMass);
	if (nJet > 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) tagfail0leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) tagfail1leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) tagfail2leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) tagfail3leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) tagfail4leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) tagfail5leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) tagfail6leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) tagfail7leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) tagfail8leadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 1) {
	  if (subleadingJet_pT >=  30.0  && subleadingJet_pT <  40.0) tagfail0subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  40.0  && subleadingJet_pT <  50.0) tagfail1subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  50.0  && subleadingJet_pT <  70.0) tagfail2subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  70.0  && subleadingJet_pT <  90.0) tagfail3subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  90.0  && subleadingJet_pT < 150.0) tagfail4subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >= 150.0) tagfail5subleadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 2) {
	  if (subsubleadingJet_pT >=  30.0  && subsubleadingJet_pT <  50.0) tagfail0subsubleadjetpt->Fill(e_ee_invMass);
	  if (subsubleadingJet_pT >=  50.0  && subsubleadingJet_pT <  150.0) tagfail1subsubleadjetpt->Fill(e_ee_invMass);
	  if (subsubleadingJet_pT >= 150.0) tagfail2subsubleadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 3) {
	  tagfail0subsubsubleadjetpt->Fill(e_ee_invMass);
	}
      }
    }
  }

  // RECO:

  if (RECO_efficiency_) {
    if (SelectionUtils::DoWP80Pf(highestptele,iEvent) && SelectionUtils::DoIso2011(highestptele, iEvent, isoVals) && TagHLTmatch_highestptele){
      probeall_pt->Fill(probe_SC->energy());
      probeall_eta->Fill(probe_SC->eta());
      probeall_mee->Fill(e_ee_invMass);
      probeall_leadjetpt->Fill(leadingJet_pT);
      if (SC_isPassingProbe){
	probepass_pt->Fill(probe_SC->energy());
	probepass_eta->Fill(probe_SC->eta());
	probepass_mee->Fill(e_ee_invMass);
	if (nJet==0) probepass0jet->Fill(e_ee_invMass);
	if (nJet==1) probepass1jet->Fill(e_ee_invMass);
	if (nJet==2) probepass2jet->Fill(e_ee_invMass);
	if (nJet==3) probepass3jet->Fill(e_ee_invMass);
	if (nJet==4) probepass4jet->Fill(e_ee_invMass);
	if (nJet>4)  probepass5jet->Fill(e_ee_invMass);
	if (nJet > 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probepass0leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probepass1leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probepass2leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probepass3leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probepass4leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probepass5leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probepass6leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probepass7leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probepass8leadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 1) {
	  if (subleadingJet_pT >=  30.0  && subleadingJet_pT <  40.0) probepass0subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  40.0  && subleadingJet_pT <  50.0) probepass1subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  50.0  && subleadingJet_pT <  70.0) probepass2subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  70.0  && subleadingJet_pT <  90.0) probepass3subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  90.0  && subleadingJet_pT < 150.0) probepass4subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >= 150.0) probepass5subleadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 2) {
	  if (subsubleadingJet_pT >=  30.0  && subsubleadingJet_pT <  50.0) probepass0subsubleadjetpt->Fill(e_ee_invMass);
	  if (subsubleadingJet_pT >=  50.0  && subsubleadingJet_pT <  150.0) probepass1subsubleadjetpt->Fill(e_ee_invMass);
	  if (subsubleadingJet_pT >= 150.0) probepass2subsubleadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 3) {
	  probepass0subsubsubleadjetpt->Fill(e_ee_invMass);
	}
      }
      else{
	probefail_pt->Fill(probe_SC->energy());
	probefail_eta->Fill(probe_SC->eta());
	probefail_mee->Fill(e_ee_invMass);
	if (nJet==0) probefail0jet->Fill(e_ee_invMass);
	if (nJet==1) probefail1jet->Fill(e_ee_invMass);
	if (nJet==2) probefail2jet->Fill(e_ee_invMass);
	if (nJet==3) probefail3jet->Fill(e_ee_invMass);
	if (nJet==4) probefail4jet->Fill(e_ee_invMass);
	if (nJet>4)  probefail5jet->Fill(e_ee_invMass);
	if (nJet > 0) {
	  if (leadingJet_pT >=  30.0  && leadingJet_pT <  40.0) probefail0leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  40.0  && leadingJet_pT <  50.0) probefail1leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  50.0  && leadingJet_pT <  70.0) probefail2leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  70.0  && leadingJet_pT <  90.0) probefail3leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >=  90.0  && leadingJet_pT < 120.0) probefail4leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 120.0  && leadingJet_pT < 150.0) probefail5leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 150.0  && leadingJet_pT < 190.0) probefail6leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 190.0  && leadingJet_pT < 230.0) probefail7leadjetpt->Fill(e_ee_invMass);
	  if (leadingJet_pT >= 230.0) probefail8leadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 1) {
	  if (subleadingJet_pT >=  30.0  && subleadingJet_pT <  40.0) probefail0subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  40.0  && subleadingJet_pT <  50.0) probefail1subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  50.0  && subleadingJet_pT <  70.0) probefail2subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  70.0  && subleadingJet_pT <  90.0) probefail3subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >=  90.0  && subleadingJet_pT < 150.0) probefail4subleadjetpt->Fill(e_ee_invMass);
	  if (subleadingJet_pT >= 150.0) probefail5subleadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 2) {
	  if (subsubleadingJet_pT >=  30.0  && subsubleadingJet_pT <  50.0) probefail0subsubleadjetpt->Fill(e_ee_invMass);
	  if (subsubleadingJet_pT >=  50.0  && subsubleadingJet_pT <  150.0) probefail1subsubleadjetpt->Fill(e_ee_invMass);
	  if (subsubleadingJet_pT >= 150.0) probefail2subsubleadjetpt->Fill(e_ee_invMass);
	}
	if (nJet > 3) {
	  probefail0subsubsubleadjetpt->Fill(e_ee_invMass);
	}
      }
    }  
  }
  
  return true;
  
}

  
// ------------ method called once each job just before starting event loop  ------------
void
EfficiencyFilter::beginJob (){

}

// ------------ method called when starting to processes a run  ------------
bool
EfficiencyFilter::beginRun(edm::Run &iRun, edm::EventSetup const& iSetup)
{
//HLT names
   std::vector<std::string>  hlNames;
   bool changed (true);
   if (hltConfig_.init(iRun,iSetup,triggerCollection_.process(),changed)) {
     if (changed) {
       hlNames = hltConfig_.triggerNames();
     }
   } else {
     edm::LogError("MyAnalyzer") << " HLT config extraction failure with process name " << triggerCollection_.process();
   }
   if (Debug) cout<<"useAllTriggers?"<<useAllTriggers_<<endl;
   if(useAllTriggers_) triggerNames_ = hlNames;
   //triggerNames_ = hlNames;
   //HLT indices
   triggerIndices_.clear();
   for(unsigned int itrig = 0; itrig < triggerNames_.size(); ++itrig) {
     if(find(hlNames.begin(),hlNames.end(),triggerNames_[itrig])!=hlNames.end())
       triggerIndices_.push_back(hltConfig_.triggerIndex(triggerNames_[itrig]));
     else
       triggerIndices_.push_back(2048);
   }

   if (Debug){
     // text (Debug) output
     int i=0;
     for(std::vector<std::string>::const_iterator it = triggerNames_.begin(); it<triggerNames_.end();++it) {
       std::cout << (i++) << " = " << (*it) << std::endl;
     } 
   }
   return true;
}



// ------------ method called once each job just after ending the event loop  ------------
void
EfficiencyFilter::endJob ()
{

//endJob

}

DEFINE_FWK_MODULE (EfficiencyFilter);

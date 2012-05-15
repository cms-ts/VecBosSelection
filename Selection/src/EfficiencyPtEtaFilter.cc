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
// Original Author:  superben
//         Created:  Wed May 11 14:53:26 CESDo2011
// $Id: EfficiencyFilter.cc,v 1.27 2012/05/12 11:34:59 schizzi Exp $



// system include files
#include <memory>

// user include files

#include "VecBosSelection/Selection/interface/EfficiencyPtEtaFilter.h"
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

bool Debug_flag=false
; //Activate with true if you wonna have verbosity for Debug

//
// member functions
//

// ------------ method called for each event  ------------
bool
EfficiencyPtEtaFilter::filter (edm::Event & iEvent, edm::EventSetup const & iSetup)
{
  if (Debug_flag) cout<<"------- NEW Event -----"<<endl;
  using namespace edm;

  if (Debug_flag) cout<<"WP80_efficiency_ = "<<WP80_efficiency_<<endl;  
  if (Debug_flag) cout<<"HLTele17_efficiency_ = "<<HLTele17_efficiency_<<endl;  
  if (Debug_flag) cout<<"HLTele8_efficiency_ = "<<HLTele8_efficiency_<<endl;  
  if (Debug_flag) cout<<"RECO_efficiency_ = "<<RECO_efficiency_<<endl;  

  // Pick up SUPERCLUSTERS:

  Handle<reco::SuperClusterCollection> superClusters_EB_h;
  iEvent.getByLabel(superClusterCollection_EB_,superClusters_EB_h );
  if ( ! superClusters_EB_h.isValid() ) return false;

  Handle<reco::SuperClusterCollection> superClusters_EE_h;
  iEvent.getByLabel(superClusterCollection_EE_,superClusters_EE_h );
  if ( ! superClusters_EE_h.isValid() ) return false;
  
  reco::SuperClusterCollection::const_iterator probe_SC;
  bool SC_isPassingProbe=false;
  
  int l=0;
  
  //EB superclusters:
  for (reco::SuperClusterCollection::const_iterator superCluster = superClusters_EB_h->begin(); superCluster != superClusters_EB_h->end(); superCluster++) {
    if (superCluster->energy()>20.0 && fabs(superCluster->eta())<=1.4442) {
      //      sc_counter++;
      if (l==0) { probe_SC=superCluster; l++; break;}
    }
  }

  //EE superclusters:  
  for (reco::SuperClusterCollection::const_iterator superCluster = superClusters_EE_h->begin(); superCluster != superClusters_EE_h->end(); superCluster++) {
    if (superCluster->energy()>20.0 && (fabs(superCluster->eta())>=1.5660 && fabs(superCluster->eta())<=2.4000)) {
      //  sc_counter++;
      if (l==0) { probe_SC=superCluster; l++; break;}
      if (l==1 && (superCluster->energy() > probe_SC->energy())) { probe_SC=superCluster; l++; break;}
    }
  }


  if (l<1) return false;


  // Pick up ELECTRONS:

  Handle < pat::ElectronCollection > electronCollection;
  iEvent.getByLabel (theElectronCollectionLabel, electronCollection);
  if (!electronCollection.isValid ()) return false;

  pat::ElectronCollection::const_iterator firstele;
  pat::ElectronCollection::const_iterator secondele;

  int i=0;
  if (electronCollection->size()<=1) return false;
  bool protection=false;
  
  for (pat::ElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
    protection=true;
    if (recoElectron->pt()>20.0 && recoElectron->superCluster()->eta()<2.4 && (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_)) {
      if (i==0 && SelectionUtils::DoWP80Pf(recoElectron,iEvent,removePU_)) { firstele=recoElectron; i++;}
      if (i==1 && !(firstele->charge() == recoElectron->charge())) { secondele=recoElectron; i++; break;}
    }
    if (recoElectron->pt()>20.0 && recoElectron->superCluster()->eta()<2.4 && RECO_efficiency_) {
      if (sqrt((probe_SC->eta()-recoElectron->superCluster()->eta())*(probe_SC->eta()-recoElectron->superCluster()->eta())+(probe_SC->phi()-recoElectron->phi())*(probe_SC->phi()-recoElectron->phi())) < 0.2) {
	SC_isPassingProbe=true;
	continue;
      }
      if (i==0 && SelectionUtils::DoWP80Pf(recoElectron,iEvent,removePU_)) { firstele=recoElectron; i++;}
    }
  }
  if (!protection) {
    cout<<"problems with PAT collection, in Efficiency.cc-->... Please check..."<<endl;    
    return false;
  }
  
  //  if (Debug_flag) cout<<"firstele->pt() = "<<firstele->pt() <<endl;  
  if (Debug_flag) cout<<"Out of the SC and RECOele loops!"<<endl;  

  if (i<1) return false;
  if (i<2 && (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_)) return false;

  if (Debug_flag) cout<<"l = "<<l<<"; i = "<<i<<"; SC_isPassingProbe  = "<<SC_isPassingProbe<<endl;  
  if (Debug_flag) cout<<"-------------" <<endl;  
  if (Debug_flag) cout<<"firstele->pt() = "<<firstele->pt() <<endl;  
  if (Debug_flag) cout<<"firstele->superCluster()->eta() = "<<firstele->superCluster()->eta() <<endl;  
  if (Debug_flag) cout<<"firstele->phi() = "<<firstele->phi() <<endl;  
  if (Debug_flag) cout<<"-------------" <<endl;  
  if (RECO_efficiency_)  {
    if (Debug_flag) cout<<"-------------" <<endl;  
    if (Debug_flag) cout<<"probe_SC->energy() = "<<probe_SC->energy() <<endl;  
    if (Debug_flag) cout<<"probe_SC->eta() = "<<probe_SC->eta() <<endl;  
    if (Debug_flag) cout<<"probe_SC->phi() = "<<probe_SC->phi() <<endl;  
    if (Debug_flag) cout<<"-------------" <<endl;  
  } else {
    if (Debug_flag) cout<<"-------------" <<endl;  
    if (Debug_flag) cout<<"secondele->pt() = "<<secondele->pt() <<endl;  
    if (Debug_flag) cout<<"secondele->superCluster()->eta() = "<<secondele->superCluster()->eta() <<endl;  
    if (Debug_flag) cout<<"secondele->phi() = "<<secondele->phi() <<endl;  
    if (Debug_flag) cout<<"-------------" <<endl;  
  }

  // TRIGGER MATCHING (set flags to TRUE if matched):

  Handle < edm::RefToBaseVector<reco::GsfElectron> > TagHLTelectronCollection;
  iEvent.getByLabel (theTagHLTElectronCollectionLabel, TagHLTelectronCollection);
  if (!TagHLTelectronCollection.isValid ()) return false;

  bool HLTmatch = false;

  double deltaReles = 999.0;

  for (edm::RefToBaseVector<reco::GsfElectron>::const_iterator TagHLTElectron = TagHLTelectronCollection->begin (); TagHLTElectron != TagHLTelectronCollection->end (); TagHLTElectron++) {
    if (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) {
      deltaReles = sqrt((firstele->superCluster()->eta()-(*TagHLTElectron)->superCluster()->eta())*
			(firstele->superCluster()->eta()-(*TagHLTElectron)->superCluster()->eta())+
			(firstele->phi()-(*TagHLTElectron)->phi())*
			(firstele->phi()-(*TagHLTElectron)->phi()));
    }
    if (RECO_efficiency_) {
      deltaReles = sqrt((probe_SC->eta()-(*TagHLTElectron)->superCluster()->eta())*
			(probe_SC->eta()-(*TagHLTElectron)->superCluster()->eta())+
			(probe_SC->phi()-(*TagHLTElectron)->phi())*
			(probe_SC->phi()-(*TagHLTElectron)->phi()));
    }
    if (deltaReles < 0.2) HLTmatch = true;
  }
  
  if (!HLTmatch) return false;

  Handle < edm::RefToBaseVector<reco::GsfElectron> > ProbeHLTelectronCollection;
  iEvent.getByLabel (theProbeHLTElectronCollectionLabel, ProbeHLTelectronCollection);
  if (!ProbeHLTelectronCollection.isValid ()) return false;
  
  if (HLTele17_efficiency_ || HLTele8_efficiency_) {
    HLTmatch = false;
    for (edm::RefToBaseVector<reco::GsfElectron>::const_iterator ProbeHLTElectron = ProbeHLTelectronCollection->begin (); ProbeHLTElectron != ProbeHLTelectronCollection->end (); ProbeHLTElectron++) {
      deltaReles = sqrt((secondele->superCluster()->eta()-(*ProbeHLTElectron)->superCluster()->eta())*
			(secondele->superCluster()->eta()-(*ProbeHLTElectron)->superCluster()->eta())+
			(secondele->phi()-(*ProbeHLTElectron)->phi())*
			(secondele->phi()-(*ProbeHLTElectron)->phi()));
      if (deltaReles < 0.2) HLTmatch = true;
    }
  }

  if (Debug_flag) cout<<"Trigger matching finished: HLTmatch = "<<HLTmatch<<endl;  

  //Calculating tag & probe stuff
  TLorentzVector tagv;
  tagv.SetPtEtaPhiM(firstele->pt(),firstele->superCluster()->eta(),firstele->phi(), 0.0);
  TLorentzVector probev;
  if (RECO_efficiency_)  {
    probev.SetXYZT(probe_SC->energy()*cos(probe_SC->phi())/cosh(probe_SC->eta()),
		   probe_SC->energy()*sin(probe_SC->phi())/cosh(probe_SC->eta()),
		   probe_SC->energy()*tanh(probe_SC->eta()),
		   probe_SC->energy());
  } else {
    probev.SetPtEtaPhiM(secondele->pt(),secondele->superCluster()->eta(),secondele->phi(), 0.0);
  }


  TLorentzVector e_pair = tagv + probev;
  double e_ee_invMass = e_pair.M ();

  //cut on the tag and probe mass...
  if (e_ee_invMass>120.0 || e_ee_invMass<60.0) return false;

  if (Debug_flag) cout<<"Start filling TAP histos..."<<endl;

  // Filling TAP distributions:

  // WP80 & HLT:
  if (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) {
    probeall_pt->Fill(secondele->pt());
    probeall_eta->Fill(secondele->superCluster()->eta());
    probeall_mee->Fill(e_ee_invMass);
    tagall_pt->Fill(firstele->pt());
    tagall_eta->Fill(firstele->superCluster()->eta());
    if ((WP80_efficiency_ && SelectionUtils::DoWP80Pf(secondele,iEvent,removePU_)) || ((HLTele17_efficiency_ || HLTele8_efficiency_) && HLTmatch)) {
      probepass_pt->Fill(secondele->pt());
      probepass_eta->Fill(secondele->superCluster()->eta());
      probepass_mee->Fill(e_ee_invMass);
      if (fabs(secondele->superCluster()->eta()) >= 0.0    && fabs(secondele->superCluster()->eta()) < 0.8) {
	if (secondele->pt() >= 20.0 && secondele->pt() <  30.0) probepass1eta1pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 30.0 && secondele->pt() <  40.0) probepass2eta2pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 40.0 && secondele->pt() <  50.0) probepass3eta3pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 50.0 && secondele->pt() < 999.0) probepass4eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(secondele->superCluster()->eta()) >= 0.8    && fabs(secondele->superCluster()->eta()) < 1.4442) {
	if (secondele->pt() >= 20.0 && secondele->pt() <  30.0) probepass2eta1pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 30.0 && secondele->pt() <  40.0) probepass2eta2pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 40.0 && secondele->pt() <  50.0) probepass2eta3pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 50.0 && secondele->pt() < 999.0) probepass2eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(secondele->superCluster()->eta()) >= 1.566  && fabs(secondele->superCluster()->eta()) < 2.0) {
	if (secondele->pt() >= 20.0 && secondele->pt() <  30.0) probepass3eta1pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 30.0 && secondele->pt() <  40.0) probepass3eta2pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 40.0 && secondele->pt() <  50.0) probepass3eta3pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 50.0 && secondele->pt() < 999.0) probepass3eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(secondele->superCluster()->eta()) >= 2.0    && fabs(secondele->superCluster()->eta()) < 2.5) {
	if (secondele->pt() >= 20.0 && secondele->pt() <  30.0) probepass4eta1pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 30.0 && secondele->pt() <  40.0) probepass4eta2pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 40.0 && secondele->pt() <  50.0) probepass4eta3pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 50.0 && secondele->pt() < 999.0) probepass4eta4pt->Fill(e_ee_invMass);
      }
    } else {
      probefail_pt->Fill(secondele->pt());
      probefail_eta->Fill(secondele->superCluster()->eta());
      probefail_mee->Fill(e_ee_invMass);
      if (fabs(secondele->superCluster()->eta()) >= 0.0    && fabs(secondele->superCluster()->eta()) < 0.8) {
	if (secondele->pt() >= 20.0 && secondele->pt() <  30.0) probefail1eta1pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 30.0 && secondele->pt() <  40.0) probefail2eta2pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 40.0 && secondele->pt() <  50.0) probefail3eta3pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 50.0 && secondele->pt() < 999.0) probefail4eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(secondele->superCluster()->eta()) >= 0.8    && fabs(secondele->superCluster()->eta()) < 1.4442) {
	if (secondele->pt() >= 20.0 && secondele->pt() <  30.0) probefail2eta1pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 30.0 && secondele->pt() <  40.0) probefail2eta2pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 40.0 && secondele->pt() <  50.0) probefail2eta3pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 50.0 && secondele->pt() < 999.0) probefail2eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(secondele->superCluster()->eta()) >= 1.566  && fabs(secondele->superCluster()->eta()) < 2.0) {
	if (secondele->pt() >= 20.0 && secondele->pt() <  30.0) probefail3eta1pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 30.0 && secondele->pt() <  40.0) probefail3eta2pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 40.0 && secondele->pt() <  50.0) probefail3eta3pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 50.0 && secondele->pt() < 999.0) probefail3eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(secondele->superCluster()->eta()) >= 2.0    && fabs(secondele->superCluster()->eta()) < 2.5) {
	if (secondele->pt() >= 20.0 && secondele->pt() <  30.0) probefail4eta1pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 30.0 && secondele->pt() <  40.0) probefail4eta2pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 40.0 && secondele->pt() <  50.0) probefail4eta3pt->Fill(e_ee_invMass);
	if (secondele->pt() >= 50.0 && secondele->pt() < 999.0) probefail4eta4pt->Fill(e_ee_invMass);
      }
    }
  }

  // RECO:
  if (RECO_efficiency_) {
    probeall_pt->Fill(probe_SC->energy());
    probeall_eta->Fill(probe_SC->eta());
    probeall_mee->Fill(e_ee_invMass);
    tagall_pt->Fill(firstele->pt());
    tagall_eta->Fill(firstele->superCluster()->eta());
    if (SC_isPassingProbe){
      probepass_pt->Fill(probe_SC->energy());
      probepass_eta->Fill(probe_SC->eta());
      probepass_mee->Fill(e_ee_invMass);
      if (fabs(probe_SC->eta()) >= 0.0    && fabs(probe_SC->eta()) < 0.8) {
	if (probe_SC->energy() >= 20.0 && probe_SC->energy() <  30.0) probepass1eta1pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 30.0 && probe_SC->energy() <  40.0) probepass2eta2pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 40.0 && probe_SC->energy() <  50.0) probepass3eta3pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 50.0 && probe_SC->energy() < 999.0) probepass4eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(probe_SC->eta()) >= 0.8    && fabs(probe_SC->eta()) < 1.4442) {
	if (probe_SC->energy() >= 20.0 && probe_SC->energy() <  30.0) probepass2eta1pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 30.0 && probe_SC->energy() <  40.0) probepass2eta2pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 40.0 && probe_SC->energy() <  50.0) probepass2eta3pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 50.0 && probe_SC->energy() < 999.0) probepass2eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(probe_SC->eta()) >= 1.566  && fabs(probe_SC->eta()) < 2.0) {
	if (probe_SC->energy() >= 20.0 && probe_SC->energy() <  30.0) probepass3eta1pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 30.0 && probe_SC->energy() <  40.0) probepass3eta2pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 40.0 && probe_SC->energy() <  50.0) probepass3eta3pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 50.0 && probe_SC->energy() < 999.0) probepass3eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(probe_SC->eta()) >= 2.0    && fabs(probe_SC->eta()) < 2.5) {
	if (probe_SC->energy() >= 20.0 && probe_SC->energy() <  30.0) probepass4eta1pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 30.0 && probe_SC->energy() <  40.0) probepass4eta2pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 40.0 && probe_SC->energy() <  50.0) probepass4eta3pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 50.0 && probe_SC->energy() < 999.0) probepass4eta4pt->Fill(e_ee_invMass);
      }
    } else {
      probefail_pt->Fill(probe_SC->energy());
      probefail_eta->Fill(probe_SC->eta());
      probefail_mee->Fill(e_ee_invMass);
      if (fabs(probe_SC->eta()) >= 0.0    && fabs(probe_SC->eta()) < 0.8) {
	if (probe_SC->energy() >= 20.0 && probe_SC->energy() <  30.0) probefail1eta1pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 30.0 && probe_SC->energy() <  40.0) probefail2eta2pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 40.0 && probe_SC->energy() <  50.0) probefail3eta3pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 50.0 && probe_SC->energy() < 999.0) probefail4eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(probe_SC->eta()) >= 0.8    && fabs(probe_SC->eta()) < 1.4442) {
	if (probe_SC->energy() >= 20.0 && probe_SC->energy() <  30.0) probefail2eta1pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 30.0 && probe_SC->energy() <  40.0) probefail2eta2pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 40.0 && probe_SC->energy() <  50.0) probefail2eta3pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 50.0 && probe_SC->energy() < 999.0) probefail2eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(probe_SC->eta()) >= 1.566  && fabs(probe_SC->eta()) < 2.0) {
	if (probe_SC->energy() >= 20.0 && probe_SC->energy() <  30.0) probefail3eta1pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 30.0 && probe_SC->energy() <  40.0) probefail3eta2pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 40.0 && probe_SC->energy() <  50.0) probefail3eta3pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 50.0 && probe_SC->energy() < 999.0) probefail3eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(probe_SC->eta()) >= 2.0    && fabs(probe_SC->eta()) < 2.5) {
	if (probe_SC->energy() >= 20.0 && probe_SC->energy() <  30.0) probefail4eta1pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 30.0 && probe_SC->energy() <  40.0) probefail4eta2pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 40.0 && probe_SC->energy() <  50.0) probefail4eta3pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 50.0 && probe_SC->energy() < 999.0) probefail4eta4pt->Fill(e_ee_invMass);
      }
    }
  }

  if (Debug_flag) cout<<"Finished event analysis!"<<endl;  
  return true;
}

  
// ------------ method called once each job just before starting event loop  ------------
void
EfficiencyPtEtaFilter::beginJob (){

}

// ------------ method called when starting to processes a run  ------------
bool
EfficiencyPtEtaFilter::beginRun(edm::Run &iRun, edm::EventSetup const& iSetup)
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
   if (Debug_flag) cout<<"useAllTriggers?"<<useAllTriggers_<<endl;
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

   if (Debug_flag){
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
EfficiencyPtEtaFilter::endJob ()
{

//endJob

}

DEFINE_FWK_MODULE (EfficiencyPtEtaFilter);

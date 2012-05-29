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
// $Id: EfficiencyPtEtaFilter.cc,v 1.5 2012/05/26 00:42:53 schizzi Exp $



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

bool Debug_flag=false;
bool Debug_flag2=false;

//
// member functions
//

// ------------ method called for each event  ------------
bool
EfficiencyPtEtaFilter::filter (edm::Event & iEvent, edm::EventSetup const & iSetup)
{
  if (Debug_flag) cout<<"------- NEW Event -----"<<endl;
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
    if (superCluster->energy()>20.0 && fabs(superCluster->eta())<=1.4442) {
      if (l==0) tag_SC=superCluster;
      if (l==1){
	if (tag_SC->energy()<superCluster->energy()){
          probe_SC=tag_SC;
          tag_SC=superCluster;
        }
        else{
          probe_SC=superCluster;
	}
      }
      if (l>1){
	if (tag_SC->energy()<superCluster->energy()){
          probe_SC=tag_SC;
          tag_SC=superCluster;
        }
        else{
          if (probe_SC->energy()<superCluster->energy()) probe_SC=superCluster;
	}
      }
      l++;
    }
  }

  //EE superclusters:  
  for (reco::SuperClusterCollection::const_iterator superCluster = superClusters_EE_h->begin(); superCluster != superClusters_EE_h->end(); superCluster++) {
    if (superCluster->energy()>20.0 && (fabs(superCluster->eta())>=1.5660 && fabs(superCluster->eta())<=2.4000)) {
      if (l==0) tag_SC=superCluster;
      if (l==1){
	if (tag_SC->energy()<superCluster->energy()){
          probe_SC=tag_SC;
          tag_SC=superCluster;
        }
        else{
          probe_SC=superCluster;
	}
      }
      if (l>1){
	if (tag_SC->energy()<superCluster->energy()){
          probe_SC=tag_SC;
          tag_SC=superCluster;
        }
        else{
          if (probe_SC->energy()<superCluster->energy()) probe_SC=superCluster;
	}
      }
      l++;
    }
  }

  if (Debug_flag && l<1) cout<<"No valid SuperCluster: exiting."<<endl;
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

  pat::ElectronCollection::const_iterator tag_ele;
  pat::ElectronCollection::const_iterator probe_ele;

  // iso deposits
  IsoDepositVals isoVals(isoValInputTags_.size());
  for (size_t j = 0; j < isoValInputTags_.size(); ++j) {
     iEvent.getByLabel(isoValInputTags_[j], isoVals[j]);
  }

  int i=0;
  if (electronCollection->size()<=1) return false;
  bool protection=false;
  
  for (pat::ElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {
    protection=true;
    if (recoElectron->pt()>20.0 && recoElectron->superCluster()->eta()<2.4 && (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_)) {
      if (i==0) tag_ele=recoElectron;
      if (i==1){
	if (tag_ele->pt()<recoElectron->pt()){
	  probe_ele=tag_ele;
	  tag_ele=recoElectron;
	}
	else{
	  probe_ele=recoElectron;
	}
      }
      if (i>1){
	if (tag_ele->pt()<recoElectron->pt()){
	  probe_ele=tag_ele;
	  tag_ele=recoElectron;
	}
	else{
	  if (probe_ele->pt()<recoElectron->pt()) probe_ele=recoElectron;
	}
      }
      i++;
    }
    if (recoElectron->pt()>20.0 && recoElectron->superCluster()->eta()<2.4 && RECO_efficiency_) {
      if (sqrt((probe_SC->eta()-recoElectron->superCluster()->eta())*(probe_SC->eta()-recoElectron->superCluster()->eta())+(probe_SC->phi()-recoElectron->phi())*(probe_SC->phi()-recoElectron->phi())) < 0.2) {
	SC_isPassingProbe=true;
	continue;
      }
      if (i==0) tag_ele=recoElectron;
      if (i>0){
	if (tag_ele->pt()<recoElectron->pt()){
	  tag_ele=recoElectron;
	}
      }
      i++;
    }
  }
  if (!protection) {
    cout<<"problems with PAT collection, in Efficiency.cc-->... Please check..."<<endl;    
    return false;
  }
  
  //  if (Debug_flag) cout<<"tag_ele->pt() = "<<tag_ele->pt() <<endl;  
  if (Debug_flag) cout<<"Out of the SC and RECOele loops!"<<endl;  

  if (i<1) return false;
  if ((WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) && (i<2 || tag_ele->charge() == probe_ele->charge())) return false;


  // Mixing TAG and PROBE electrons

  if (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) {
    int tagchoice = rand()%2;
    if (tagchoice == 1) {
      pat::ElectronCollection::const_iterator temp_ele;
      temp_ele = tag_ele;
      tag_ele = probe_ele;
      probe_ele = temp_ele;
    }
  }

  if (Debug_flag) cout << "Finished mixing Tag and Probe eles." << endl;

  // TRIGGER MATCHING (set flags to TRUE if matched):

  Handle < edm::RefToBaseVector<reco::GsfElectron> > TagHLTelectronCollection;
  iEvent.getByLabel (theTagHLTElectronCollectionLabel, TagHLTelectronCollection);
  if (!TagHLTelectronCollection.isValid ()) return false;

  bool HLTmatch = false;

  double deltaReles = 999.0;

  for (edm::RefToBaseVector<reco::GsfElectron>::const_iterator TagHLTElectron = TagHLTelectronCollection->begin (); TagHLTElectron != TagHLTelectronCollection->end (); TagHLTElectron++) {
    deltaReles = sqrt((tag_ele->superCluster()->eta()-(*TagHLTElectron)->superCluster()->eta())*
		      (tag_ele->superCluster()->eta()-(*TagHLTElectron)->superCluster()->eta())+
		      (tag_ele->phi()-(*TagHLTElectron)->phi())*
		      (tag_ele->phi()-(*TagHLTElectron)->phi()));
    if (deltaReles < 0.2) HLTmatch = true;
  }
  
  if (!(HLTmatch && SelectionUtils::DoWP80Pf(tag_ele,iEvent) && SelectionUtils::DoIso2011(tag_ele, iEvent, isoVals))) return false;

  Handle < edm::RefToBaseVector<reco::GsfElectron> > ProbeHLTelectronCollection;
  iEvent.getByLabel (theProbeHLTElectronCollectionLabel, ProbeHLTelectronCollection);
  if (!ProbeHLTelectronCollection.isValid ()) return false;
  
  if (HLTele17_efficiency_ || HLTele8_efficiency_) {
    HLTmatch = false;
    for (edm::RefToBaseVector<reco::GsfElectron>::const_iterator ProbeHLTElectron = ProbeHLTelectronCollection->begin (); ProbeHLTElectron != ProbeHLTelectronCollection->end (); ProbeHLTElectron++) {
      deltaReles = sqrt((probe_ele->superCluster()->eta()-(*ProbeHLTElectron)->superCluster()->eta())*
			(probe_ele->superCluster()->eta()-(*ProbeHLTElectron)->superCluster()->eta())+
			(probe_ele->phi()-(*ProbeHLTElectron)->phi())*
			(probe_ele->phi()-(*ProbeHLTElectron)->phi()));
      if (deltaReles < 0.2) HLTmatch = true;
    }
  }

  if (Debug_flag) cout<<"Trigger matching finished: HLTmatch = "<<HLTmatch<<endl;  

  //Calculating tag & probe stuff
  TLorentzVector tagv;
  tagv.SetPtEtaPhiM(tag_ele->pt(),tag_ele->superCluster()->eta(),tag_ele->phi(), 0.0);
  TLorentzVector probev;
  if (RECO_efficiency_)  {
    probev.SetXYZT(probe_SC->energy()*cos(probe_SC->phi())/cosh(probe_SC->eta()),
		   probe_SC->energy()*sin(probe_SC->phi())/cosh(probe_SC->eta()),
		   probe_SC->energy()*tanh(probe_SC->eta()),
		   probe_SC->energy());
  } else {
    probev.SetPtEtaPhiM(probe_ele->pt(),probe_ele->superCluster()->eta(),probe_ele->phi(), 0.0);
  }


  TLorentzVector e_pair = tagv + probev;
  double e_ee_invMass = e_pair.M ();

  //cut on the tag and probe mass...
  if (e_ee_invMass>120.0 || e_ee_invMass<60.0) return false;

  if (Debug_flag) cout<<"Start filling TAP histos..."<<endl;

  if (Debug_flag2) {
    if (RECO_efficiency_) {
      if (probe_SC->energy() > tag_ele->pt()) { cout << "PROBE wins!" << endl;} else { cout << "TAG wins!" << endl;}
    } else {
      if (probe_ele->pt() > tag_ele->pt()) { cout << "PROBE wins!" << endl;} else { cout << "TAG wins!" << endl;}
    }
  }
  
  if (Debug_flag) cout<<"l = "<<l<<"; i = "<<i<<"; SC_isPassingProbe  = "<<SC_isPassingProbe<<endl;  
  if (Debug_flag) cout<<"-------------" <<endl;  
  if (Debug_flag) cout<<"tag_ele->pt() = "<<tag_ele->pt() <<endl;  
  if (Debug_flag) cout<<"tag_ele->superCluster()->eta() = "<<tag_ele->superCluster()->eta() <<endl;  
  if (Debug_flag) cout<<"tag_ele->phi() = "<<tag_ele->phi() <<endl;  
  if (Debug_flag) cout<<"-------------" <<endl;  
  if (RECO_efficiency_)  {
    if (Debug_flag) cout<<"-------------" <<endl;  
    if (Debug_flag) cout<<"probe_SC->energy() = "<<probe_SC->energy() <<endl;  
    if (Debug_flag) cout<<"probe_SC->eta() = "<<probe_SC->eta() <<endl;  
    if (Debug_flag) cout<<"probe_SC->phi() = "<<probe_SC->phi() <<endl;  
    if (Debug_flag) cout<<"-------------" <<endl;  
  } else {
    if (Debug_flag) cout<<"-------------" <<endl;  
    if (Debug_flag) cout<<"probe_ele->pt() = "<<probe_ele->pt() <<endl;  
    if (Debug_flag) cout<<"probe_ele->superCluster()->eta() = "<<probe_ele->superCluster()->eta() <<endl;  
    if (Debug_flag) cout<<"probe_ele->phi() = "<<probe_ele->phi() <<endl;  
    if (Debug_flag) cout<<"-------------" <<endl;  
  }
  if (Debug_flag) cout<<"-------------" <<endl;  
  if (Debug_flag) cout<<"ee Invariant mass = "<< e_ee_invMass <<endl;    
  if (Debug_flag) cout<<"-------------" <<endl;
  
  // Filling TAP distributions:

  // WP80 & HLT:
  if (WP80_efficiency_ || HLTele17_efficiency_ || HLTele8_efficiency_) {
    probeall_pt->Fill(probe_ele->pt());
    probeall_eta->Fill(probe_ele->superCluster()->eta());
    probeall_mee->Fill(e_ee_invMass);
    tagall_pt->Fill(tag_ele->pt());
    tagall_eta->Fill(tag_ele->superCluster()->eta());
    if ((WP80_efficiency_ && SelectionUtils::DoWP90Pf(probe_ele,iEvent) && SelectionUtils::DoIso2011(probe_ele, iEvent, isoVals)) || ((HLTele17_efficiency_ || HLTele8_efficiency_) && HLTmatch)) {
      probepass_pt->Fill(probe_ele->pt());
      probepass_eta->Fill(probe_ele->superCluster()->eta());
      probepass_mee->Fill(e_ee_invMass);
      if (fabs(probe_ele->superCluster()->eta()) >= 0.0    && fabs(probe_ele->superCluster()->eta()) < 0.8) {
	if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probepass1eta1pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probepass1eta2pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probepass1eta3pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probepass1eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(probe_ele->superCluster()->eta()) >= 0.8    && fabs(probe_ele->superCluster()->eta()) < 1.4442) {
	if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probepass2eta1pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probepass2eta2pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probepass2eta3pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probepass2eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(probe_ele->superCluster()->eta()) >= 1.566  && fabs(probe_ele->superCluster()->eta()) < 2.0) {
	if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probepass3eta1pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probepass3eta2pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probepass3eta3pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probepass3eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(probe_ele->superCluster()->eta()) >= 2.0    && fabs(probe_ele->superCluster()->eta()) < 2.5) {
	if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probepass4eta1pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probepass4eta2pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probepass4eta3pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probepass4eta4pt->Fill(e_ee_invMass);
      }
    } else {
      probefail_pt->Fill(probe_ele->pt());
      probefail_eta->Fill(probe_ele->superCluster()->eta());
      probefail_mee->Fill(e_ee_invMass);
      if (fabs(probe_ele->superCluster()->eta()) >= 0.0    && fabs(probe_ele->superCluster()->eta()) < 0.8) {
	if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probefail1eta1pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probefail1eta2pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probefail1eta3pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probefail1eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(probe_ele->superCluster()->eta()) >= 0.8    && fabs(probe_ele->superCluster()->eta()) < 1.4442) {
	if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probefail2eta1pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probefail2eta2pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probefail2eta3pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probefail2eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(probe_ele->superCluster()->eta()) >= 1.566  && fabs(probe_ele->superCluster()->eta()) < 2.0) {
	if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probefail3eta1pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probefail3eta2pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probefail3eta3pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probefail3eta4pt->Fill(e_ee_invMass);
      }
      if (fabs(probe_ele->superCluster()->eta()) >= 2.0    && fabs(probe_ele->superCluster()->eta()) < 2.5) {
	if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probefail4eta1pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probefail4eta2pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probefail4eta3pt->Fill(e_ee_invMass);
	if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probefail4eta4pt->Fill(e_ee_invMass);
      }
    }
  }

  // RECO:
  if (RECO_efficiency_) {
    probeall_pt->Fill(probe_SC->energy());
    probeall_eta->Fill(probe_SC->eta());
    probeall_mee->Fill(e_ee_invMass);
    tagall_pt->Fill(tag_ele->pt());
    tagall_eta->Fill(tag_ele->superCluster()->eta());
    if (SC_isPassingProbe){
      probepass_pt->Fill(probe_SC->energy());
      probepass_eta->Fill(probe_SC->eta());
      probepass_mee->Fill(e_ee_invMass);
      if (fabs(probe_SC->eta()) >= 0.0    && fabs(probe_SC->eta()) < 0.8) {
	if (probe_SC->energy() >= 20.0 && probe_SC->energy() <  30.0) probepass1eta1pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 30.0 && probe_SC->energy() <  40.0) probepass1eta2pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 40.0 && probe_SC->energy() <  50.0) probepass1eta3pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 50.0 && probe_SC->energy() < 999.0) probepass1eta4pt->Fill(e_ee_invMass);
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
	if (probe_SC->energy() >= 30.0 && probe_SC->energy() <  40.0) probefail1eta2pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 40.0 && probe_SC->energy() <  50.0) probefail1eta3pt->Fill(e_ee_invMass);
	if (probe_SC->energy() >= 50.0 && probe_SC->energy() < 999.0) probefail1eta4pt->Fill(e_ee_invMass);
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

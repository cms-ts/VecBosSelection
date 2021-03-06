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



// system include files
#include <memory>

// user include files

#include "VecBosSelection/Selection/interface/EfficiencyGSFtoPFfilter.h"
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

bool debugging_flag=false;

//
// member functions
//

// ------------ method called for each event  ------------
bool
EfficiencyGSFtoPFfilter::filter (edm::Event & iEvent, edm::EventSetup const & iSetup)
{
  if (debugging_flag) cout<<"------- NEW Event -----"<<endl;
  using namespace edm;

  Handle < reco::GsfElectronCollection > GSFelectronCollection;
  iEvent.getByLabel (theGSFElectronCollectionLabel, GSFelectronCollection);
  if (!GSFelectronCollection.isValid()) return false;

  Handle < edm::RefToBaseVector<reco::GsfElectron> > TagHLTelectronCollection;
  iEvent.getByLabel (theTagHLTElectronCollectionLabel, TagHLTelectronCollection);
  if (!TagHLTelectronCollection.isValid()) return false;

  Handle < pat::ElectronCollection > electronPFCollection;
  iEvent.getByLabel (thePFElectronCollectionLabel, electronPFCollection);
  if (!electronPFCollection.isValid ()) return false;

  reco::GsfElectronCollection::const_iterator tag_ele;
  reco::GsfElectronCollection::const_iterator probe_ele;

  // iso deposits
  IsoDepositVals isoVals(isoValInputTags_.size());
  for (size_t j = 0; j < isoValInputTags_.size(); ++j) {
     iEvent.getByLabel(isoValInputTags_[j], isoVals[j]);
  }

  int i=0;
  int gsfEleNumber =-1;
  bool protection=false;
  bool HLTmatch=false;
  double deltaPhi=0.0;
  if (GSFelectronCollection->size()<1) {
    if (debugging_flag) cout << "Few electrons for a Z... Stop here." << endl; 
    return false;
  }

  if (debugging_flag) cout<<"Begin main loop..."<<endl;
  for (reco::GsfElectronCollection::const_iterator recoElectron = GSFelectronCollection->begin ();
       recoElectron != GSFelectronCollection->end (); recoElectron++) {
    protection=true;
    gsfEleNumber++; // Needed for the gsf ref (Iso computation)
    reco::GsfElectronRef myElectronRef(GSFelectronCollection,gsfEleNumber);
    if (((recoElectron->pt()-myElectronRef->pt())>1.0) && debugging_flag) 
      cout << "Warning: WRONG ref to GSF, please check!" << endl;
    if (recoElectron->pt()>=20.0 && fabs(recoElectron->eta())<=2.4
	&& !(fabs(recoElectron->eta())>1.4442 && fabs(recoElectron->eta())<1.566)
	&& SelectionUtils::DoWP90PfGSF(recoElectron,iEvent)
	&& SelectionUtils::DoIso2011GSF(myElectronRef,iEvent,isoVals)) {
      for (edm::RefToBaseVector<reco::GsfElectron>::const_iterator TagHLTelectron = TagHLTelectronCollection->begin (); TagHLTelectron != TagHLTelectronCollection->end (); TagHLTelectron++) {
	if ((*TagHLTelectron)->pt()>15.0) {
	  if (fabs(recoElectron->phi()-(*TagHLTelectron)->phi()) < 3.1416) {
	    deltaPhi = recoElectron->phi()-(*TagHLTelectron)->phi();
	  } else {
	    deltaPhi = fabs(recoElectron->phi()-(*TagHLTelectron)->phi()) - 6.2832;
	  }
	  if (fabs(deltaPhi)<0.1 &&
	      fabs(recoElectron->eta()-(*TagHLTelectron)->eta())<0.1 &&
	      fabs(recoElectron->pt()-(*TagHLTelectron)->pt())<(0.1 * recoElectron->pt())) {
	    HLTmatch=true;
	  }
	} 
      }
      if ((HLTmatch && !HLT_efficiency_) || (HLT_efficiency_)) {
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
    }
  }
  if (!protection) {
    return false;
  }
  
  if (debugging_flag) cout<<"Out of main loop!"<<endl;  
  if (i<2) return false;
  if (tag_ele->charge() == probe_ele->charge()) return false;

  // Mixing TAG and PROBE electrons:
  int tagchoice = rand()%2;
  if (tagchoice == 1) {
    reco::GsfElectronCollection::const_iterator temp_ele;
    temp_ele = tag_ele;
    tag_ele = probe_ele;
    probe_ele = temp_ele;
  }
  
  if (debugging_flag) cout << "Finished mixing Tag and Probe eles." << endl;

  // Tight cut on Tag ele:
  if (!(SelectionUtils::DoWP80PfGSF(tag_ele,iEvent))) {
    if (debugging_flag) cout << "Tag elle not satisfying ID requirements, exit." << endl;
    return false;
  }
  //  if(HLT_efficiency_) {
  //  HLTmatch=false;
  //    for (edm::RefToBaseVector<reco::GsfElectron>::const_iterator TagHLTSinglElectron = TagHLTelectronCollection->begin (); TagHLTSinglElectron != TagHLTelectronCollection->end (); TagHLTSinglElectron++) {
  //      if ((*TagHLTSinglElectron)->pt()>15.0) {
  //	if (fabs(tag_ele->phi()-(*TagHLTSinglElectron)->phi()) < 3.1416) {
  //	  deltaPhi = tag_ele->phi()-(*TagHLTSinglElectron)->phi();
  //	} else {
  //	  deltaPhi = fabs(tag_ele->phi()-(*TagHLTSinglElectron)->phi()) - 6.2832;
  //	}
  //	if (fabs(deltaPhi)<0.1 &&
  //	    fabs(tag_ele->eta()-(*TagHLTSinglElectron)->eta())<0.1 &&
  //	    fabs(tag_ele->pt()-(*TagHLTSinglElectron)->pt())<(0.1 * tag_ele->pt())) {
  //	  HLTmatch=true;
  //	}
  //      } 
  //    }
  //    if (!HLTmatch) {
  //      if (debugging_flag) cout << "Tag elle not satisfying HLT requirements, exit." << endl;
  //      return false;
  //    }
//  }

  /////////////////////////
  // Matching the probe  //
  /////////////////////////

  Handle < edm::RefToBaseVector<reco::GsfElectron> > ProbeHLTelectronCollection;
  iEvent.getByLabel (theProbeHLTElectronCollectionLabel, ProbeHLTelectronCollection);
  if (!ProbeHLTelectronCollection.isValid()) return false;

  bool PROBEmatch = false;
  if (!HLT_efficiency_) {
    for (pat::ElectronCollection::const_iterator PFElectron = electronPFCollection->begin (); PFElectron != electronPFCollection->end (); PFElectron++) {
      if (PFElectron->pt()>15.0) {
	//      if (PFElectron->gsfTrackRef()==probe_ele->gsfTrack()) PROBEmatch=true;
	if (fabs(probe_ele->phi()-PFElectron->phi()) < 3.1416) {
	  deltaPhi = probe_ele->phi()-PFElectron->phi();
	} else {
	  deltaPhi = fabs(probe_ele->phi()-PFElectron->phi()) - 6.2832;
	}
	if (fabs(deltaPhi)<0.1 &&
	    fabs(probe_ele->eta()-PFElectron->eta())<0.1 &&
	    fabs(probe_ele->pt()-PFElectron->pt())<(0.1*probe_ele->pt())) {
	  PROBEmatch=true;
	}
      } 
    }
  } else {
    for (edm::RefToBaseVector<reco::GsfElectron>::const_iterator ProbeHLTelectron = ProbeHLTelectronCollection->begin (); ProbeHLTelectron != ProbeHLTelectronCollection->end (); ProbeHLTelectron++) {
      if ((*ProbeHLTelectron)->pt()>15.0) {
	if (fabs(probe_ele->phi()-(*ProbeHLTelectron)->phi()) < 3.1416) {
	  deltaPhi = probe_ele->phi()-(*ProbeHLTelectron)->phi();
	} else {
	  deltaPhi = fabs(probe_ele->phi()-(*ProbeHLTelectron)->phi()) - 6.2832;
	}
	if (fabs(deltaPhi)<0.1 &&
	    fabs(probe_ele->eta()-(*ProbeHLTelectron)->eta())<0.1 &&
	    fabs(probe_ele->pt()-(*ProbeHLTelectron)->pt())<(0.1 * probe_ele->pt())) {
	  PROBEmatch=true;
	}
      } 
    }
  }

  //////////////
  // MC match //
  //////////////
  if (matchMC_) {
    Handle<reco::GenParticleCollection> genPart;
    iEvent.getByLabel(genParticleCollection_, genPart);
    if (genPart.isValid()) {
      int nMatch = 0;
      for  (reco::GenParticleCollection::const_iterator itgen=genPart->begin(); itgen!=genPart->end(); itgen++) {
        if (itgen->status() == 3 && fabs(itgen->pdgId()) == 11 && itgen->mother()->pdgId() == 23) {
          double deltaPhi1, deltaR1;
          if (fabs(itgen->phi()-tag_ele->phi()) < 3.1416) {
            deltaPhi1 = itgen->phi()-tag_ele->phi();
          } else {
            deltaPhi1 = fabs(itgen->phi()-tag_ele->phi()) - 6.2832;
          }
          deltaR1 = sqrt( pow(itgen->eta()-tag_ele->eta(), 2) + deltaPhi1*deltaPhi1 );
          double deltaPhi2, deltaR2;
	  if (fabs(itgen->phi()-probe_ele->phi()) < 3.1416) {
	    deltaPhi2 = itgen->phi()-probe_ele->phi();
	  } else {
	    deltaPhi2 = fabs(itgen->phi()-probe_ele->phi()) - 6.2832;
	  }
	  deltaR2 = sqrt( pow(itgen->eta()-probe_ele->eta(), 2) + deltaPhi2*deltaPhi2 );
          if ( deltaR1 < 0.1 && tag_ele->charge() == itgen->charge() ) {
            nMatch++;
            if (debugging_flag) std::cout << "match 1" << std::endl;
          } else if ( deltaR2 < 0.1 && probe_ele->charge() == itgen->charge() ) {
            nMatch++;
            if (debugging_flag) std::cout << "match 2" << std::endl;
          } else {
            if (debugging_flag) std::cout << "no match" << std::endl;
          }
        }
      }
      if (debugging_flag) std::cout << "# match = " << nMatch << std::endl;
      if (nMatch > 2) std::cout << "TOO MANY MATCHED ELECTRONS" << std::endl;
      if (nMatch != 2) return false;
    }
  }

  ///////////////////////////////////
  // Calculating tag & probe stuff //
  ///////////////////////////////////

  TLorentzVector tagv;
  TLorentzVector probev;

  tagv.SetPtEtaPhiM(tag_ele->pt(),tag_ele->eta(),tag_ele->phi(), 0.0);
  probev.SetPtEtaPhiM(probe_ele->pt(),probe_ele->eta(),probe_ele->phi(), 0.0);

  TLorentzVector e_pair = tagv + probev;
  double e_ee_invMass = e_pair.M ();

  //cut on the tag and probe mass...
  if (e_ee_invMass>120.0 || e_ee_invMass<60.0) return false;

  if (debugging_flag) cout<<"Start filling TAP histos..."<<endl;
  

  ////////////////////////////////
  // Filling TAP distributions: //
  ////////////////////////////////

  probeall_pt->Fill(probe_ele->pt());
  probeall_eta->Fill(probe_ele->eta());
  probeall_mee->Fill(e_ee_invMass);
  tagall_pt->Fill(tag_ele->pt());
  tagall_eta->Fill(tag_ele->eta());
  if (PROBEmatch) {
    probepass_pt->Fill(probe_ele->pt());
    probepass_eta->Fill(probe_ele->eta());
    probepass_mee->Fill(e_ee_invMass);
    if (fabs(probe_ele->eta()) >= 0.0    && fabs(probe_ele->eta()) < 0.8) {
      if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probepass1eta1pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probepass1eta2pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probepass1eta3pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probepass1eta4pt->Fill(e_ee_invMass);
    }
    if (fabs(probe_ele->eta()) >= 0.8    && fabs(probe_ele->eta()) < 1.4442) {
      if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probepass2eta1pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probepass2eta2pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probepass2eta3pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probepass2eta4pt->Fill(e_ee_invMass);
    }
    if (fabs(probe_ele->eta()) >= 1.566  && fabs(probe_ele->eta()) < 2.0) {
      if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probepass3eta1pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probepass3eta2pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probepass3eta3pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probepass3eta4pt->Fill(e_ee_invMass);
    }
    if (fabs(probe_ele->eta()) >= 2.0    && fabs(probe_ele->eta()) < 2.4) {
      if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probepass4eta1pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probepass4eta2pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probepass4eta3pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probepass4eta4pt->Fill(e_ee_invMass);
    }
  } else {
    probefail_pt->Fill(probe_ele->pt());
    probefail_eta->Fill(probe_ele->eta());
    probefail_mee->Fill(e_ee_invMass);
    if (fabs(probe_ele->eta()) >= 0.0    && fabs(probe_ele->eta()) < 0.8) {
      if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probefail1eta1pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probefail1eta2pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probefail1eta3pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probefail1eta4pt->Fill(e_ee_invMass);
    }
    if (fabs(probe_ele->eta()) >= 0.8    && fabs(probe_ele->eta()) < 1.4442) {
      if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probefail2eta1pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probefail2eta2pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probefail2eta3pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probefail2eta4pt->Fill(e_ee_invMass);
    }
    if (fabs(probe_ele->eta()) >= 1.566  && fabs(probe_ele->eta()) < 2.0) {
      if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probefail3eta1pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probefail3eta2pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probefail3eta3pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probefail3eta4pt->Fill(e_ee_invMass);
    }
    if (fabs(probe_ele->eta()) >= 2.0    && fabs(probe_ele->eta()) < 2.4) {
      if (probe_ele->pt() >= 20.0 && probe_ele->pt() <  30.0) probefail4eta1pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 30.0 && probe_ele->pt() <  40.0) probefail4eta2pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 40.0 && probe_ele->pt() <  50.0) probefail4eta3pt->Fill(e_ee_invMass);
      if (probe_ele->pt() >= 50.0 && probe_ele->pt() < 999.0) probefail4eta4pt->Fill(e_ee_invMass);
    }
  }

  if (debugging_flag) cout<<"Finished event analysis!"<<endl;  
  return true;
}

  
// ------------ method called once each job just before starting event loop  ------------
void
EfficiencyGSFtoPFfilter::beginJob (){

}

// ------------ method called when starting to processes a run  ------------
bool
EfficiencyGSFtoPFfilter::beginRun(edm::Run &iRun, edm::EventSetup const& iSetup)
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
   if (debugging_flag) cout<<"useAllTriggers?"<<useAllTriggers_<<endl;
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

   if (debugging_flag){
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
EfficiencyGSFtoPFfilter::endJob ()
{

//endJob

}

DEFINE_FWK_MODULE (EfficiencyGSFtoPFfilter);

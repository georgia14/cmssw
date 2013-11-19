#ifndef METTESTER_H
#define METTESTER_H

// author: Mike Schmitt (The University of Florida)
// date: 8/24/2006
// modification: Bobby Scurlock 
// date: 03.11.2006
// note: added RMS(METx) vs SumET capability 
// modification: Rick Cavanaugh
// date: 05.11.2006 
// note: added configuration parameters 
// modification: Mike Schmitt
// date: 02.28.2007
// note: code rewrite

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonMETCorrectionData.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
//#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "TMath.h"


//#include "FWCore/PluginManager/interface/ModuleDef.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
//
//

//

//#include "FWCore/Framework/interface/ESHandle.h"


//

//  
//
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
//
//

//#include <vector>
//#include <utility>
//#include <ostream>
//#include <fstream>
//#include <iostream>
//#include <algorithm>
//#include <cmath>
//#include <memory>
//
//#include "TMath.h"



//#include <string>
//#include <map>
//





class METTester: public edm::EDAnalyzer {
public:

  explicit METTester(const edm::ParameterSet&);

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //virtual void beginJob() ;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&) ;
  //virtual void beginJob() ;
  virtual void endJob() ;
//  virtual void endRun(const edm::Run&, const edm::EventSetup&);
	void FillpfMETRes();


 private:



  // DAQ Tools
  //DQMStore* dbe_;
  std::map<std::string, MonitorElement*> me;

  // Inputs from Configuration File

  std::string METType_;
  edm::InputTag mInputCollection_;
  std::string     mOutputFile;
  //std::string FolderName_;
  std::string sample_;
  edm::InputTag inputMETLabel_;
 
  edm::InputTag inputCaloMETLabel_;
  edm::InputTag inputTrackLabel_;
  edm::InputTag inputMuonLabel_;
  edm::InputTag inputElectronLabel_;
  edm::InputTag inputBeamSpotLabel_;
 

  bool isGoodTrack( const reco::TrackRef, float d0corr );

  int minhits_;
  double maxd0_;
  double maxchi2_;
  double maxeta_;
  double maxpt_;
  double maxPtErr_;
  std::vector<int> trkQuality_;
  std::vector<int> trkAlgos_;

 // Events variables
  MonitorElement* mNvertex;

 // Common variables
  MonitorElement* mMEx;
  MonitorElement* mMEy;
  MonitorElement* mMETSig;
  MonitorElement* mMET;
  MonitorElement* mMETFine;
  MonitorElement* mMET_Nvtx;
  MonitorElement* mMETPhi;
  MonitorElement* mSumET;
  MonitorElement* mMETResolution_GenMETTrue;
  MonitorElement* mMETPhiResolution_GenMETTrue;
  MonitorElement* mMETResolution_GenMETCalo;
  MonitorElement* mMETPhiResolution_GenMETCalo;
  //CaloMET variables

  MonitorElement* mCaloMaxEtInEmTowers;
  MonitorElement* mCaloMaxEtInHadTowers;
  MonitorElement* mCaloEtFractionHadronic;
  MonitorElement* mCaloEmEtFraction;
  MonitorElement* mCaloHadEtInHB;
  MonitorElement* mCaloHadEtInHO;
  MonitorElement* mCaloHadEtInHE;
  MonitorElement* mCaloHadEtInHF;
  MonitorElement* mCaloHadEtInEB;
  MonitorElement* mCaloHadEtInEE;
  MonitorElement* mCaloEmEtInHF;
  MonitorElement* mCaloSETInpHF;
  MonitorElement* mCaloSETInmHF;
  MonitorElement* mCaloEmEtInEE;
  MonitorElement* mCaloEmEtInEB;

  //GenMET variables
  MonitorElement* mNeutralEMEtFraction;
  MonitorElement* mNeutralHadEtFraction;
  MonitorElement* mChargedEMEtFraction;
  MonitorElement* mChargedHadEtFraction;
  MonitorElement* mMuonEtFraction; 
  MonitorElement* mInvisibleEtFraction;

  //MET variables

  //PFMET variables
  MonitorElement* mMETResolution_GenMETTrue_MET0to20;
  MonitorElement* mMETResolution_GenMETTrue_MET20to40;
  MonitorElement* mMETResolution_GenMETTrue_MET40to60;
  MonitorElement* mMETResolution_GenMETTrue_MET60to80;
  MonitorElement* mMETResolution_GenMETTrue_MET80to100;
  MonitorElement* mMETResolution_GenMETTrue_MET100to150;
  MonitorElement* mMETResolution_GenMETTrue_MET150to200;
  MonitorElement* mMETResolution_GenMETTrue_MET200to300;
  MonitorElement* mMETResolution_GenMETTrue_MET300to400;
  MonitorElement* mMETResolution_GenMETTrue_MET400to500;
  MonitorElement* mMETResolution_GenMETTrue_METResolution;

  
  //TCMET specific variables  
  MonitorElement* mMExCorrection;
  MonitorElement* mMEyCorrection;
  MonitorElement* mMuonCorrectionFlag;
  MonitorElement* mtrkPt;
  MonitorElement* mtrkEta;
  MonitorElement* mtrkNhits;
  MonitorElement* mtrkChi2;
  MonitorElement* mtrkD0;
  MonitorElement* mtrkQuality;
  MonitorElement* mtrkAlgo;
  MonitorElement* mtrkPtErr;
  MonitorElement* melePt;
  MonitorElement* meleEta;
  MonitorElement* meleHoE;

  MonitorElement* mmuPt;
  MonitorElement* mmuEta;
  MonitorElement* mmuNhits;
  MonitorElement* mmuChi2;
  MonitorElement* mmuD0;
  MonitorElement* mnMus;
  MonitorElement* mnMusPis;
  MonitorElement* mmuSAhits;
  MonitorElement* mmuTesthits;
  MonitorElement* mnEls;
  MonitorElement* mfracTrks;
  MonitorElement* mdMET;
  MonitorElement* mdMETx;
  MonitorElement* mdMETy;
  MonitorElement* mdMEy;
  MonitorElement* mdMUx;
  MonitorElement* mdMUy;

};

#endif // METTESTER_H



// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/RecoPhotonFilter
// Class:      RecoPhotonFilter
// 
/**\class RecoPhotonFilter RecoPhotonFilter.cc PhaseTwoAnalysis/RecoPhotonFilter/plugins/RecoPhotonFilter.cc

Description: adds a vector of reco Photons

Implementation:
*/
//
// Original Author:  Elvire Bouvier
//         Created:  Sun, 02 Jul 2017 20:54:15 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <vector>
#include "Math/GenVector/VectorUtil.h"
#include "TMVA/Reader.h" 
#include "TMVA/MethodBDT.h" 

//
// class declaration
//

class RecoPhotonFilter : public edm::stream::EDProducer<> {
  public:
    explicit RecoPhotonFilter(const edm::ParameterSet&);
    ~RecoPhotonFilter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    enum PhotonMatchType {UNMATCHED = 0,
      TRUE_PROMPT_Photon,
      TRUE_Photon_FROM_TAU,
      TRUE_NON_PROMPT_Photon};      

  private:
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    bool isLoosePho(const reco::Photon & recopho, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal);
    bool isMediumPho(const reco::Photon & recopho, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal);
    bool isTightPho(const reco::Photon & recopho, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal);

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------

    edm::EDGetTokenT<std::vector<reco::Photon>> phosToken_;
    edm::EDGetTokenT<reco::BeamSpot> bsToken_;
    edm::EDGetTokenT<std::vector<reco::Conversion>> convToken_;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartsToken_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;    

    //PUPPI isolation tokens
    edm::EDGetTokenT<edm::ValueMap<float> > PUPPIIsolation_charged_hadrons_;
    edm::EDGetTokenT<edm::ValueMap<float> > PUPPIIsolation_neutral_hadrons_;
    edm::EDGetTokenT<edm::ValueMap<float> > PUPPIIsolation_photons_;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
RecoPhotonFilter::RecoPhotonFilter(const edm::ParameterSet& iConfig):
  phosToken_(consumes<std::vector<reco::Photon>>(iConfig.getParameter<edm::InputTag>("photons"))),
  bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
  convToken_(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions"))),
  genPartsToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
  verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
  PUPPIIsolation_charged_hadrons_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiIsolationChargedHadrons"))),
  PUPPIIsolation_neutral_hadrons_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiIsolationNeutralHadrons"))),
  PUPPIIsolation_photons_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppiIsolationPhotons")))
{
  produces<std::vector<reco::Photon>>("LoosePhotons");
  produces<std::vector<double>>("LoosePhotonRelIso");
  produces<std::vector<reco::Photon>>("MediumPhotons");
  produces<std::vector<double>>("MediumPhotonRelIso");
  produces<std::vector<reco::Photon>>("TightPhotons");
  produces<std::vector<double>>("TightPhotonRelIso");

}


RecoPhotonFilter::~RecoPhotonFilter()
{

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
  void
RecoPhotonFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;


  Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(verticesToken_, vertices);
  // Vertices
  int prVtx = -1;
  for (size_t i = 0; i < vertices->size(); i++) {
    if (vertices->at(i).isFake()) continue;
    if (vertices->at(i).ndof() <= 4) continue;
    if (prVtx < 0) prVtx = i;
  }

  Handle<std::vector<reco::Photon>> phos;
  iEvent.getByToken(phosToken_, phos);
  Handle<reco::ConversionCollection> conversions;
  iEvent.getByToken(convToken_, conversions);
  Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(bsToken_, bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();
  Handle<std::vector<reco::GenParticle>> genParts;
  iEvent.getByToken(genPartsToken_, genParts);

  edm::Handle<edm::ValueMap<float>> PUPPIIsolation_charged_hadrons;
  edm::Handle<edm::ValueMap<float>> PUPPIIsolation_neutral_hadrons;
  edm::Handle<edm::ValueMap<float>> PUPPIIsolation_photons;
  iEvent.getByToken(PUPPIIsolation_charged_hadrons_, PUPPIIsolation_charged_hadrons);
  iEvent.getByToken(PUPPIIsolation_neutral_hadrons_, PUPPIIsolation_neutral_hadrons);
  iEvent.getByToken(PUPPIIsolation_photons_, PUPPIIsolation_photons);  

  std::unique_ptr<std::vector<reco::Photon>> filteredLoosePhotons;
  std::unique_ptr<std::vector<double>> filteredLoosePhotonRelIso;
  std::unique_ptr<std::vector<reco::Photon>> filteredMediumPhotons;
  std::unique_ptr<std::vector<double>> filteredMediumPhotonRelIso;
  std::unique_ptr<std::vector<reco::Photon>> filteredTightPhotons;
  std::unique_ptr<std::vector<double>> filteredTightPhotonRelIso;
  std::vector<reco::Photon> looseVec;
  std::vector<double> looseIsoVec;
  std::vector<reco::Photon> mediumVec;
  std::vector<double> mediumIsoVec;
  std::vector<reco::Photon> tightVec;
  std::vector<double> tightIsoVec;

  for(size_t i = 0; i < phos->size(); i++) { 
    if (phos->at(i).pt() < 10.) continue;
    if (fabs(phos->at(i).eta()) > 3.) continue;
    Ptr<const reco::Photon> el_edmPtr(phos,i);
    double elpt = phos->at(i).pt();

    double relIso = 0.;
    relIso += (*PUPPIIsolation_charged_hadrons)[el_edmPtr];
    relIso += (*PUPPIIsolation_neutral_hadrons)[el_edmPtr];
    relIso += (*PUPPIIsolation_photons)[el_edmPtr];
    relIso /= elpt;

    double elMVAVal = -1.;
    bool isLoose  = isLoosePho(phos->at(i),conversions,beamspot,elMVAVal);    
    bool isMedium = isMediumPho(phos->at(i),conversions,beamspot,elMVAVal);    
    bool isTight  = isTightPho(phos->at(i),conversions,beamspot,elMVAVal);    

    if (!isLoose) continue;
    looseVec.push_back(phos->at(i));
    looseIsoVec.push_back(relIso);

    if (!isMedium) continue;
    mediumVec.push_back(phos->at(i));
    mediumIsoVec.push_back(relIso);

    if (!isTight) continue;
    tightVec.push_back(phos->at(i));
    tightIsoVec.push_back(relIso);

  }

  filteredLoosePhotons.reset(new std::vector<reco::Photon>(looseVec));
  filteredLoosePhotonRelIso.reset(new std::vector<double>(looseIsoVec));
  filteredMediumPhotons.reset(new std::vector<reco::Photon>(mediumVec));
  filteredMediumPhotonRelIso.reset(new std::vector<double>(mediumIsoVec));
  filteredTightPhotons.reset(new std::vector<reco::Photon>(tightVec));
  filteredTightPhotonRelIso.reset(new std::vector<double>(tightIsoVec));

  iEvent.put(std::move(filteredLoosePhotons), "LoosePhotons");
  iEvent.put(std::move(filteredLoosePhotonRelIso), "LoosePhotonRelIso");
  iEvent.put(std::move(filteredMediumPhotons), "MediumPhotons");
  iEvent.put(std::move(filteredMediumPhotonRelIso), "MediumPhotonRelIso");
  iEvent.put(std::move(filteredTightPhotons), "TightPhotons");
  iEvent.put(std::move(filteredTightPhotonRelIso), "TightPhotonRelIso");

  return;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
  void
RecoPhotonFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
RecoPhotonFilter::endStream() {
}


// ------------ loose pho ID -----------
bool 
RecoPhotonFilter::isLoosePho(const reco::Photon & recopho, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal) {
  return true;
}

// ------------ medium pho ID -----------
bool 
RecoPhotonFilter::isMediumPho(const reco::Photon & recopho, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal) {
  return true;  
}

// ------------ tight pho ID -----------
bool 
RecoPhotonFilter::isTightPho(const reco::Photon & recopho, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, double MVAVal) {
  return true; 
}

// ------------ method called when starting to processes a run  ------------
/*
   void
   RecoPhotonFilter::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   RecoPhotonFilter::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   RecoPhotonFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   RecoPhotonFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecoPhotonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecoPhotonFilter);

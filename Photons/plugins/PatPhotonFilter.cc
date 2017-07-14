// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/PatPhotonFilter
// Class:      PatPhotonFilter
// 
/**\class PatPhotonFilter PatPhotonFilter.cc PhaseTwoAnalysis/PatPhotonFilter/plugins/PatPhotonFilter.cc

Description: adds a vector of pat Photons

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

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include <vector>

//
// class declaration
//

class PatPhotonFilter : public edm::stream::EDProducer<> {
    public:
        explicit PatPhotonFilter(const edm::ParameterSet&);
        ~PatPhotonFilter();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginStream(edm::StreamID) override;
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        virtual void endStream() override;

        bool isLoosePho(const pat::Photon & patPho, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
        bool isMediumPho(const pat::Photon & patPho, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
        bool isTightPho(const pat::Photon & patPho, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 

        //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

        // ----------member data ---------------------------
        edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
        edm::EDGetTokenT<std::vector<pat::Photon>> phosToken_;
        edm::EDGetTokenT<reco::BeamSpot> bsToken_;
        edm::EDGetTokenT<std::vector<reco::Conversion>> convToken_;
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
PatPhotonFilter::PatPhotonFilter(const edm::ParameterSet& iConfig):
    phosToken_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons"))),
    bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
    convToken_(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions")))  
{
    produces<std::vector<pat::Photon>>("LoosePhotons");
    produces<std::vector<double>>("LoosePhotonRelIso");
    produces<std::vector<pat::Photon>>("MediumPhotons");
    produces<std::vector<double>>("MediumPhotonRelIso");
    produces<std::vector<pat::Photon>>("TightPhotons");
    produces<std::vector<double>>("TightPhotonRelIso");

}


PatPhotonFilter::~PatPhotonFilter()
{

    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
    void
PatPhotonFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    Handle<std::vector<pat::Photon>> phos;
    iEvent.getByToken(phosToken_, phos);
    Handle<reco::ConversionCollection> conversions;
    iEvent.getByToken(convToken_, conversions);
    Handle<reco::BeamSpot> bsHandle;
    iEvent.getByToken(bsToken_, bsHandle);
    const reco::BeamSpot &beamspot = *bsHandle.product();  
    std::unique_ptr<std::vector<pat::Photon>> filteredLoosePhotons;
    std::unique_ptr<std::vector<double>> filteredLoosePhotonRelIso;
    std::unique_ptr<std::vector<pat::Photon>> filteredMediumPhotons;
    std::unique_ptr<std::vector<double>> filteredMediumPhotonRelIso;
    std::unique_ptr<std::vector<pat::Photon>> filteredTightPhotons;
    std::unique_ptr<std::vector<double>> filteredTightPhotonRelIso;
    std::vector<pat::Photon> looseVec;
    std::vector<double> looseIsoVec;
    std::vector<pat::Photon> mediumVec;
    std::vector<double> mediumIsoVec;
    std::vector<pat::Photon> tightVec;
    std::vector<double> tightIsoVec;
    for (size_t i = 0; i < phos->size(); i++) {
        if (phos->at(i).pt() < 10.) continue;
        if (fabs(phos->at(i).eta()) > 3.) continue;

        bool isLoose = isLoosePho(phos->at(i),conversions,beamspot);    
        bool isMedium = isMediumPho(phos->at(i),conversions,beamspot);    
        bool isTight = isTightPho(phos->at(i),conversions,beamspot);    

        double relIso = (phos->at(i).puppiChargedHadronIso() + phos->at(i).puppiNeutralHadronIso() + phos->at(i).puppiPhotonIso()) / phos->at(i).pt();

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

    filteredLoosePhotons.reset(new std::vector<pat::Photon>(looseVec));
    filteredLoosePhotonRelIso.reset(new std::vector<double>(looseIsoVec));
    filteredMediumPhotons.reset(new std::vector<pat::Photon>(mediumVec));
    filteredMediumPhotonRelIso.reset(new std::vector<double>(mediumIsoVec));
    filteredTightPhotons.reset(new std::vector<pat::Photon>(tightVec));
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
PatPhotonFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
PatPhotonFilter::endStream() {
}

// ------------ method check that an e passes loose ID ----------------------------------
    bool
PatPhotonFilter::isLoosePho(const pat::Photon & patPho, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
    return true;
}

// ------------ method check that an e passes medium ID ----------------------------------
    bool
PatPhotonFilter::isMediumPho(const pat::Photon & patPho, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
    return true;
}

// ------------ method check that an e passes tight ID ----------------------------------
    bool
PatPhotonFilter::isTightPho(const pat::Photon & patPho, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
    return true;
}

// ------------ method called when starting to processes a run  ------------
/*
   void
   PatPhotonFilter::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   PatPhotonFilter::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   PatPhotonFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   PatPhotonFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PatPhotonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PatPhotonFilter);

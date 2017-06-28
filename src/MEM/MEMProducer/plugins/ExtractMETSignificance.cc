/*
** class  : ExtractMETSignificance
** author : L. Cadamuro (LLR)
** date   : 4 November 2016
** brief  : takes pat::MET in input and produces a standalone significance and covariance collection from it
**          used as a replacement of the previous standalone significance producer
*/

#include "MEM/MEMProducer/plugins/ExtractMETSignificance.h"

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include <FWCore/Framework/interface/EDProducer.h>
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/METReco/interface/CommonMETData.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

ExtractMETSignificance::ExtractMETSignificance(const edm::ParameterSet& iConfig)
{
    theMETTag      = consumes<pat::METCollection>     (iConfig.getParameter<edm::InputTag>("mets2"));

    theSigTag      = consumes<double>                 (iConfig.getParameter<edm::InputTag>("srcSig"));
    theCovTag      = consumes<math::Error<2>::type>   (iConfig.getParameter<edm::InputTag>("srcCov"));
    produces<pat::METCollection>();
    produces<double>("METSignificance");
    produces<math::Error<2>::type>("METCovariance");
}

ExtractMETSignificance::~ExtractMETSignificance()
{
}

void ExtractMETSignificance::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::cout << "produce ExtractMETSignificance" << std::endl;
    edm::Handle<pat::METCollection> METHandle;
    iEvent.getByToken(theMETTag, METHandle);

    const pat::MET& patMET = (*METHandle)[0];//
    
    const reco::METCovMatrix cov = patMET.getSignificanceMatrix();
    double sig = patMET.significance();

    std::auto_ptr<double> significance (new double);
   (*significance) = sig;

    std::cout << "nouvelle methode" << std::endl;
    std::auto_ptr<math::Error<2>::type> covPtr(new math::Error<2>::type());
    (*covPtr)(0,0) = cov(0,0);
    (*covPtr)(1,0) = cov(1,0);
    (*covPtr)(1,1) = cov(1,1);
    std::cout << "(*covPtr)(0,0) = " << (*covPtr)(0,0) << std::endl;
    std::cout << "(*covPtr)(1,0) = " << (*covPtr)(1,0) << std::endl;
    std::cout << "(*covPtr)(1,1) = " << (*covPtr)(1,1) << std::endl;

    std::cout << "methode classique" << std::endl;
    Handle<math::Error<2>::type> covHandle;
    iEvent.getByToken (theCovTag, covHandle);
    cout <<"covHandle ok" << endl;

    cout << "put significance" << endl;
    iEvent.put( significance, "METSignificance" );/**/
    cout << "put covPtr" << endl;
    iEvent.put( covPtr, "METCovariance" );
    
    cout << "fin producer ExtractMETSignificance" << endl;
}

void
ExtractMETSignificance::beginStream(edm::StreamID)
{
  std::cout << "begin ExtractMETSignificance" << std::endl;
}

void
ExtractMETSignificance::endStream() {
  std::cout << "\nJe sors de ExtractMETSignificance" << std::endl;
}

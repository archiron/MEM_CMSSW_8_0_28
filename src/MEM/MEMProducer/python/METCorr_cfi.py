import FWCore.ParameterSet.Config as cms

corrMETSignificance = cms.EDProducer ("ExtractMETSignificance",
    mets2=cms.InputTag("slimmedMETs"), # mets2=cms.InputTag("slimmedMETs","","OWNPARTICLES"),

    srcSig = cms.InputTag("METSignificance", "METSignificance"),
    srcCov = cms.InputTag("METSignificance", "METCovariance"),

    configName = cms.string('ExtractMETSignificance.py'),

)

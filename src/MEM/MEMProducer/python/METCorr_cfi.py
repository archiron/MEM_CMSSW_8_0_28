import FWCore.ParameterSet.Config as cms

corrMETSignificance = cms.EDProducer ("ExtractMETSignificance",
    mets2=cms.InputTag("slimmedMETs"),

    srcSig = cms.InputTag("METSignificance", "METSignificance"),
    srcCov = cms.InputTag("METSignificance", "METCovariance"),

    configName = cms.string('ExtractMETSignificance.py'),

)

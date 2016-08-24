from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = 'ST_tW_antitop_5f_weight25-v'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.inputFiles = ['Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt','Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt','Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt','Spring16_25nsV6_MC_L1FastJet_AK8PFchs.txt','Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt','Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt','Spring16_25nsV6_MC_L1FastJet_AK8PFPuppi.txt','Spring16_25nsV6_MC_L2Relative_AK8PFPuppi.txt','Spring16_25nsV6_MC_L3Absolute_AK8PFPuppi.txt']
# Name of the CMSSW configuration file
#config.JobType.psetName    = 'bkg_ana.py'
config.JobType.psetName    = 'analysis.py'
#config.JobType.allowUndistributedCMSSW = True
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
#config.Data.inputDataset = '/WJetsToLNu_13TeV-madgraph-pythia8-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDataset = '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'
#config.Data.inputDataset = '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.totalUnits = -1
config.Data.publication = False

# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'ST_tW_antitop_5f_weight25'

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'

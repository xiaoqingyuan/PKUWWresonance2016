
hadd -f mu_PKUTree_SingleTop_xww.root mu_out_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root        mu_out_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root mu_out_ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root  mu_out_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root mu_out_ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root

#ln -s mu_out_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root mu_PKUTree_TTBARpowheg_xww.root 
cp mu_out_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root mu_PKUTree_TTBARpowheg_xww.root 

hadd -f mu_PKUTree_VV_xww.root mu_out_WWToLNuQQ_13TeV-powheg.root mu_out_WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.root mu_out_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root

hadd -f  mu_PKUTree_WJetsPt180_xww.root mu_out_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root

hadd -f mu_PKUTree_allBkg_xww.root  mu_PKUTree_SingleTop_xww.root mu_PKUTree_TTBARpowheg_xww.root mu_PKUTree_VV_xww.root mu_PKUTree_WJetsPt180_xww.root

#ln -s mu_out_BulkGravWW600.root mu_PKUTree_BulkGravWW600.root
#ln -s mu_out_BulkGravWW700.root mu_PKUTree_BulkGravWW700.root
#ln -s mu_out_BulkGravWW750.root mu_PKUTree_BulkGravWW750.root
#ln -s mu_out_BulkGravWW800.root mu_PKUTree_BulkGravWW800.root
#ln -s mu_out_BulkGravWW900.root mu_PKUTree_BulkGravWW900.root
#ln -s mu_out_BulkGravWW1000.root mu_PKUTree_BulkGravWW1000.root
cp mu_out_BulkGravWW600.root mu_PKUTree_BulkGravWW600.root
cp mu_out_BulkGravWW800.root mu_PKUTree_BulkGravWW800.root
cp mu_out_BulkGravWW1000.root mu_PKUTree_BulkGravWW1000.root

root MakePseudoData.C\(\"mu\"\) -q
hadd -f mu_PKUTree_16B.root mu_out_singleMuon16BMay-v2.root


hadd -f el_PKUTree_SingleTop_xww.root el_out_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root        el_out_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root el_out_ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root  el_out_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root el_out_ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root

#ln -s el_out_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root el_PKUTree_TTBARpowheg_xww.root 
cp el_out_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root el_PKUTree_TTBARpowheg_xww.root 

hadd -f el_PKUTree_VV_xww.root el_out_WWToLNuQQ_13TeV-powheg.root el_out_WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.root el_out_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root

hadd -f  el_PKUTree_WJetsPt180_xww.root el_out_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root

hadd -f el_PKUTree_allBkg_xww.root  el_PKUTree_SingleTop_xww.root el_PKUTree_TTBARpowheg_xww.root el_PKUTree_VV_xww.root el_PKUTree_WJetsPt180_xww.root

#ln -s el_out_BulkGravWW600.root el_PKUTree_BulkGravWW600.root
#ln -s el_out_BulkGravWW700.root el_PKUTree_BulkGravWW700.root
#ln -s el_out_BulkGravWW750.root el_PKUTree_BulkGravWW750.root
#ln -s el_out_BulkGravWW800.root el_PKUTree_BulkGravWW800.root
#ln -s el_out_BulkGravWW900.root el_PKUTree_BulkGravWW900.root
#ln -s el_out_BulkGravWW1000.root el_PKUTree_BulkGravWW1000.root
cp el_out_BulkGravWW600.root el_PKUTree_BulkGravWW600.root
cp el_out_BulkGravWW800.root el_PKUTree_BulkGravWW800.root
cp el_out_BulkGravWW1000.root el_PKUTree_BulkGravWW1000.root

root MakePseudoData.C\(\"el\"\) -q
hadd -f el_PKUTree_16B.root el_out_singleEl15BMay-v2.root

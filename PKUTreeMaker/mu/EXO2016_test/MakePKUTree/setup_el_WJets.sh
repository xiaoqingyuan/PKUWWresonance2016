
#ln -s el_out_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root el_PKUTree_TTBARpowheg_xww.root 
cp el_out_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root el_PKUTree_TTBARpowheg_xww.root 


hadd -f  el_PKUTree_WJetsPt180_xww.root  el_out_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root el_out_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root el_out_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root el_out_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root el_out_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root el_out_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root el_out_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root

hadd -f el_PKUTree_allBkg_xww.root  el_PKUTree_SingleTop_xww.root el_PKUTree_TTBARpowheg_xww.root el_PKUTree_VV_xww.root el_PKUTree_WJetsPt180_xww.root

#ln -s el_out_BulkGravWW600.root el_PKUTree_BulkGravWW600.root
#ln -s el_out_BulkGravWW700.root el_PKUTree_BulkGravWW700.root
#ln -s el_out_BulkGravWW750.root el_PKUTree_BulkGravWW750.root
#ln -s el_out_BulkGravWW800.root el_PKUTree_BulkGravWW800.root
#ln -s el_out_BulkGravWW900.root el_PKUTree_BulkGravWW900.root
#ln -s el_out_BulkGravWW1000.root el_PKUTree_BulkGravWW1000.root

root MakePseudoData.C\(\"el\"\) -q

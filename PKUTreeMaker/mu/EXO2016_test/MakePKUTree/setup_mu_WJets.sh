
#ln -s mu_out_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root mu_PKUTree_TTBARpowheg_xww.root 
cp mu_out_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root mu_PKUTree_TTBARpowheg_xww.root 


hadd -f  mu_PKUTree_WJetsPt180_xww.root  mu_out_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root mu_out_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root

hadd -f mu_PKUTree_allBkg_xww.root  mu_PKUTree_SingleTop_xww.root mu_PKUTree_TTBARpowheg_xww.root mu_PKUTree_VV_xww.root mu_PKUTree_WJetsPt180_xww.root

#ln -s mu_out_BulkGravWW600.root mu_PKUTree_BulkGravWW600.root
#ln -s mu_out_BulkGravWW700.root mu_PKUTree_BulkGravWW700.root
#ln -s mu_out_BulkGravWW750.root mu_PKUTree_BulkGravWW750.root
#ln -s mu_out_BulkGravWW800.root mu_PKUTree_BulkGravWW800.root
#ln -s mu_out_BulkGravWW900.root mu_PKUTree_BulkGravWW900.root
#ln -s mu_out_BulkGravWW1000.root mu_PKUTree_BulkGravWW1000.root

root MakePseudoData.C\(\"mu\"\) -q

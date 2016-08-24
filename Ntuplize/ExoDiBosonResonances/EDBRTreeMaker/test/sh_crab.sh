#!/bin/bash
crab submit -c crab3_analysisST_s.py
crab submit -c crab3_analysisST_t.py
crab submit -c crab3_analysisST_tW.py
crab submit -c crab3_analysisSTbar_t.py
crab submit -c crab3_analysisSTbar_tW.py
#crab submit -c crab3_analysisTT.py
#crab submit -c crab3_analysisTTMLM.py
crab submit -c crab3_analysisTTpowheg.py
crab submit -c crab3_analysisvv_w.py
crab submit -c crab3_analysisvv_wz.py
crab submit -c crab3_analysisvv_z.py
crab submit -c crab3_analysisw100.py
crab submit -c crab3_analysisw200.py
crab submit -c crab3_analysisw400.py
#crab submit -c crab3_analysisw600.py
crab submit -c crab3_analysisw600-800.py
crab submit -c crab3_analysisw800.py
crab submit -c crab3_analysisw1200.py
crab submit -c crab3_analysisw2500.py
#crab submit -c crab3_analysiswjets25.py



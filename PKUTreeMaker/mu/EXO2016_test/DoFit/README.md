Control Plots:

python g1_exo_doFit_class.py --control -c el

python g1_exo_doFit_class.py --control -c mu

Simple Try Analysis:

python g1_exo_doFit_class.py  -b -c el

python g1_exo_doFit_class.py  -b -c mu



Do Analysis:

python runLimitsEXO750_mu.py --channel mu --makeCards -b > mu.log

python runLimitsEXO750_mu.py --channel el --makeCards -b > el.log

mv cards_EXO_mu_HP750_g1/ cards_allCats

python runLimitsEXO750_mu.py --channel mu --computeLimits

python runLimitsEXO750_mu.py --channel el --computeLimits

python runLimitsEXO750_mu.py --channel mu --plotLimits  

python runLimitsEXO750_mu.py --channel el --plotLimits  


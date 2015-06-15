################################ 
#      _                   _   #
#     (_)                 | |  #
#  ___ _  __ _ _ __   __ _| |  #
# / __| |/ _` | '_ \ / _` | |  #
# \__ \ | (_| | | | | (_| | |  #
# |___/_|\__, |_| |_|\__,_|_|  #
#         __/ |                #
#        |___/                 #
#			       #
################################

python createAndSubmitMC.py -d Outputs_Spring15_AK4cors_v2 -v ggRSGgg_25nsAsympt -i Inputs_Spring15/InputList_ggRSGgg_25nsAsympt_Spring15.txt -t Inputs_Spring15/crab3_template_juska.py -c ../flat-signal-cfg_miniAOD.py -n $USER -p ggRSGgg_25nsAsympt --submit &&

python createAndSubmitMC.py -d Outputs_Spring15_AK4cors_v2 -v ggRSGgg_50nsAsympt -i Inputs_Spring15/InputList_ggRSGgg_50nsAsympt_Spring15.txt -t Inputs_Spring15/crab3_template_juska.py -c ../flat-signal-cfg_miniAOD.py -n $USER -p ggRSGgg_50nsAsympt --submit &&

python createAndSubmitMC.py -d Outputs_Spring15_AK4cors_v2 -v qgQstarqg_25nsAsympt -i Inputs_Spring15/InputList_qgQstarqg_25nsAsympt_Spring15.txt -t Inputs_Spring15/crab3_template_juska.py -c ../flat-signal-cfg_miniAOD.py -n $USER -p qgQstarqg_25nsAsympt --submit &&

python createAndSubmitMC.py -d Outputs_Spring15_AK4cors_v2 -v qgQstarqg_50nsAsympt -i Inputs_Spring15/InputList_qgQstarqg_50nsAsympt_Spring15.txt -t Inputs_Spring15/crab3_template_juska.py -c ../flat-signal-cfg_miniAOD.py -n $USER -p qgQstarqg_50nsAsympt --submit &&

python createAndSubmitMC.py -d Outputs_Spring15_AK4cors_v2 -v qqRSGqq_25nsAsympt -i Inputs_Spring15/InputList_qqRSGqq_25nsAsympt_Spring15.txt -t Inputs_Spring15/crab3_template_juska.py -c ../flat-signal-cfg_miniAOD.py -n $USER -p qqRSGqq_25nsAsympt --submit &&

python createAndSubmitMC.py -d Outputs_Spring15_AK4cors_v2 -v qqRSGqq_50nsAsympt -i Inputs_Spring15/InputList_qqRSGqq_50nsAsympt_Spring15.txt -t Inputs_Spring15/crab3_template_juska.py -c ../flat-signal-cfg_miniAOD.py -n $USER -p qqRSGqq_50nsAsympt --submit &&

###############################################################
#  _                _                                   _     #
# | |              | |                                 | |    #
# | |__   __ _  ___| | ____ _ _ __ ___  _   _ _ __   __| |    #
# | '_ \ / _` |/ __| |/ / _` | '__/ _ \| | | | '_ \ / _` |    #
# | |_) | (_| | (__|   < (_| | | | (_) | |_| | | | | (_| |    #
# |_.__/ \__,_|\___|_|\_\__, |_|  \___/ \__,_|_| |_|\__,_|    #
#                        __/ |                                #
#                       |___/                                 #
#							      #
###############################################################

python createAndSubmitMC.py -d Outputs_Spring15_AK4cors_v2 -v QCD_pthatBins_Pythia8_25nsAsympt -i Inputs_Spring15/InputList_QCD_pthatBins_Pythia8_25nsAsympt_Spring15.txt -t Inputs_Spring15/crab3_template_juska.py -c ../flat-signal-cfg_miniAOD.py -n $USER -p QCD_pthatBins_Pythia8_25nsAsympt --submit &&

python createAndSubmitMC.py -d Outputs_Spring15_AK4cors_v2 -v QCD_pthatBins_Pythia8_50nsAsympt -i Inputs_Spring15/InputList_QCD_pthatBins_Pythia8_50nsAsympt_Spring15.txt -t Inputs_Spring15/crab3_template_juska.py -c ../flat-signal-cfg_miniAOD.py -n $USER -p QCD_pthatBins_Pythia8_50nsAsympt --submit

#!/usr/bin/env python
# coding: utf-8

import os
import re
import time
import subprocess
import glob
import tarfile
import shutil
import getpass

DelExe    = '../tupleTest'
OutDir = '/store/user/%s/StopStudy' %  getpass.getuser()
tempdir = '/uscms_data/d3/%s/condor_temp/' % getpass.getuser()
ProjectName = 'signal_and_BG'
#ProjectName = 'T2tt_singal_scan'

Process = {
     # "ProcessName" : ['Path to filelist', split lines]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SM ~~~~~
     "WJetsToLNu_HT_70to100"    : ['', 10],
     #"WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt"    : ['', 10],
     "WJetsToLNu_HT_100to200"   : ['', 50],
     #"WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt"   : ['', 50],
     "WJetsToLNu_HT_200to400"   : ['', 30],
     #"WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt"   : ['', 30],
     "WJetsToLNu_HT_400to600"   : ['', 50],
     #"WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt"   : ['', 50],
     "WJetsToLNu_HT_600to800"   : ['', 20],
     #"WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt"   : ['', 20],
     "WJetsToLNu_HT_800to1200"  : ['', 7],
     #"WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt"  : ['', 7],
     "WJetsToLNu_HT_1200to2500" : ['', 4],
     #"WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt" : ['', 6],
     "WJetsToLNu_HT_2500toInf"  : ['', 4],
     #"WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt"  : ['', 4],
#
     "TTbarDiLep"                         : ['', 10],
     #"TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt"                         : ['', 20],
     "TTbarSingleLepTbar"                 : ['', 80],
     #"TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt"                 : ['', 80],
     "TTbarSingleLepT"                    : ['', 80],
     #"TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt"                    : ['', 80],

     "ZJetsToNuNu_HT_100to200"   : ['', 13],
     #"ZJetsToNuNu_HT-100To200_13TeV-madgraph.txt"   : ['', 13],
     "ZJetsToNuNu_HT_200to400"   : ['', 15],
     #"ZJetsToNuNu_HT-200To400_13TeV-madgraph.txt"   : ['', 15],
     "ZJetsToNuNu_HT_400to600"   : ['', 7],
     #"ZJetsToNuNu_HT-400To600_13TeV-madgraph.txt"   : ['', 7],
     "ZJetsToNuNu_HT_600to800"   : ['', 7],
     #"ZJetsToNuNu_HT-600To800_13TeV-madgraph.txt"   : ['', 7],
     "ZJetsToNuNu_HT_800to1200"  : ['', 10],
     #"ZJetsToNuNu_HT-800To1200_13TeV-madgraph.txt"  : ['', 10],
     "ZJetsToNuNu_HT_1200to2500" : ['', 10],
     #"ZJetsToNuNu_HT-1200To2500_13TeV-madgraph.txt" : ['', 5],
     "ZJetsToNuNu_HT_2500toInf"  : ['', 10],
     #"ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph.txt"  : ['', 6],
#
#    "QCD_HT100to200"              : ['', 40],
     "QCD_HT200to300"              : ['', 25],
     "QCD_HT300to500"              : ['', 32],
     "QCD_HT500to700"              : ['', 50],
     "QCD_HT700to1000"             : ['', 50],
     "QCD_HT1000to1500"            : ['', 13],
     "QCD_HT1500to2000"            : ['', 11],
     "QCD_HT2000toInf"             : ['', 10],

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SingleTop ~~~~~
     "ST_s"             : ['', 4],
     #"ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.txt"             : ['', 4],
     "ST_t_antitop"     : ['', 40],
     #"ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.txt"     : ['', 40],
     "ST_t_top"         : ['', 40],
     #"ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.txt"         : ['', 80],
     "tW_antitop_incl"  : ['', 60],
     "tW_top_incl"      : ['', 60],

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Top Associated Production ~~~~~
###    "ST_tWll"      : ['', 2],
###    "ST_tWnunu"    : ['', 2],
###    "TTGJets"      : ['', 16],
###    "ttHTobb"      : ['', 12],
###    "ttHToNonbb"   : ['', 10],
###    "TTTT"         : ['', 2],
###    "TTWJetsToLNu" : ['', 16],
###    "TTWJetsToQQ"  : ['', 2],
    "TTZToLLNuNu"  : ['', 10],
    #"TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt"  : ['', 10],
    "TTZToQQ"      : ['', 2],
    #"TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt"      : ['', 2],
##    "tZq_ll"       : ['', 30],
##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Diboson ~~~~~
    "ZZTo2L2Nu"   : ['', 22],
    "ZZTo2Q2Nu"   : ['', 50],
##    "ZZTo4L"      : ['', 10],
##    "ZZTo4Q"      : ['', 50],
    "WZ"          : ['', 8],
    "WWTo2L2Nu"   : ['', 2],
##    "WWTo4Q"      : ['', 2],
    "WWToLNuQQ"   : ['', 12],
##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TriBoson ~~~~~
##    "WWG" : ['', 2],
##    "WWW" : ['', 2],
##    "WWZ" : ['', 2],
##    "WZG" : ['', 2],
##    "WZZ" : ['', 2],
##    "ZZZ" : ['', 2],
#
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Signal ~~~~~
     "Signal_fastsim_T1tttt_mGluino1200_mLSP800"    : ['',1],
     "Signal_fastsim_T1tttt_mGluino1500_mLSP100"    : ['',1],
     "Signal_fastsim_T1tttt_mGluino2000_mLSP100"    : ['',1],
    # "Signal_fastsim_T2tt_425_325"       : ['',1],
     "Signal_fastsim_T2tt_mStop500_mLSP325"       : ['',1],
    # "Signal_fastsim_T2tt_650_350"       : ['',1],
     "Signal_fastsim_T2tt_mStop850_mLSP100"       : ['',1],
     "Signal_fastsim_T2tt_mStop1000_mLSP1"      : ['',1],
     "Signal_fastsim_T2fbd_mStop500_mLSP490"      : ['',1],
     "Signal_fastsim_T2fbd_mStop600_mLSP520"      : ['',1],
     "Signal_fastsim_T2cc_mStop500_mLSP420"      : ['',1],
     "Signal_fastsim_T2cc_mStop500_mLSP490"      : ['',1],
     "Signal_fastsim_T2bw_mStop500_mLSP325"      : ['',1],
     "Signal_fastsim_T2bw_mStop850_mLSP100"      : ['',1],
     "Signal_fastsim_T2bwC_mStop500_mLSP490"      : ['',1],
     "Signal_fastsim_T2bwC_mStop600_mLSP520"      : ['',1],
     "Signal_fastsim_T5ttcc_mGluino1000_mLSP800"      : ['',1],
     "Signal_fastsim_T5ttcc_mGluino1500_mLSP100"      : ['',1],

    # "Signal_T2tt_mStop250_mLSP150"      : ['',1],
    # "Signal_T2tt_mStop425_mLSP325"      : ['',1],
    # "Signal_T2tt_mStop500_mLSP325"      : ['',1],
    # "Signal_T2tt_mStop850_mLSP100"      : ['',1],
    # "Signal_T1tttt_mGluino1200_mLSP800" : ['',1],
    # "Signal_T1tttt_mGluino1500_mLSP100" : ['',1],
    # "Signal_T1tttt_mGluino2000_mLSP100" : ['',1],

###===========================T2tt singal scan======================
#"Signal_fastsim_T2tt_mStop1000_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop1000_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop1000_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop1000_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop1000_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop1000_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop1000_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop1000_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop1000_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop1000_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop1000_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop1000_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop1000_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop1000_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop1050_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop1050_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop1050_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop1050_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop1050_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop1050_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop1050_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop1050_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop1050_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop1050_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop1050_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop1050_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop1050_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop1050_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop1100_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop1100_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop1100_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop1100_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop1100_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop1100_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop1100_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop1100_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop1100_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop1100_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop1100_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop1100_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop1100_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop1100_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop1150_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop1150_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop1150_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop1150_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop1150_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop1150_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop1150_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop1150_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop1150_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop1150_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop1150_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop1150_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop1150_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop1150_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop1200_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop1200_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop1200_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop1200_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop1200_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop1200_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop1200_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop1200_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop1200_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop1200_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop1200_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop1200_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop1200_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop1200_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop150_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop150_mLSP25"    : ['',1],
#"Signal_fastsim_T2tt_mStop150_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop150_mLSP63"    : ['',1],
#"Signal_fastsim_T2tt_mStop167_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop175_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop175_mLSP25"    : ['',1],
#"Signal_fastsim_T2tt_mStop175_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop175_mLSP75"    : ['',1],
#"Signal_fastsim_T2tt_mStop175_mLSP88"    : ['',1],
#"Signal_fastsim_T2tt_mStop183_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop192_mLSP25"    : ['',1],
#"Signal_fastsim_T2tt_mStop200_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop200_mLSP113"    : ['',1],
#"Signal_fastsim_T2tt_mStop200_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop200_mLSP25"    : ['',1],
#"Signal_fastsim_T2tt_mStop200_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop200_mLSP75"    : ['',1],
#"Signal_fastsim_T2tt_mStop208_mLSP25"    : ['',1],
#"Signal_fastsim_T2tt_mStop217_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop225_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop225_mLSP125"    : ['',1],
#"Signal_fastsim_T2tt_mStop225_mLSP138"    : ['',1],
#"Signal_fastsim_T2tt_mStop225_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop225_mLSP25"    : ['',1],
#"Signal_fastsim_T2tt_mStop225_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop225_mLSP75"    : ['',1],
#"Signal_fastsim_T2tt_mStop233_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop242_mLSP75"    : ['',1],
#"Signal_fastsim_T2tt_mStop250_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop250_mLSP125"    : ['',1],
#"Signal_fastsim_T2tt_mStop250_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop250_mLSP163"    : ['',1],
#"Signal_fastsim_T2tt_mStop250_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop250_mLSP25"    : ['',1],
#"Signal_fastsim_T2tt_mStop250_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop250_mLSP75"    : ['',1],
#"Signal_fastsim_T2tt_mStop258_mLSP75"    : ['',1],
#"Signal_fastsim_T2tt_mStop267_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop275_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop275_mLSP125"    : ['',1],
#"Signal_fastsim_T2tt_mStop275_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop275_mLSP175"    : ['',1],
#"Signal_fastsim_T2tt_mStop275_mLSP188"    : ['',1],
#"Signal_fastsim_T2tt_mStop275_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop275_mLSP25"    : ['',1],
#"Signal_fastsim_T2tt_mStop275_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop275_mLSP75"    : ['',1],
#"Signal_fastsim_T2tt_mStop283_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop292_mLSP125"    : ['',1],
#"Signal_fastsim_T2tt_mStop300_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop300_mLSP125"    : ['',1],
#"Signal_fastsim_T2tt_mStop300_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop300_mLSP175"    : ['',1],
#"Signal_fastsim_T2tt_mStop300_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop300_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop300_mLSP213"    : ['',1],
#"Signal_fastsim_T2tt_mStop300_mLSP25"    : ['',1],
#"Signal_fastsim_T2tt_mStop300_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop300_mLSP75"    : ['',1],
#"Signal_fastsim_T2tt_mStop308_mLSP125"    : ['',1],
#"Signal_fastsim_T2tt_mStop317_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop325_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop325_mLSP125"    : ['',1],
#"Signal_fastsim_T2tt_mStop325_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop325_mLSP175"    : ['',1],
#"Signal_fastsim_T2tt_mStop325_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop325_mLSP225"    : ['',1],
#"Signal_fastsim_T2tt_mStop325_mLSP238"    : ['',1],
#"Signal_fastsim_T2tt_mStop325_mLSP25"    : ['',1],
#"Signal_fastsim_T2tt_mStop325_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop325_mLSP75"    : ['',1],
#"Signal_fastsim_T2tt_mStop333_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop342_mLSP175"    : ['',1],
#"Signal_fastsim_T2tt_mStop350_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop350_mLSP125"    : ['',1],
#"Signal_fastsim_T2tt_mStop350_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop350_mLSP175"    : ['',1],
#"Signal_fastsim_T2tt_mStop350_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop350_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop350_mLSP225"    : ['',1],
#"Signal_fastsim_T2tt_mStop350_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop350_mLSP263"    : ['',1],
#"Signal_fastsim_T2tt_mStop350_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop350_mLSP75"    : ['',1],
#"Signal_fastsim_T2tt_mStop358_mLSP175"    : ['',1],
#"Signal_fastsim_T2tt_mStop367_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop375_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop375_mLSP125"    : ['',1],
#"Signal_fastsim_T2tt_mStop375_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop375_mLSP175"    : ['',1],
#"Signal_fastsim_T2tt_mStop375_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop375_mLSP225"    : ['',1],
#"Signal_fastsim_T2tt_mStop375_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop375_mLSP275"    : ['',1],
#"Signal_fastsim_T2tt_mStop375_mLSP288"    : ['',1],
#"Signal_fastsim_T2tt_mStop375_mLSP75"    : ['',1],
#"Signal_fastsim_T2tt_mStop383_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop392_mLSP225"    : ['',1],
#"Signal_fastsim_T2tt_mStop400_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop400_mLSP125"    : ['',1],
#"Signal_fastsim_T2tt_mStop400_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop400_mLSP175"    : ['',1],
#"Signal_fastsim_T2tt_mStop400_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop400_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop400_mLSP225"    : ['',1],
#"Signal_fastsim_T2tt_mStop400_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop400_mLSP275"    : ['',1],
#"Signal_fastsim_T2tt_mStop400_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop400_mLSP313"    : ['',1],
#"Signal_fastsim_T2tt_mStop400_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop408_mLSP225"    : ['',1],
#"Signal_fastsim_T2tt_mStop417_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop425_mLSP125"    : ['',1],
#"Signal_fastsim_T2tt_mStop425_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop425_mLSP175"    : ['',1],
#"Signal_fastsim_T2tt_mStop425_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop425_mLSP225"    : ['',1],
#"Signal_fastsim_T2tt_mStop425_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop425_mLSP275"    : ['',1],
#"Signal_fastsim_T2tt_mStop425_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop425_mLSP325"    : ['',1],
#"Signal_fastsim_T2tt_mStop425_mLSP338"    : ['',1],
#"Signal_fastsim_T2tt_mStop433_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop442_mLSP275"    : ['',1],
#"Signal_fastsim_T2tt_mStop450_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop450_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop450_mLSP175"    : ['',1],
#"Signal_fastsim_T2tt_mStop450_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop450_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop450_mLSP225"    : ['',1],
#"Signal_fastsim_T2tt_mStop450_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop450_mLSP275"    : ['',1],
#"Signal_fastsim_T2tt_mStop450_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop450_mLSP325"    : ['',1],
#"Signal_fastsim_T2tt_mStop450_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop450_mLSP363"    : ['',1],
#"Signal_fastsim_T2tt_mStop450_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop458_mLSP275"    : ['',1],
#"Signal_fastsim_T2tt_mStop467_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop475_mLSP175"    : ['',1],
#"Signal_fastsim_T2tt_mStop475_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop475_mLSP225"    : ['',1],
#"Signal_fastsim_T2tt_mStop475_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop475_mLSP275"    : ['',1],
#"Signal_fastsim_T2tt_mStop475_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop475_mLSP325"    : ['',1],
#"Signal_fastsim_T2tt_mStop475_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop475_mLSP375"    : ['',1],
#"Signal_fastsim_T2tt_mStop475_mLSP388"    : ['',1],
#"Signal_fastsim_T2tt_mStop483_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop492_mLSP325"    : ['',1],
#"Signal_fastsim_T2tt_mStop500_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop500_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop500_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop500_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop500_mLSP225"    : ['',1],
#"Signal_fastsim_T2tt_mStop500_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop500_mLSP275"    : ['',1],
#"Signal_fastsim_T2tt_mStop500_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop500_mLSP325"    : ['',1],
#"Signal_fastsim_T2tt_mStop500_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop500_mLSP375"    : ['',1],
#"Signal_fastsim_T2tt_mStop500_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop500_mLSP413"    : ['',1],
#"Signal_fastsim_T2tt_mStop500_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop508_mLSP325"    : ['',1],
#"Signal_fastsim_T2tt_mStop517_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop525_mLSP225"    : ['',1],
#"Signal_fastsim_T2tt_mStop525_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop525_mLSP275"    : ['',1],
#"Signal_fastsim_T2tt_mStop525_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop525_mLSP325"    : ['',1],
#"Signal_fastsim_T2tt_mStop525_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop525_mLSP375"    : ['',1],
#"Signal_fastsim_T2tt_mStop525_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop525_mLSP425"    : ['',1],
#"Signal_fastsim_T2tt_mStop525_mLSP438"    : ['',1],
#"Signal_fastsim_T2tt_mStop533_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop542_mLSP375"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP275"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP325"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP375"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP425"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP463"    : ['',1],
#"Signal_fastsim_T2tt_mStop550_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop558_mLSP375"    : ['',1],
#"Signal_fastsim_T2tt_mStop567_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop575_mLSP275"    : ['',1],
#"Signal_fastsim_T2tt_mStop575_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop575_mLSP325"    : ['',1],
#"Signal_fastsim_T2tt_mStop575_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop575_mLSP375"    : ['',1],
#"Signal_fastsim_T2tt_mStop575_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop575_mLSP425"    : ['',1],
#"Signal_fastsim_T2tt_mStop575_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop575_mLSP475"    : ['',1],
#"Signal_fastsim_T2tt_mStop575_mLSP488"    : ['',1],
#"Signal_fastsim_T2tt_mStop583_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop592_mLSP425"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP325"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP375"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP425"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP475"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop600_mLSP513"    : ['',1],
#"Signal_fastsim_T2tt_mStop608_mLSP425"    : ['',1],
#"Signal_fastsim_T2tt_mStop617_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop625_mLSP325"    : ['',1],
#"Signal_fastsim_T2tt_mStop625_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop625_mLSP375"    : ['',1],
#"Signal_fastsim_T2tt_mStop625_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop625_mLSP425"    : ['',1],
#"Signal_fastsim_T2tt_mStop625_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop625_mLSP475"    : ['',1],
#"Signal_fastsim_T2tt_mStop625_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop625_mLSP525"    : ['',1],
#"Signal_fastsim_T2tt_mStop625_mLSP538"    : ['',1],
#"Signal_fastsim_T2tt_mStop633_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop642_mLSP475"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP375"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP425"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP475"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP525"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop650_mLSP563"    : ['',1],
#"Signal_fastsim_T2tt_mStop658_mLSP475"    : ['',1],
#"Signal_fastsim_T2tt_mStop667_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop675_mLSP375"    : ['',1],
#"Signal_fastsim_T2tt_mStop675_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop675_mLSP425"    : ['',1],
#"Signal_fastsim_T2tt_mStop675_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop675_mLSP475"    : ['',1],
#"Signal_fastsim_T2tt_mStop675_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop675_mLSP525"    : ['',1],
#"Signal_fastsim_T2tt_mStop675_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop675_mLSP575"    : ['',1],
#"Signal_fastsim_T2tt_mStop675_mLSP588"    : ['',1],
#"Signal_fastsim_T2tt_mStop683_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop692_mLSP525"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP425"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP475"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP525"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP575"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop700_mLSP613"    : ['',1],
#"Signal_fastsim_T2tt_mStop708_mLSP525"    : ['',1],
#"Signal_fastsim_T2tt_mStop717_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop725_mLSP425"    : ['',1],
#"Signal_fastsim_T2tt_mStop725_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop725_mLSP475"    : ['',1],
#"Signal_fastsim_T2tt_mStop725_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop725_mLSP525"    : ['',1],
#"Signal_fastsim_T2tt_mStop725_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop725_mLSP575"    : ['',1],
#"Signal_fastsim_T2tt_mStop725_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop725_mLSP625"    : ['',1],
#"Signal_fastsim_T2tt_mStop725_mLSP638"    : ['',1],
#"Signal_fastsim_T2tt_mStop733_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop742_mLSP575"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP475"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP525"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP575"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP625"    : ['',1],
#"Signal_fastsim_T2tt_mStop750_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop758_mLSP575"    : ['',1],
#"Signal_fastsim_T2tt_mStop767_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop775_mLSP475"    : ['',1],
#"Signal_fastsim_T2tt_mStop775_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop775_mLSP525"    : ['',1],
#"Signal_fastsim_T2tt_mStop775_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop775_mLSP575"    : ['',1],
#"Signal_fastsim_T2tt_mStop775_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop775_mLSP625"    : ['',1],
#"Signal_fastsim_T2tt_mStop775_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop783_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop792_mLSP625"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP525"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP575"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP625"    : ['',1],
#"Signal_fastsim_T2tt_mStop800_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop808_mLSP625"    : ['',1],
#"Signal_fastsim_T2tt_mStop817_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop825_mLSP525"    : ['',1],
#"Signal_fastsim_T2tt_mStop825_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop825_mLSP575"    : ['',1],
#"Signal_fastsim_T2tt_mStop825_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop825_mLSP625"    : ['',1],
#"Signal_fastsim_T2tt_mStop825_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop833_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP575"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP625"    : ['',1],
#"Signal_fastsim_T2tt_mStop850_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop875_mLSP575"    : ['',1],
#"Signal_fastsim_T2tt_mStop875_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop875_mLSP625"    : ['',1],
#"Signal_fastsim_T2tt_mStop875_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP625"    : ['',1],
#"Signal_fastsim_T2tt_mStop900_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop925_mLSP625"    : ['',1],
#"Signal_fastsim_T2tt_mStop925_mLSP650"    : ['',1],
#"Signal_fastsim_T2tt_mStop950_mLSP100"    : ['',1],
#"Signal_fastsim_T2tt_mStop950_mLSP150"    : ['',1],
#"Signal_fastsim_T2tt_mStop950_mLSP1"    : ['',1],
#"Signal_fastsim_T2tt_mStop950_mLSP200"    : ['',1],
#"Signal_fastsim_T2tt_mStop950_mLSP250"    : ['',1],
#"Signal_fastsim_T2tt_mStop950_mLSP300"    : ['',1],
#"Signal_fastsim_T2tt_mStop950_mLSP350"    : ['',1],
#"Signal_fastsim_T2tt_mStop950_mLSP400"    : ['',1],
#"Signal_fastsim_T2tt_mStop950_mLSP450"    : ['',1],
#"Signal_fastsim_T2tt_mStop950_mLSP500"    : ['',1],
#"Signal_fastsim_T2tt_mStop950_mLSP50"    : ['',1],
#"Signal_fastsim_T2tt_mStop950_mLSP550"    : ['',1],
#"Signal_fastsim_T2tt_mStop950_mLSP600"    : ['',1],
#"Signal_fastsim_T2tt_mStop950_mLSP650"    : ['',1],

}

Mergeblock = """#!/usr/bin/env python

# File        : merge.py
# Author      : Ben Wu
# Contact     : benwu@fnal.gov
# Date        : 2015 Jul 20
#
# Description : Code to merge output hists

import re
import glob
import os
import subprocess
import multiprocessing

def MergeFile(prod):
    print "Processing %s" % prod
    g = glob.glob("%s*.root" % prod)
    logfile = open("%s.log" % prod, 'w')
    sub = re.compile(r'^%s_\d+\.root$' % prod)
    allfile = set()
    goodfile = set()
    for f in g:
        if sub.match(f) is not None:
            allfile.add(f)
            if os.path.getsize(f) != 0:
                goodfile.add(f)
    run = "hadd -f merged/%s.root " % prod
    run += " ".join(goodfile)
    process = subprocess.Popen(run, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    logfile.write(out)
    logfile.write(err)
    logfile.close()


if __name__ == "__main__":
    cmd_exists = lambda x: any(os.access(os.path.join(path, x), os.X_OK)
                               for path in os.environ["PATH"].split(os.pathsep))
    if cmd_exists('hadd'):
        if not os.path.isdir("merged"):
            os.mkdir("merged")
    else:
        HEADER = '[95m'
        OKBLUE = '[94m'
        OKGREEN = '[92m'
        WARNING = '[93m'
        FAIL = '[91m'
        ENDC = '[0m'
        BOLD = '[1m'
        UNDERLINE = '[4m'
        print(FAIL + "Warning: no hadd available! Please setup ROOT!!" + ENDC)
        exit()

    pattern = re.compile(r'^(.*)_\d+\.root$')
    g = glob.glob("*.root")

    ## Get all the process
    process = set()
    for files in g:
        match = pattern.match(files)
        if match is not None:
            process.add(match.group(1))
        else:
            print files
            cmd = "cp %s merged/" % files
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()

    print process
    ## Run with multiprocessing Pool
    pool = multiprocessing.Pool(processes = multiprocessing.cpu_count()/3)
    pool.map(MergeFile, process)
"""

def Condor_Sub(condor_file):
    curdir = os.path.abspath(os.path.curdir)
    os.chdir(os.path.dirname(condor_file))
    print "To submit condor with " + condor_file
    os.system("condor_submit " + condor_file)
    os.chdir(curdir)


def SplitPro(key, file, fraction):
    splitedfiles = []
    filelistdir = tempdir + '/' + "FileList"
    try:
        os.makedirs(filelistdir)
    except OSError:
        pass


    filename = os.path.abspath(file)
    if fraction == 1:
        splitedfiles.append(os.path.abspath(filename))
        shutil.copy2(os.path.abspath(filename), "%s/%s" % (filelistdir, os.path.basename(filename)))
        return splitedfiles

    f = open(filename, 'r')
    lines = f.readlines()
    if len(lines) <= fraction:
        lineperfile = 1
        fraction = len(lines)
    else:
        lineperfile = len(lines) / fraction
        if len(lines) % fraction > 0:
            lineperfile += 1


    for i in range(0, fraction):
        wlines = []
        if i == fraction - 1 :
            wlines = lines[lineperfile*i :]
        else:
            wlines = lines[lineperfile*i : lineperfile*(i+1)]
        if len(wlines) > 0:
            outf = open("%s/%s.%d.list" % (filelistdir, key, i), 'w')
            outf.writelines(wlines)
            splitedfiles.append(os.path.abspath("%s/%s.%d.list" % (filelistdir, key, i)))
        outf.close()

    return splitedfiles

def my_process():
    ## temp dir for submit
    global tempdir
    global Mergeblock
    global ProjectName
    ProjectName = time.strftime('%b%d') + ProjectName
    tempdir = tempdir + os.getlogin() + "/" + ProjectName +  "/"
    try:
        os.makedirs(tempdir)
    except OSError:
        pass

    ## Create the output directory
    outdir = OutDir +  "/" + ProjectName + "/"
    try:
        os.makedirs("/eos/uscms/%s" % outdir)
    except OSError:
        pass

    ## Update RunHT.csh with DelDir and pileups
    RunHTFile = tempdir + "/" + "RunExe.csh"
    with open(RunHTFile, "wt") as outfile:
        for line in open("RunExe.csh", "r"):
            line = line.replace("DELSCR", os.environ['SCRAM_ARCH'])
            line = line.replace("DELDIR", os.environ['CMSSW_VERSION'])
            line = line.replace("DELEXE", DelExe.split('/')[-1])
            line = line.replace("OUTDIR", outdir)
            outfile.write(line)

    ## Script for merging output histograms
    MergeFile = tempdir + "/" + "merge.py"
    f = open("%s/merge.py" % tempdir, 'wt')
    f.writelines(Mergeblock)
    f.close()
    shutil.copy2("%s/merge.py" % tempdir, "/eos/uscms/%s/merge.py" % outdir)
    ### Keeping track of running script
    # shutil.copy2("../src/testMain.cc", "/eos/uscms/%s/testMain.cc" % outdir)

    ### Create Tarball
    NewNpro = {}
    Tarfiles = []
    for key, value in Process.items():
        if value[0] == "":
            value[0] = "../FileList/"+key+".list"
            #value[0] = "../FileList/"+key
        if not os.path.isfile(value[0]):
            continue
        npro = GetProcess(key, value)
        Tarfiles+=npro
        NewNpro[key] = len(npro)

    Tarfiles.append(os.path.abspath(DelExe))
    Tarfiles += GetNeededFileList(key)
    tarballname ="%s/%s.tar.gz" % (tempdir, ProjectName)
    with tarfile.open(tarballname, "w:gz", dereference=True) as tar:
        [tar.add(f, arcname=f.split('/')[-1]) for f in Tarfiles]
        tar.close()

    ### Update condor files
    for key, value in Process.items():
        if NewNpro[key] > 1:
            arg = "\nArguments = %s.$(Process).list %s_$(Process).root \nQueue %d \n" % (key, key, NewNpro[key])
        else:
            arg = "\nArguments = %s.list %s.root \n Queue\n" % (key, key)

        ## Prepare the condor file
        condorfile = tempdir + "/" + "condor_" + ProjectName +"_" + key
        with open(condorfile, "wt") as outfile:
            for line in open("condor_template", "r"):
                line = line.replace("EXECUTABLE", os.path.abspath(RunHTFile))
                line = line.replace("TARFILES", tarballname)
                line = line.replace("TEMPDIR", tempdir)
                line = line.replace("PROJECTNAME", ProjectName)
                line = line.replace("ARGUMENTS", arg)
                outfile.write(line)

        Condor_Sub(condorfile)

def GetProcess(key, value):
    if len(value) == 1:
        return SplitPro(key, value[0], 1)
    else :
        return SplitPro(key, value[0], value[1])

def GetNeededFileList(key):
    relist = []
    g = glob.glob("../FileList/*.tar.gz")
    relist += [os.path.abspath(h) for h in g]
    g = glob.glob("../*root")
    relist += [os.path.abspath(h) for h in g]
    g = glob.glob("../*csv")
    relist += [os.path.abspath(h) for h in g]
    g = glob.glob("../*cfg")
    relist += [os.path.abspath(h) for h in g]
    g = glob.glob("../*model")
    relist += [os.path.abspath(h) for h in g]
    process = subprocess.Popen( "ldd %s " % os.path.abspath(DelExe) , shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for l in process.stdout:
        if os.getenv('USER') in l:
            relist.append(l.strip().split(' ')[2])
    return relist

if __name__ == "__main__":
    my_process()

%% MO-ASMO-II :: scriptExtractDataStructure script
% 1. Extract variables from data table structure
% Usage:
%  scriptExtractDataStructure
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

% Setup data variables (94 total) from data structure
c01_smpX = data.c01_smpX{1,1};
c02_hffF = data.c02_hffF{1,1};
c02_hffC = data.c02_hffC{1,1};
c02_hffCEQ = data.c02_hffCEQ{1,1};

c03_smpX_valid = data.c03_smpX_valid{1,1};
c04_hffF_valid = data.c04_hffF_valid{1,1};
c04_hffC_valid = data.c04_hffC_valid{1,1};
c04_hffCEQ_valid = data.c04_hffCEQ_valid{1,1};
c05_smpX_invalid = data.c05_smpX_invalid{1,1};
c06_hffF_invalid = data.c06_hffF_invalid{1,1};
c06_hffC_invalid = data.c06_hffC_invalid{1,1};
c06_hffCEQ_invalid = data.c06_hffCEQ_invalid{1,1};

c07_poolX_valid = data.c07_poolX_valid{1,1};
c08_poolHffF_valid = data.c08_poolHffF_valid{1,1};
c08_poolHffC_valid = data.c08_poolHffC_valid{1,1};
c08_poolHffCEQ_valid = data.c08_poolHffCEQ_valid{1,1};
c09_poolX_invalid = data.c09_poolX_invalid{1,1};
c10_poolHffF_invalid = data.c10_poolHffF_invalid{1,1};
c10_poolHffC_invalid = data.c10_poolHffC_invalid{1,1};
c10_poolHffCEQ_invalid = data.c10_poolHffCEQ_invalid{1,1};

c11_surrogateF = data.c11_surrogateF{1,1};
c11_surrogateC = data.c11_surrogateC{1,1};
c11_surrogateCEQ = data.c11_surrogateCEQ{1,1};
c12_invalidRegion = data.c12_invalidRegion{1,1};
c13_startingPoints = data.c13_startingPoints{1,1};

c14_parX = data.c14_parX{1,1};
c15_parSurF = data.c15_parSurF{1,1};
c15_parSurC = data.c15_parSurC{1,1};
c15_parSurCEQ = data.c15_parSurCEQ{1,1};
c16_parSurOut = data.c16_parSurOut{1,1};
c17_parHffF = data.c17_parHffF{1,1};
c17_parHffC = data.c17_parHffC{1,1};
c17_parHffCEQ = data.c17_parHffCEQ{1,1};

c18_parX_valid = data.c18_parX_valid{1,1};
c19_parSurF_valid = data.c19_parSurF_valid{1,1};
c19_parSurC_valid = data.c19_parSurC_valid{1,1};
c19_parSurCEQ_valid = data.c19_parSurCEQ_valid{1,1};
c20_parHffF_valid = data.c20_parHffF_valid{1,1};
c20_parHffC_valid = data.c20_parHffC_valid{1,1};
c20_parHffCEQ_valid = data.c20_parHffCEQ_valid{1,1};
c21_parX_invalid = data.c21_parX_invalid{1,1};
c22_parSurF_invalid = data.c22_parSurF_invalid{1,1};
c22_parSurC_invalid = data.c22_parSurC_invalid{1,1};
c22_parSurCEQ_invalid = data.c22_parSurCEQ_invalid{1,1};
c23_parHffF_invalid = data.c23_parHffF_invalid{1,1};
c23_parHffC_invalid = data.c23_parHffC_invalid{1,1};
c23_parHffCEQ_invalid = data.c23_parHffCEQ_invalid{1,1};

c24_valX = data.c24_valX{1,1};
c25_valSurF = data.c25_valSurF{1,1};
c25_valSurC = data.c25_valSurC{1,1};
c25_valSurCEQ = data.c25_valSurCEQ{1,1};
c26_valHffF = data.c26_valHffF{1,1};
c26_valHffC = data.c26_valHffC{1,1};
c26_valHffCEQ = data.c26_valHffCEQ{1,1};

c27_valX_valid = data.c27_valX_valid{1,1};
c28_valSurF_valid = data.c28_valSurF_valid{1,1};
c28_valSurC_valid = data.c28_valSurC_valid{1,1};
c28_valSurCEQ_valid = data.c28_valSurCEQ_valid{1,1};
c29_valHffF_valid = data.c29_valHffF_valid{1,1};
c29_valHffC_valid = data.c29_valHffC_valid{1,1};
c29_valHffCEQ_valid = data.c29_valHffCEQ_valid{1,1};
c30_valX_invalid = data.c30_valX_invalid{1,1};
c31_valSurF_invalid = data.c31_valSurF_invalid{1,1};
c31_valSurC_invalid = data.c31_valSurC_invalid{1,1};
c31_valSurCEQ_invalid = data.c31_valSurCEQ_invalid{1,1};
c32_valHffF_invalid = data.c32_valHffF_invalid{1,1};
c32_valHffC_invalid = data.c32_valHffC_invalid{1,1};
c32_valHffCEQ_invalid = data.c32_valHffCEQ_invalid{1,1};

c33_EDvecF = data.c33_EDvecF{1,1};
c33_EDavgF = data.c33_EDavgF(1,1);
c33_prevEDarrF = data.c33_prevEDarrF{1,1};

c34_HVpred = data.c34_HVpred(1,1);
c34_HVpredHistory = data.c34_HVpredHistory{1,1};
c34_HVRpred = data.c34_HVRpred(1,1);
c34_HVRpredHistory = data.c34_HVRpredHistory{1,1};

c35_resHVpred = data.c35_resHVpred(1,1);
c35_resHVpredHistory = data.c35_resHVpredHistory{1,1};
c35_minHVpred = data.c35_minHVpred(1,1);
c35_maxHVpred = data.c35_maxHVpred(1,1);
c35_resHVRpred = data.c35_resHVRpred(1,1);
c35_resHVRpredHistory = data.c35_resHVRpredHistory{1,1};
c35_minHVRpred = data.c35_minHVRpred(1,1);
c35_maxHVRpred = data.c35_maxHVRpred(1,1);

c36_HVhff = data.c36_HVhff(1,1);
c36_HVhffHistory = data.c36_HVhffHistory{1,1};
c36_HVRhff = data.c36_HVRhff(1,1);
c36_HVRhffHistory = data.c36_HVRhffHistory{1,1};

c37_resHVhff = data.c37_resHVhff(1,1);
c37_resHVhffHistory = data.c37_resHVhffHistory{1,1};
c37_minHVhff = data.c37_minHVhff(1,1);
c37_maxHVhff = data.c37_maxHVhff(1,1);
c37_resHVRhff = data.c37_resHVRhff(1,1);
c37_resHVRhffHistory = data.c37_resHVRhffHistory{1,1};
c37_minHVRhff = data.c37_minHVRhff(1,1);
c37_maxHVRhff = data.c37_maxHVRhff(1,1);

% Setup global variables from problem structure
verbose = problem.control.verbose;
randomseed = problem.control.randomseed;
randomseedRTincr = problem.control.randomseedRTincr;
randomgenerator = problem.control.randomgenerator;

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

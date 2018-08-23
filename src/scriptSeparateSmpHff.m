%% MO-ASMO-II :: scriptSeparateSmpHff script
% 1. Separate valid and invalid data for high fidelity results of training samples
% Usage:
%  scriptSeparateSmpHff
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

c03_smpX_valid = c01_smpX(idxValid,:);                      % x         valid
c05_smpX_invalid = c01_smpX(idxInvalid,:);                  % x         invalid
c04_hffF_valid = c02_hffF(idxValid,:);                      % f         valid
c06_hffF_invalid = c02_hffF(idxInvalid,:);                  % f         invalid
if (size(c02_hffC, 1) ~= 0)
    c04_hffC_valid = c02_hffC(idxValid,:);                  % c         valid
    c06_hffC_invalid = c02_hffC(idxInvalid,:);              % c         invalid
else
    c04_hffC_valid = [];
    c06_hffC_invalid = [];
end
if (size(c02_hffCEQ, 1) ~= 0)
    c04_hffCEQ_valid = c02_hffCEQ(idxValid,:);              % ceq       valid
    c06_hffCEQ_invalid = c02_hffCEQ(idxInvalid,:);          % ceq       invalid
else
    c04_hffCEQ_valid = [];
    c06_hffCEQ_invalid = [];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

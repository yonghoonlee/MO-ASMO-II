%% MO-ASMO-II :: scriptSeparateParSurHff script
% 1. Separate valid and invalid data for surrogate results and high fidelity results of Pareto set
% Usage:
%  scriptSeparateParSurHff
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

c18_parX_valid = c14_parX(idxValid,:);                      % x         valid
c21_parX_invalid = c14_parX(idxInvalid,:);                  % x         invalid
c19_parSurF_valid = c15_parSurF(idxValid,:);                % f(sur)    valid
c22_parSurF_invalid = c15_parSurF(idxInvalid,:);            % f(sur)    invalid
if (size(c15_parSurC, 1) ~= 0)
    c19_parSurC_valid = c15_parSurC(idxValid,:);            % c(sur)    valid
    c22_parSurC_invalid = c15_parSurC(idxInvalid,:);        % c(sur)    invalid
else
    c19_parSurC_valid = [];
    c22_parSurC_invalid = [];
end
if (size(c15_parSurCEQ, 1) ~= 0)
    c19_parSurCEQ_valid = c15_parSurCEQ(idxValid,:);        % ceq(sur)  valid
    c22_parSurCEQ_invalid = c15_parSurCEQ(idxInvalid,:);    % ceq(sur)  invalid
else
    c19_parSurCEQ_valid = [];
    c22_parSurCEQ_invalid = [];
end
c20_parHffF_valid = c17_parHffF(idxValid,:);                % f(hff)    valid
c23_parHffF_invalid = c17_parHffF(idxInvalid,:);            % f(hff)    invalid
if (size(c17_parHffC, 1) ~= 0)
    c20_parHffC_valid = c17_parHffC(idxValid,:);            % c(hff)    valid
    c23_parHffC_invalid = c17_parHffC(idxInvalid,:);        % c(hff)    invalid
else
    c20_parHffC_valid = [];
    c23_parHffC_invalid = [];
end
if (size(c17_parHffCEQ, 1) ~= 0)
    c20_parHffCEQ_valid = c17_parHffCEQ(idxValid,:);        % ceq(hff)  valid
    c23_parHffCEQ_invalid = c17_parHffCEQ(idxInvalid,:);    % ceq(hff)  invalid
else
    c20_parHffCEQ_valid = [];
    c23_parHffCEQ_invalid = [];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

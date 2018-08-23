%% MO-ASMO-II :: scriptSeparateValSurHff script
% 1. Separate valid and invalid data for surrogate results and high fidelity results of val. pts
% Usage:
%  scriptSeparateValSurHff
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

c27_valX_valid = c24_valX(idxValid,:);                      % x         valid
c30_valX_invalid = c24_valX(idxInvalid,:);                  % x         invalid
c28_valSurF_valid = c25_valSurF(idxValid,:);                % f(sur)    valid
c31_valSurF_invalid = c25_valSurF(idxInvalid,:);            % f(sur)    invalid
if (size(c25_valSurC, 1) ~= 0)
    c28_valSurC_valid = c25_valSurC(idxValid,:);            % c(sur)    valid
    c31_valSurC_invalid = c25_valSurC(idxInvalid,:);        % c(sur)    invalid
else
    c28_valSurC_valid = [];
    c31_valSurC_invalid = [];
end
if (size(c25_valSurCEQ, 1) ~= 0)
    c28_valSurCEQ_valid = c25_valSurCEQ(idxValid,:);        % ceq(sur)  valid
    c31_valSurCEQ_invalid = c25_valSurCEQ(idxInvalid,:);    % ceq(sur)  invalid
else
    c28_valSurCEQ_valid = [];
    c31_valSurCEQ_invalid = [];
end
c29_valHffF_valid = c26_valHffF(idxValid,:);                % f(hff)    valid
c32_valHffF_invalid = c26_valHffF(idxInvalid,:);            % f(hff)    invalid
if (size(c26_valHffC, 1) ~= 0)
    c29_valHffC_valid = c26_valHffC(idxValid,:);            % c(hff)    valid
    c32_valHffC_invalid = c26_valHffC(idxInvalid,:);        % c(hff)    invalid
else
    c29_valHffC_valid = [];
    c32_valHffC_invalid = [];
end
if (size(c26_valHffCEQ, 1) ~= 0)
    c29_valHffCEQ_valid = c26_valHffCEQ(idxValid,:);        % ceq(hff)  valid
    c32_valHffCEQ_invalid = c26_valHffCEQ(idxInvalid,:);    % ceq(hff)  invalid
else
    c29_valHffCEQ_valid = [];
    c32_valHffCEQ_invalid = [];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

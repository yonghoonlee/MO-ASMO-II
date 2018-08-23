%% MO-ASMO-II :: scriptMergePrevAndCurrentPool script
% 1. Merge (concatenate) and store previous and current high fidelity function results
% Usage:
%  scriptMergePrevAndCurrentPool
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

if k == 1
    % valid data
    c07_poolX_valid = c03_smpX_valid;
    c08_poolHffF_valid = c04_hffF_valid;
    c08_poolHffC_valid = c04_hffC_valid;
    c08_poolHffCEQ_valid = c04_hffCEQ_valid;
    % invalid data
    c09_poolX_invalid = c05_smpX_invalid;
    c10_poolHffF_invalid = c06_hffF_invalid;
    c10_poolHffC_invalid = c06_hffC_invalid;
    c10_poolHffCEQ_invalid = c06_hffCEQ_invalid;
else
    % valid data
    c07_poolX_valid = [c07_poolX_valid; c27_valX_valid; c03_smpX_valid];
    c08_poolHffF_valid = [c08_poolHffF_valid; c29_valHffF_valid; c04_hffF_valid];
    c08_poolHffC_valid = [c08_poolHffC_valid; c29_valHffC_valid; c04_hffC_valid];
    c08_poolHffCEQ_valid = [c08_poolHffCEQ_valid; c29_valHffCEQ_valid; c04_hffCEQ_valid];
    % invalid data
    c09_poolX_invalid ...
        = [c09_poolX_invalid; c30_valX_invalid; c05_smpX_invalid];
    c10_poolHffF_invalid ...
        = [c10_poolHffF_invalid; c32_valHffF_invalid; c06_hffF_invalid];
    c10_poolHffC_invalid ...
        = [c10_poolHffC_invalid; c32_valHffC_invalid; c06_hffC_invalid];
    c10_poolHffCEQ_invalid ...
        = [c10_poolHffCEQ_invalid; c32_valHffCEQ_invalid; c06_hffCEQ_invalid];
end

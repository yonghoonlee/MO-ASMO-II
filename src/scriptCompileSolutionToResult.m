%% MO-ASMO-II :: scriptCompileSolutionToResult script
% 1. Compile solution data (xopt and fopt) to result structure
% Usage:
%  scriptCompileSolutionToResult
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

xpoolall = [c07_poolX_valid; c27_valX_valid];
fpoolall = [c08_poolHffF_valid; c29_valHffF_valid];
cpoolall = [c08_poolHffC_valid; c29_valHffC_valid];
ceqpoolall = [c08_poolHffCEQ_valid; c29_valHffCEQ_valid];
ipoolall = ones(size(fpoolall,1), 1);

% constraint enforcement
if size(cpoolall,2) > 0
    cpoolall = max(cpoolall, [], 2);
    icpool = zeros(size(ipoolall));
    icpool((cpoolall - problem.control.tolC) <= 0) = 1;
    ipoolall = ipoolall & icpool;
end
if size(ceqpoolall,2) > 0
    ceqpoolall = max(sqrt(ceqpoolall.^2), [], 2);
    iceqpool = zeros(size(ipoolall));
    iceqpool((ceqpoolall - problem.control.tolCEQ) <= 0) = 1;
    ipoolall = ipoolall & iceqpool;
end
ipoolall = logical(ipoolall);
xpool = xpoolall(ipoolall, :);
fpool = fpoolall(ipoolall, :);

% non-dominated sorting for obtaining solution
[xnds, fnds, inds] = ndSort(xpool, fpool); % non-dominated sorting
xnds = xnds(inds==1, :);
fnds = fnds(inds==1, :);
if size(fnds, 1) == 0
    result.xopt = [];
    result.fopt = [];
else
    [~, inds] = sortrows(fnds, 1);
    result.xopt = xnds(inds, :);
    result.fopt = fnds(inds, :);
end

% obtaining number of high fidelity function evaluations
result.n_hff.valid_feasible = size(xpool, 1);
result.n_hff.valid_infeasible = size(xpoolall, 1) - size(xpool, 1);
result.n_hff.invalid = size(c09_poolX_invalid, 1) + size(c30_valX_invalid, 1);
result.n_hff.total ...
    = result.n_hff.valid_feasible + result.n_hff.valid_infeasible + result.n_hff.invalid;
%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

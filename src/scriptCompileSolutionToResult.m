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

xpool = [c07_poolX_valid; c27_valX_valid];
fpool = [c08_poolHffF_valid; c29_valHffF_valid];
cpool = [c08_poolHffC_valid; c29_valHffC_valid];
ceqpool = [c08_poolHffCEQ_valid; c29_valHffCEQ_valid];
ipool = ones(size(fpool,1), 1);
% constraint enforcement
if size(cpool,2) > 0
    cpool = max(cpool, [], 2);
    icpool = zeros(size(ipool));
    icpool((cpool - problem.control.tolC) <= 0) = 1;
    ipool = ipool & icpool;
end
if size(ceqpool,2) > 0
    ceqpool = max(sqrt(ceqpool.^2), [], 2);
    iceqpool = zeros(size(ipool));
    iceqpool((ceqpool - problem.control.tolCEQ) <= 0) = 1;
    ipool = ipool & iceqpool;
end
xpool = xpool(ipool, :);
fpool = fpool(ipool, :);
% non-dominated sorting
[xpool, fpool, ipool] = ndSort(xpool, fpool);
xpool = xpool(ipool==1, :);
fpool = fpool(ipool==1, :);
[~, ipool] = sortrows(fpool, 1);
result.xopt = xpool(ipool, :);
result.fopt = fpool(ipool, :);

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

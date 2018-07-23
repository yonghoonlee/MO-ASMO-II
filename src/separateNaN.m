%% MO-ASMO-II :: separateNaN function
% 1. Identify rows without and with NaN in any variable provided as varargin array
% 2. Output indices of rows with and without NaN
% Usage:
%  [idxValid, idxInvalid] = separateNaN(varargin)
% Arguments:
%  {var1, var2, ...}
%  all var#s should have the same number of rows
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [idxValid, idxInvalid] = separateNaN(varargin)
    declareGlobalVariables;
    if verbose, disp('Separating results including invalid values...'); end

    catmatrix = varargin{1};
    number = size(catmatrix, 1);
    if nargin >= 2
        for idx = 2:nargin
            catmatrix = [catmatrix, varargin{idx}];
        end
    end

    idxInvalid = find(any((isnan(catmatrix) + isinf(catmatrix)),2));
    idxValid = find(~any((isnan(catmatrix) + isinf(catmatrix)),2));
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

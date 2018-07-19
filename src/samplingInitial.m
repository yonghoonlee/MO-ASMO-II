%% MO-ASMO-II :: samplingInitial function
% 1. 
% Usage:
%  x = samplingInitial(problem, k)
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function x = samplingInitial(problem, k)
    method = problem.sampling.initial.method;
    number = problem.sampling.initial.number;
    force_number = problem.sampling.initial.force_number;
    dimension = problem.bound.num_x;
    xlb = problem.bound.xlb;
    xub = problem.bound.xub;

    switch lower(method)
    case 'lhs'
        % random permutation to create Latin Hypercube
        xt = zeros(number, dimension);
        for idx = 1:dimension
            xt(:,idx) = randperm(number)';
        end
        % randomize within Latin Hypercube
        xt = (xt - 0.5)/number;
        rn = rand(size(xt));
        rn = (rn - 0.5)/number;
        xt = xt + rn;
    case 'random'
        xt = rand(number, dimension);
    case 'cci'
        if (dimension < 2), error(['Variables in x need to be 2 or larger for ', method]); end
        if (dimension > 10), error(['Too many variables in x for ', method]); end
        xt = (ccdesign(dimension, 'Type', 'Inscribed') + 1.0)./2;
        xt = unique(xt, 'rows');
    case 'bbd'
        if (dimension < 3), error(['Variables in x need to be 3 or larger for ', method]); end
        if (dimension > 10), error(['Too many variables in x for ', method]); end
        xt = (bbdesign(dimension) + 1.0)./2;
        xt = unique(xt, 'rows');
    otherwise
        error([method,' not supported']);
    end

    % Force number of samples
    if (force_number == true) && (number < size(xt,1))
        xt = datasample(xt, number, 'Replace', false);
    end

    % Descale
    x = samplingScale(xt, xlb, xub, 'descale');
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

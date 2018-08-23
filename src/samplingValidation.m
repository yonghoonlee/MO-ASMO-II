%% MO-ASMO-II :: samplingValidation function
% 1. Sample validation points from predicted Pareto set
% Usage:
%  x = samplingValidation(problem, k)
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [valX, valF, valC, valCEQ] = samplingValidation(problem, parX, parF, parC, parCEQ)
    declareGlobalVariables;
    if verbose, disp('Sampling for validation...'); end

    num_x = size(parX, 2);
    num_f = size(parF, 2);
    num_c = size(parC, 2);
    num_ceq = size(parCEQ, 2);
    
    xnumber = size(parX, 1);
    Vnumber = problem.sampling.validation.number;
    Vmethod = problem.sampling.validation.method;

    choose_idx = zeros(xnumber, 1);
    
    if (xnumber <= Vnumber)
        choose_idx = ones(xnumber, 1);
    else
        switch lower(Vmethod)
        case 'uniform'
            % Choose anchor points
            for idx1 = 1:num_f
                [~, idx2] = min(parF(:, idx1), [], 1);
                choose_idx(idx2, 1) = 1;
            end
            % Choose end points of discontinued set
            for idx1 = 1:num_f
                [~, idx2] = sortrows(parF, idx1);
                d = zeros(length(idx2) - 1, 1);
                for idx3 = 1:length(idx2)
                    if (idx3 == 1), continue; end
                    d(idx3 - 1) = sqrt(sum((parF(idx2(idx3), :) - parF(idx2(idx3 - 1), :)).^2));
                end
                davg = mean(d);
                for idx3 = 1:length(idx2)
                    if (idx3 == 1), continue; end
                    if (d(idx3 - 1) > 3*davg)
                        choose_idx(idx2(idx3 - 1), 1) = 1;
                        choose_idx(idx2(idx3), 1) = 1;
                    end
                end
            end
            % Count how many to smaple further
            if Vnumber >= 2*sum(choose_idx)
                Vnumber_remaining = Vnumber - sum(choose_idx);
            elseif Vnumber > sum(choose_idx)
                Vnumber_remaining = Vnumber - ceil(sum(choose_idx)/2);
            else
                Vnumber_remaining = ceil(Vnumber / 2);
            end
            % Choose from remaining points
            for idx1 = 1:Vnumber_remaining
                d = zeros(xnumber, 1);
                for idx2 = 1:xnumber
                    selectedParF = parF(choose_idx == 1, :);
                    dnew = 1e9;
                    for idx3 = 1:size(selectedParF, 1)
                        newd = sqrt(sum((parF(idx2, :) - selectedParF(idx3, :)).^2));
                        dnew = min(dnew, newd);
                    end
                    d(idx2, 1) = dnew;
                end
                [~, idx4] = max(d);
                choose_idx(idx4, 1) = 1;
            end
        case 'random'
            error([Vmethod, ' not implemented']);
        otherwise
            error([Vmethod, ' not implemented']);
        end
    end
    % Return selected values
    valX = parX(choose_idx == 1, :);
    valF = parF(choose_idx == 1, :);
    if (size(parC, 1) == size(choose_idx, 1))
        valC = parC(choose_idx == 1, :);
    else
        valC = [];
    end
    if (size(parCEQ, 1) == size(choose_idx, 1))
        valCEQ = parCEQ(choose_idx == 1, :);
    else
        valCEQ = [];
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

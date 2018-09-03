%% MO-ASMO-II :: evaluateHff function
% 1. Run high fidelity function
% 2. Obtain objective function and constraint function values
% Usage:
%  [f, c, ceq, t_elapsed] = evaluateHff(problem, x)
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [f, c, ceq, t_elapsed] = evaluateHff(problem, x)
    declareGlobalVariables;
    
    t_elapsed = 0;
    tic;
    
    nx = size(x,1);
    mf = problem.bound.num_f;

    f = [];
    c = [];
    ceq = [];
    c_cheap = [];
    ceq_cheap = [];
    if (nx == 0), return; end

    npool = parallelPoolSize();
    if ((problem.functions.hifi_parallel == true) && (npool > 1)) % Parallel
        currentpool = gcp('nocreate');
        % Expensive functions evaluation
        if (size(problem.functions.hifi_combined_exp, 1) ~= 0) % F,C,CEQ defined together
            for idx = 1:nx
                hFuture(idx) = parfeval( ...
                    currentpool, problem.functions.hifi_combined_exp, 3, ...
                    x(idx,:), problem.parameter);
            end
            for idx = 1:nx
                [completedIdx, Fvalue, Cvalue, CEQvalue] = fetchNext(hFuture);
                f(completedIdx, :) = reshape(Fvalue, 1, mf);
                c(completedIdx, :) = reshape(Cvalue, 1, numel(Cvalue));
                ceq(completedIdx, :) = reshape(CEQvalue, 1, numel(CEQvalue));
            end
        else % F,C,CEQ defined separately
            if (size(problem.functions.hifi_obj_exp, 1) ~= 0)
                for idx = 1:nx
                    hFuture(idx) = parfeval( ...
                        currentpool, problem.functions.hifi_obj_exp, 1, ...
                        x(idx,:), problem.parameter);
                end
                for idx = 1:nx
                    [completedIdx, Fvalue] = fetchNext(hFuture);
                    f(completedIdx, :) = reshape(Fvalue, 1, mf);
                end
            else
                error('No objective function found.');
            end
            clear hFuture;
            if (size(problem.functions.hifi_nonlcon_exp, 1) ~= 0)
                for idx = 1:nx
                    hFuture(idx) = parfeval( ...
                        currentpool, problem.functions.hifi_nonlcon_exp, 2, ...
                        x(idx,:), problem.parameter);
                end
                for idx = 1:nx
                    [completedIdx, Cvalue, CEQvalue] = fetchNext(hFuture);
                    c(completedIdx, :) = reshape(Cvalue, 1, numel(Cvalue));
                    ceq(completedIdx, :) = reshape(CEQvalue, 1, numel(CEQvalue));
                end
            else
                c = zeros(nx,0);
                ceq = zeros(nx,0);
            end
        end
        % Cheap constraint function evaluation
        if (size(problem.functions.hifi_nonlcon_cheap, 1) ~= 0)
            for idx = 1:nx
                hFuture(idx) = parfeval( ...
                    currentpool, problem.functions.hifi_nonlcon_cheap, 2, ...
                    x(idx,:), problem.parameter);
            end
            for idx = 1:nx
                [completedIdx, Cvalue, CEQvalue] = fetchNext(hFuture);
                c_cheap(completedIdx, :) = reshape(Cvalue, 1, numel(Cvalue));
                ceq_cheap(completedIdx, :) = reshape(CEQvalue, 1, numel(CEQvalue));
            end
            c = [c, c_cheap];
            ceq = [ceq, ceq_cheap];
        end
    else % Serial
        if (problem.functions.hifi_vectorized == true) % Vectorized functions
            % Expensive functions evaluation
            if (size(problem.functions.hifi_combined_exp, 1) ~= 0) % F,C,CEQ defined together
                [f, c, ceq] = feval(problem.functions.hifi_combined_exp, x, problem.parameter);
            else % F,C,CEQ defined separately
                if (size(problem.functions.hifi_obj_exp, 1) ~= 0)
                    f = feval(problem.functions.hifi_obj_exp, x, problem.parameter);
                else
                    error('No objective function found.');
                end
                if (size(problem.functions.hifi_nonlcon_exp, 1) ~= 0)
                    [c, ceq] = feval(problem.functions.hifi_nonlcon_exp, x, problem.parameter);
                else
                    c = zeros(nx,0);
                    ceq = zeros(nx,0);
                end
            end
            % Cheap constraint function evaluation
            if (size(problem.functions.hifi_nonlcon_cheap, 1) ~= 0)
                [c_cheap, ceq_cheap] = feval( ...
                    problem.functions.hifi_nonlcon_cheap, x, problem.parameter);
                c = [c, c_cheap];
                ceq = [ceq, ceq_cheap];
            end
        else % Evaluate individually and sequentially
            % Expensive functions evaluation
            if (size(problem.functions.hifi_combined_exp, 1) ~= 0) % F,C,CEQ defined together
                for idx = 1:nx
                    [Fvalue, Cvalue, CEQvalue] = feval( ...
                        problem.functions.hifi_combined_exp, x(idx,:), problem.parameter);
                    f(idx,:) = reshape(Fvalue, 1, mf);
                    c(idx,:) = reshape(Cvalue, 1, numel(Cvalue));
                    ceq(idx,:) = reshape(CEQvalue, 1, numel(CEQvalue));
                end
            else % F,C,CEQ defined separately
                if (size(problem.functions.hifi_obj_exp, 1) ~= 0)
                    for idx = 1:nx
                        Fvalue = feval(problem.functions.hifi_obj_exp, x(idx,:), problem.parameter);
                        f(idx,:) = reshape(Fvalue, 1, mf);
                    end
                else
                    error('No objective function found.');
                end
                if (size(problem.functions.hifi_nonlcon_exp, 1) ~= 0)
                    for idx = 1:nx
                        [Cvalue, CEQvalue] = feval( ...
                            problem.functions.hifi_nonlcon_exp, x(idx,:), problem.parameter);
                        c(idx,:) = reshape(Cvalue, 1, numel(Cvalue));
                        ceq(idx,:) = reshape(CEQvalue, 1, numel(CEQvalue));
                    end
                else
                    c = zeros(nx,0);
                    ceq = zeros(nx,0);
                end
            end
            % Cheap constraint function evaluation
            if (size(problem.functions.hifi_nonlcon_cheap, 1) ~= 0)
                for idx = 1:nx
                   [Cvalue, CEQvalue] = feval( ...
                        problem.functions.hifi_nonlcon_cheap, x(idx,:), problem.parameter);
                   c_cheap(idx,:) = reshape(Cvalue, 1, numel(Cvalue));
                   ceq_cheap(idx,:) = reshape(CEQvalue, 1, numel(CEQvalue));
                end
                c = [c, c_cheap];
                ceq = [ceq, ceq_cheap];
            end
        end
    end
    
    t_elapsed = toc;
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

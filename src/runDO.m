%% MO-ASMO-II :: runDO function
% 1. Read problem structure
% 2. Run direct optimization algorithm
% 3. Return result structure
% Usage:
%  result = runMOASMO(problem, solver)
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function result = runDO(problem, solver, varargin)
    declareGlobalVariables;

    switch lower(solver)
    case 'nsga-ii'
        if verbose, disp('Direct optimization using NSGA-II'); end
        % Get initial populations
        if (nargin >= 3)
            initpop = varargin{1};
        else
            initpop = [];
        end
        % Combined function vs separate functions
        if (size(problem.functions.hifi_combined_exp, 1) ~= 0)
            % Combined function
            [xopt, fopt, out] = runCombinedObjNSGA2(problem, initpop);
        elseif (size(problem.functions.hifi_obj_exp, 1) ~= 0)
            % Separate function
            [xopt, fopt, out] = runSeparateObjNSGA2(problem, initpop);
        else
            error('No objective function provided');
        end
    case 'epsilon-constraints'
        if verbose, disp('Direct optimization using Epsilon-Constraints'); end
        % Counter
        FuncEval = 0;
        % Get previous solutions
        if (nargin >= 6)
            xprev = varargin{1};
            fprev = varargin{2};
            cprev = varargin{3};
            ceqprev = varargin{4};
        else
            xprev = (problem.bound.xlb + problem.bound.xub)./2;
            fprev = ones(1, problem.bound.num_f);
            cprev = [];
            ceqprev = [];
        end
        
        % constraint enforcement
        iprev = ones(size(fprev,1), 1);
        if size(cprev,2) > 0
            cprev = max(cprev, [], 2);
            icprev = zeros(size(iprev));
            icprev((cprev - problem.control.tolC) <= 0) = 1;
            iprev = iprev & icprev;
        end
        if size(ceqprev,2) > 0
            ceqprev = max(sqrt(ceqprev.^2), [], 2);
            iceqprev = zeros(size(iprev));
            iceqprev((ceqprev - problem.control.tolCEQ) <= 0) = 1;
            iprev = iprev & iceqprev;
        end
        xprev = xprev(iprev, :);
        fprev = fprev(iprev, :);
        
        % function structure
        if (size(problem.functions.hifi_combined_exp, 1) ~= 0)
            fn.combined = problem.functions.hifi_combined_exp;
            fn.nonlcon_cheap = problem.functions.hifi_nonlcon_cheap;
        elseif (size(problem.functions.hifi_obj_exp, 1) ~= 0)
            fn.obj = problem.functions.hifi_obj_exp;
            fn.nonlcon_exp = problem.functions.hifi_nonlcon_exp;
            fn.nonlcon_cheap = problem.functions.hifi_nonlcon_cheap;
        else
            error('No objective function provided');
        end
        
        % Find anchor points
        mx = problem.bound.num_x;
        mf = problem.bound.num_f;
        xopt_anchor = zeros(mf, mx);
        fopt_anchor = zeros(mf, mf);
        npool = parallelPoolSize();
        if ((npool > 1) && problem.functions.hifi_parallel)
            if (size(problem.functions.hifi_combined_exp, 1) ~= 0)
                funcCount = [];
                parfor idx = 1:mf
                    [~,idxsrt] = sortrows(fprev, idx);
                    xprevtmp = xprev(idxsrt, :);
                    x0 = xprevtmp(1, :);
                    epscon = [];
                    [xoptAP, ~, o3] = runECs(problem, idx, x0, epscon, 'Combined'); % run ECs
                    [foptAP, ~, ~] = computeallCombinedObj(xoptAP, problem.parameter, fn);
                    xopt_anchor(idx, :) = xoptAP;
                    fopt_anchor(idx, :) = foptAP;
                    funcCount(idx) = o3.output.funcCount + 1;
                end
                FuncEval = FuncEval + sum(funcCount);
            elseif (size(problem.functions.hifi_obj_exp, 1) ~= 0)
                funcCount = [];
                parfor idx = 1:mf
                    [~,idxsrt] = sortrows(fprev, idx);
                    xprevtmp = xprev(idxsrt, :);
                    x0 = xprevtmp(1, :);
                    epscon = [];
                    [xoptAP, ~, o3] = runECs(problem, idx, x0, epscon, 'Separated'); % run ECs
                    [foptAP, ~, ~] = computeallSeparatedObj(xoptAP, problem.parameter, fn);
                    xopt_anchor(idx, :) = xoptAP;
                    fopt_anchor(idx, :) = foptAP;
                    funcCount(idx) = o3.output.funcCount + 1;
                end
                FuncEval = FuncEval + sum(funcCount);
            end
        else
            if (size(problem.functions.hifi_combined_exp, 1) ~= 0)
                for idx = 1:mf
                    [~,idxsrt] = sortrows(fprev, idx);
                    xprevtmp = xprev(idxsrt, :);
                    x0 = xprevtmp(1, :);
                    epscon = [];
                    [xoptAP, ~, o3] = runECs(problem, idx, x0, epscon, 'Combined'); % run ECs
                    [foptAP, ~, ~] = computeallCombinedObj(xoptAP, problem.parameter, fn);
                    xopt_anchor(idx, :) = xoptAP;
                    fopt_anchor(idx, :) = foptAP;
                    FuncEval = FuncEval + o3.output.funcCount + 1;
                end
            elseif (size(problem.functions.hifi_obj_exp, 1) ~= 0)
                for idx = 1:mf
                    [~,idxsrt] = sortrows(fprev, idx);
                    xprevtmp = xprev(idxsrt, :);
                    x0 = xprevtmp(1, :);
                    epscon = [];
                    [xoptAP, ~, o3] = runECs(problem, idx, x0, epscon, 'Separated'); % run ECs
                    [foptAP, ~, ~] = computeallSeparatedObj(xoptAP, problem.parameter, fn);
                    xopt_anchor(idx, :) = xoptAP;
                    fopt_anchor(idx, :) = foptAP;
                    FuncEval = FuncEval + o3.output.funcCount + 1;
                end
            end
        end
        
        % Find non-dominated solutions
        nECs = problem.optimization.EC.num_per_dim;
        iECs = problem.optimization.EC.obj_num;
        fopt_min = min(fopt_anchor, [], 1);
        fopt_max = max(fopt_anchor, [], 1);
        % constraint combinations
        idxcon = ones(1, size(fopt_min, 2));
        idxcon(iECs) = 0;
        idxcon = logical(idxcon);
        epscon_arr = [];
        for idx = 1:size(fopt_min, 2)
            if idx ~= iECs
                tmp = transpose(linspace(fopt_min(idx), fopt_max(idx), nECs));
                epscon_arr{1, idx} = tmp(2:(end-1),:);
            else
                epscon_arr{1, idx} = fopt_max(1, idx);
            end
        end
        combinations = cell(1,numel(epscon_arr));
        [combinations{end:-1:1}] = ndgrid(epscon_arr{end:-1:1});
        combinations = cat(numel(epscon_arr)+1, combinations{:});
        combinations = reshape(combinations,[],numel(epscon_arr));
        ncomb = size(combinations, 1);
        % run optimization
        xopt_solution = zeros(ncomb, size(xopt_anchor, 2));
        fopt_solution = zeros(ncomb, size(fopt_anchor, 2));
        if ((npool > 1) && problem.functions.hifi_parallel)
            currentpool = gcp('nocreate');
            if (size(problem.functions.hifi_combined_exp, 1) ~= 0)
                hFuture = [];
                funcCount = [];
                for idx = 1:ncomb
                    epscon = combinations(idx, :);
                    epscon_compact = repmat(epscon(1, idxcon), size(fprev, 1), 1);
                    fprev_compact = fprev(:, idxcon);
                    [~, idxclose] = min(sum((epscon_compact - fprev_compact).^2, 2), [], 1);
                    x0 = xprev(idxclose, :);
                    hFuture(idx) = parfeval( ...
                        currentpool, runECs, 3, problem, iECs, x0, epscon, 'Combined');
                end
                for idx = 1:ncomb
                    [completedIdx, o1, ~, o3] = fetchNext(hFuture);
                    xopt_solution(completedIdx, :) = reshape(o1, 1, numel(o1));
                    funcCount(completedIdx) = o3.output.funcCount + 1;
                end
                hFuture = [];
                for idx = 1:ncomb
                    hFuture(idx) = parfeval( ...
                        currentpool, computeallCombinedObj, 3, ...
                        xopt_solution(idx, :), problem.parameter, fn);
                end
                for idx = 1:ncomb
                    [completedIdx, o1, ~, ~] = fetchNext(hFuture);
                    fopt_solution(completedIdx, :) = reshape(o1, 1, numel(o1));
                end
                FuncEval = FuncEval + sum(funcCount);
            elseif (size(problem.functions.hifi_obj_exp, 1) ~= 0)
                hFuture = [];
                funcCount = [];
                for idx = 1:ncomb
                    epscon = combinations(idx, :);
                    epscon_compact = repmat(epscon(1, idxcon), size(fprev, 1), 1);
                    fprev_compact = fprev(:, idxcon);
                    [~, idxclose] = min(sum((epscon_compact - fprev_compact).^2, 2), [], 1);
                    x0 = xprev(idxclose, :);
                    hFuture(idx) = parfeval( ...
                        currentpool, runECs, 3, problem, iECs, x0, epscon, 'Separated');
                end
                for idx = 1:ncomb
                    [completedIdx, o1, ~, o3] = fetchNext(hFuture);
                    xopt_solution(completedIdx, :) = reshape(o1, 1, numel(o1));
                    funcCount(completedIdx) = o3.output.funcCount + 1;
                end
                hFuture = [];
                for idx = 1:ncomb
                    hFuture(idx) = parfeval( ...
                        currentpool, computeallSeparatedObj, 3, ...
                        xopt_solution(idx, :), problem.parameter, fn);
                end
                for idx = 1:ncomb
                    [completedIdx, o1, ~, ~] = fetchNext(hFuture);
                    fopt_solution(completedIdx, :) = reshape(o1, 1, numel(o1));
                end
                FuncEval = FuncEval + sum(funcCount);
            end
        else
            if (size(problem.functions.hifi_combined_exp, 1) ~= 0)
                for idx = 1:ncomb
                    epscon = combinations(idx, :);
                    epscon_compact = repmat(epscon(1, idxcon), size(fprev, 1), 1);
                    fprev_compact = fprev(:, idxcon);
                    [~, idxclose] = min(sum((epscon_compact - fprev_compact).^2, 2), [], 1);
                    x0 = xprev(idxclose, :);
                    [o1, ~, o3] = runECs(problem, iECs, x0, epscon, 'Combined');
                    xopt_solution(idx, :) = reshape(o1, 1, numel(o1));
                    FuncEval = FuncEval + o3.output.funcCount + 1;
                end
                for idx = 1:ncomb
                    [o1, ~, ~] = computeallCombinedObj( ...
                        xopt_solution(idx, :), problem.parameter, fn);
                    fopt_solution(idx, :) = reshape(o1, 1, numel(o1));
                end
            elseif (size(problem.functions.hifi_obj_exp, 1) ~= 0)
                for idx = 1:ncomb
                    epscon = combinations(idx, :);
                    epscon_compact = repmat(epscon(1, idxcon), size(fprev, 1), 1);
                    fprev_compact = fprev(:, idxcon);
                    [~, idxclose] = min(sum((epscon_compact - fprev_compact).^2, 2), [], 1);
                    x0 = xprev(idxclose, :);
                    [o1, ~, o3] = runECs(problem, iECs, x0, epscon, 'Separated');
                    xopt_solution(idx, :) = reshape(o1, 1, numel(o1));
                    FuncEval = FuncEval + o3.output.funcCount + 1;
                end
                for idx = 1:ncomb
                    [o1, ~, ~] = computeallSeparatedObj( ...
                        xopt_solution(idx, :), problem.parameter, fn);
                    fopt_solution(idx, :) = reshape(o1, 1, numel(o1));
                end
            end
        end
        xopt = [xopt_anchor; xopt_solution];
        fopt = [fopt_anchor; fopt_solution];
        out.FuncEval = FuncEval;
    otherwise
        error([solver, ' not supported']);
    end

    result.problem = problem;
    result.DO.xopt = xopt;
    result.DO.fopt = fopt;
    result.DO.out = out;

end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [xopt, fopt, out] = runECs(problem, nobj, x0, epscon, fntype)
    xLast = [];
    myf = [];
    myc = [];
    myceq = [];
    switch lower(fntype)
    case 'combined'
        fun = @objfnCombinedObj;
        cfun = @constrfnCombinedObj;
        fn.combined = problem.functions.hifi_combined_exp;
        fn.nonlcon_cheap = problem.functions.hifi_nonlcon_cheap;
    case 'separated'
        fun = @objfnSeparatedObj;
        cfun = @constrfnSeparatedObj;
        fn.obj = problem.functions.hifi_obj_exp;
        fn.nonlcon_exp = problem.functions.hifi_nonlcon_exp;
        fn.nonlcon_cheap = problem.functions.hifi_nonlcon_cheap;
    end
    A = problem.lincon.A;
    b = problem.lincon.b;
    Aeq = problem.lincon.Aeq;
    beq = problem.lincon.beq;
    lb = problem.bound.xlb;
    ub = problem.bound.xub;
    
    npool = parallelPoolSize();
    opt = optimoptions('fmincon');
    opt.Algorithm = problem.optimization.fmincon.solver;
    if (problem.control.verbose > 1), opt.Display = 'iter-detailed'; end

    if ((npool > 1) && problem.functions.hifi_parallel)
        opt.UseParallel = true;
    else
        opt.UseParallel = false;
    end
    
    [xopt, fopt, exitflag, output] = fmincon( ...
        @(x) fun(x, problem.parameter, fn, nobj), x0, A, b, Aeq, beq, lb, ub, ...
        @(x) cfun(x, problem.parameter, fn, nobj, epscon), opt);
    out.exitflag = exitflag;
    out.output = output;
    
    %----1---------2---------3---------4---------5---------6---------7---------8---------9---------0
    
    function f = objfnCombinedObj(x, p, fn, nobj)
        if ~isequal(x, xLast)
            [myf, myc, myceq] = computeallCombinedObj(x, p, fn);
            xLast = x;
        end
        f = myf(:,nobj);
    end
    
    %----1---------2---------3---------4---------5---------6---------7---------8---------9---------0
    
    function f = objfnSeparatedObj(x, p, fn, nobj)
        if ~isequal(x, xLast)
            [myf, myc, myceq] = computeallSeparatedObj(x, p, fn);
            xLast = x;
        end
        f = myf(:,nobj);
    end

    %----1---------2---------3---------4---------5---------6---------7---------8---------9---------0
    
    function [c, ceq] = constrfnCombinedObj(x, p, fn, nobj, epscon)
        if ~isequal(x, xLast)
            [myf, myc, myceq] = computeallCombinedObj(x, p, fn);
            xLast = x;
        end
        if size(epscon,1) == 0
            c = myc;
        else
            epsconmat = repmat(epscon, size(myf,1), 1);
            ceps = myf - epsconmat;
            ceps(:, nobj) = 0;
            if size(myc, 1) == 0
                c = ceps;
            else
                c = [myc, ceps];
            end
        end
        ceq = myceq;
    end

    %----1---------2---------3---------4---------5---------6---------7---------8---------9---------0
    
    function [c, ceq] = constrfnSeparatedObj(x, p, fn, nobj, epscon)
        if ~isequal(x, xLast)
            [myf, myc, myceq] = computeallSeparatedObj(x, p, fn);
            xLast = x;
        end
        if size(epscon,1) == 0
            c = myc;
        else
            epsconmat = repmat(epscon, size(myf,1), 1);
            ceps = myf - epsconmat;
            ceps(:, nobj) = 0;
            if size(myc, 1) == 0
                c = ceps;
            else
                c = [myc, ceps];
            end
        end
        ceq = myceq;
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [xopt, fopt, out] = runCombinedObjNSGA2(problem, initpop)
    xLast = [];
    myf = [];
    myc = [];
    myceq = [];
    fun = @objfnCombinedObj;
    cfun = @constrfnCombinedObj;
    fn.combined = problem.functions.hifi_combined_exp;
    fn.nonlcon_cheap = problem.functions.hifi_nonlcon_cheap;
    A = problem.lincon.A;
    b = problem.lincon.b;
    Aeq = problem.lincon.Aeq;
    beq = problem.lincon.beq;
    lb = problem.bound.xlb;
    ub = problem.bound.xub;
    num_x = problem.bound.num_x;
    
    npool = parallelPoolSize();
    opt = gaoptimset(@gamultiobj);
    if ((npool > 1) && problem.functions.hifi_parallel)
        opt.Vectorized = 'off';
        opt.UseParallel = true;
    else
        opt.UseParallel = false;
        if problem.functions.hifi_vectorized,
            opt.Vectorized = 'on';
        else
            opt.Vectorized = 'off';
        end
    end
    popsize = 10*problem.optimization.nsga2.popsize;
    if size(initpop, 1) == 0
        initpopx = [];
        initpopf = [];
    else
        [initpopx, initpopf, ~] = ndSort(initpop.x, initpop.f);
        if (size(initpopx, 1) > popsize)
            initpopx = initpopx(1:popsize, :);
            initpopf = initpopf(1:popsize, :);
        end
    end
    opt.PopulationSize = popsize;
    opt.InitialPopulation = initpopx;
    opt.InitialScores = initpopf;
    opt.ParetoFraction = problem.optimization.nsga2.paretofrac;
    opt.StallGenLimit = problem.optimization.nsga2.stallgenlimit;
    opt.PlotFcns = @gaplotpareto;

    % Call gamultiobj (NSGA-II) solver
    [xopt, fopt, exitflag, output] = gamultiobj( ...
        @(x) fun(x, problem.parameter, fn), num_x, A, b, Aeq, beq, lb, ub, ...
        @(x) cfun(x, problem.parameter, fn), opt);
    out.exitflag = exitflag;
    out.output = output;

    %----1---------2---------3---------4---------5---------6---------7---------8---------9---------0

    function f = objfnCombinedObj(x, p, fn)
        if ~isequal(x, xLast)
            [myf, myc, myceq] = computeallCombinedObj(x, p, fn);
            xLast = x;
        end
        f = myf;
    end

    %----1---------2---------3---------4---------5---------6---------7---------8---------9---------0

    function [c, ceq] = constrfnCombinedObj(x, p, fn)
        if ~isequal(x, xLast)
            [myf, myc, myceq] = computeallCombinedObj(x, p, fn);
            xLast = x;
        end
        c = myc;
        ceq = myceq;
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [xopt, fopt, out] = runSeparateObjNSGA2(problem, initpop)
    objfn = problem.functions.hifi_obj_exp;
    fn.nonlcon1 = problem.functions.hifi_nonlcon_exp;
    fn.nonlcon2 = problem.functions.hifi_nonlcon_cheap;
    A = problem.lincon.A;
    b = problem.lincon.b;
    Aeq = problem.lincon.Aeq;
    beq = problem.lincon.beq;
    lb = problem.bound.xlb;
    ub = problem.bound.xub;
    num_x = problem.bound.num_x;
    
    npool = parallelPoolSize();
    opt = gaoptimset(@gamultiobj);
    if ((npool > 1) && problem.functions.hifi_parallel)
        opt.Vectorized = 'off';
        opt.UseParallel = true;
    else
        opt.UseParallel = false;
        if problem.functions.hifi_vectorized,
            opt.Vectorized = 'on';
        else
            opt.Vectorized = 'off';
        end
    end
    popsize = problem.optimization.nsga2.popsize;
    if size(initpop, 1) == 0
        initpopx = [];
        initpopf = [];
    else
        [initpopx, initpopf, ~] = ndSort(initpop.x, initpop.f);
        if (size(initpopx, 1) > popsize)
            initpopx = initpopx(1:popsize, :);
            initpopf = initpopf(1:popsize, :);
        end
    end
    opt.PopulationSize = popsize;
    opt.InitialPopulation = initpopx;
    opt.InitialScores = initpopf;
    opt.ParetoFraction = problem.optimization.nsga2.paretofrac;
    opt.StallGenLimit = 10*problem.optimization.nsga2.stallgenlimit;
    opt.PlotFcns = @gaplotpareto;

    % Call gamultiobj (NSGA-II) solver
    if (size(fn.nonlcon1, 1) ~= 0) && (size(fn.nonlcon2, 1) ~= 0)
        [xopt, fopt, exitflag, output] = gamultiobj( ...
            @(x) objfn(x, problem.parameter), num_x, A, b, Aeq, beq, lb, ub, ...
            @(x) constrfnSeparateObj(x, problem.parameter, fn), opt);
    elseif (size(fn.nonlcon1, 1) ~= 0) && (size(fn.nonlcon2, 1) == 0)
        nonlconfun = fn.nonlcon1;
        [xopt, fopt, exitflag, output] = gamultiobj( ...
            @(x) objfn(x, problem.parameter), num_x, A, b, Aeq, beq, lb, ub, ...
            @(x) nonlconfun(x, problem.parameter), opt);
    elseif (size(fn.nonlcon1, 1) == 0) && (size(fn.nonlcon2, 1) ~= 0)
        nonlconfun = fn.nonlcon2;
        [xopt, fopt, exitflag, output] = gamultiobj( ...
            @(x) objfn(x, problem.parameter), num_x, A, b, Aeq, beq, lb, ub, ...
            @(x) nonlconfun(x, problem.parameter), opt);
    elseif (size(fn.nonlcon1, 1) == 0) && (size(fn.nonlcon2, 1) == 0)
        [xopt, fopt, exitflag, output] = gamultiobj( ...
            @(x) objfn(x, problem.parameter), num_x, A, b, Aeq, beq, lb, ub, [], opt);
    end
    out.exitflag = exitflag;
    out.output = output;
    
    %----1---------2---------3---------4---------5---------6---------7---------8---------9---------0
    
    function [c, ceq] = constrfnSeparateObj(x, p, fn)
        [c1, ceq1] = feval(fn.nonlcon1, x, p);
        [c2, ceq2] = feval(fn.nonlcon2, x, p);
        c = [c1, c2];
        ceq = [ceq1, ceq2];
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [f, c, ceq] = computeallCombinedObj(x, p, fn)
    [f, c1, ceq1] = feval(fn.combined, x, p);
    if size(fn.nonlcon_cheap, 1) == 0
        c2 = []; ceq2 = [];
    else
        [c2, ceq2] = feval(fn.nonlcon_cheap, x, p);
    end
    c = [c1, c2];
    ceq = [ceq1, ceq2];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [f, c, ceq] = computeallSeparatedObj(x, p, fn)
    f = feval(fn.obj, x, p);
    if size(fn.nonlcon_exp, 1) == 0
        c1 = []; ceq1 = [];
    else
        [c1, ceq1] = feval(fn.nonlcon_exp, x, p);
    end
    if size(fn.nonlcon_cheap, 1) == 0
        c2 = []; ceq2 = [];
    else
        [c2, ceq2] = feval(fn.nonlcon_cheap, x, p);
    end
    c = [c1, c2];
    ceq = [ceq1, ceq2];
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

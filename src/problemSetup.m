%% MO-ASMO-II :: problemSetup function
% 1. Check problem structure.
% 2. Run settings function.
% 3. Save settings to the problem structure and return.
% Usage:
%  pout = problemSetup(pinp)
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function pout = problemSetup(pinp)
    declareGlobalVariables;

    % Default values
    pout = defaultProblemStructure();

    % Parsed arguments
    if isfield(pinp.functions,'settingsfun') ...
            && isa(pinp.functions.settingsfun,'function_handle')
        pout.functions.settingsfun = pinp.functions.settingsfun; end
    if isfield(pinp.functions,'hifi_obj_exp') ...
            && isa(pinp.functions.hifi_obj_exp,'function_handle')
        pout.functions.hifi_obj_exp = pinp.functions.hifi_obj_exp; end
    if isfield(pinp.functions,'hifi_nonlcon_exp') ...
            && isa(pinp.functions.hifi_nonlcon_exp,'function_handle')
        pout.functions.hifi_nonlcon_exp = pinp.functions.hifi_nonlcon_exp; end
    if isfield(pinp.functions,'hifi_combined_exp') ...
            && isa(pinp.functions.hifi_combined_exp,'function_handle')
        pout.functions.hifi_combined_exp = pinp.functions.hifi_combined_exp; end
    if isfield(pinp,'control') ...
            && isfield(pinp.control,'casefile') ...
            && (exist(pinp.control.casefile,'file') == 2)
        pout.control.casefile = pinp.control.casefile;
    else
        error('Provide full path of the casefile in the problem definition.');
    end
    if isfield(pinp,'nested')
        if isfield(pinp.nested,'outeriter') && isfield(pinp.nested,'innercase')
            pout.nested.outeriter = pinp.nested.outeriter;
            pout.nested.innercase = pinp.nested.innercase;
        else
            pout.nested.outeriter = 0;
            pout.nested.innercase = 0;
        end
    end

    % Other provided values
    % problem.bound
    if isfield(pinp,'bound')
        if isfield(pinp.bound,'xlb'), pout.bound.xlb = pinp.bound.xlb; end
        if isfield(pinp.bound,'xub'), pout.bound.xub = pinp.bound.xub; end
        if isfield(pinp.bound,'num_x'), pout.bound.num_x = pinp.bound.num_x; end
        if (pout.bound.num_x == 0)
            pout.bound.num_x = max(numel(pout.bound.xlb), numel(pout.bound.xub));
            if (pout.bound.num_x == 0), error('num_x, xlb, xub are not defined.'); end
        end
        xlb = pout.bound.xlb;
        if (min(xlb) < -1e6), error('Please scale the problem: x in [-1e6,+1e6]'); end
        if (numel(xlb) == 0), xlb = -1e6*ones(pout.bound.num_x, 1); end
        xlb(xlb<-1e6) = -1e6;
        pout.bound.xlb = xlb;
        xub = pout.bound.xub;
        if (max(xub) > 1e6), error('Please scale the problem: x in [-1e6,+1e6]'); end
        if (numel(xub) == 0), xub = 1e6*ones(pout.bound.num_x, 1); end
        xub(xub>1e6) = 1e6;
        pout.bound.xub = xub;
        if ~((pout.bound.num_x == numel(pout.bound.xlb)) ...
                && (pout.bound.num_x == numel(pout.bound.xub)))
            error('value of num_x, dimensions of xlb, and xub are not matching.');
        end
        if isfield(pinp.bound,'flb'), pout.bound.flb = pinp.bound.flb; end
        if isfield(pinp.bound,'fub'), pout.bound.fub = pinp.bound.fub; end
        if isfield(pinp.bound,'num_f'), pout.bound.num_f = pinp.bound.num_f; end
        if (pout.bound.num_f == 0)
            pout.bound.num_f = max(numel(pout.bound.flb), numel(pout.bound.fub));
            if (pout.bound.num_f == 0), error('num_f, flb, fub are not defined.'); end
        end
        flb = pout.bound.flb;
        fub = pout.bound.fub;
        if (min(flb) < -1e6), error('Please scale the problem: f in [-1e6,+1e6]'); end
        if (max(fub) > 1e6), error('Please scale the problem: f in [-1e6,+1e6]'); end
        if ((numel(flb) == 0) && (numel(fub) == 0))
            flb = zeros(1, pout.bound.num_f);
            fub = ones(1, pout.bound.num_f);
        else
            if (numel(flb) == 0), flb = -1e6*ones(1, pout.bound.num_f); end
            if (numel(fub) == 0), fub = 1e6*ones(1, pout.bound.num_f); end
        end
        pout.bound.flb = flb;
        pout.bound.fub = fub;
        if ~((pout.bound.num_f == numel(pout.bound.flb)) ...
                && (pout.bound.num_f == numel(pout.bound.fub)))
            error('value of num_f, dimensions of flb, and fub are not matching.');
        end
        if isfield(pinp.bound,'adaptive')
            pout.bound.adaptive = pinp.bound.adaptive; end
    end
    % problem.control
    [casepath, casefile] = fileparts(which(pinp.control.casefile));
    pout.control.casefile = pinp.control.casefile;
    pout.control.case = casefile;
    pout.control.path = casepath;
    pout.control.solpath = fullfile(casepath,'solution');
    pout.control.plotpath = fullfile(casepath,'plot');
    if isfield(pinp.control,'plotexport')
        pout.control.plotexport = pinp.control.plotexport; end
    if isfield(pinp.control,'verbose')
        pout.control.verbose = pinp.control.verbose; end
    if isfield(pinp.control,'randomseed')
        pout.control.randomseed = pinp.control.randomseed; end
    if isfield(pinp.control,'randomseedRTincr')
        pout.control.randomseedRTincr = pinp.control.randomseedRTincr; end
    if isfield(pinp.control,'randomgenerator')
        pout.control.randomgenerator = pinp.control.randomgenerator; end
    verbose = pout.control.verbose;
    randomseed = pout.control.randomseed;
    randomseedRTincr = pout.control.randomseedRTincr;
    randomgenerator = pout.control.randomgenerator;
    if isfield(pinp.control,'tolC')
        pout.control.tolC = pinp.control.tolC; end
    if isfield(pinp.control,'tolCEQ')
        pout.control.tolCEQ = pinp.control.tolCEQ; end
    % problem.stop
    if isfield(pinp,'stop')
        if isfield(pinp.stop,'maxiter')
            pout.stop.maxiter = pinp.stop.maxiter; end
        if isfield(pinp.stop,'residual')
            if isfield(pinp.stop.residual,'satisfaction_method')
                pout.stop.residual.satisfaction_method ...
                    = pinp.stop.residual.satisfaction_method; end
            if isfield(pinp.stop.residual,'satisfaction_continuous')
                pout.stop.residual.satisfaction_continuous ...
                    = pinp.stop.residual.satisfaction_continuous; end
            if isfield(pinp.stop.residual,'ED_max')
                pout.stop.residual.ED_max = pinp.stop.residual.ED_max; end
            if isfield(pinp.stop.residual,'ED_avg')
                pout.stop.residual.ED_avg = pinp.stop.residual.ED_avg; end
            if isfield(pinp.stop.residual,'HV_size')
                pout.stop.residual.HV_size = pinp.stop.residual.HV_size; end
            if isfield(pinp.stop.residual,'HV_ratio')
                pout.stop.residual.HV_ratio = pinp.stop.residual.HV_ratio; end
            if isfield(pinp.stop.residual,'HV_data')
                pout.stop.residual.HV_data = pinp.stop.residual.HV_data; end
            if isfield(pinp.stop.residual,'HV_eval')
                pout.stop.residual.HV_eval = pinp.stop.residual.HV_eval; end
        end
    end
    % problem.functions
    if isfield(pinp.functions,'hifi_nonlcon_cheap')
        pout.functions.hifi_nonlcon_cheap = pinp.functions.hifi_nonlcon_cheap;
    end
    if isfield(pinp.functions,'hifi_expensive')
        pout.functions.hifi_expensive = pinp.functions.hifi_expensive; end
    if isfield(pinp.functions,'hifi_vectorized')
        pout.functions.hifi_vectorized = pinp.functions.hifi_vectorized; end
    if isfield(pinp.functions,'hifi_parallel')
        pout.functions.hifi_parallel = pinp.functions.hifi_parallel; end
    % problem.lincon
    if isfield(pinp,'lincon')
        if (isfield(pinp.lincon,'A') && isfield(pinp.lincon,'b'))
            pout.lincon.A = pinp.lincon.A;
            pout.lincon.b = pinp.lincon.b;
        end
        if (isfield(pinp.lincon,'Aeq') && isfield(pinp.lincon,'beq'))
            pout.lincon.Aeq = pinp.lincon.Aeq;
            pout.lincon.beq = pinp.lincon.beq;
        end
    end
    % problem.parameter
    if isfield(pinp,'parameter')
        pout.parameter = pinp.parameter;
    end
    % problem.sampling
    if isfield(pinp,'sampling')
        if isfield(pinp.sampling,'tolerance')
            if isfield(pinp.sampling.tolerance,'ineq')
                pout.sampling.tolerance.ineq = pinp.sampling.tolerance.ineq; end
            if isfield(pinp.sampling.tolerance,'eq')
                pout.sampling.tolerance.eq = pinp.sampling.tolerance.eq; end
        end
        if isfield(pinp.sampling,'initial')
            if isfield(pinp.sampling.initial,'method')
                pout.sampling.initial.method = pinp.sampling.initial.method; end
            if isfield(pinp.sampling.initial,'number')
                pout.sampling.initial.number = pinp.sampling.initial.number; end
            if isfield(pinp.sampling.initial,'force_number')
                pout.sampling.initial.force_number = pinp.sampling.initial.force_number; end
        end
        if isfield(pinp.sampling,'update')
            if isfield(pinp.sampling.update,'explore')
                if isfield(pinp.sampling.update.explore,'method')
                    pout.sampling.update.explore.method = pinp.sampling.update.explore.method; end
                if isfield(pinp.sampling.update.explore,'number')
                    pout.sampling.update.explore.number = pinp.sampling.update.explore.number; end
                if isfield(pinp.sampling.update.explore,'w1_distance')
                    pout.sampling.update.explore.w1_distance ...
                    = pinp.sampling.update.explore.w1_distance; end
                if isfield(pinp.sampling.update.explore,'w2_disperse')
                    pout.sampling.update.explore.w2_disperse ...
                    = pinp.sampling.update.explore.w2_disperse; end
            end
            if isfield(pinp.sampling.update,'exploit')
                if isfield(pinp.sampling.update.exploit,'method')
                    pout.sampling.update.exploit.method = pinp.sampling.update.exploit.method; end
                if isfield(pinp.sampling.update.exploit,'number')
                    pout.sampling.update.exploit.number = pinp.sampling.update.exploit.number; end
            end
        end
        if isfield(pinp.sampling,'validation')
            if isfield(pinp.sampling.validation,'number')
                pout.sampling.validation.number = pinp.sampling.validation.number; end
            if isfield(pinp.sampling.validation,'method')
                pout.sampling.validation.method = pinp.sampling.validation.method; end
        end
    end
    % problem.surrogate
    if isfield(pinp,'surrogate')
        if isfield(pinp.surrogate,'method')
            pout.surrogate.method = pinp.surrogate.method; end
        if isfield(pinp.surrogate,'scale')
            pout.surrogate.scale = pinp.surrogate.scale; end
        if isfield(pinp.surrogate,'gpr')
            if isfield(pinp.surrogate.gpr,'sigma')
                pout.surrogate.gpr.sigma = pinp.surrogate.gpr.sigma; end
            if isfield(pinp.surrogate.gpr,'constantsigma')
                pout.surrogate.gpr.constantsigma = pinp.surrogate.gpr.constantsigma; end
        end
        if isfield(pinp.surrogate,'rbf')
            if isfield(pinp.surrogate.rbf,'basisfn')
                pout.surrogate.rbf.basisfn = pinp.surrogate.rbf.basisfn; end
            if isfield(pinp.surrogate.rbf,'epsilon')
                pout.surrogate.rbf.epsilon = pinp.surrogate.rbf.epsilon; end
        end
        if isfield(pinp.surrogate,'snn')
            if isfield(pinp.surrogate.snn,'trainratio')
                pout.surrogate.snn.trainratio = pinp.surrogate.snn.trainratio; end
            if isfield(pinp.surrogate.snn,'valratio')
                pout.surrogate.snn.valratio = pinp.surrogate.snn.valratio; end
            if isfield(pinp.surrogate.snn,'testratio')
                pout.surrogate.snn.testratio = pinp.surrogate.snn.testratio; end
            if isfield(pinp.surrogate.snn,'trainfcn')
                pout.surrogate.snn.trainfcn = pinp.surrogate.snn.trainfcn; end
            if isfield(pinp.surrogate.snn,'maxfail')
                pout.surrogate.snn.maxfail = pinp.surrogate.snn.maxfail; end
        end
        if isfield(pinp.surrogate,'dace')
            if isfield(pinp.surrogate.dace,'regfn')
                pout.surrogate.dace.regfn = pinp.surrogate.dace.regfn; end
            if isfield(pinp.surrogate.dace,'corrfn')
                pout.surrogate.dace.corrfn = pinp.surrogate.dace.corrfn; end
            if isfield(pinp.surrogate.dace,'theta_guess')
                pout.surrogate.dace.theta_guess = pinp.surrogate.dace.theta_guess; end
            if isfield(pinp.surrogate.dace,'theta_lb')
                pout.surrogate.dace.theta_lb = pinp.surrogate.dace.theta_lb; end
            if isfield(pinp.surrogate.dace,'theta_ub')
                pout.surrogate.dace.theta_ub = pinp.surrogate.dace.theta_ub; end
        end
    end
    % problem.invalidregion
    if isfield(pinp,'invalidregion')
        if isfield(pinp.invalidregion,'use')
            pout.invalidregion.use = pinp.invalidregion.use; end
        if isfield(pinp.invalidregion,'C')
            pout.invalidregion.C = pinp.invalidregion.C; end
        if isfield(pinp.invalidregion,'q')
            pout.invalidregion.q = pinp.invalidregion.q; end
        if isfield(pinp.invalidregion,'epsilon')
            pout.invalidregion.epsilon = pinp.invalidregion.epsilon; end
    end
    % problem.optimization
    if isfield(pinp,'optimization')
        if isfield(pinp.optimization,'solver')
            pout.optimization.solver = pinp.optimization.solver; end
        if isfield(pinp.optimization,'nsga2')
            if isfield(pinp.optimization.nsga2,'popsize')
                pout.optimization.nsga2.popsize = pinp.optimization.nsga2.popsize; end
            if isfield(pinp.optimization.nsga2,'paretofrac')
                pout.optimization.nsga2.paretofrac = pinp.optimization.nsga2.paretofrac; end
            if isfield(pinp.optimization.nsga2,'stallgenlimit')
                pout.optimization.nsga2.stallgenlimit = pinp.optimization.nsga2.stallgenlimit; end
            if isfield(pinp.optimization.nsga2,'maxgen')
                pout.optimization.nsga2.maxgen = pinp.optimization.nsga2.maxgen; end
        end
    end

    % Run settings function
    settingsfun = pout.functions.settingsfun;
    if isa(settingsfun, 'function_handle')
        try
            pout = feval(settingsfun,pout);
        catch
            error('@settings function not available.');
        end
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function p = defaultProblemStructure()
    % Problem bounds
    p.bound.xlb = [];
    p.bound.xub = [];
    p.bound.num_x = 0;
    p.bound.flb = [];
    p.bound.fub = [];
    p.bound.num_f = 0;
    p.bound.adaptive = true;

    % Control variables
    p.control.casefile = [];
    p.control.case = [];
    p.control.path = [];
    p.control.solpath = [];
    p.control.plotpath = [];
    p.control.plotexport = true;
    p.control.verbose = 1;
    p.control.randomseed = 0;
    p.control.randomseedRTincr = true;
    p.control.randomgenerator = 'twister';
    p.control.tolC = 1e-3;
    p.control.tolCEQ = 1e-3;

    % Stopping criteria
    p.stop.maxiter = 30;
    p.stop.residual.satisfaction_method = 'AND'; % ['AND'] or 'OR'
    p.stop.residual.satisfaction_continuous = 1;
    p.stop.residual.ED_max = 1e-3;
    p.stop.residual.ED_avg = 1e-3;
    p.stop.residual.HV_size = 1e-2;
    p.stop.residual.HV_ratio = 1e-2;
    p.stop.residual.HV_eval = 1000000;
    p.stop.residual.HV_data = 'highfidelity'; % ['highfidelity'], 'predicted', or 'both'

    % Function handles
    p.functions.settingsfun = [];
    p.functions.hifi_obj_exp = [];
    p.functions.hifi_nonlcon_exp = [];
    p.functions.hifi_combined_exp = [];
    p.functions.hifi_nonlcon_cheap = [];
    p.functions.hifi_expensive = true;
    p.functions.hifi_vectorized = false;
    p.functions.hifi_parallel = false;
    
    % Nested case
    p.nested.outeriter = 0;
    p.nested.innercase = 0;

    % Parameters
    p.parameter = [];

    % Linear constraints
    p.lincon.A = [];
    p.lincon.b = [];
    p.lincon.Aeq = [];
    p.lincon.beq = [];

    % Sampling
    p.sampling.tolerance.ineq = 1e-3;
    p.sampling.tolerance.eq = 1e-3;
    p.sampling.initial.method = 'LHS'; % ['LHS'], 'RANDOM', ('CCI', 'BBD' for num_x <= 10)
    p.sampling.initial.number = 10;
    p.sampling.initial.force_number = false;
    p.sampling.update.explore.method = 'LHS'; % ['LHS'], 'RANDOM'
    p.sampling.update.explore.number = 5;
    p.sampling.update.explore.w1_distance = 1e-4;
    p.sampling.update.explore.w2_disperse = 1e-8;
    p.sampling.update.exploit.method = 'FDL'; % ['FDL']
    p.sampling.update.exploit.number = 5;
    p.sampling.validation.number = 10;
    p.sampling.validation.method = 'uniform'; % ['uniform'], random

    % Surrogate model
    p.surrogate.method = 'GPR'; % ['GPR'], 'RBF', 'RBN', 'SNN', 'DACE'
    p.surrogate.scale = true;
    % GPR
    p.surrogate.gpr.sigma = 1e-3;
    p.surrogate.gpr.constantsigma = true;
    % RBF
    p.surrogate.rbf.basisfn = 'TPS'; % ['TPS'], 'Linear', 'Cubic', 'Gaussian', 'MQ', 'InvMQ'
    p.surrogate.rbf.epsilon = 1;
    % RBN
    % SNN
    p.surrogate.snn.trainratio = 8.0;
    p.surrogate.snn.valratio = 2.0;
    p.surrogate.snn.testratio = 0.0;
    p.surrogate.snn.trainfcn = 'trainlm'; % ['trainlm'], 'trainbr'
    p.surrogate.snn.maxfail = 30;
    % DACE
    p.surrogate.dace.regfn = @regpoly1;
    p.surrogate.dace.corrfn = @corrspherical;
    p.surrogate.dace.theta_guess= 1.0;
    p.surrogate.dace.theta_lb = 0.1;
    p.surrogate.dace.theta_ub = 20;

    % Invalid region model
    p.invalidregion.use = true;
    p.invalidregion.C = 0.4;
    p.invalidregion.q = 0.4;
    p.invalidregion.epsilon = 1e-6;

    % Optimization algorithm
    p.optimization.solver = 'NSGA-II'; % ['NSGA-II'], 'epsilon-Constraints'
    p.optimization.nsga2.popsize = 1000;
    p.optimization.nsga2.paretofrac = 0.15;
    p.optimization.nsga2.stallgenlimit = 20;
    p.optimization.nsga2.maxgen = 100;
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

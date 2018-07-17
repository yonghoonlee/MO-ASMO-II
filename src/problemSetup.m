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
        if isfield(pinp.bound,'flb'), pout.bound.flb = pinp.bound.flb; end
        if isfield(pinp.bound,'fub'), pout.bound.fub = pinp.bound.fub; end
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
            if isfeidl(pinp.stop.residual,'HV_eval')
                pout.stop.residual.HV_eval = pinp.stop.residual.HV_eval; end
        end
    end
    % problem.functions
    if isfield(pinp.functions,'hifi_nonlcon_cheap')
        pout.functions.hifi_nonlcon_cheap = pinp.functions.hifi_nonlcon_cheap;
    end
    if isfield(pinp.functions,'hifi_expensive')
        pout.functions.hifi_expensive = pinp.functions.hifi_expensive; end
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
    % problem.surrogate
    if isfield(pinp,'surrogate')
        if isfield(pinp.surrogate,'method')
            pout.surrogate.method = pinp.surrogate.method; end
        if isfield(pinp.surrogate,'scale')
            pout.surrogate.scale = pinp.surrogate.scale; end
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
    p.bound.flb = [];
    p.bound.fub = [];
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
    p.stop.residual.satisfaction_continuous = 3;
    p.stop.residual.ED_max = 1e-3;
    p.stop.residual.ED_avg = 1e-3;
    p.stop.residual.HV_size = 1e-3;
    p.stop.residual.HV_ratio = 1e-3;
    p.stop.residual.HV_eval = 1000000;
    p.stop.residual.HV_data = 'highfidelity'; % ['highfidelity'] or 'predicted'

    % Function handles
    p.functions.settingsfun = [];
    p.functions.hifi_obj_exp = [];
    p.functions.hifi_nonlcon_exp = [];
    p.functions.hifi_combined_exp = [];
    p.functions.hifi_nonlcon_cheap = [];
    p.functions.hifi_expensive = true;
    
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

    % Surrogate model
    p.surrogate.method = 'GPR';
    p.surrogate.scale = true;

end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

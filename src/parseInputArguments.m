%% MO-ASMO-II :: parseInputArguments function
% 1. Parse input arguments, put into problem structure and return
% Usage:
%  problem = parseInputArguments(varargin)
% Arguments:
%  * Recommended:
%   - {problem}
%   - {result} -- for restart problem. problem structure should be a field in the result structure
%  * Obsolete (but supported):
%   - {(@settings), @hifi_combined, casefile}
%   - {(@settings), @hifi_obj, (@hifi_nonlcon), casefile}
%   - {(@settings), @hifi_combined, casefile, outeriter, innercase}
%   - {(@settings), @hifi_obj, (@hifi_nonlcon), casefile, outeriter, innercase}
% Note:
%  Functions @hifi_obj, @hifi_nonlcon, @hifi_combined are assumed expensive.
%  To include cheap constraints that can be evaluated prior to the objective function evaluation,
%  please use problem structure as the only input argument.
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [problem, result, restart] = parseInputArguments(varargin)
    result = [];
    restart = false;
    
    initializeEnvironment();
    
    switch nargin
    case 1  % {problem} or {result}
        if isstruct(varargin{1})
            tmp = varargin{1};
            if isfield(tmp,'problem') && isfield(tmp,'data')
                result = tmp;
                restart = true;
                problem = tmp.problem;
            else
                problem = tmp;
            end
            if ~isfield(problem,'functions')
                error('No function provided.');
            elseif ~(isfield(problem.functions,'hifi_obj_exp') ...
                    || isfield(problem.functions,'hifi_combined_exp'))
                error('No objective function provided.');
            end
        else
            error('Wrong argument type. See Usage of runMOASMO.');
        end
    case 3  % {(@settings), @highfidelity, casefile}
        if ((isa(varargin{1}, 'function_handle') || sum(size(varargin{1})) == 0) && ...
                isa(varargin{2}, 'function_handle') && ...
                (ischar(varargin{3}) && size(varargin{3},2)))
            problem.functions.settingsfun = varargin{1};
            problem.functions.hifi_obj_exp = [];
            problem.functions.hifi_nonlcon_exp = [];
            problem.functions.hifi_combined_exp = varargin{2};
            problem.control.casefile = varargin{3};
            problem.nested.outeriter = 0;
            problem.nested.innercase = 0;
        else
            error('Wrong argument type. See Usage of runMOASMO.');
        end
    case 4  % {(@settings), @obj, (@nonlcon), casefile}
        if ((isa(varargin{1}, 'function_handle') || sum(size(varargin{1})) == 0) && ...
                isa(varargin{2}, 'function_handle') && ...
                (isa(varargin{3}, 'function_handle') || sum(size(varargin{3})) == 0) && ...
                (ischar(varargin{4}) && size(varargin{4},2)))
            problem.functions.settingsfun = varargin{1};
            problem.functions.hifi_obj_exp = varargin{2};
            problem.functions.hifi_nonlcon_exp = varargin{3};
            problem.functions.hifi_combined_exp = [];
            problem.control.casefile = varargin{4};
            problem.nested.outeriter = 0;
            problem.nested.innercase = 0;
        else
            error('Wrong argument type. See Usage of runMOASMO.');
        end
    case 5  % {(@settings), @highfidelity, casefile, outeriter, innercase}
        if ((isa(varargin{1}, 'function_handle') || sum(size(varargin{1})) == 0) && ...
                isa(varargin{2}, 'function_handle') && ...
                (ischar(varargin{3}) && size(varargin{3},2)) && ...
                isnumeric(varargin{4}) && ...
                isnumeric(varargin{5}))
            problem.functions.settingsfun = varargin{1};
            problem.functions.hifi_obj_exp = [];
            problem.functions.hifi_nonlcon_exp = [];
            problem.functions.hifi_combined_exp = varargin{2};
            problem.control.casefile = varargin{3};
            problem.nested.outeriter = varargin{4};
            problem.nested.innercase = varargin{5};
        else
            error('Wrong argument type. See Usage of runMOASMO.');
        end
    case 6  % {(@settings), @obj, (@nonlcon), casefile, outeriter, innercase}
        if ((isa(varargin{1}, 'function_handle') || sum(size(varargin{1})) == 0) && ...
                isa(varargin{2}, 'function_handle') && ...
                (isa(varargin{3}, 'function_handle') || sum(size(varargin{3})) == 0) && ...
                (ischar(varargin{4}) && size(varargin{4},2)) && ...
                isnumeric(varargin{5}) && ...
                isnumeric(varargin{6}))
            problem.functions.settingsfun = varargin{1};
            problem.functions.hifi_obj_exp = varargin{2};
            problem.functions.hifi_nonlcon_exp = varargin{3};
            problem.functions.hifi_combined_exp = [];
            problem.control.casefile = varargin{4};
            problem.nested.outeriter = varargin{5};
            problem.nested.innercase = varargin{6};
        else
            error('Wrong argument type. See Usage of runMOASMO.');
        end
    otherwise
        error('Wrong argument type. See Usage of runMOASMO.');
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function initializeEnvironment()   
    declareDebugVariables;
    
    set(0,'DefaultAxesTickLabelInterpreter','latex');
    set(0,'DefaultColorbarTickLabelInterpreter','latex');
    set(0,'DefaultLegendInterpreter','latex');
    set(0,'DefaultPolaraxesTickLabelInterpreter','latex');
    set(0,'DefaultTextInterpreter','latex');
    set(0,'DefaultTextarrowshapeInterpreter','latex');
    set(0,'DefaultTextboxshapeInterpreter','latex');
    set(0,'DefaultAxesFontSize', 14);
    set(0,'DefaultColorbarFontSize', 12);
    set(0,'DefaultLegendFontSize', 12);
    set(0,'DefaultPolaraxesFontSize', 14);
    set(0,'DefaultTextFontSize', 14);
    set(0,'DefaultTextarrowshapeFontSize', 14);
    set(0,'DefaultTextboxshapeFontSize', 14);
    set(0,'DefaultUibuttongroupFontSize', 14);
    set(0,'DefaultUicontrolFontSize', 14);
    set(0,'DefaultUipanelFontSize', 14);
    set(0,'DefaultUitableFontSize', 14);
    set(0,'defaultFigurePosition', [680 620 560 320]);
    
    fg_debug1 = figure('Color',[1 1 1],'Position',[10 720 560 320],'Visible','off');
    fg_debug2 = figure('Color',[1 1 1],'Position',[460 710 560 320],'Visible','off');
    fg_debug3 = figure('Color',[1 1 1],'Position',[910 700 560 320],'Visible','off');
    fg_debug4 = figure('Color',[1 1 1],'Position',[1360 690 560 320],'Visible','off');
    fg_debug5 = figure('Color',[1 1 1],'Position',[30 360 560 320],'Visible','off');
    fg_debug6 = figure('Color',[1 1 1],'Position',[480 350 560 320],'Visible','off');
    fg_debug7 = figure('Color',[1 1 1],'Position',[930 340 560 320],'Visible','off');
    fg_debug8 = figure('Color',[1 1 1],'Position',[1380 330 560 320],'Visible','off');
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

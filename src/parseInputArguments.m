%% MO-ASMO-II :: parseInputArguments function
% 1. Parse input arguments, put into problem structure and return
% Usage:
%  problem = parseInputArguments(varargin)
% Arguments:
%  {problem}
%  {(@settings), @hifi_combined, casefile}
%  {(@settings), @hifi_obj, (@hifi_nonlcon), casefile}
%  {(@settings), @hifi_combined, casefile, outeriter, innercase}
%  {(@settings), @hifi_obj, (@hifi_nonlcon), casefile, outeriter, innercase}
% Note:
%  Functions @hifi_obj, @hifi_nonlcon, @hifi_combined are assumed expensive.
%  To include cheap constraints that can be evaluated prior to the objective function evaluation,
%  please use problem structure as the only input argument.
%
% Multiobjective Adaptive Surrogate Modeling-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function [problem, result, restart] = parseInputArguments(varargin)
    result = [];
    restart = false;
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

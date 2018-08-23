%% MO-ASMO-II :: prepRandomSeed function
% 1. Prepare random number generation
% 2. If needed, use incremental random seed
% Usage:
%  prepRandomSeed();
%
% Multi-Objective Adaptive Surrogate Model-based Optimization (MO-ASMO) Code :: version II
% Link: https://github.com/yonghoonlee/MO-ASMO-II
% Contact: ylee196@illinois.edu, yonghoonlee@outlook.com
% Copyright (c) 2018, Yong Hoon Lee. All rights reserved. (See the LICENSE file)

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

function prepRandomSeed()
    declareGlobalVariables;

    if randomseedRTincr
        rng(randomseed, randomgenerator);
        randomseed = randomseed + 1;
        if (verbose == 2) % debug level verbosity
            disp(['Restart random generator...Seed:', ...
            num2str(randomseed),',Type:',randomgenerator,'.']); end
    else
        test = rng;
        if (test.Seed ~= randomseed) || ~strcmp(test.Type,randomgenerator)
            rng(randomseed, randomgenerator);
            if (verbose == 2) % debug level verbosity
                disp(['Restart random generator...Seed:', ...
                num2str(randomseed),',Type:',randomgenerator,'.']); end
        end
    end
end

%--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0

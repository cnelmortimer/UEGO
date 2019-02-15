function [sass_optimizer, sass_config] = SASS()
%SASS A stochastic hillclimber for real spaces
%  This module is expected to be used within UEGO
    sass_optimizer = @SASS_optimizer;
    sass_config = buildConfig();
end

% Internal auxiliary functions:
function config = buildConfig()
    config.name = 'SASS';
end

function [new_point, new_value] = SASS_optimizer(start_point, radius, initial_val, bounds, func, max_steps, optimizer_config)
    new_point = [];
    new_value = 0;
    disp(optimizer_config.name);
end

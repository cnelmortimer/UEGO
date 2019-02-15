function [local_optimizer, optim_config] = LoadUEGOOptimizer(name)
% LoadOptimizer  Returns a valid local optimizer for UEGO. It also loads
% the default configuration that the optimizer requires (when necessary).
    addpath('./LocalOptimizers/')
    if exist('name', 'var')
        switch name
            case 'SASS'
                [local_optimizer, optim_config] = SASS();
            otherwise
                disp('Warning: Unknown local optimizer. Loading SASS by default');
                [local_optimizer, optim_config] = SASS();
        end
    else
        disp('Warning: No local optimizer requested. Loading SASS by default');
        [local_optimizer, optim_config] = SASS();
    end
end

function [fmc_optimizer, fmc_config] = FMC()
%FMC Is a wrapper over the Matlab's FMinCon in its interior-point flavour.
%  It is prepared to be used within UEGO. The Optimization Toolbox is required
    fmc_optimizer = @FMC_optimizer;
    fmc_config = buildConfig();
end

% Internal auxiliary functions:
function config = buildConfig()
    config.name = 'FMC'; % Default configuration to work with
    config.lazyness = 0.1; % Do not try to restart fmincon if the number of allowed evaluations is <= this
    config.maxNonImprovingRestarts = 20; % After restarting fmincon 20 times without improving the best known, surrender!
    config.equalityThr = 1.0000e-05; % When do you think that two real are too similar?
    config.options = optimoptions('fmincon','Algorithm','interior-point');
    config.options.Display = 'off';
end

function [c,ceq] = MyNonLnrConstr(x, myStartPoint, radius) % See: https://www.mathworks.com/help/optim/ug/nonlinear-constraints.html
    c = norm(x-myStartPoint) - radius;
    ceq = [];
end

function newStart = GetNewStartingPoint(center, radius, bounds)
    compass = (bounds(:,2) - bounds(:, 1)).*rand(size(bounds, 1), 1) + bounds(:, 1); % (max-min)*[0,1] + min; Random point in the search space
    dir = (compass - center); % Vector from the current point to the compass one (the window is moving!)
    dist = norm(dir); % Distance between the current point and the compass one
    if(dist<=radius) % The new point is already in range
        newStart = compass;
    else
        newStart = center + radius*rand()*(dir/dist); % Move from the center
        newStart(newStart < bounds(:, 1)) = bounds(newStart < bounds(:, 1), 1); % If too small: saturate
        newStart(newStart > bounds(:, 2)) = bounds(newStart > bounds(:, 2), 2); % If too big: saturate
    end
end

function [x, current_value] = FMC_optimizer(start_point, radius, initial_val, bounds, func, max_evals, config)
    current_value = initial_val;
    x = start_point;
    x0 = start_point;
    fails = 0;
    config.options.MaxFunctionEvaluations = max_evals;
    
    while config.options.MaxFunctionEvaluations > (config.lazyness * max_evals) % Do not restart fmincon if allowed evals <= lazyness% of MaxEv.
        [sol, val, flag, output] = fmincon(func, x0, [], [], [], [], bounds(:,1), bounds(:,2), @(x)MyNonLnrConstr(x, x0, radius), config.options); % Do not go too far for now
        if val<current_value
            x = sol;
            current_value = val;
            fails = 0; % We found a better solution than the previous one!
        else
            fails = fails + 1; % We could not... another fail :-(
        end
        if fails < config.maxNonImprovingRestarts && flag ~= 0 && config.options.MaxFunctionEvaluations > output.funcCount % As long as I can keep trying
            config.options.MaxFunctionEvaluations = config.options.MaxFunctionEvaluations - output.funcCount; % Updating the counter!
            if norm(sol-x0)>config.equalityThr && abs(norm(sol - x0) - radius) <= config.equalityThr
                x0 = sol; % If you moved && you finished at the border of the radius: Restart from the same point (move the window!)
            else % Otherwise, look for a different start point
                x0 = GetNewStartingPoint(x, radius, bounds); % The center is the current point whatever it is: the window is moving!
                if config.options.MaxFunctionEvaluations > 0 % Let us check if we randomly found a better point within the region!
                    val = func(x0);
                    if val<current_value
                        x = x0;
                        current_value = val; % fails = 0; % FMinCon did not affect here, so do not change the fail count 
                    end
                    config.options.MaxFunctionEvaluations = config.options.MaxFunctionEvaluations - 1; % Another evaluation has been done :-/
                end
            end
        else
            break;
        end
    end
end

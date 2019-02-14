function [config, spec_list, spec_radii, spec_values] = UEGO(evals, levels, max_spec_num, min_r, bounds, func)
% UEGO  Implementation of Universal Evolutionary Global Optimizer (UEGO), 
% which is a stochastic population-based (memetic) optimization algorithm.
% [config, spec_list, spec_radii, spec_values] = UEGO(evals, levels, 
%   max_spec_num, min_r, bounds, func)
% Implementation by: N.C. Cruz (2019). Supervision by: J.L. Redondo, 
%   J.D. Alvarez, M. Berenguel and P.M. Ortigosa (UAL, Spain)
% See [M. Jelasity, P.M. Ortigosa & I. Garcia. UEGO, an abstract clustering 
% technique for multimodal global optimization. Journal of Heuristics, 7(3), 
% 215-233, 2001.] for further information about the method.
    config = ComputeConfig(evals, levels, max_spec_num, min_r, bounds);
    if ~isempty(config)
        [spec_list, spec_radii, spec_values] = Init_Spec_List(config.radii(1), bounds, func); % Level 1
        if config.n(1)>0
           disp('Local optimization pending'); % Optimize...
        end
        for i=2:1:levels
            [spec_list, spec_radii, spec_values] = Create_Species(spec_list, spec_radii, spec_values, bounds, func, config.new(i), config.radii(i));
            [spec_list, spec_radii, spec_values] = Fuse_Species(spec_list, spec_radii, spec_values, config.radii(i));
            [spec_list, spec_radii, spec_values] = Shorten_Spec_List(spec_list, spec_radii, spec_values, max_spec_num);
            disp(['Level: ' num2str(i) ' -> Local optimization pending']); % Optimize...
            [spec_list, spec_radii, spec_values] = Fuse_Species(spec_list, spec_radii, spec_values, config.radii(i));
        end
        [spec_values, indices] = sort(spec_values, 'ascend'); % Sort the final list -> Lower is better
        spec_radii = spec_radii(indices);
        spec_list = spec_list(:, indices);
    end
end

% Internal auxiliary functions:
function val=speed(dim, r)
    V = [0.0, 0.25, 2/(3*pi), 3/16, 8/(15*pi), 10/64, 16/(35*pi), 35/256, 128/(315*pi), 63/512, 256/(693*pi), 231/2048, 1024/(3003*pi), 858/8192]; % Speed constants
    if dim < length(V)
        val = V(dim+1)*r; % WARNING: To access the same element as the C version
    else
       val = r / sqrt( pi * (2 * (dim+1)) );
    end
end

function config = ComputeConfig(evals, levels, max_spec_num, min_r, bounds)
    R_1 = norm(bounds(:,1) - bounds(:,2)); % Diameter of the search space
    config.radii = R_1 * ( (min_r / R_1) .^ ((0:1:(levels-1))/(levels-1))); % The radius of each level [1, L]
    
    config.new = 3*max_spec_num*ones(1, levels); % Func. evaluations for spec creation
    config.new(1) = 0; % (No new specs at the first level)
    
    evals1 = evals - (levels-1)*3*max_spec_num;
    dim = size(bounds, 1); % Number of dimensions
    if evals1>=0
        config.n = zeros(1, levels); % Func. evaluations for local optimization
        if levels<2
            config.n(1) = evals; % UEGO will not affect
        else
            config.n(1) = 0;
            accum = 0;
            for i=2:1:levels
                val = max_spec_num*R_1/speed(dim, config.radii(i));
                accum = accum + fix(val); % Warning: Discarding decimals as in the C version!!
                config.n(i) = val;
            end
            config.n(2:levels) = fix( config.n(2:levels)*(evals1/accum) ); % threshold (evals1/accum) % Warning: Discarding decimals as in the C version!!
        end
    else
        disp('Configuration error: All evaluations (N) too small, parameters not set.');
        config = []; 
    end
end

function [spec_list, spec_radii, spec_values] = Init_Spec_List(R_1, bounds, func)
    spec_list = (bounds(:,2)-bounds(:,1)) .* rand(size(bounds, 1), 1) + bounds(:,1); % Set of species (matrix). Each column is a solution point
    spec_radii = R_1; % The radius of the first species is equal to the diameter of the search space
    spec_values = func(spec_list(:, 1)); % Call the objective function
end

function new_point = RandomPoint(center, radius, bounds) % Function to randomly generate a new point from a species
    compass = (bounds(:,2) - bounds(:, 1)).*rand(size(bounds, 1), 1) + bounds(:, 1); % (max-min)*[0,1] + min; Random point in the search space
    dir = (compass - center);
    dist = norm(dir);
    if(dist<=radius) % The new point is already in range
        new_point = compass;
    else
        new_point = center + radius*rand()*(dir/dist); % Move from the center
        new_point(new_point < bounds(:, 1)) = bounds(new_point < bounds(:, 1), 1); % If too small: saturate
        new_point(new_point > bounds(:, 2)) = bounds(new_point > bounds(:, 2), 2); % If too big: saturate
    end
end

function [spec_list, spec_radii, spec_values] = Create_Species(spec_list, spec_radii, spec_values, bounds, func, evals, current_rad)
    num_species = size(spec_list, 2); % The number of columns is equal to the length of the list of species
    num_points = fix( (sqrt(1+8*(evals / num_species)) - 1)/2 ); % The number of points to generate per species
    
    pool = zeros(size(bounds, 1), num_points); % Where to save the new points of each species
    pool_values = zeros(1, num_points);
    
    for i=1:1:num_species % For each species:
        for j=1:1:num_points
            pool(:, j) = RandomPoint(spec_list(:, i), spec_radii(i), bounds); % New point
            pool_values(j) = func(pool(:, j));
        end
        inserted = []; % Have any of the new points been already inserted ?
        for j=1:1:num_points % Making pairs
            if pool_values(j) < spec_values(i) % Replace the center of species 'i' (but not the radius)
                spec_list(:, i) = pool(:, j);
                spec_values(i) = pool_values(j);
            end
            for k=j+1:1:num_points
                middle = (pool(:, j) + pool(:, k))/2;
                middle_val = func(middle);
                if middle_val < spec_values(i) % Replace the center of species 'i' (but not the radius)
                    spec_list(:, i) = middle;
                    spec_values(i) = middle_val;
                end
                if middle_val > pool_values(j) && middle_val > pool_values(k) % The middle is worse than the values of the pair
                    if ~ismember(j, inserted)
                        spec_list = [spec_list pool(:, j)];
                        spec_radii = [spec_radii current_rad]; % The radius of the new species is the current one, not the original one
                        spec_values = [spec_values pool_values(j)];
                        inserted = [inserted j]; % Do not add duplicates
                    end
                    if ~ismember(k, inserted)
                        spec_list = [spec_list pool(:, k)];
                        spec_radii = [spec_radii current_rad];
                        spec_values = [spec_values pool_values(k)];
                        inserted = [inserted k];
                    end
                end
            end
        end
    end
end

function [spec_list, spec_radii, spec_values] = Fuse_Species(spec_list, spec_radii, spec_values, current_rad)
    num_species = size(spec_list, 2); % The number of columns is equal to the length of the list of species
    for i=1:1:(num_species - 1)
       for j=i+1:1:num_species
           distance = norm(spec_list(:, i) - spec_list(:, j));
           if distance < current_rad % Fuse both species
               if spec_values(i) < spec_values(j) % i is better -> new center
                  kill_focus = j;
                  keep_focus = i;
               else % j is better (or equal) -> new center
                   kill_focus = i;
                   keep_focus = j;
               end
               spec_list(:, kill_focus) = []; % Remove the discarded center
               spec_values(kill_focus) = []; % Remove the discarded value
               if(spec_radii(i)>spec_radii(j)) % The radius will be the larger one
                   spec_radii(keep_focus) = spec_radii(i);
               else
                   spec_radii(keep_focus) = spec_radii(j);
               end
               spec_radii(kill_focus) = [];
           end
       end
    end
end

function [spec_list, spec_radii, spec_values] = Shorten_Spec_List(spec_list, spec_radii, spec_values, max_spec_num)
    num_species = size(spec_list, 2); % The number of columns is equal to the length of the list of species
    if num_species > max_spec_num % We have too many species:
        num_to_kill = num_species - max_spec_num;
        [~ , indices] = sort(spec_radii, 'ascend'); % Higher level species (i.e., LOWER RADIUS) are deleted first
        for i=1:1:num_to_kill
            spec_list(:, indices(i)) = [];
            spec_radii(indices(i)) = [];
            spec_values(indices(i)) = [];
        end
    end
end

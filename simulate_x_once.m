%
% By Junior R. Ribeiro, Oct 16, 2021, jrodrib@usp.br
%
% This function simulates the state variable once.
%
% struct2 = simulate_x_once(Struct, x0, get_only_x)
%
% The input parameters are 
% (A) Struct: from "parse_mmjls" or "validate_mmjls",
% (B) x0(or empty[]): (nx1) the initial value of x. If empty, it is
%      generated randomly.
% (C) get_only_x: (1x1) If true, it returns a structure with only the "x"
%      field. It is useful for doing a Monte-Carlo simulation.
%      Otherwise, other informations are returned too in their respective
%      fields, regarding the continuous/discrete domain and Markovian states,
%      as well as the norm of x through time.
%
%
% The Struct returned is explained below. Some NaN values are introduced to
% make plotting easy. You can type "plot(struct.c_domain, struct.norm_x)",
% for instance, to get only the c-time path.
%     struct.x(:, time)
%     struct.c_domain(time); they are NaN where the system goes in d-time.
%     struct.d_domain(time); they are NaN where the system goes in c-time.
%     struct.c_states(time); they are NaN where the system goes in d-time.
%     struct.d_states(time); they are NaN where the system goes in c-time.
%     struct.norm_x(time)
%



function ANS = simulate_x_once(Struct, x0, get_only_x)
% input validations
assert(nargin >= 2,'Please, give at least 2 parameters.')
assert(isfield(Struct,'time_domain'),'You must to "define_domain(...)" first.');
n = size(Struct.Ac, 1);
if isempty(x0)
    x0 = randn(n,1);
end
if nargin < 3
    get_only_x = false;
end
%

dt = Struct.time_domain.c_sampling_time;
discr_dt = Struct.time_domain.d_sampling_time;

if get_only_x && discr_dt > dt
    discr_dt = dt;
end
t_max = Struct.time_domain.max_time;

% pre-allocating domain
domain = 0:dt:t_max;
states = domain * 0;

% initial condition;
x = zeros(n, numel(domain));
x(:,1) = x0;

% first Markov state
global_counter = 1; % counter to pass through all the domain.
states(global_counter) = find(rand < cumsum(Struct.init_distrib), 1);

THEEND = false;

while ~THEEND
    %%% CONTINUOUS-TIME BRANCH
    
    % exponential arrival time (length of stay in continuous-time domain)
    exp_time = exprnd(1/Struct.sigma);
    
    % continuous-time Markov chain
    ctime_mkch = c_markovch(Struct.RateMatrix, exp_time, ...
        abs(states(global_counter)), [], true); % struct
    
    % sub-domains, one for each jump of the Markov chain
    branches_domain = get_branch_domain(ctime_mkch, dt); % cell array
    
    % for each sub-domain
    for branch_counter = 1:numel(branches_domain)
        branch = branches_domain{branch_counter}; % back to array
        
        % for each point in the sub-domain
        for k = 1:numel(branch)
            global_counter = global_counter + 1;
            if global_counter > numel(domain) || domain(global_counter) >= t_max
                THEEND = true;
                break;
            end
            
            % updating Markov states and x
            states(global_counter) = ctime_mkch.states(1 + branch_counter);
            x(:,global_counter) = Struct.Ac_d(:,:,abs(states(global_counter-1))) * ...
                x(:, global_counter-1);
        end
        if THEEND
            break;
        end
    end
    if THEEND
        break;
    end
    
    %%% DISCRETE-TIME BRANCH
    
    % the "mu" distributed discrete random variable
    % (length of stay in discrete domain)
    n_discrete_steps = find(rand < cumsum(Struct.mu), 1) - 1;
    
    % the Markov states of the discrete-time Markov chain
    dtime_mkch = d_markovch(Struct.ProbMatrix, n_discrete_steps,...
        states(global_counter)); % array
    
    for k = 1:n_discrete_steps
        global_counter = global_counter + 1;
        if global_counter > numel(domain) || domain(global_counter) >= t_max
            THEEND = true;
            break;
        end
        
        % moving forward the remaining domain as the d-spacing orders
        domain(global_counter:end) = domain(global_counter:end) + discr_dt;
        
        % signalizing the discrete domain with a minus sign
        domain(global_counter) = - domain(global_counter);
        
        % updating Markov states and x
        states(global_counter) = - dtime_mkch(k);
        x(:,global_counter) = Struct.Ad(:,:,abs(states(global_counter-1))) * ...
            x(:, global_counter-1);
    end
    
end


% retrieving data
if ~get_only_x
    x(:, domain > t_max) = [];
    states(domain > t_max) = [];
    domain(domain > t_max) = [];
end

ANS.x = x;
if get_only_x
    return
end

c_domain = domain;
d_domain = -domain;
c_domain(c_domain<0) = nan;
d_domain(d_domain<0) = nan;

c_states = states;
d_states = -states;
c_states(c_states<0) = nan;
d_states(d_states<0) = nan;

ANS.c_domain = c_domain;
ANS.d_domain = d_domain;
ANS.c_states = c_states;
ANS.d_states = d_states;

% computing the norm of x
ANS.norm_x = 1:numel(domain);
for k = 1:numel(domain)
    ANS.norm_x(k) = norm(x(:,k));
end

% to get the initial point of the discrete domain
% that was ignored until now
for k = 2:numel(domain)-1
    if all(isnan(d_domain(k:k+1))==[1,0]) % ex: [nan, 3.2]
        % get the last point of the continuous domain
        ANS.d_domain(k) = ANS.c_domain(k);
    end
    if all(isnan(c_domain(k:k+1))==[1,0]) % ex: [nan, 3.2]
        % get the last point of the discrete domain
        ANS.c_domain(k) = ANS.d_domain(k);
    end
end
ANS.d_domain(1) = nan;
end

function domain = get_branch_domain(ctime_mkch, dt)
domain = {};
for branch = 1:numel(ctime_mkch.cumulative_times) - 1
    domain{branch} = ctime_mkch.cumulative_times(branch):dt:ctime_mkch.cumulative_times(branch+1);
    current_domain = domain{branch};
    ctime_mkch.cumulative_times(branch+1) = current_domain(end);
end
end
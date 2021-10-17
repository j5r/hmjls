%
% By Junior R. Ribeiro, Oct 16, 2021, jrodrib@usp.br
%
% This function do a Monte Carlo simulation for the state variable "x'. 
% 
% !! This function just makes calls to the function "simulate_x_once". For
% more details, please see the help of <a
% href="matlab:web('simulate_x_once.m')">simulate_x_once</a>.
%
% struct2 = simulate_x_montecarlo(Struct, x0, MC)
%
% The input parameters are
% (A) Struct: from "parse_mmjls" or "validate_mmjls",
% (B) x0: the initial state (you can let it empty x0=[] to be generated
%        automatically),
% (D) MC: the number of Monte Carlo simulations to be done.
%
% The struct returned follows the structure:
%      struct.MC = MC
%      struct.x(:, time)
%      struct.norm_x(time)
%      struct.domain(time)      
%

function ANS = simulate_x_montecarlo(Struct, x0, MC)
% input validations
assert(nargin==3,'Please, give the 3 parameters needed.');
assert(MC > 2,'"MC" must be > 2. Usually it is a big number');
n = size(Struct.Ac, 1);
if isempty(x0)
    x0 = randn(n,1);
end

% pre-allocating 
sum_x = simulate_x_once(Struct, x0, true);
sum_x.x = sum_x.x / MC;

w_bar = waitbar(0,{'Monte Carlo simulation is running...',...
    'Do not close this wait bar!'});
% running Monte Carlo
for mc = 2:MC
    if ~mod(mc,17)
        percent = mc / MC;
        waitbar(percent, w_bar);
    end
    new_x = simulate_x_once(Struct, x0, true);
    sum_x.x = sum_x.x +  new_x.x / MC;    
end
close(w_bar);

% retrieving data
ANS.MC = MC;
ANS.x = sum_x.x;
ANS.norm_x = 1:size(sum_x.x,2);
for k = 1:size(sum_x.x, 2)
    ANS.norm_x(k) = norm(sum_x.x(:,k));
end
ANS.domain = 0:Struct.time_domain.c_sampling_time:Struct.time_domain.max_time;
end
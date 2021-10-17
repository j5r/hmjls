%
% By Junior R. Ribeiro, Oct 8, 2021, jrodrib@usp.br
%
% This function do a Monte Carlo simulation for the 2nd moment X of the MMJLS 
% for a set of arrival times t_r. This process is referred as Y. 
% 
% !! This function just makes calls to the function "simulate_Yr_once". For
% more details, please see the help of <a
% href="matlab:web('simulate_Yr_once.m')">simulate_Yr_once</a>.
%
% struct2 = simulate_Yr_montecarlo(Struct, x0, n_switch_times, MC)
%
% The input parameters are
% (A) Struct: from "parse_mmjls" or "validate_mmjls",
% (B) x0: the initial state (you can let it empty x0=[] to be generated
%        automatically),
% (C) n_switch_times: the number of t_r for which Y shall be computed, i.e.,
%        Y(:,:,nmarkov, 1:n_switch_times+1)
% (D) MC: the number of Monte Carlo simulations to be done.
%
% The struct returned follows the same structure as explained at <a
% href="matlab:web('simulate_Yr_once.m')">simulate_Yr_once</a>.
%

function ANS = simulate_Yr_montecarlo(Struct, x0, n_switch_times, MC)
if nargin == 0
    ANS = [];
    return
end
n = size(Struct.Ac,1); % the length of the vector x0
if isempty(x0)
    x0 = randn(n, 1);
end

ANS = simulate_Yr_once(Struct, x0, n_switch_times);
ANS.values = ANS.values / MC;

w_bar = waitbar(0,{'Monte Carlo simulation is running...',...
    'Do not close this wait bar!'});
% running Monte Carlo
for mc = 2:MC
    if ~mod(mc,17)
        percent = mc / MC;
        waitbar(percent, w_bar);
    end
Y = simulate_Yr_once(Struct, x0, n_switch_times);
ANS.values = ANS.values + Y.values / MC;
end
close(w_bar);

ANS.MC = MC;
end
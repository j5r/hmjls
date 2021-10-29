%
% By Junior R. Ribeiro, Oct 8, 2021, jrodrib@usp.br
%
% This function simulates a continuous-time Markov Chain.
%
% struct = c_markovc(RateMatrix, t_max, init_state, init_distrib, truncate)
% 
% The input parameters are
% (A) RateMatrix. It is a rate matrix, an infinitesimal generator of probability 
%       matrices. 
% (B) t_max that is the time until which you  want to simulate the chain,
% (C) init_state(optional) is the first state of the chain. If init_state is
%       not given, it will be generated randomly.
% (D) init_distrib(optional or empty []) is the initial distribution of the 
%       chain. If you give an initial distribution, init_state will be ignored.
% (E) truncate(optional). If true, it truncates the last time to be equal to t_max.
%
% The return is a struct with the fields explained below.
%   struct.states: (1x?) is the states themselves
%   struct.times: (1x?) are the staying times of each state
%   struct.cumulative_times: (1x?) are the cumulative sum of struct.times.
%

function ANS = c_markovch(RateMatrix, t_max, init_state, init_distrib, truncate)
if nargin == 0
    ANS = [];
    fprintf('\tc_markovc(RateMatrix, t_max, initial_state, init_distrib, truncate)');
    fprintf('\n\tsimulates a Markov chain for time in [0, t_max]');
    fprintf('\n\n\tIf you give an initial_distrib, initial_state will be ignored.');
    fprintf('\n\tIf truncate==true, the last cumulative time will be t_max.\n');
    fprintf('\n\tType "help c_markovch" for more info.\n')
    return
end

if nargin < 3
    init_state = randi([1, size(RateMatrix,1)]);
end

if nargin > 3 && ~isempty(init_distrib)
    assert( size(RateMatrix, 1)==numel(init_distrib),...
        'init_distrib size and RateMatrix size do not match.')
    init_state = find(rand < cumsum(init_distrib),1);
end

if nargin < 5
    truncate = false;
end

msg = sprintf('init_state must be between 1 and %d',size(RateMatrix,1));
assert(init_state >= 1, msg)
assert(init_state <= size(RateMatrix,1), msg)

states = init_state;

% the 1st arrival time
arrival_times = [0, exprnd(-1/RateMatrix(init_state,init_state))];
% the probability matrix: it is proportional to the RateMatrix
probMtrix = abs(RateMatrix);
for k=1:size(probMtrix, 1)
    probMtrix(k,:) = probMtrix(k,:) / probMtrix(k,k);
    probMtrix(k,k) = 0;
end
%
next_state = find(rand < cumsum(probMtrix(init_state,:)), 1);
states = [states, next_state];
%
cum_sum_times = sum(arrival_times);
while cum_sum_times < t_max    
    init_state = next_state;
    new_time = exprnd(-1/RateMatrix(init_state,init_state));
    arrival_times = [arrival_times, new_time];
    cum_sum_times = cum_sum_times + new_time;    
    next_state = find(rand < cumsum(probMtrix(init_state,:)), 1);
    states = [states, next_state];
end

ANS.states = states;
ANS.times = arrival_times;
ANS.cumulative_times = cumsum(arrival_times);
if truncate% truncating the last time
    difference = ANS.cumulative_times(end) - t_max;
    ANS.cumulative_times(end) = t_max;
    ANS.times(end) = ANS.times(end) - difference;
end
end

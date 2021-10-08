%
% By Junior R. Ribeiro, Oct 8, 2021, jrodrib@usp.br
%
% This function simulates a discrete-time Markov Chain.
%
% MarkovStates = d_markovch(ProbMatrix, t_max, init_state, init_distrib)
%
% The input parameters are
% (A) ProbMatrix, a row-wise stochastic matrix. It is the probability matrix.
% (B) t_max that is the time until which you want to simulate the chain,
% (C) init_state(optional) is the first state of the chain. If it is not
%         given, it will be generated randomly.
% (D) init_distrib(optional) is the initial distribution of the chain. 
%         If you give an initial distribution, init_state will be ignored. 
%
% The return is a vector of states.
%

function MarkovStates = d_markovch(ProbMatrix, t_max, init_state, init_distrib)
if nargin == 0
    MarkovStates = [];
    fprintf('\tMarkovStates = d_markovch(ProbMatrix, t_max, init_state, init_distrib)');
    fprintf('\n\tsimulates a Markov chain for time = [0,...,t_max]\n');
    fprintf('\n\tType "help d_markovch" for more info.\n')
    return
end

if nargin < 2
    init_state = randi([1, size(ProbMatrix,1)]);
end

if nargin > 3 && ~isempty(init_distrib)
    assert( size(RateMatrix, 1) == numel(init_distrib),...
        'init_distrib size and ProbMatrix size do not match.')
    init_state = find(rand < cumsum(init_distrib), 1);
end


msg = sprintf('init_state must be between 1 and %d', size(ProbMatrix,1));
assert(init_state >= 1, msg)
assert(init_state <= size(ProbMatrix,1), msg)

states = zeros(1,t_max+1);
states(1) = init_state;
for k=1:t_max
    states(k+1) = find(rand < cumsum(ProbMatrix(states(k),:)), 1);
end
MarkovStates = states;
end
%
% By Junior R. Ribeiro, Oct 8, 2021, jrodrib@usp.br
%
% This function generates randomly a set of parameters for the MMJLS.
%
% Struct = generate_parameters(n_markov, sigma)
%
% The input parameters are
% (A) n_markov(optional): the number of Markov states.
% (B) sigma(optional): exponential rate for time in continuous domain.
%
% The return is a struct with all parameters generated and others computed.
% This function uses the "parse_mmjls" function to build the struct.
%
% For more details, see the documentation of <a href="matlab:web('parse_mmjls.m')">parse_mmjls</a>.
%

function Struct = generate_parameters(n_markov, sigma)
if nargin == 1
    % from 0.2 to 10
    sigma = randi([2,100]) / 10;
elseif nargin == 0 || isempty(n_markov)
    % nmarkov from 2 to 10
    n_markov = randi([2,10]);
    % from 0.2 to 10
    sigma = randi([2,100]) / 10;
end

% dimension of x
n = randi([2,10]);

Ad = randn(n, n, n_markov) * randn;
Ac = randn(n, n, n_markov) * randn;

init_distrib = rand(n_markov, 1);
init_distrib = init_distrib/sum(init_distrib);

ProbMatrix = rand(n_markov, n_markov);
RateMatrix = randi([0, randi([3,100])], n_markov, n_markov) / randi([1,100]);
RateMatrix = RateMatrix - diag(diag(RateMatrix));

% adjusting ProbMatrix and RateMatrix
for i = 1:n_markov
    ProbMatrix(i,:) = ProbMatrix(i,:) / sum(ProbMatrix(i,:));
    RateMatrix(i,i) = -sum(RateMatrix(i,:));
end

n_mu = randi([2,100]);
mu = abs( randi([1,50], n_mu, 1) .* randn(n_mu,1) );
mu = mu/sum(mu);

Struct = parse_mmjls(Ac, Ad, RateMatrix, ProbMatrix, init_distrib, mu, sigma);
end
%
% By Junior R. Ribeiro, Oct 8, 2021, jrodrib@usp.br
%
% This function generates randomly a set of parameters for the MMJLS. It
% returns a struct with the parameters of a stable MMJLS.
%
% Struct = generate_parameters_with_T_stable(n_markov, sigma)
%
% The input parameters are
% (A) n_markov(optional): the number of Markov states.
% (B) sigma(optional): exponential rate for time in continuous domain.
%
% The return is a struct with all parameters generated and others computed.
% This function uses the "generate_parameters" function to generate the
% parameters.
% This function uses the "parse_mmjls" function to build the struct.
%
% For more details, see the documentation of <a href="matlab:web('parse_mmjls.m')">parse_mmjls</a>.
%

function Struct = generate_parameters_with_T_stable(n_markov, sigma)
if nargin == 1
    % from 0.2 to 10
    sigma = randi([2,100]) / 10;
elseif nargin == 0 || isempty(n_markov)
    % nmarkov from 2 to 10
    n_markov = randi([2,10]);
    % from 0.2 to 10
    sigma = randi([2,100]) / 10;
end

while true
    Struct = generate_parameters(n_markov, sigma);
    if Struct.opT_is_stable
        break
    end
end
end
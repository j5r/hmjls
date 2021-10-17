%
% By Junior R. Ribeiro, Oct 17, 2021, jrodrib@usp.br
%
% This function recomputes the parameters that depend on
%      Struct.sigma
%      Struct.Ac
%      Struct.Ad
%      Struct.RateMatrix
%      Struct.ProbMatrix
%      Struct.init_distrib
%      Struct.mu
%
% Run this function whenever you change some of these parameters. If you
% want to do a validation raising errors, set "needs_validation" to true,
% otherwise, you can skip this parameter.
%
% Struct = recompute_parameters(Struct)
% Struct = recompute_parameters(Struct, needs_validation)
% 

function Struct = recompute_parameters(Struct, needs_validation)
if nargin == 1 || ~needs_validation
    needs_validation = false;
end

if needs_validation
    Struct = validate_mmjls(Struct.Ac, Struct.Ad, Struct.RateMatrix,...
        Struct.ProbMatrix, Struct.init_distrib, Struct.mu, Struct.sigma);
else
    Struct = parse_mmjls(Struct.Ac, Struct.Ad, Struct.RateMatrix,...
        Struct.ProbMatrix, Struct.init_distrib, Struct.mu, Struct.sigma);
end
end
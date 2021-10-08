%
% By Junior R. Ribeiro, Oct 8, 2021, jrodrib@usp.br
%
% This function appends information about the domain to the struct given.
% The new field is "time_domain" with the subfields "max_time" and
% "sampling_time".
% It also discretizes the continuous-time parameters: 
%      Ac -> Ac_d
%    opAc -> opAc_d
% according the sampling time given. If you want to change previous
% definition, just run this function again, with the new infos. All those
% data will be replaced.
%
% Struct = define_domain(Struct, t_max, sampling_time)
%
% The input parameters are 
% (A) Struct: from "parse_mmjls" or "validate_mmjls",
% (B) t_max: (1x1) the time for which you want to simulate,
% (C) sampling_time: (1x1) the sampling time to discretize the domain.
%       Usually it is a small positive number.
%
% The Struct returned is the same Struct given as input, but added by the
% fields:
%     struct.time_domain.max_time = t_max
%     struct.time_domain.sampling_time = sampling_time
%     struct.opAc_d = opAc discretized
%     struct.Ac_d = Ac matrices discretized
%

function Struct = define_domain(Struct, t_max, sampling_time)
assert(t_max > 0, '"t_max" must be > 0.');
assert(sampling_time > 0, '"sampling_time" must be > 0.')
assert(t_max > sampling_time, '"t_max" must be > "sampling_time".');

domain.max_time = t_max;
domain.sampling_time = sampling_time;
Struct.time_domain = domain;

% discretizing continuous-time operator: opAc_d
Struct.opAc_d =  expm(sampling_time * Struct.opAc);
% discretizing continuous-time matrices:   Ac_d
for theta=1:Struct.nmarkov
    Struct.Ac_d(:,:,theta) = expm(sampling_time * Struct.Ac(:,:,theta));
end
end
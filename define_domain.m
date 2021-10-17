%
% By Junior R. Ribeiro, Oct 8, 2021, jrodrib@usp.br
%
% This function appends information about the domain to the struct given.
% The new field is "time_domain" with the subfields "max_time", 
% "c_sampling_time", and "d_sampling_time".
% It also discretizes the continuous-time parameters: 
%      Ac -> Ac_d
%    opAc -> opAc_d
% according the sampling time given. If you want to change previous
% definition, just run this function again, with the new infos. All those
% data will be replaced.
%
% Struct = define_domain(Struct, t_max, c_sampling_time, d_sampling_time)
%
% The input parameters are 
% (A) Struct: from "parse_mmjls" or "validate_mmjls",
% (B) t_max: (1x1) the time for which you want to simulate,
% (C) c_sampling_time: (1x1) the sampling time to discretize the 
%       continuous-time domain. Usually it is a small positive number.
% (D) d_sampling_time(optional): (1x1) the sampling time to space 
%       the discrete-time domain. Must be greater or equal to c_sampling_time.
%       If it is not given, it will be equal to c_sampling_time.
%
% The (C) and (D) parameters regards to the spacing between the points,
% that can be different. One restriction is required, that (D) must be not
% less than (C). Besides that, (D) is compelled to be a multiple of (C).
%
% The Struct returned is the same Struct given as input, but added by the
% fields:
%     struct.time_domain.max_time = t_max
%     struct.time_domain.c_sampling_time = c_sampling_time
%     struct.time_domain.d_sampling_time = d_sampling_time
%     struct.opAc_d = opAc discretized
%     struct.Ac_d = Ac matrices discretized
%

function Struct = define_domain(Struct, t_max, c_sampling_time, d_sampling_time)
% validating input
assert(t_max > 0, '"t_max" must be > 0.');
assert(c_sampling_time > 0, '"sampling_time" must be > 0.')
assert(t_max > c_sampling_time, '"t_max" must be > "sampling_time".');
if nargin < 4
    d_sampling_time = c_sampling_time;
end
assert(d_sampling_time >= c_sampling_time, ...
    '"d_sampling_time" must be >= "c_sampling_time".');

%
domain.max_time = t_max;
domain.c_sampling_time = c_sampling_time;

% enforcing d-spacing to be multiple of c-spacing.
domain.d_sampling_time = floor(d_sampling_time/c_sampling_time)*c_sampling_time;
Struct.time_domain = domain;

% discretizing continuous-time operator: opAc_d
Struct.opAc_d =  expm(c_sampling_time * Struct.opAc);
% discretizing continuous-time matrices:   Ac_d
for theta=1:Struct.nmarkov
    Struct.Ac_d(:,:,theta) = expm(c_sampling_time * Struct.Ac(:,:,theta));
end
end
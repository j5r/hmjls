%
% By Junior R. Ribeiro, Oct 8, 2021, jrodrib@usp.br
%
% This function simulates the 2nd moment X of the MMJLS for a set of
% arrival times t_r. This process is referred as Y.
%
% struct2 = simulate_Yr_once(Struct, x0, n_switch_times)
%
% The input parameters are
% (A) Struct: from "parse_mmjls" or "validate_mmjls",
% (B) x0: the initial state (you can let it empty x0=[] to be generated
%        automatically),
% (C) n_switch_times: the number of t_r for which Y shall be computed, i.e.,
%        Y(:,:,nmarkov, 1:n_switch_times+1)
%
% The struct2 returned has the fields explained below.
%   struct2.values: (nnN x n_switch_times+1) the values of Y themselves, as
%        vectorized matrices
%   struct2.shape_vec: (1x2) the shape of the vectorized Y,
%            shape_vec = [nnN, n_switch_times+1]
%   struct2.shape_full: (1x4) the shape of the 4D array Y(:,:,markov,time)
%            shape_full = [n, n, N, n_switch_times+1]
%        You can get Y back as a 4D array by doing
%            reshape(struct2.values, struct2.shape_full)
%   struct2.trace_indexes: (nxN) the indexes of the trace of Y, for time=1.
%        You can get the values of the trace of Y by doing
%        |  for time = 1:n_switch_times+1
%        |      Y_ = Y_full(:,:,:,time);
%        |      for i = 1:N % N is the number of Markov states
%        |          Y_trace(i,time) = sum( Y_(struct2.trace_indexes(:,i)) );
%        |      end
%        |  end
%   struct2.x0: (nx1) the x(0) given or the x(0) generated.
%

function ANS = simulate_Yr_once(Struct, x0, n_switch_times)
if nargin == 0
    ANS = [];
    return
end
n = size(Struct.Ac,1); % the length of the vector x0
if isempty(x0)
    x0 = randn(n, 1);
end

% Y pre-allocation
shape_vec = [n * n * Struct.nmarkov, n_switch_times + 1];
shape_full = [n, n, Struct.nmarkov, n_switch_times + 1];

%
Y_vec = zeros(shape_vec);

% computing the initial value of Y, i.e., Y(0)
Y0 = zeros(n, n, Struct.nmarkov);
for i = 1:Struct.nmarkov
    Y0(:,:,i) = Struct.init_distrib(i) * (x0*x0'); % X_i(0) = x0*x0' * [pi0(i)]
end
% assigning the initial value of Y to Y_vec
Y_vec(:,1) = Y0(:);

% getting the indexes of the trace of Y
Y_ = reshape( 1:n*n*Struct.nmarkov, n, n, Struct.nmarkov);
trace_indexes = zeros(n, Struct.nmarkov);
for i = 1:Struct.nmarkov
    trace_indexes(:,i) = diag(Y_(:,:,i));
end
%

for k = 2:n_switch_times
    % applying ctime operators
    ctime = exprnd(1/Struct.sigma);
    Y_vec(:,k) = expm(ctime * Struct.opAc) * Y_vec(:,k-1);
    % applying dtime operators
    dtime = find(rand < cumsum(Struct.mu), 1) - 1;
    Y_vec(:,k) = (Struct.opAd^dtime) * Y_vec(:,k);
end

ANS.values = Y_vec;
ANS.shape_vec = shape_vec;
ANS.shape_full = shape_full;
ANS.trace_indexes = trace_indexes;
ANS.x0 = x0;
end
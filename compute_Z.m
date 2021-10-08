%
% By Junior R. Ribeiro, Oct 8, 2021, jrodrib@usp.br
%
% This function computes the 2nd moment of X irrespectively to the arrival
% times of the continuous domain, variable called Z. This variable performs
% a discrete process using the operator opT.
%
% struct2 = compute_Z(Struct, x0, n_steps)
% 
% The parameters given are 
% (A) Struct: from "parse_mmjls" or "validate_mmjls",
% (B) x0: the initial state (you can let it empty x0=[] to be generated
%        automatically), 
% (C) "n_steps" is the number of the branchs continuous+discrete domains, 
%        i.e., for computing Z(:,:,r) for r = 1:n_steps + 1
%        (plus one because the initial state is done in r = 1, not r = 0).
%
% The struct2 returned has the fields explained below.
%   struct2.values: (nnN x n_steps+1) the values of Z themselves, as vectorized 
%        matrices
%   struct2.shape_vec: (1x2) the shape of the vectorized Z, 
%            shape_vec = [nnN, n_steps+1]
%   struct2.shape_full: (1x4) the shape of the 4D array Z(:,:,markov,time)
%            shape_full = [n, n, N, n_steps+1]
%        You can get Z back as a 4D array by doing 
%            reshape(struct2.values, struct2.shape_full)
%   struct2.trace_indexes: (nxN) the indexes of the trace of Z, for time=1.
%        You can get the values of the trace of Z by doing
%        |  for time = 1:n_steps+1
%        |      Z_ = Z(:,:,:,time);
%        |      for i = 1:N % N is the number of Markov states
%        |          Z_trace(i,time) = sum( Z_(struct2.trace_indexes(:,i)) );
%        |      end
%        |  end
%   struct2.x0: (nx1) the x(0) given or the x(0) generated.
%

function ANS = compute_Z(Struct, x0, n_steps)
n = size(Struct.Ac, 1);

if isempty(x0)
    disp('As you do not gave x0, it will be generated randomly.');
    x0 = randn(n,1);
end

shape_full = [n, n, Struct.nmarkov, n_steps + 1];
shape_vec  = [n * n * Struct.nmarkov, n_steps + 1];

Z0 = zeros(n, n, Struct.nmarkov);
for i = 1:Struct.nmarkov
    % X_i(0) = x0*x0' * [init_distrib(i)]
    Z0(:,:,i) = Struct.init_distrib(i) * (x0*x0');
end

% initializing Z
Z = zeros(shape_vec);
Z(:,1) = Z0(:);

% getting the indexes of the trace of Y
Z_ = reshape(1:n*n*Struct.nmarkov,n,n,Struct.nmarkov);
trace_indexes = zeros(n,Struct.nmarkov);
for i = 1:Struct.nmarkov
    trace_indexes(:,i) = diag(Z_(:,:,i));
end

% computing Z
for k = 2:n_steps + 1
    Z(:,k) = Struct.opT * Z(:, k-1);
end
ANS.values = Z;
ANS.shape_vec = shape_vec;
ANS.shape_full = shape_full;
ANS.trace_indexes = trace_indexes;
ANS.x0 = x0;
end
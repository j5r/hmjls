%
% By Junior R. Ribeiro, Oct 7, 2021, jrodrib@usp.br
%
% This function validates the parameters given, raising errors if they do
% not make sense. After all validations, it uses the "parse_mmjls" to get
% a struct with all the data given and others computed.
%
% For more details, see the documentation of <a href="matlab:web('parse_mmjls.m')">parse_mmjls</a>.

function Struct = validate_mmjls(Ac, Ad, RateMatrix, ProbMatrix, init_distrib, mu, sigma)
if nargin == 0
    disp('@ function Struct = validate(Ac, Ad, RateMatrix, ProbMatrix, init_distrib, mu, sigma)');
    fprintf('\tAc: matrices of the continuous-time process\n')
    fprintf('\tAd: matrices of the discrete-time process\n')
    fprintf('\tRateMatrix: rate matrix for the continuous-time process\n')
    fprintf('\tProbMatrix: probability matrix for the discrete-time process\n')
    fprintf('\tinit_distrib: initial distribution for Markov states\n')
    fprintf('\tmu: distribution for the length of stay in discrete domain\n')
    fprintf('\tsigma: parameter of the exponential distribuion for the length')
    fprintf('\n\t\t\tof stay in continuous domain\n');
    disp('@ prompt "fieldnames(Struct)" to see the fields of the struct obtained.');
    Struct = [];
    return
end

% assert evaluates only scalar boolean values, not arrays. So I will use all()
% For numeric error issues, I will consider eps(100) as zero.
assert (sigma > 0,'"sigma" must be > 0.')
sizeA = size(Ac);
nmarkov = sizeA(3);

%% ASSERTS
% [Ac] is square
assert ( sizeA(1) == sizeA(2), ...
    '"Ac" must be square like (n x n x N) - look at n x n.');

% [Ac] and [Ad] have the same shape
assert ( numel(sizeA) == numel(size(Ad)) && all(sizeA == size(Ad)),...
    '"Ac" and "Ad" must have same shape.');

% [Ac] and [Ad] are three dimensional array
assert ( numel(sizeA) == 3 ,'"Ac" must be 3-dimensional array.');

% [RateMatrix] is square with nmarkov rows/cols
assert ( all(size(RateMatrix) == [nmarkov,nmarkov]),...
    ['RateMatrix be square with ', num2str(nmarkov), ' rows/cols.']);

% [RateMatrix] and [ProbMatrix] have same shape
assert ( all(size(RateMatrix) == size(ProbMatrix)), ...
    '"ProbMatrix" must have the same shape as "RateMatrix".');

% [RateMatrix] rows add up to zero
assert ( all(abs(sum(RateMatrix,2)) < eps(100) ), ...
    '"RateMatrix" rows must add up to zero.');

% [ProbMatrix] rows add up to one
assert ( all(abs(sum(ProbMatrix,2) - 1) < eps(100)) ,...
    '"ProbMatrix" rows must add up to one.');

% [ProbMatrix] is >= 0
assert ( all(all(ProbMatrix>=0)), ...
    '"ProbMatrix" elements must be >= 0.');

% Distribution [init_distrib] has nmarkov dimensions
assert ( nmarkov == numel(init_distrib), ...
    ['"init_distrib" must have ', num2str(nmarkov), ' elements.']);

% Distribution [init_distrib] add up to one
assert ( abs(sum(init_distrib) - 1) < eps(100),...
    '"init_distrib" must add up to one.');

% Distribution [init_distrib] is >=0
assert ( all(init_distrib >=0 ) ,...
    '"init_distrib" elements must be >= 0.');

% Distribution [mu] adds up to one
assert ( abs(sum(mu(:)) - 1) < eps(100), '"mu" must add up to one.');

% Distribution [mu] is >=0
assert ( all(mu(:) >=0 ), '"mu" elements must be >= 0.');

Struct = parse_mmjls(Ac, Ad, RateMatrix, ProbMatrix, init_distrib, mu, sigma);
end

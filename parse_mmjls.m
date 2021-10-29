%
% By Junior R. Ribeiro, Oct 7, 2021, jrodrib@usp.br
%
% This function parses the parameters given into a structure.
%
%  struct = parse_mmjls(Ac, Ad, RateMatrix, ProbMatrix, initial_distrib, mu, sigma)
% 
% !! If you are sure that the parameters make sense, just go ahead and use "parse_mmjls".
% !! If you want to ensure that the parameters make sense, use "validate_mmjls"
% instead. It will raise the correspondent errors.
%
% "parse_mmjls" returns a structure with all the parameters given and some other 
% parameters computed, as follows.
%
% GIVEN PARAMETERS
%   struct.sigma := (1x1) exponential rate for time in continuous domain.
%   struct.Ac := (nxnxN) parameter of the continuous process.
%   struct.Ad := (nxnxN) parameter of the discrete process.
%   struct.nmarkov := (1x1), N - the number of Markov states. This
%        parameter is not given, but it is inferred by the 3rd dimension of Ac.
%   struct.RateMatrix := (NxN) the rate transition matrix of the
%        continuous Markov chain.
%   struct.ProbMatrix := (NxN) the probability transition matrix of the
%        discrete Markov chain.
%   struct.init_distrib := (Nx1) initial distribution of the Markov states.
%   struct.mu := (rx1) distribution for time in discrete domain.
%
% COMPUTED PARAMETERS
%   struct.opAd := (nnNxnnN) operator for 2nd moment X vectorized of the 
%        discrete process (used to evaluate stability in discrete-only
%        process, when all its eigenvalues are inside the unit circle).
%        || See operator A_1 at eq. 3.12d of this book
%        || https://dx.doi.org//10.1007/b138575
%   struct.opAc := (nnNxnnN) operator for 2nd moment X vectorized of the 
%        continuous process (used to evaluate stability in continuous-only
%        process, when all its eigenvalues are left-side the complex axis).
%        || See operator A at eq. 3.25 of this book
%        || https://dx.doi.org/10.1007/978-3-642-34100-7
%   struct.opL := (nnNxnnN) integral operator playing the role of the action 
%        of the <opAc> as its expected value.
%   struct.opBoldA := (nnNxnnN) the expected value of the action of the <opAd>.
%   struct.opT_is_well_defined := (1x1) boolean:  max(real(eig(opAc))) < sigma.
%   struct.opT := (nnNxnnN) operator for 2nd moment X vectorized of the 
%        continuous-discrete process (the discrete 2nd moment process).
%   struct.radiuses := (struct) computes the main information about eigenvalues
%        of all the operators: the greatest eigenvalue module in case of 
%        discrete-time operator, or the greatest real part of the eigenvalues 
%        in case of continuous-time operator.
%   struct.opT_is_stable := (1x1) boolean: max(abs(eig(opT))) < 1.
%

function Struct = parse_mmjls(Ac, Ad, RateMatrix, ProbMatrix, init_distrib, mu, sigma)
if nargin ~= 7
    assert(false, 'u:stuffed:it', ['Incorrect syntax. All the parameters ',...
        'are mandatory. Please, use\n>>\t',...
        'parse_mmjls(Ac, Ad, RateMatrix, ProbMatrix, init_distrib, mu, sigma)']);
end
sizeA = size(Ac);
nmarkov = sizeA(3);

%% Computing operators
diagD = []; % diagonal block matrix related to discrete-time domain
diagC = []; % diagonal block matrix related to continuous-time domain
for i = 1:nmarkov
    diagD = blkdiag( diagD, ...
        kron(Ad(:,:,i), Ad(:,:,i))  );
    
    diagC = blkdiag( diagC, ...
        kron(Ac(:,:,i), eye(sizeA(1))   ) + ...
        kron(eye(sizeA(1)), Ac(:,:,i))   );
end

% calligraph operators [Ac] and [Ad]
opAc  = kron(RateMatrix, eye(sizeA(1)^2) ) + diagC;
opAd = kron(ProbMatrix', eye(sizeA(1)^2) ) * diagD;

% boldface operator [A]: expected action of the calAd
boldA = opAd *0;
powerOpAd = eye(size(opAd));
for i = 1:numel(mu)    
    boldA = boldA + powerOpAd * mu(i);
    powerOpAd = powerOpAd * opAd;
end

%% condition for definiteness of calligraph L (and hence T)
max_re_eig_opAc = max(real(eig(opAc)));
if max_re_eig_opAc < sigma
    opT_is_well_defined = true;
    fprintf(['\t@ Max(Real(Eigen(opAc))) < sigma : TRUE;',...
        '\n\t@ opL and hence opT IS well defined.',...
        '\n\t@ see: Struct.radiuses.ctime.opAc\n']);
    opL = - (opAc - sigma*eye(size(opAc)) ) \ eye(size(opAc));
    opT = sigma*boldA*opL;
else
    opT_is_well_defined = false;
    fprintf(['\t@ Max(Real(Eigen(opAc))) < sigma : FALSE;',...
        '\n\t@ opL and hence opT IS NOT well defined.',...
        '\n\t@ see: Struct.radiuses.ctime.opAc\n']);
    opL = []; % empty if 
    opT = []; % empty
end

%% Computing convergence radiuses for all parameters
% for the matrices
ctime.Ac = zeros(1,nmarkov); % regarding continuous-time
dtime.Ad = zeros(1,nmarkov); % regarding discrete-time
for i=1:nmarkov
    dtime.Ad(i) =  max(abs(eig(Ad(:,:,i))));  % unit circle
    ctime.Ac(i) =  max(real(eig(Ac(:,:,i)))); % left-side complex plane
end
% for the operators
ctime.opAc = max(real(eig(opAc)));      % left-side complex plane
ctime.opL = max(real(eig(opL)));        % left-side complex plane
dtime.opAd = max(abs(eig(opAd)));       % unit circle
dtime.opBoldA = max(abs(eig(boldA)));   % unit circle
dtime.opT = max(abs(eig(opT)));         % unit circle
% some help
ctime.help = sprintf('max(real(eig(opAc))) < sigma means CTIME-stability (sigma = %.2f)',sigma);
dtime.help = {'max(abs(eig(opAd))) < 1 means DTIME-stability';
    'max(abs(eig(opT))) < 1 means MMJLS-stability'};
radius.ctime = ctime;
radius.dtime = dtime;

%% saving data into struct
Struct.sigma = sigma;
Struct.Ac = Ac;
Struct.Ad = Ad;
Struct.nmarkov = nmarkov; % number of states of the Markov chain
Struct.RateMatrix = RateMatrix;
Struct.ProbMatrix = ProbMatrix;
Struct.init_distrib = init_distrib;
Struct.mu = mu;
Struct.opAd = opAd;
Struct.opAc = opAc;
Struct.opL = opL;
Struct.opBoldA = boldA;
Struct.opT_is_well_defined = opT_is_well_defined;
Struct.opT = opT;
Struct.radiuses = radius;
Struct.opT_is_stable = Struct.radiuses.dtime.opT < 1;

end

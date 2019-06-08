%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB code used to compute the moment bound in Section V.C of 
% https://arxiv.org/abs/1702.05468
% Prerequisites:
% - YALMIP (c.f. https://yalmip.github.io/)
% - mpYALMIP (c.f. https://github.com/aeroimperial-optimization/mpYALMIP‎)
% - SDPA-GMP (c.f. http://sdpa.sourceforge.net/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Declare model parameters

k1 = 50; k2 = 1/2; % Rate constants

%% Declare optimisation parameters

d = 7; % Approximation order
db = 2; % Maximum degree of the propensities
objective = 7; % Order of moment to be bound 
opts = sdpsettings('solver','sdpa_gmp','sdpa_gmp.precision',300); % Solver 
                                                                  % options

%% Construct SDPs

% Declare optimisation variables

y = sdpvar(d+1,1); % In 1-D, #d = d + 1

% Declare mass constraint

Constraints = y(1) == 1;

% Declare stationary moment equations constraints

for n = 0:d-db+1
    temp = 0;
    for m = 0:n-1
        temp = temp + nchoosek(n,m)*2^(n-m)*(k1*y(m+1)+k2*(-1)^(n-m)*...
        (y(m+3)-y(m+2)));
    end
    Constraints = [Constraints; temp == 0];
end

% Declare moment matrix constraint

for n = 0:floor(d/2)
    for m = 0:floor(d/2)
        M0(n+1,m+1) = y(n+m+1);
    end
end

Constraints = [Constraints; M0 >= 0];

% Declare localising matrix constraint

for n = 0:floor((d-1)/2)
    for m = 0:floor((d-1)/2)
        M1(n+1,m+1) = y(n+m+2);
    end
end

Constraints = [Constraints; M1 >= 0];

%% Solve SDPs

% Lower bound 

diag{1} = optimize(Constraints,y(objective+1),opts); % Solve SDP
lowerbound = double(y(objective+1)); % Extract optimal value (i.e., bound)

% Upper bound

diag{2} = optimize(Constraints,-y(objective+1),opts); % Solve SDP
upperbound = double(y(objective+1)); % Extract optimal value (i.e., bound)



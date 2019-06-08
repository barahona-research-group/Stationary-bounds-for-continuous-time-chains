%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB code used to generate Figure 8 in Section V.C of 
% https://arxiv.org/abs/1702.05468
% Prerequisites:
% - YALMIP (c.f. https://yalmip.github.io/)
% - CPLEX (c.f. https://www.ibm.com/uk-en/products/ilog-cplex-optimization-
% studio‎)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Declare model parameters

k1 = 50; k2 = 1/2; % Rate constant
d = 7; % Order of moment
c = 4.8814e07; % Moment bound (<x^d> <= c) computed using the SDP approach


%% Declare optimisation parameters

ts = 25; % Truncation size
r = ts^d; % Truncation parameter
opts = sdpsettings('solver','cplex','cplex.lpmethod',2,... % Solver options
'cplex.emphasis.numerical',1,'cplex.preprocessing.presolve',0); 


%% Construct LPs

% Declare decision variables

pi = sdpvar(1,ts);

% Declare stationary equations constraints

Q = sparse(ts,ts-2); % Initialise truncated rate matrix

Q(1,1) = -k1; Q(3,1) = 2*k2; % Populate first row
Q(2,2) = -k1; Q(4,2) = 6*k2; % Populate second row

for x = 2:ts-3 % Populate remaining rows
    Q(x-1,x+1) = k1;
    Q(x+1,x+1) = -k1-k2*x*(x-1);
    Q(x+3,x+1) = k2*(x+2)*(x+1);
end

Constraints = pi*Q == 0; % Declare constraints

% Declare non-negativity constraints

Constraints = [Constraints; pi >= 0];

% Declare mass constraints

Constraints = [Constraints; 1 - c/r <= sum(pi); sum(pi) <= 1]; 

% Declare moment constraint

Constraints = [Constraints; pi*(0:ts-1)'.^d <= c];



%% Compute bounds on the set of stationary distributions

for x = 0:ts-1
    % Upper bounds
    diagu{x+1} = optimize(Constraints,-pi(x+1),opts);
    piupper(x+1) = double(pi(x+1));
    % Lower bounds
    diagl{x+1} = optimize(Constraints,pi(x+1),opts);
    pilower(x+1) = double(pi(x+1));
end

%% Compute approximation for the ergodic distribution with support on the 
%% even numbers

diag0 = optimize(Constraints,-pi(1),opts);
pieven = double(pi);


%% Compute approximation for the ergodic distribution with support on the 
%% odd numbers

diag1 = optimize(Constraints,-pi(2),opts);
piodd = double(pi);


%% Compute bounds on the ergodic distribtuion with support on the even 
%% numbers

% Restrict the support

Constraintseven = [Constraints; pi(2:2:2*floor((ts-2)/2)+2) == 0]; 

% Compute bounds

for x = 0:2:2*floor((ts-1)/2)
    % Upper bounds
    diague{x/2+1} = optimize(Constraintseven,-pi(x+1),opts);
    piuppereven(x+1) = double(pi(x+1));
    % Lower bounds
    diagle{x/2+1} = optimize(Constraintseven,pi(x+1),opts);
    pilowereven(x+1) = double(pi(x+1));
end

% Populate odd entries with NaNs

for x = 1:2:2*floor((ts-2)/2)+1
    piuppereven(x+1) = NaN;
    pilowereven(x+1) = NaN;
end


%% Compute bounds on the ergodic distribtuion with support on the odd 
%% numbers

% Restrict the support

Constraintsodd = [Constraints; pi(1:2:2*floor((ts-1)/2)+1) == 0]; 

% Compute bounds

for x = 1:2:2*floor((ts-2)/2)+1
    % Upper bounds
    diaguo{(x+1)/2} = optimize(Constraintsodd,-pi(x+1),opts);
    piupperodd(x+1) = double(pi(x+1));
    % Lower bounds
    diaglo{(x+1)/2} = optimize(Constraintsodd,pi(x+1),opts);
    pilowerodd(x+1) = double(pi(x+1));
end

% Populate even entries with NaNs

for x = 0:2:2*floor((ts-1)/2)
    piupperodd(x+1) = NaN;
    pilowerodd(x+1) = NaN;
end

%% Generate figures

% First figure: Bounds on the set of stationary distributions

figure(1)
hold on
plot(0:ts-1,piupper,'.k','MarkerSize',25)
plot(0:ts-1,pilower,'ok','MarkerSize',6.5)
hold off
ylim([0,0.27])
xlabel('Molecule count')
ylabel('Probability')
title('Bounds on the entire set of stationary distributions (dots upper, circles lower)')

% Second figure: Approximations of the ergodic distributions

figure(2)
hold on
plot(0:ts-1,pieven,'.b','MarkerSize',25)
plot(0:ts-1,piodd,'.r','MarkerSize',25)
hold off
ylim([0,0.27])
xlabel('Molecule count')
ylabel('Probability')
title('Approximations of ergodic distributions (blue even, red odd)')

% Third figure: Bounds on the ergodic distributions (dots upper, circles
% lower, blue even, red odd)

figure(3)
hold on
plot(0:ts-1,piuppereven,'.b','MarkerSize',25)
plot(0:ts-1,pilowereven,'ob','MarkerSize',6.5)
plot(0:ts-1,piupperodd,'.r','MarkerSize',25)
plot(0:ts-1,pilowerodd,'or','MarkerSize',6.5)
hold off
ylim([0,0.27])
xlabel('Molecule count')
ylabel('Probability')
title('Bounds on the ergodic distributions (dots upper, circles lower, blue even, red odd)')

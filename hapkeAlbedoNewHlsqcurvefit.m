% author: Jennifer Harris

% a routine for calculating the SSA given the reflectance spectra using the Hapke model.
% this will assume non-isotropic scattering taken care of in called function (hapke_reflectanceLegrendeP), therefore only solving for w here

% based on Sklute's MasterPhase1_PP.m from paper https://doi.org/10.2138/am-2015-4824


N = length(X);
M = length(Rc); % the length of the input reflectance spectrum

% define angle parameters
inc = 30;
emi = 0;
g = 30;

mu = cosd(emi);
mu0 = cosd(inc);
mug = cosd(g);

%DEFINE MINIMIZATION PARAMETERS
%Maximum function evaluations
maxfun=1000000000000000; %search MaxFunEvals for more info
%Function Tolerance (smaller is closer match)
funtol=0.000000000000000000000001; %search FunTol for more info
xtol=0.0000000000000000000000001; %serach XTol for more info
%Maximum iterations
maxit=2000;
%number of start points
spts=2;

% bounds for the calculated SSA
lb = 0.0;

% bounds for the calculated SSA
ub = 1.0;

% need to provide an initial guess as to the value of w, this is a single
% value and it doesn't affect the final result so only needs to be vaguely
% correct

w = wguess;

% Will set up a loop to calculate w for each value of Rc

% first set up an empty array to hold the results, the same length as the input reflectance spectrum
W_m = zeros(M,1);

% set up index for saving the results into the results array
index = 0;

for i = 1:1:M
    %  Create an anonymous function in order to pass g, inc, and emi as 
    %  extra parameters to the objective function.

    myObjFcn = @(w,X)hapke_reflectanceLegrendeP(w,X); % this function requires estimates of b and c to account for the non-isotropy, alt function is hapke_reflectanceSimple
    % USE LSQCURVEFIT TO PERFORM LEAST SQUARES MINIMIZATION between data 
    % and model. Format is [output]=lsqcurvefit(function,x0,xdata,ydata,lb,ub)
    options=optimoptions(@lsqcurvefit,'Algorithm','trust-region-reflective',...
       'Display','iter','MaxIter',maxit,'FinDiffType','central',...
       'MaxFunEvals',maxfun, 'TolFun',funtol,'TolX', xtol);
    problem = createOptimProblem('lsqcurvefit','x0',w,'objective',myObjFcn,...
       'xdata',X(i),'ydata',Rc(i),'lb',lb,'ub',ub,'options',options);
    ms = MultiStart('StartPointsToRun','bounds'); 
    [xmulti,errormulti]=run(ms,problem,spts);
    %  xmulti = hapke_reflectance(coefg,X);
    index = index + 1;
    W_m(index,:) = xmulti;
end
% SSAnonIsoDLCMP091AC1DL91ACPX = W_m;
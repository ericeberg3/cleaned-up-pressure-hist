function optParams = optimize_SC_bayes(mSCguess, lb, ub, xopt, xtilt, yopt, ytilt, zopt, u1d, ...
    invStdPWRL, tiltstd, tiltreduced, nanstatbeginning)

f = @(m)green_residuals(m, [xopt, xtilt], [yopt, ytilt], [zopt, 0], u1d, ...
    [invStdPWRL(1), invStdPWRL(2), invStdPWRL(3),1/((tiltstd)), 1/((tiltstd))], tiltreduced(1:2), nanstatbeginning);

params = [
    optimizableVariable('dpHMM', [lb(1), ub(1)], 'Type', 'real');
    optimizableVariable('dip', [lb(2), ub(2)], 'Type', 'real');
    optimizableVariable('strike', [lb(3), ub(3)], 'Type', 'real');
    optimizableVariable('x1', [lb(4), ub(4)], 'Type', 'real');
    optimizableVariable('x2', [lb(5), ub(5)], 'Type', 'real');
    optimizableVariable('x3', [lb(6), ub(6)], 'Type', 'real');
    optimizableVariable('dpSC', [lb(7), ub(7)], 'Type', 'real');
];

% Perform the optimization
results = bayesopt(@(m) f([m.dpHMM, m.dip, m.strike, m.x1, m.x2, m.x3, m.dpSC]), ...
                   params,'AcquisitionFunctionName', 'expected-improvement-plus', ...
                   'MaxObjectiveEvaluations', 100);

optParams = results.XAtMinObjective;
disp('Optimal Parameters:');
disp(optParams);
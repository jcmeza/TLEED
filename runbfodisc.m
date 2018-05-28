%% Run BFO with combined discrete/continuous variables

for run=1:10
    disp(run);
    
%%  Set initial iterate value
    init_x0 = tleed_x0;
    x0 = init_x0(run).x; 
    p0  = cell2mat(init_x0(run).p);
    xp0 = [ p0 x0']; % optimize with respect to both p and (z,x,y)
    disp(xp0'); % display initial guess

%% First 14 variables are discrete; all others continuous
    xdisc  = repmat('i',14,1)';
    xcont = repmat('c',42,1)';
    xtype = strcat(xdisc,xcont);

%% Other options
    plower = ones(14,1)'; % Lower/Upper bound for discrete variables
    pupper = 2*ones(14,1)';
    delta  = 1.0;         % possibly use smaller initial box for constraints
    xlower = x0 - delta*max(1.0,abs(x0));
    xupper = x0 + delta*max(1.0,abs(x0));
    xplower = [plower xlower']; % combine all bounds
    xpupper = [pupper xupper'];
    fbound = 1.0; % Not used for now

%% Call BFO
    [ x, fx, msg, wrn, neval, f_hist] = bfo(@(x)tleedfcn2(x), xp0, ...
            'epsilon', 0.0001, 'maxeval', 3000, ... 
            'xlower', xplower, 'xupper', xpupper, ...
            'save-freq', 10, 'restart-file', 'bfodisc.restart', ...
            'verbosity', 'low'); %, ...
            %'f-call-type', 'with-bound', ...
            %'f-bound', 1.0)
 
%% Save output
    xstar_bfodisc(run).x = x;
    xstar_bfodisc(run).fx = fx;
    xstar_bfodisc(run).neval = neval;
end
 
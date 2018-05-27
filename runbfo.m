% set numerical run
for run=1:10
    disp(run);
    
%%  Set initial iterate value
    iterate = tleed_nomadm_x0;
    init_x0 = tleed_x0;
    x0 = init_x0(run).x;
    p  = cell2mat(init_x0(run).p);

%% Set optimization parameters

    delta  = 1.0; % possibly use smaller initial box for constraints
    xlower = x0 - delta*max(1.0,abs(x0));
    xupper = x0 + delta*max(1.0,abs(x0));
    fbound = 1.0;
%% Call BFO

    [x, fx, msg, wrn, neval, f_hist] = bfo(@(x)tleedfcn(x,p), x0, ...
            'epsilon', 0.0001, 'maxeval', 2000, ... 
            'xlower', xlower, 'xupper', xupper, ...
            'save-freq', 10, 'restart-file', 'bfo.restart');
    %    'f-call-type', 'with-bound', ...
    %    'f-bound', fbound);
%% Store results
    xstar_bfo(run).x  = x;
    xstar_bfo(run).fx = fx;
    xstar_bfo(run).neval = neval;
end
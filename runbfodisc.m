%
% set numerical run
run = 1;
%
% set initial iterate value
init_x0 = tleed_x0;
x0 = init_x0(run).x;
p0  = cell2mat(init_x0(run).p);
xp = [ p0 x0']; % optimize with respect to both p and (z,x,y)
xp0 = xp;
disp(xp0); % display initial guess

% first 14 variables are discrete; all others continuous
xdisc  = repmat('i',14,1)';
xcont = repmat('c',42,1)';
xtype = strcat(xdisc,xcont);
%
% other options
plower = ones(14,1)';
pupper = 2*ones(14,1)';
delta  = 1.0; % possibly use smaller initial box for constraints
xlower = x0 - delta*max(1.0,abs(x0));
xupper = x0 + delta*max(1.0,abs(x0));
xplower = [plower xlower'];
xpupper = [pupper xupper'];

fbound = 1.0; % Not used for now

[ x, fx, msg, wrn, neval, f_hist] = bfo(@(x)tleedfcn2(x), xp0, ...
        'epsilon', 0.0001, 'maxeval', 3000, ... 
        'xlower', xplower, 'xupper', xpupper, ...
        'save-freq', 10, 'restart-file', 'bfodisc.restart', ...
        'verbosity', 'low'); %, ...
        %'f-call-type', 'with-bound', ...
        %'f-bound', 1.0)
 xstar(run).x = x;
 xstar(run).fx = fx;
 
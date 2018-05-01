%
% set initial iterate value
iterate = tleed_nomadm_x0;
x0 = iterate(1).x;
p  = cell2mat(iterate(1).p);
xp = [ p x0']; % optimize with respect to both p and (z,x,y)
xp0 = xp;

% first 14 variables are discrete; all others continuous
xdisc  = repmat('i',14,1)';
xcont = repmat('c',42,1)';
xtype = strcat(xdisc,xcont);
%
% other options
plower = ones(14,1)';
pupper = 2*ones(14,1)';
xlower = x0 - 0.5*max(1.0,abs(x0));
xupper = x0 + 0.5*max(1.0,abs(x0));
xplower = [plower xlower'];
xpupper = [pupper xupper'];

fbound = 1.0;

[ x, fx, msg, wrn, neval, f_hist] = bfo(@(x)tleedfcn2(x), xp0, ...
        'epsilon', 0.0001, 'maxeval', 1000, ... 
        'xlower', xplower, 'xupper', xpupper, ...
        'save-freq', 10, 'restart-file', 'bfodisc.restart', ...
        'verbosity', 'low');
    %    'f-call-type', 'with-bound', ...
    %    'f-bound', fbound);
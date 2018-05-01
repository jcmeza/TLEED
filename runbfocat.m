%
% set initial iterate value
iterate = tleed_nomadm_x0;
x0 = iterate(1).x;
p  = cell2mat(iterate(1).p);

% first 14 variables are categorical; all others continuous
xp = {{'Ni', 'Ni', 'Ni', 'Ni', 'Ni', 'Li', 'Li', 'Li', 'Li', ...
    'Li', 'Li', 'Li', 'Li', 'Li', x0'}};

xcat  = repmat('s',14,1)';
xcont = repmat('c',42,1)';
xtype = strcat(xcat,xcont);
%
% other options
xlower = x0 - 0.5*max(1.0,abs(x0));
xupper = x0 + 0.5*max(1.0,abs(x0));
fbound = 1.0;

[ x, fx, msg, wrn, neval, f_hist] = bfo(@(x)kleedfcnc(x), x0, ...
        'epsilon', 0.0001, 'maxeval', 5000, ... 
        'xlower', xlower, 'xupper', xupper, ...
        'save-freq', 10, 'restart-file', 'bfo.restart');
        'f-call-type', 'with-bound', ...
        'f-bound', fbound);
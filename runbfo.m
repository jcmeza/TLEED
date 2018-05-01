%
% set initial iterate value
iterate = tleed_nomadm_x0;
x0 = iterate(1).x;
p  = cell2mat(iterate(1).p);

%
% other options
xlower = x0 - 0.5*max(1.0,abs(x0));
xupper = x0 + 0.5*max(1.0,abs(x0));
fbound = 1.0;

[ x, fx, msg, wrn, neval, f_hist] = bfo(@(x)kleedfcn(x,p), x0, ...
        'epsilon', 0.0001, 'maxeval', 100, ... 
        'xlower', xlower, 'xupper', xupper, ...
        'save-freq', 10, 'restart-file', 'bfo.restart');
    %    'f-call-type', 'with-bound', ...
    %    'f-bound', fbound);
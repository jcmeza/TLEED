function [fx,cx] = tleed_nomadm_Sur(x,p);
surrogate = getappdata(0,'SUR');
isLocal = 0;
n = length(x);
minDist = Inf;
for k = 1:size(surrogate.X,1)
   normterm = norm(transpose(x) - surrogate.X(k,1:n));
   minDist  = min(minDist, normterm);
   if normterm <= surrogate.trust
      isLocal = 1;
      break
   end
end
if isLocal
   y = x;
   if nargin > 1
      for k = 1:length(p)
         if ischar(p{k})
            y = [y; find(strcmp(surrogate.plist{k},p{k})) ];
         else
            y = [y; p{k}];
         end
      end
   end
   fx = feval(surrogate.evaluator,y,surrogate.f)
   fx = fx + feval(surrogate.searchFile,x,p)
   fx = fx - minDist*surrogate.dist;
end
if ~isLocal || ~isfinite(fx)
   fx = 1/eps;
end
cx = [];
return

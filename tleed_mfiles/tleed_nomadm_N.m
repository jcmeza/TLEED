
%*******************************************************************************
% tleed_nomadm_N:  User-supplied function defining set of neighbors for a
%                a given vector of categorical variables p.
% ------------------------------------------------------------------------------
%   Variables:
%     iterate = iterate for whom discrete neighbor point will be found.
%     plist   = cell array of possible p values: p{i}{j}, i = variable, j = list
%     Problem = structure holding Omega(p)
%     N       = vector of iterates who are neighbors of p.
%     delta   = mesh size parameter (used to make sure neighbors stay on mesh)
%*******************************************************************************
function N = tleed_nomadm_N(Problem,iterate,plist,delta);

N = [];
return

Param = getappdata(0,'PARAM');
%Mesh = inline('delta*round(x/delta)','x','delta');

Nump = length(iterate.p);
N = [];

% Include neighbors in which 1 atom changes its identity
for k = 1:Nump
   for j = 1:length(plist{k})
      if iterate.p{k} ~= plist{k}{j}
         neighbor.x = iterate.x;
         neighbor.p = iterate.p;
         neighbor.p{k} = plist{k}{j};
         N = [N neighbor];
      end
   end
end

% Include random neighbors without changing .x 
neighbor.x=iterate.x;
%nt=mod(floor(rand(Nump*10,Nump)*10), 2)+1;
nt=mod(floor(rand(Nump,Nump)*10), 2)+1;
%for II=1:Nump*10
for II=1:Nump
   for k = 1:Nump
       neighbor.p{k}=nt(II,k);
   end
   N = [N neighbor];
end

%%flip all identities
%for k = 1:Nump
%   for j = 1:length(plist{k})
%      if iterate.p{k} ~= plist{k}{j}
%         neighbor.p{k} = plist{k}{j};
%      end
%   end
%end
%neighbor.x=iterate.x;
%N = [N neighbor];

%%in favor of {1 1 1 ... 1 2 2 2 ...2}
%y=floor(rand(1,14)*10);
%ind=mod(sum(y),10);
%for i=1:ind
%       neighbor.p{i}=1;
%end
%for i=1+ind:Nump
%       neighbor.p{i}=2;
%end
%neighbor.x=iterate.x;
%N = [N neighbor];
%
% Include neighbors in which 2 atoms change their identities
%for k = 1:Nump-1
%	for l = k+1:Nump
%   		for j = 1:length(plist{k})
%      			if iterate.p{k} ~= plist{k}{j}  
%         			neighbor.x = iterate.x;
%         			neighbor.p = iterate.p;
%         			neighbor.p{k} = plist{k}{j};
%				for i=1:length(plist{l})
%                        		if iterate.p{l} ~= plist{l}{i}    
%                                		neighbor.p{l} = plist{l}{i};
%         					N = [N neighbor];
%					end
%				end
%      			end
%   		end
%	end
%end

%tleed:make point set physically feasible by going through the eval.f
%call:	iterates_filtered=eval_filter(iterates)
function iterates_filtered=eval_filter(iterates)
i=0;
for k = 1:length(iterates)
   [xtmp,ptmp,isValid]= deal_xp(iterates(k).x,iterates(k).p);
   if isValid
      i=i+1;
      iterates_filtered(i).x=xtmp;
      iterates_filtered(i).p=ptmp;
   end
end



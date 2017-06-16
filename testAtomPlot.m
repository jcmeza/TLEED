%===============================================================================
% testAtomPlot:  Plot molecular configurations from of LEED data
% ------------------------------------------------------------------------------
% Calls:  plotAtoms
% VARIABLES:
%   Param      = name of parameter file
%     iterate0 =   vector of initial iterates
%       .x     =     continuous  variable values
%       .p     =     categorical variable values
%   xVar       = matrix of atom xyz coordinates for current iterate
%===============================================================================
function testAtomPlot
 Param = tleed_Param;
 for k = 1:length(Param.iterate0)
    xVar = reshape(Param.iterate0(k).x,Param.nAtoms,Param.nDim);
    plotAtoms(k,xVar(:,[2,3,1]),Param.iterate0(k).p);
    saveas(gcf,['atoms',int2str(k),'.pdf'],'pdf');
 end
return

%===============================================================================
% plotAtoms:  Plot molecular configurations from of TLEED/KLEED data
% ------------------------------------------------------------------------------
% Called by:  testAtomPlot
% VARIABLES:
%   plotID = plot ID number for plot
%   xVar   = matrix of plot data, each row being xyz coordinates of an atom
%   p      = cell array containing type of each atom
%   scale  = scale factor for increasing the size of plot symbols
%   xVar1  = slice of xVar matrix with categorical variable value of 1
%   xVar2  = slice of xVar matrix with categorical variable value of 2
%   h      = vector of plot handles
%===============================================================================
function plotAtoms(plotID,xVar,p)
 if size(xVar,1) ~= length(p)
    error('xVar matrix must have same number of rows as elements in vector p');
 end
 scale = 20;
 figure;
 xVar1 = xVar([p{:}]==1,:);
 xVar2 = xVar([p{:}]==2,:);
 h(1) = scatter3(xVar1(:,1),xVar1(:,2),xVar1(:,3),'ko');
 hold on
 h(2) = scatter3(xVar2(:,1),xVar2(:,2),xVar2(:,3),'ko');
 axis([-8, 8, -8, 8, -8, 8]);
 set(h(1),'MarkerEdgeColor','k','MarkerFaceColor',[1.,.65,0.]);
 set(h(2),'MarkerEdgeColor','k','MarkerFaceColor','g');
 set(h(1),'SizeData',scale*get(h(1),'SizeData'));
 set(h(2),'SizeData',scale*get(h(2),'SizeData'));
 title(['Molecular Configuration:  Point',int2str(plotID)]);
 xlabel('x'); ylabel('y'); zlabel('z');
 legend('Ni','Li');
 hold off
return

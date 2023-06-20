function [numElemWithHeatGen, avTempAnnulus, globalAvTemp] = ...
    poisson(tempIntBd, tempExtBd, heatGen, iplots)
%close all
%clc

kc = 1.0; % W/K/cm

rHeatGen = 124.9;
rMin = 10.1;
rMax = 149.9;
R1 = 39.9;
R2 = 60.1;

eval('meshCercleForatQuad')

numNodes = size(nodes,1);
numElem = size(elem,1);

nodes = nodes(:,1:2);

if iplots == 1
    numbering = 0;
    plotElementsOld(nodes, elem, numbering)
    hold on
end

xx = nodes(:,1);
yy = nodes(:,2);
rads = sqrt(xx.^2 + yy.^2);

nodsIntBd = find(rads < rMin);
nodsExtBd = find(rads > rMax);
nodsAnnulus = find(rads > R1 & rads < R2);


a11 = kc;    %Coefficients for Exercise 2
a12 = 0.0;
a21 = a12;
a22 = a11;
a00 = 0.0;

coeff = [a11,a12,a21,a22,a00,0];
K = zeros(numNodes);
Q = zeros(numNodes,1);
F=zeros(numNodes,1);

numElemWithHeatGen = 0;

for e=1:numElem
    nodsElem = elem(e,:);
    vertexs = nodes(nodsElem,:);
    if min( sqrt(vertexs(:,1).^2 + vertexs(:,2).^2) > rHeatGen )
        numElemWithHeatGen = numElemWithHeatGen + 1;
        coeff(6) = heatGen;
        if iplots == 1
            fill(vertexs(:,1),vertexs(:,2),'yellow')
        end
    else
        coeff(6) = 0;
    end
    [Ke,Fe]=bilinearQuadElement(coeff,nodes,elem,e);
    %
    % Assemble the elements
    %
    rows=[elem(e,1); elem(e,2); elem(e,3); elem(e,4)];
    colums= rows;
    K(rows,colums)=K(rows,colums)+Ke; %assembly
    if (coeff(6) ~= 0)
        F(rows)=F(rows)+Fe;
    end
end 

if iplots == 1
    plot(nodes(nodsIntBd,1),nodes(nodsIntBd,2),'o',...
     'MarkerFaceColor','red','MarkerEdgeColor','black','MarkerSize',8)
    plot(nodes(nodsExtBd,1),nodes(nodsExtBd,2),'o',...
        'MarkerFaceColor','green','MarkerEdgeColor','black','MarkerSize',8)
    plot(nodes(nodsAnnulus,1),nodes(nodsAnnulus,2),'o',...
      'MarkerFaceColor','magenta','MarkerEdgeColor','black','MarkerSize',8)
    hold off
end

fixedNodes = [nodsIntBd',nodsExtBd']; %Fixed Nodes (global num.)
freeNodes = setdiff(1:numNodes,fixedNodes); %Complementary of fixed nodes

% ------------ Natural BC
Q(freeNodes)=0; %all of them are zero

% ------------ Essential BC
u=zeros(numNodes,1); %initialize u vector
u(nodsIntBd) = tempIntBd; 
u(nodsExtBd) = tempExtBd;
Fm=F(freeNodes)-K(freeNodes,fixedNodes)*u(fixedNodes);%here u can be 
                                                      %different from zero 
                                                      % modify the linear system, this is only valid if BC = 0.
%------------- Reduced system
Km=K(freeNodes,freeNodes);
Fm=Fm+Q(freeNodes);%only for fixed nodes

%Compute the solution 
%solve the System
um=Km\Fm;
u(freeNodes)=um;
%u

if iplots == 1
    %PostProcess: plot results
    titol='Equation solution';
    colorScale='jet';
    plotContourSolution(nodes,elem,u,titol,colorScale)
end

% Solutions
avTempAnnulus = sum(u(nodsAnnulus))/length(nodsAnnulus);
globalAvTemp = sum(u)/numNodes;

end
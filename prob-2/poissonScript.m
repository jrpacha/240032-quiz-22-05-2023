clearvars
close all
%clc

kc = 1.0;     
tempIntBd = 65.0;
tempExtBd = 95.0;
heatGen = 5.0;
elem_hint = 121;
rHeatGen = 124.9;
rMin = 10.1;
rMax = 149.9;
R1 = 39.9;
R2 = 60.1;

fprintf('Problem-2\n')
fprintf('Temperature at internal boundary, tempIntBd = %e\n', tempIntBd)
fprintf('Temperature at external boundary, tempExtBd = %e\n', tempExtBd)
fprintf('Heat generation at R > %f, f = %e\n\n', rHeatGen, heatGen)

eval('meshCercleForatQuad')

numNodes = size(nodes,1);
numElem = size(elem,1);

nodes = nodes(:,1:2);

numbering = 0;
plotElementsOld(nodes, elem, numbering)
hold on

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
        fill(vertexs(:,1),vertexs(:,2),'yellow')
    else
        coeff(6) = 0;
    end
    [Ke,Fe]=bilinearQuadElement(coeff,nodes,elem,e);
    if e == elem_hint
        fprintf(['Part (a)\n',...
            '   for element e = %d:\n',...
            '           F(%d) = %22.15e\n',...
            'Hint-1. K(%d,%d) = %22.15e\n\n'],...
            elem_hint,4,Fe(4),3,3,Ke(3,3))
    end
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

plot(nodes(nodsIntBd,1),nodes(nodsIntBd,2),'o',...
    'MarkerFaceColor','red','MarkerEdgeColor','black','MarkerSize',8)
plot(nodes(nodsExtBd,1),nodes(nodsExtBd,2),'o',...
    'MarkerFaceColor','green','MarkerEdgeColor','black','MarkerSize',8)
plot(nodes(nodsAnnulus,1),nodes(nodsAnnulus,2),'o',...
    'MarkerFaceColor','magenta','MarkerEdgeColor','black','MarkerSize',8)
hold off

%Boundary Conditions (BC)

fixedNodes = [nodsIntBd', nodsExtBd']; %Fixed Nodes (global num.)
freeNodes = setdiff(1:numNodes,fixedNodes); %Complementary of fixed nodes

% ------------ Natural BC
Q(freeNodes)=0; %all of them are zero

% ------------ Essential BC
u = zeros(numNodes,1); %initialize u vector
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
format short e; %just to a better view of the numbers
format compact;
um=Km\Fm;
u(freeNodes)=um;

%PostProcess: compute Q and plot results
Q = K*u - F;
titol='Equation solution';
colorScale='jet';
plotContourSolution(nodes,elem,u,titol,colorScale)

% Solutions
avQ = sum(Q(nodsExtBd))/length(nodsExtBd);
avTempAnnulus = sum(u(nodsAnnulus))/length(nodsAnnulus);
globalAvTemp = sum(u)/numNodes;


%fprintf('Number of elements with heat generation: %d\n', numElemWithHeatGen)
fprintf('Part (b)\n')
fprintf('Average of the Q values at the outer border (R = %f), <Q(outer)> = %.4e\n',...
    R2, avQ)
fprintf('Hint-2. The number of elements in C2 is: %d\n\n',...
    numElem - numElemWithHeatGen)
fprintf('Part (c)\n')
fprintf('Averaged temperature in the annulus (%f < R < %f), <T(annulus)> = %.4e\n',...
    R1, R2, avTempAnnulus)
fprintf('Hint3. The number of nodes of this intermediate annulus is: %d\n',...
    length(nodsAnnulus))
%fprintf('Hint: global Averaged temperature, <T> = %.4e\n', globalAvTemp)

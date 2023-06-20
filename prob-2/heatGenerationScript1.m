clearvars
close all
%clc

kc = 1.0;            % W/K/cm
tempIntBd = 55.0;
tempExtBd = 85.0;
heatGen = 1.0;
rHeatGen = 124.9;
rMin = 10.1;
rMax = 149.9;
R = 40.0;
R1 = 39.9;
R2 = 60.1;
s = 2.0;

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
K = zeros(numNodes+1);
Q = zeros(numNodes+1,1);
F=zeros(numNodes+1,1);

numElemWithHeatGen = 0;

for e=1:numElem
    nodsElem = elem(e,:);
    vertexs = nodes(nodsElem,:);
    if min( sqrt(vertexs(:,1).^2 + vertexs(:,2).^2) > rHeatGen )
        numElemWithHeatGen = numElemWithHeatGen + 1;
        coeff(6) = heatGen;
        fill(vertexs(:,1),vertexs(:,2),'green')
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

plot(nodes(nodsIntBd,1),nodes(nodsIntBd,2),'o',...
    'MarkerFaceColor','red','MarkerEdgeColor','black','MarkerSize',8)
plot(nodes(nodsExtBd,1),nodes(nodsExtBd,2),'o',...
    'MarkerFaceColor','green','MarkerEdgeColor','black','MarkerSize',8)
plot(nodes(nodsAnnulus,1),nodes(nodsAnnulus,2),'o',...
    'MarkerFaceColor','magenta','MarkerEdgeColor','black','MarkerSize',8)
hold off

K(1:end,end) = -F;
K(end,nodsAnnulus) = 1;

%Boundary Conditions (BC)
fixedNodes = [nodsIntBd', nodsExtBd']; %Fixed Nodes (global num.)
freeNodes = setdiff(1:numNodes+1,fixedNodes); %Complementary of fixed nodes

% ------------ Natural BC
Q(freeNodes)=0; %all of them are zero
Q(end) = tempIntBd*s*length(nodsAnnulus);

% ------------ Essential BC
u = zeros(numNodes,1); %initialize u vector
u(nodsIntBd) = tempIntBd; 
u(nodsExtBd) = tempExtBd;
Fm=Q(freeNodes)-K(freeNodes,fixedNodes)*u(fixedNodes);%here u can be 
                                                      %different from zero 
                                                      % modify the linear system, this is only valid if BC = 0.

%------------- Reduced system
Km=K(freeNodes,freeNodes);

%Compute the solution 
%solve the System
format long e; %just to a better view of the numbers
format compact;
um=Km\Fm;
u(freeNodes)=um;

%PostProcess: plot results
titol='Equation solution';
colorScale='jet';
plotContourSolution(nodes,elem,u(1:end-1),titol,colorScale)

% Solutions
hetGenLinearSystem = u(end);
fprintf('Temperature at internal boundary, tempIntBd = %e\n', tempIntBd)
fprintf('Temperature at external boundary, tempExtBd = %e\n', tempExtBd)
fprintf('Heat generation at R > %f, f = %e\n', rHeatGen, hetGenLinearSystem)
fprintf('Number of elements with heat generation: %d\n', numElemWithHeatGen)

% Check Solution
[numElemWithHeatGen, avTempAnnulus, globalAvTemp] = ...
    poisson(tempIntBd, tempExtBd, hetGenLinearSystem, 1);

fprintf('Averaged temprature in the annulus %f < R < %f, <T> = %.4e\n',...
    R1, R2, avTempAnnulus)
fprintf('Err:= |<T at R> - %d x intBdTemp| = |<T at R> - %f | = %.4e\n',...
    s, s*tempIntBd, abs(avTempAnnulus-s*tempIntBd))
fprintf('Hint: global Averaged temperature, <T> = %.4e\n', globalAvTemp)
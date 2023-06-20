clearvars
close all

Y = 150.0;         %GPa
areaSec = 3400.0;  %mm^2

%Geometry
nodes = [0, 0;...
    1250.0, 1000.0;...
    1750.0, 1000.0;...
    3000.0, 0;...
    1500.0, 2000.0];

%nodes = 1.0e3*nodes;

elem = [1, 2;...
    2, 3;...
    3, 4;...
    1, 5;...
    2, 5;...
    3, 5;...
    4, 5];

numNodes = size(nodes,1);
numElem = size(elem,1);
ndim = size(nodes,2);

numbering = 1;
plotElementsOld(nodes, elem, numbering)

%% Real constants
E = Y*ones(1,numElem);
A = areaSec*ones(1,numElem);

K = zeros(ndim*numNodes);
F = zeros(ndim*numNodes,1);
Q = zeros(ndim*numNodes,1);

for e = 1:numElem
    Ke=planeLinkStiffMatrix(nodes,elem,e,E,A);
    rows = [ndim*elem(e,1)-1, ndim*elem(e,1),...
        ndim*elem(e,2)-1, ndim*elem(e,2)];
    cols = rows;
    K(rows,cols) = K(rows,cols) + Ke;
end

% Node 1 is attached while node 4 only 
% displacement in the x-direction is allowed
fixedNodes = [1, 2, 8];
freeNodes = setdiff(1:ndim*numNodes,fixedNodes);

%Natural B.C.
%In node a force of 200kN is loaded in the 
%x-direction, and a force of the same 
%magnitude is loaded in node 5 but in the 
%direction of the bar 7 from towards node 4.
nod = 3;
Q(ndim*nod-1) = 200; %kN

nod = 5;
cosAlpha = 0.8;
sinAlpha = 0.6;

Q(ndim*nod-1) = sinAlpha*200; %kN
Q(ndim*nod) = -cosAlpha*200;  %kN

%Esential B.C.
u = zeros(ndim*numNodes,1);
u(fixedNodes) = 0; %not necessary

%Reduced system 
Im = Q(freeNodes);
Km = K(freeNodes,freeNodes);

%Solve the reduced system
um = Km\Im;
u(freeNodes) = um;

desplX = u(1:ndim:end); %displacements in x-directon
desplY = u(2:ndim:end); %displacements in y-directon

format long e
format compact
[desplX, desplY]

%Final lenght of bar e = 7
e = 7;
n1 = elem(e,1); n2 = elem(e,2);
z1 = nodes(n1,:) + [desplX(n1), desplY(n1)]; 
z2 = nodes(n2,:) + [desplX(n2), desplY(n2)]; 

finalLengthElem = norm(z2-z1) 

%Final position of y component of node 5
nod = 5;
finalPositionNod = nodes(nod,2) + desplY(nod) 




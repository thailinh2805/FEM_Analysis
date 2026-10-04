clc; clear all;
%% input data
E=2e11;
nu=0.3;
thickness=1; %thickness plane strain = 1 
r = 0.5;
R = 1;
p_in = 2e7;
p_out = -1e7;

Db=((1-nu)*E/((1+nu)*(1-2*nu)))*[1, nu/(1-nu), 0; 
                                nu/(1-nu), 1, 0; 
                                0, 0, (1-2*nu)/(2*(1-nu))];

%% mesh
nnr = 20;
nnphi = 20;
elemType = 'Q4';  % Switch between 'Q4' and 'T3'
switch elemType
    case 'Q4'
        nnode = 4;
        quadorder = 3;
    case 'T3'
        nnode = 3;
        quadorder = 3;  % or 1 for fewer points
end

[node, element,numnode,numElem] = mesh2D(elemType,r,R,nnr,nnphi);

numelem = size(element,1);
% Vẽ mesh
figure; hold on;

for e = 1:numelem
    nodes_e = element(e,:); % lay data node cua tung element
    coord = node(nodes_e,:); % toa do tuong ung cua node
    patch(coord(:,1), coord(:,2), 'w','edgecolor','k');
end

plot(node(:,1), node(:,2), 'r.');

axis equal;
xlim([0 R]);
ylim([0 R]);
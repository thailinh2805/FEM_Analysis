function [node, element,numnode, numelem] = mesh2D(type,r,R,nnr,nnphi)

rnode = linspace(r, R, nnr);
phinode = linspace(0, pi/2, nnphi);

node = [];

% tao node + lay gtri (x,y)
for i = 1:nnr
    for j = 1:nnphi
        x = rnode(i)*cos(phinode(j));
        y = rnode(i)*sin(phinode(j));
        node = [node; 
                x y];
    end
end

numnode = size(node,1);

% xếp lại sao cho các nút cùng dr nằm trên 1 hàng, cùng dp trên 1 cột
node_id = reshape(1:numnode, nnphi, nnr)';

n1 = node_id(1:end-1, 2:end);
n2 = node_id(1:end-1, 1:end-1);
n3 = node_id(2:end, 1:end-1);
n4 = node_id(2:end, 2:end);

e1 = reshape(n1', [], 1); 
e2 = reshape(n2', [], 1);
e3 = reshape(n3', [], 1);
e4 = reshape(n4', [], 1);

switch type
     case 'Q4'
         element = [e1, e2, e3, e4];
 
     case 'T3'
        elem1 = [e1, e2, e4];
        elem2 = [e2, e3, e4];
        element = [elem1; elem2];
 end

numelem = size(element,1);

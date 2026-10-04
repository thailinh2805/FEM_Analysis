function [Nv, dNdxi] = Larange2D(type, pt)
xi = pt(1); 
eta = pt(2);

    switch type
        case 'T3'
            N = [1 - xi - eta;
                 xi;
                 eta];

            dNdxi = [-1, -1;
                      1,  0;
                      0,  1];

        case 'Q4'
            N = 1/4 * [(1 - xi)*(1 - eta);  % Nút 1: (-1, -1)
                       (1 + xi)*(1 - eta);  % Nút 2: ( 1, -1)
                       (1 + xi)*(1 + eta);  % Nút 3: ( 1,  1)
                       (1 - xi)*(1 + eta)]; % Nút 4: (-1,  1)

            dNdxi = 1/4 * [-(1 - eta), -(1 - xi);  % dN1/dxi, dN1/deta
                            (1 - eta), -(1 + xi);  % dN2/dxi, dN2/deta
                            (1 + eta),  (1 + xi);  % dN3/dxi, dN3/deta
                           -(1 + eta),  (1 - xi)]; % dN4/dxi, dN4/deta
    end

    I = eye(2);
    Nv = [];
    for a = 1:size(N,1)
        Nv = [Nv; I*N(a)];
    end
end
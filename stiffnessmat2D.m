function [K] = stiffnessmat2D(elemType, quadorder, numElem, numnode, node, element, Db, thickness)
    ndof = 2; 
    rho = 7830; % Khối lượng riêng (kg/m3)
    tdof = ndof * numnode;
    K = sparse(tdof, tdof);
    M = sparse(tdof, tdof);
    nnode = size(element, 2); 
    
    for iel = 1:numElem
        conn = element(iel,:);    
        xe_coord = node(conn, 1); 
        ye_coord = node(conn, 2); 
        
        entry = zeros(1, ndof * nnode);
        for i = 1:nnode
            entry(ndof*i-1) = ndof*conn(i)-1;
            entry(ndof*i)   = ndof*conn(i);
        end

        switch elemType
            case 'T3'
                Ae = 0.5 * abs(xe_coord(1)*(ye_coord(2) - ye_coord(3)) + ...
                               xe_coord(2)*(ye_coord(3) - ye_coord(1)) + ...
                               xe_coord(3)*(ye_coord(1) - ye_coord(2)));

                b = [ye_coord(2)-ye_coord(3), ye_coord(3)-ye_coord(1), ye_coord(1)-ye_coord(2)] / (2*Ae);
                c = [xe_coord(3)-xe_coord(2), xe_coord(1)-xe_coord(3), xe_coord(2)-xe_coord(1)] / (2*Ae);
                
                B = zeros(3, ndof * 3);    
                    for m = 1:3
                        B(:,(2*m-1):(2*m)) = [b(m), 0;
                                                0, c(m);
                                              c(m), b(m)];
                    end
                    Ke = (B' * Db * B) * Ae * thickness;
                    K(entry, entry) = K(entry, entry) + Ke;

            case 'Q4'
                [w, q] = quadrature2D(elemType, quadorder); 
                xe_mat = node(conn, :); 
                
                for ig = 1:size(q,1)
                    pt = q(ig,:);
                    weight = w(ig);
                    [~, dNdxi] = Larange2D(elemType, pt);
                    
                    J0 = dNdxi' * xe_mat; 
                    jac = det(J0);
                    dNdx = J0 \ dNdxi'; 
                    
                    B = zeros(3, ndof * 4);
                    for m = 1:4
                        B(:,(2*m-1):(2*m)) = [dNdx(1,m), 0;
                                                0, dNdx(2,m);
                                              dNdx(2,m), dNdx(1,m)];
                    end
                    Ke = (B' * Db * B) * weight * jac * thickness;
                    K(entry, entry) = K(entry, entry) + Ke;
                end
        end
    end
end
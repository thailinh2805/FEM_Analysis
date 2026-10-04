function [vm, sigma_r, sigma_theta, sigma_z, sigma_x, sigma_y, tau_xy] = stress_cal(elemType, quadorder, numElem, node, element, Db, nu, u)
    ndof = 2;
    vm = zeros(numElem, 1);
    sigma_r = zeros(numElem, 1);
    sigma_theta = zeros(numElem, 1);
    sigma_z_val = zeros(numElem, 1);
    sigma_x     = zeros(numElem, 1);
    sigma_y     = zeros(numElem, 1);
    tau_xy      = zeros(numElem, 1);

    for iel = 1:numElem
        conn = element(iel,:);    
        xe_coord = node(conn, 1); 
        ye_coord = node(conn, 2); 
        
        u_element = [];
        for i = 1:length(conn)
            u_element = [u_element; u(2*conn(i)-1); u(2*conn(i))];
        end
        
        switch elemType
            case 'T3'
                Ae = 0.5 * abs(xe_coord(1)*(ye_coord(2) - ye_coord(3)) + ...
                               xe_coord(2)*(ye_coord(3) - ye_coord(1)) + ...
                               xe_coord(3)*(ye_coord(1) - ye_coord(2)));
                
                b_val = [ye_coord(2)-ye_coord(3), ye_coord(3)-ye_coord(1), ye_coord(1)-ye_coord(2)] / (2*Ae);
                c_val = [xe_coord(3)-xe_coord(2), xe_coord(1)-xe_coord(3), xe_coord(2)-xe_coord(1)] / (2*Ae);
                
                B = zeros(3, ndof*3);    
                for m = 1:3
                    B(:,(2*m-1):(2*m)) = [b_val(m), 0;
                                            0, c_val(m);
                                          c_val(m), b_val(m)];
                end
                stress = Db * B * u_element;
                
            case 'Q4'
                [w, q] = quadrature2D(elemType, quadorder); 
                xe_mat = node(conn, :); 
                for ig = 1:size(q,1)
                    pt = q(ig,:);
                    [~, dNdxi] = Larange2D(elemType, pt); 
                    
                    J0 = dNdxi' * xe_mat; 
                    dNdx = J0 \ dNdxi'; 
                
                B = zeros(3, ndof*4);
                for m = 1:4
                    B(:,(2*m-1):(2*m)) = [dNdx(1,m), 0;
                                            0, dNdx(2,m);
                                          dNdx(2,m), dNdx(1,m)];
                end
                stress = Db * B * u_element;
                end
                
        end 

        % dùng biến tạm scalar
        sx  = stress(1);
        sy  = stress(2);
        txy = stress(3);
    

        % dùng biến tạm scalar
        sx  = stress(1);
        sy  = stress(2);
        txy = stress(3);
        sz  = nu * (sx + sy);   % plane strain

        sigma_x(iel) = sx;
        sigma_y(iel) = sy;
        tau_xy(iel)  = txy;
        sigma_z(iel) = sz;


        vm(iel) = sqrt(0.5 * ( ...
            (sx - sy)^2 + ...
            (sy - sz)^2 + ...
            (sz - sx)^2 + ...
            6*txy^2 ));
        
        center = mean(node(conn,:));
        theta = atan2(center(2), center(1));
        c = cos(theta); s = sin(theta);

        sigma_r(iel)     = sx*c^2 + sy*s^2 + 2*txy*s*c;
        sigma_theta(iel) = sx*s^2 + sy*c^2 - 2*txy*s*c;
        
    end 
end
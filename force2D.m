function [f, nodes_in, nodes_out] = force2D(node, element, thickness, p_in, r_in, p_out, r_out)
    numnode = size(node, 1);
    tdof = 2 * numnode;
    f = zeros(tdof, 1); 
        
    dist = sqrt(node(:,1).^2 + node(:,2).^2);
    nodes_in = find(abs(dist - r_in) < 1e-6);
    nodes_out = find(abs(dist - r_out) < 1e-6);
    
    coord_in = node(nodes_in, :);
    coord_out = node(nodes_out, :);
    
    % inner boundary (p_in)
    num_nodes_in = length(nodes_in);
    for i = 1 : (num_nodes_in - 1)
        node_a = nodes_in(i);
        node_b = nodes_in(i+1); 
        
        xa = node(node_a, 1); ya = node(node_a, 2);
        xb = node(node_b, 1); yb = node(node_b, 2);
        
        L = sqrt((xb - xa)^2 + (yb - ya)^2);
        
        % phap tuyen trung diem
        mid_x = (xa + xb) / 2;
        mid_y = (ya + yb) / 2;
        r_mid = sqrt(mid_x^2 + mid_y^2);
        
        nx = mid_x / r_mid; 
        ny = mid_y / r_mid;
        
        % do lon luc phanbo  p*L*h/2
        mag_in = (p_in * L * thickness) / 2;
        
        % cong vao vector tong
        f(2*node_a-1) = f(2*node_a-1) + mag_in * nx;
        f(2*node_a)   = f(2*node_a)   + mag_in * ny;
        
        f(2*node_b-1) = f(2*node_b-1) + mag_in * nx;
        f(2*node_b)   = f(2*node_b)   + mag_in * ny;
    end

    % outter boundary 
    num_nodes_out = length(nodes_out);
    for i = 1 : (num_nodes_out - 1)
        node_c = nodes_out(i);
        node_d = nodes_out(i+1); 
        
        xc = node(node_c, 1); yc = node(node_c, 2);
        xd = node(node_d, 1); yd = node(node_d, 2);
        
        L = sqrt((xd - xc)^2 + (yd - yc)^2);
        
        mid_x = (xc + xd) / 2;
        mid_y = (yc + yd) / 2;
        r_mid = sqrt(mid_x^2 + mid_y^2);
        
        nx = mid_x / r_mid; 
        ny = mid_y / r_mid;
        
        mag_out = (p_out * L * thickness) / 2;
        
        f(2*node_c-1) = f(2*node_c-1) + mag_out * nx;
        f(2*node_c)   = f(2*node_c)   + mag_out * ny;
        f(2*node_d-1) = f(2*node_d-1) + mag_out * nx;
        f(2*node_d)   = f(2*node_d)   + mag_out * ny;
    end

end

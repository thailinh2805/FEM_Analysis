function [K, M] = stiffnessmat2D_modal(elemType, quadorder, numElem, numnode, node, element, Db, thickness)
    ndof = 2; 
    rho = 7830; % Khối lượng riêng (kg/m3)
    tdof = ndof * numnode;
    
    % Sử dụng ma trận thưa để tối ưu bộ nhớ
    K = sparse(tdof, tdof);
    M = sparse(tdof, tdof);
    
    nnode_per_elem = size(element, 2); 
    
    for iel = 1:numElem
        conn = element(iel, :); % Danh sách các nút của phần tử hiện tại
        xe_mat = node(conn, :); % Tọa độ (x, y) của 4 nút
        
        % Thiết lập vector chỉ số bậc tự do (DOFs) cho phần tử
        % entry sẽ có dạng: [2*n1-1, 2*n1, 2*n2-1, 2*n2, ...]
        entry = zeros(1, ndof * nnode_per_elem);
        for i = 1:nnode_per_elem
            entry(2*i-1) = 2*conn(i) - 1;
            entry(2*i)   = 2*conn(i);
        end
        
        switch elemType
            case 'Q4'

                % Các điểm tích phân Gauss 2x2 (Trọng số w = 1 nên không cần nhân)
                gp = 1 / sqrt(3);
                Q = [-gp, -gp; 
                      gp, -gp; 
                      gp,  gp; 
                     -gp,  gp];
                
                Ke = zeros(8, 8);
                Me = zeros(8, 8);
                
                for g = 1:size(Q, 1)
                    xi = Q(g, 1);
                    eta = Q(g, 2);
                    
                    % Ma trận hàm dạng N và Nm
                    N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
                    Nm = zeros(2, 8);
                    Nm(1, 1:2:end) = N;
                    Nm(2, 2:2:end) = N;
                    
                    % Đạo hàm hàm dạng trong hệ tọa độ tự nhiên
                    dxi = 0.25 * [-(1-eta), 1-eta, 1+eta, -(1+eta)];
                    deta = 0.25 * [-(1-xi), -(1+xi), 1+xi, 1-xi];
                    
                    % TÍNH JACOBIAN BẰNG XE_MAT (Tối ưu hóa phép nhân ma trận)
                    J = [dxi; deta] * xe_mat; 
                    detJ = det(J);
                    
                    % Đạo hàm hàm dạng trong hệ Đề-các
                    dxy = J \ [dxi; deta];
                    
                    % Ma trận biến dạng - chuyển vị B
                    B = zeros(3, 8);
                    B(1, 1:2:end) = dxy(1, :);
                    B(2, 2:2:end) = dxy(2, :);
                    B(3, 1:2:end) = dxy(2, :); B(3, 2:2:end) = dxy(1, :);
                    
                    % Tích phân cộng dồn ma trận phần tử
                    Ke = Ke + (B' * Db * B) * thickness * detJ;
                    Me = Me + rho * thickness * (Nm' * Nm) * detJ;
                end
                
                % Lắp ghép cục bộ vào tổng thể
                K(entry, entry) = K(entry, entry) + Ke;
                M(entry, entry) = M(entry, entry) + Me;

            case 'T3'
                % Tính diện tích Ae
                x = xe_mat(:,1); y = xe_mat(:,2);
                Ae = 0.5 * abs(x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2)));
                
                % Đạo hàm hàm dạng cho T3 là hằng số
                b = [y(2)-y(3); y(3)-y(1); y(1)-y(2)] / (2*Ae);
                c = [x(3)-x(2); x(1)-x(3); x(2)-x(1)] / (2*Ae);
                
                B = zeros(3, 6);
                for m = 1:3
                    B(:, 2*m-1:2*m) = [b(m), 0; 0, c(m); c(m), b(m)];
                end
                
                K(entry, entry) = K(entry, entry) + (B' * Db * B) * Ae * thickness;
                
                % Ma trận khối lượng nhất quán cho T3
                Me_const = (rho * Ae * thickness / 12) * ...
                           [2, 0, 1, 0, 1, 0;
                            0, 2, 0, 1, 0, 1;
                            1, 0, 2, 0, 1, 0;
                            0, 1, 0, 2, 0, 1;
                            1, 0, 1, 0, 2, 0;
                            0, 1, 0, 2, 0, 2];
                M(entry, entry) = M(entry, entry) + Me_const;
        end
    end
end
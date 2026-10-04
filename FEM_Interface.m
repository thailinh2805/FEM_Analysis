classdef FEM_Interface < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                     matlab.ui.Figure
        CONVERGCEButton              matlab.ui.control.Button
        TYPEDropDown_2               matlab.ui.control.DropDown
        TYPEDropDown_2Label          matlab.ui.control.Label
        MODESpinner                  matlab.ui.control.Spinner
        MODESpinnerLabel             matlab.ui.control.Label
        MODALButton                  matlab.ui.control.Button
        RESETButton                  matlab.ui.control.Button
        TextArea                     matlab.ui.control.TextArea
        TextAreaLabel                matlab.ui.control.Label
        FEMPROGRAMINTERFACELabel     matlab.ui.control.Label
        APPLIEDFORCEButton           matlab.ui.control.Button
        MESHINGButton                matlab.ui.control.Button
        TYPEDropDown                 matlab.ui.control.DropDown
        TYPEDropDownLabel            matlab.ui.control.Label
        StatusLabel                  matlab.ui.control.Label
        RUNANALYSISButton            matlab.ui.control.Button
        ElementTypeSwitch            matlab.ui.control.Switch
        ElementTypeSwitchLabel       matlab.ui.control.Label
        NodesAngularEditField        matlab.ui.control.NumericEditField
        NodesAngularEditFieldLabel   matlab.ui.control.Label
        NodesRadialEditField         matlab.ui.control.NumericEditField
        NodesRadialEditFieldLabel    matlab.ui.control.Label
        OuterPressureEditField       matlab.ui.control.NumericEditField
        OuterPressureEditFieldLabel  matlab.ui.control.Label
        InnerPressureEditField       matlab.ui.control.NumericEditField
        InnerPressureEditFieldLabel  matlab.ui.control.Label
        OuterRadiusREditField        matlab.ui.control.NumericEditField
        OuterRadiusREditFieldLabel   matlab.ui.control.Label
        InnerRadiusrEditField        matlab.ui.control.NumericEditField
        InnerRadiusrEditFieldLabel   matlab.ui.control.Label
        ThicknessEditField           matlab.ui.control.NumericEditField
        ThicknessEditFieldLabel      matlab.ui.control.Label
        PoisonEditField              matlab.ui.control.NumericEditField
        PoisonEditFieldLabel         matlab.ui.control.Label
        YoungModulusEditField        matlab.ui.control.NumericEditField
        YoungModulusEditFieldLabel   matlab.ui.control.Label
        UIAxes_2                     matlab.ui.control.UIAxes
        UIAxes                       matlab.ui.control.UIAxes
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: YoungModulusEditField
        function YoungModulusEditFieldValueChanged(app, event)
            value = app.YoungModulusEditField.Value;
            
        end

        % Value changed function: PoisonEditField
        function PoisonEditFieldValueChanged(app, event)
            value = app.PoisonEditField.Value;
            
        end

        % Value changed function: ThicknessEditField
        function ThicknessEditFieldValueChanged(app, event)
            value = app.ThicknessEditField.Value;
            
        end

        % Value changed function: InnerRadiusrEditField
        function InnerRadiusrEditFieldValueChanged(app, event)
         value = app.InnerRadiusrEditField.Value;
            
        end

        % Value changed function: OuterRadiusREditField
        function OuterRadiusREditFieldValueChanged(app, event)
            value = app.OuterRadiusREditField.Value;
            
        end

        % Value changed function: InnerPressureEditField
        function InnerPressureEditFieldValueChanged(app, event)
          value = app.InnerPressureEditField.Value;
            
        end

        % Value changed function: OuterPressureEditField
        function OuterPressureEditFieldValueChanged(app, event)
           value = app.OuterPressureEditField.Value;
            
        end

        % Value changed function: NodesRadialEditField
        function NodesRadialEditFieldValueChanged(app, event)
            value = app.NodesRadialEditField.Value;
           
        end

        % Value changed function: NodesAngularEditField
        function NodesAngularEditFieldValueChanged(app, event)
          value = app.NodesAngularEditField.Value;
            
        end

        % Value changed function: ElementTypeSwitch
        function ElementTypeSwitchValueChanged(app, event)
            value = app.ElementTypeSwitch.Value;
                switch value
                    case 'Q4'
                        nnode = 4;
                        quadorder = 3;
                    case 'T3'
                        nnode = 3;
                        quadorder = 3;
                end
         
        end

        % Button pushed function: RUNANALYSISButton
        function RUNANALYSISButtonPushed(app, event)
             try
                %% input
                E = app.YoungModulusEditField.Value;
                nu = app.PoisonEditField.Value;
                thickness = app.ThicknessEditField.Value;
                r = app.InnerRadiusrEditField.Value;
                R = app.OuterRadiusREditField.Value;
                p_in = app.InnerPressureEditField.Value;
                p_out = app.OuterPressureEditField.Value;
            
                Db=((1-nu)*E/((1+nu)*(1-2*nu)))*[1, nu/(1-nu), 0; 
                                                 nu/(1-nu), 1, 0; 
                                                 0, 0, (1-2*nu)/(2*(1-nu))];
            
                nnr = app.NodesRadialEditField.Value;
                nnphi = app.NodesAngularEditField.Value;
                elemType = app.ElementTypeSwitch.Value;
            
                switch elemType
                    case 'Q4'
                        nnode = 4;
                        quadorder = 3;
                    case 'T3'
                        nnode = 3;
                        quadorder = 3;
                end
            
                [node, element, numnode, numElem] = mesh2D(elemType, r, R, nnr, nnphi);
            
                ndof = 2;
                tdof = ndof * numnode;
                K = stiffnessmat2D(elemType, quadorder, numElem, numnode, node, element, Db, thickness);
            
                f = force2D(node, element, thickness, p_in, r, p_out, R);
            
                bottom = find(abs(node(:,2)) < 1e-6);
                left = find(abs(node(:,1)) < 1e-6);
            
                DOFRoller = unique([bottom*2; left*2-1]);
            
                sdof = (1:tdof)';
                UnknownDof = sdof;
                UnknownDof(DOFRoller) = [];
            
                nUnknownDof = length(UnknownDof);
            
                Tbd = zeros(tdof, nUnknownDof);
                for i = 1:nUnknownDof
                    Tbd(UnknownDof(i), i) = 1;
                end
            
                Kred = Tbd' * K * Tbd;
                fred = Tbd' * f;
            
                displacement = Kred \ fred;
                u = Tbd * displacement;
            
                ux = u(1:2:end);
                uy = u(2:2:end);
            
                x = node(:,1);
                y = node(:,2);
            
                theta = atan2(y, x);
            
                ur = ux .* cos(theta) + uy .* sin(theta);
                utheta = -ux .* sin(theta) + uy .* cos(theta);
                
    
                [vm, sigma_r, sigma_theta] = stress_cal(elemType, quadorder, numElem, node, element, Db, nu, u);
                app.TextArea.Value = {
                    '--- KET QUA CHUYEN VI ---'
                    ['Ux = ' num2str(max(abs(ux)), '%.4e') ...
                     '   |   Uy = ' num2str(max(abs(uy)), '%.4e') ...
                     '   |   Ur = ' num2str(max(abs(ur)), '%.4e') ...
                     '   |   Utheta = ' num2str(max(abs(utheta)), '%.4e')]
                '-================================================='
                    '--- KET QUA UNG SUAT ---'
                    ['Sigma_r = ' num2str(max(sigma_r), '%.4e') ...
                     '   |   Sigma_theta = ' num2str(max(sigma_theta), '%.4e') ...
                     '   |   VonMises max = ' num2str(max(vm), '%.4e') ...
                     '   |   VonMises min = ' num2str(min(vm), '%.4e')]
                };             
                                
                
                % SUCCESS
                app.StatusLabel.Text = '✔ Analysis completed successfully';
                app.StatusLabel.FontColor = [0 0.6 0]; % xanh lá

            catch ME
        %ERROR
                app.StatusLabel.Text = ['✖ Error: ' ME.message];
                app.StatusLabel.FontColor = [1 0 0]; % đỏ
            end
        end

        % Button pushed function: MESHINGButton
        function MESHINGButtonPushed(app, event)
            E = app.YoungModulusEditField.Value;
            nu = app.PoisonEditField.Value;
            thickness = app.ThicknessEditField.Value;
            r = app.InnerRadiusrEditField.Value;
            R = app.OuterRadiusREditField.Value;
            p_in = app.InnerPressureEditField.Value;
            p_out = app.OuterPressureEditField.Value;
        
            Db=((1-nu)*E/((1+nu)*(1-2*nu)))*[1, nu/(1-nu), 0; 
                                             nu/(1-nu), 1, 0; 
                                             0, 0, (1-2*nu)/(2*(1-nu))];
        
            nnr = app.NodesRadialEditField.Value;
            nnphi = app.NodesAngularEditField.Value;
            elemType = app.ElementTypeSwitch.Value;
        
            switch elemType
                case 'Q4'
                    nnode = 4;
                    quadorder = 3;
                case 'T3'
                    nnode = 3;
                    quadorder = 3;
            end
        
            [node, element,~,numElem] = mesh2D(elemType,r,R,nnr,nnphi);
            
            cla(app.UIAxes)
            hold(app.UIAxes, 'on')
            numelem = size(element,1);
            for e = 1:numElem
                nodes_e = element(e,:);
                coord = node(nodes_e,:);
                coord = [coord; coord(1,:)];
                patch(app.UIAxes, coord(:,1), coord(:,2), 'w', 'EdgeColor', 'k');
            end
    
            plot(app.UIAxes, node(:,1), node(:,2), 'r.')
            axis(app.UIAxes, 'equal')
            xlim(app.UIAxes, [0 R])
            ylim(app.UIAxes, [0 R])
    
            hold(app.UIAxes, 'off')
                title(app.UIAxes, ...
                'Finite Element Mesh Generation')
                app.StatusLabel.Text = 'Applied meshing plotted successfully';
                app.StatusLabel.FontColor = [0 0.6 0];
        end

        % Button down function: UIAxes
        function UIAxesButtonDown(app, event)
     
                node = app.node;
                element = app.element;
                numElem = app.numElem;
                R = app.R;
        
                cla(app.UIAxes)
                hold(app.UIAxes, 'on')
            
                for e = 1:numElem
                    nodes_e = element(e,:);
                    coord = node(nodes_e,:);
        
                    coord = [coord; coord(1,:)];
        
                    patch(app.UIAxes, coord(:,1), coord(:,2), ...
                        'w', 'EdgeColor', 'k');
                end
        
                plot(app.UIAxes, node(:,1), node(:,2), 'r.')
        
                axis(app.UIAxes, 'equal')
                xlim(app.UIAxes, [0 R])
                ylim(app.UIAxes, [0 R])
        
                hold(app.UIAxes, 'off')

        end

        % Button pushed function: APPLIEDFORCEButton
        function APPLIEDFORCEButtonPushed(app, event)
        %% Input
            E = app.YoungModulusEditField.Value;
            nu = app.PoisonEditField.Value;
            thickness = app.ThicknessEditField.Value;
            r = app.InnerRadiusrEditField.Value;
            R = app.OuterRadiusREditField.Value;
            p_in = app.InnerPressureEditField.Value;
            p_out = app.OuterPressureEditField.Value;
    
            nnr = app.NodesRadialEditField.Value;
            nnphi = app.NodesAngularEditField.Value;
            elemType = app.ElementTypeSwitch.Value;
    
            Db = ((1-nu)*E/((1+nu)*(1-2*nu))) * ...
                [1, nu/(1-nu), 0;
                 nu/(1-nu), 1, 0;
                 0, 0, (1-2*nu)/(2*(1-nu))];
    
            switch elemType
                case 'Q4'
                    nnode = 4;
                    quadorder = 3;
                case 'T3'
                    nnode = 3;
                    quadorder = 3;
            end
    
            %% Mesh + force
            [node, element] = mesh2D(elemType, r, R, nnr, nnphi);
            f = force2D(node, element, thickness, p_in, r, p_out, R);
    
            %% Plot on UIAxes
            cla(app.UIAxes)
            hold(app.UIAxes, 'on')
    
            patch(app.UIAxes, ...
                'Faces', element, ...
                'Vertices', node, ...
                'FaceColor', 'none', ...
                'EdgeColor', [0.8 0.8 0.8]);
    
            quiver(app.UIAxes, ...
                node(:,1), node(:,2), ...
                f(1:2:end), f(2:2:end), ...
                1.5, 'r');
    
            axis(app.UIAxes, 'equal')
            xlim(app.UIAxes, [0 R])
            ylim(app.UIAxes, [0 R])
    
            title(app.UIAxes, ...
                ['Force distribution diagram p_{in} = ', num2str(p_in),'N', ...
                 ' and p_{out} = ', num2str(p_out),'N'])
    
            hold(app.UIAxes, 'off')
    
            app.StatusLabel.Text = 'Applied force plotted successfully';
            app.StatusLabel.FontColor = [0 0.6 0];
        end

        % Value changed function: TextArea
        function TextAreaValueChanged(app, event)
            value = app.TextArea.Value;
            
        end

        % Clicked callback: TYPEDropDown
        function TYPEDropDownClicked(app, event)
                E = app.YoungModulusEditField.Value;
                nu = app.PoisonEditField.Value;
                thickness = app.ThicknessEditField.Value;
                r = app.InnerRadiusrEditField.Value;
                R = app.OuterRadiusREditField.Value;
                p_in = app.InnerPressureEditField.Value;
                p_out = app.OuterPressureEditField.Value;
        
                Db = ((1-nu)*E/((1+nu)*(1-2*nu))) * ...
                    [1, nu/(1-nu), 0;
                     nu/(1-nu), 1, 0;
                     0, 0, (1-2*nu)/(2*(1-nu))];
        
                nnr = app.NodesRadialEditField.Value;
                nnphi = app.NodesAngularEditField.Value;
                elemType = app.ElementTypeSwitch.Value;
        
                switch elemType
                    case 'Q4'
                        nnode = 4;
                        quadorder = 3;
                    case 'T3'
                        nnode = 3;
                        quadorder = 3;
                end
        
                %% MESH
                [node, element, numnode, numElem] = mesh2D(elemType, r, R, nnr, nnphi);
        
                %% STIFFNESS
                ndof = 2;
                tdof = ndof * numnode;
                K = stiffnessmat2D(elemType, quadorder, numElem, numnode, node, element, Db, thickness);
        
                %% FORCE
                f = force2D(node, element, thickness, p_in, r, p_out, R);
        
                %% BC
                bottom = find(abs(node(:,2)) < 1e-6);
                left = find(abs(node(:,1)) < 1e-6);
        
                DOFRoller = unique([bottom*2; left*2-1]);
        
                sdof = (1:tdof)';
                UnknownDof = sdof;
                UnknownDof(DOFRoller) = [];
        
                nUnknownDof = length(UnknownDof);
        
                Tbd = zeros(tdof, nUnknownDof);
                for i = 1:nUnknownDof
                    Tbd(UnknownDof(i), i) = 1;
                end
        
                Kred = Tbd' * K * Tbd;
                fred = Tbd' * f;
        
                %% SOLVE
                displacement = Kred \ fred;
                u = Tbd * displacement;
        
                %% STRESS
                
                ux = u(1:2:end); 
                uy = u(2:2:end);
                x = node(:,1);
                y = node(:,2);
                
                theta = atan2(y, x);
                
                ur     = ux .* cos(theta) + uy .* sin(theta);
                utheta = -ux .* sin(theta) + uy .* cos(theta);
              
                [vm, vr, vh, sz, sx] = ...
                    stress_cal(elemType, quadorder, numElem, node, element, Db, nu, u);
                switch elemType
                    case 'T3'
                        tri = element;
                
                    case 'Q4'
                        % Mỗi Q4 tách thành 2 tam giác:
                        % [1 2 3] và [1 3 4]
                        tri = [element(:, [1 2 3]);
                               element(:, [1 3 4])];
                
                end    
                % ====== PLOT TRÊN UIAxes ======
                        
                % Chọn dữ liệu cần vẽ
                switch app.TYPEDropDown.Value
                    case 'None'
                app.StatusLabel.Text = 'Ban chua chon loai phan tu';
                app.StatusLabel.FontColor = [1 0 0];    
                    case 'σ_x (Pa)'
                        cla(app.UIAxes)
                        hold(app.UIAxes, 'on')
                        sigma_x_node = zeros(numnode,1);
                        node_count = zeros(numnode,1);
                        
                        for iel = 1:numElem
                            nodes_iel = element(iel,:);
                            sigma_x_node(nodes_iel) = sigma_x_node(nodes_iel) + sx(iel);
                            node_count(nodes_iel) = node_count(nodes_iel) + 1;
                        end
                        switch elemType
                            case 'Q4'
                                patch(app.UIAxes,'Faces', element, ...
                                              'Vertices', node, ...
                                              'FaceVertexCData',  sigma_x_node, ...
                                              'FaceColor', 'int', ...
                                              'EdgeColor', [0.3 0.3 0.3], ...
                                              'LineWidth', 0.5);
                            case 'T3'
                                patch(app.UIAxes,'Faces', element, ...
                                              'Vertices', node, ...
                                              'FaceVertexCData',  sigma_x_node, ...
                                              'FaceColor', 'flat', ...
                                              'EdgeColor', [0.3 0.3 0.3], ...
                                              'LineWidth', 0.5);
                        end
                                  
                         % =========================================================
                        % 3. TÌM VÀ GẮN TAG ĐIỂM MAX, MIN
                        % =========================================================
                        
                        % Tìm giá trị và vị trí (index) của điểm Max, Min
                        [max_val, max_idx] = max(sigma_x_node);
                        [min_val, min_idx] = min(sigma_x_node);
                        
                        % Lấy tọa độ (x, y) của nút Max, Min từ ma trận 'node'
                        x_max = node(max_idx, 1);
                        y_max = node(max_idx, 2);
                        
                        x_min = node(min_idx, 1);
                        y_min = node(min_idx, 2);
                        % --- ĐIỂM MAX (TAM GIÁC ĐỎ + HỘP CHỮ ĐỎ) ---
                        plot(app.UIAxes, x_max, y_max, '^', 'MarkerSize', 10, ...
                             'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r'); % Tam giác: Viền đen, Nền đỏ
                        
                        text(app.UIAxes, x_max, y_max, sprintf('Max: %.2e ', max_val), ...
                             'Color', 'w', ...              % Chữ màu trắng (để nổi trên nền đỏ)
                             'FontWeight', 'bold', ...
                             'FontSize', 10, ...
                             'BackgroundColor', 'r', ...    % NỀN ĐỎ (Khớp với tam giác)
                             'EdgeColor', 'k', ...          % VIỀN ĐEN (Khớp với tam giác)
                             'LineWidth', 0.1, ...            % Độ dày viền đen
                             'Margin', 0.4, ...               % Khoảng cách giữa chữ và viền
                             'VerticalAlignment', 'baseline'); % Đẩy chữ lên trên một chút để không đè vào tam giác
                        
                        % --- ĐIỂM MIN (TAM GIÁC XANH + HỘP CHỮ XANH) ---
                        plot(app.UIAxes, x_min, y_min, 'v', 'MarkerSize', 10, ...
                             'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b'); % Tam giác: Viền đen, Nền xanh
                        
                        text(app.UIAxes, x_min, y_min, sprintf('Min: %.2e ', min_val), ...
                             'Color', 'w', ...              % Chữ màu trắng
                             'FontWeight', 'bold', ...
                             'FontSize', 10, ...
                             'BackgroundColor', 'b', ...    % NỀN XANH (Khớp với tam giác)
                             'EdgeColor', 'k', ...          % VIỀN ĐEN (Khớp với tam giác)
                             'LineWidth', 0.1, ...
                             'Margin', 0.4, ...
                             'VerticalAlignment', 'cap');             
                        colorbar(app.UIAxes)
                        title(app.UIAxes, '')
                        xlabel(app.UIAxes, 'x')
                        ylabel(app.UIAxes, 'y')
                        axis(app.UIAxes, 'equal')
                        colormap(app.UIAxes, 'jet')
                        view(app.UIAxes, 2)
                        hold(app.UIAxes, 'off')
                    case 'σ_y (Pa)'
                        cla(app.UIAxes)
                        hold(app.UIAxes, 'on')
                        sigma_y_node = zeros(numnode,1);
                        node_count = zeros(numnode,1);
                        
                        for iel = 1:numElem
                            nodes_iel = element(iel,:);
                            sigma_y_node(nodes_iel) = sigma_y_node(nodes_iel) + sigma_y_elem(iel);
                            node_count(nodes_iel) = node_count(nodes_iel) + 1;
                        end
                        switch elemType
                            case 'Q4'
                                patch(app.UIAxes,'Faces', element, ...
                                              'Vertices', node, ...
                                              'FaceVertexCData',  sigma_y_node, ...
                                              'FaceColor', 'int', ...
                                              'EdgeColor', [0.3 0.3 0.3], ...
                                              'LineWidth', 0.5);
                            case 'T3'
                                patch(app.UIAxes,'Faces', element, ...
                                              'Vertices', node, ...
                                              'FaceVertexCData',  sigma_y_node, ...
                                              'FaceColor', 'flat', ...
                                              'EdgeColor', [0.3 0.3 0.3], ...
                                              'LineWidth', 0.5);
                        end
                                  
                         % =========================================================
                        % 3. TÌM VÀ GẮN TAG ĐIỂM MAX, MIN
                        % =========================================================
                        
                        % Tìm giá trị và vị trí (index) của điểm Max, Min
                        [max_val, max_idx] = max(sigma_y_node);
                        [min_val, min_idx] = min(sigma_y_node);
                        
                        % Lấy tọa độ (x, y) của nút Max, Min từ ma trận 'node'
                        x_max = node(max_idx, 1);
                        y_max = node(max_idx, 2);
                        
                        x_min = node(min_idx, 1);
                        y_min = node(min_idx, 2);
                        % --- ĐIỂM MAX (TAM GIÁC ĐỎ + HỘP CHỮ ĐỎ) ---
                        plot(app.UIAxes, x_max, y_max, '^', 'MarkerSize', 10, ...
                             'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r'); % Tam giác: Viền đen, Nền đỏ
                        
                        text(app.UIAxes, x_max, y_max, sprintf('Max: %.2e ', max_val), ...
                             'Color', 'w', ...              % Chữ màu trắng (để nổi trên nền đỏ)
                             'FontWeight', 'bold', ...
                             'FontSize', 10, ...
                             'BackgroundColor', 'r', ...    % NỀN ĐỎ (Khớp với tam giác)
                             'EdgeColor', 'k', ...          % VIỀN ĐEN (Khớp với tam giác)
                             'LineWidth', 0.1, ...            % Độ dày viền đen
                             'Margin', 0.4, ...               % Khoảng cách giữa chữ và viền
                             'VerticalAlignment', 'baseline'); % Đẩy chữ lên trên một chút để không đè vào tam giác
                        
                        % --- ĐIỂM MIN (TAM GIÁC XANH + HỘP CHỮ XANH) ---
                        plot(app.UIAxes, x_min, y_min, 'v', 'MarkerSize', 10, ...
                             'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b'); % Tam giác: Viền đen, Nền xanh
                        
                        text(app.UIAxes, x_min, y_min, sprintf('Min: %.2e ', min_val), ...
                             'Color', 'w', ...              % Chữ màu trắng
                             'FontWeight', 'bold', ...
                             'FontSize', 10, ...
                             'BackgroundColor', 'b', ...    % NỀN XANH (Khớp với tam giác)
                             'EdgeColor', 'k', ...          % VIỀN ĐEN (Khớp với tam giác)
                             'LineWidth', 0.1, ...
                             'Margin', 0.4, ...
                             'VerticalAlignment', 'cap');             
                        colorbar(app.UIAxes)
                        title(app.UIAxes, '')
                        xlabel(app.UIAxes, 'x')
                        ylabel(app.UIAxes, 'y')
                        axis(app.UIAxes, 'equal')
                        colormap(app.UIAxes, 'jet')
                        view(app.UIAxes, 2)
                        hold(app.UIAxes, 'off')
                    case 'σ_r (Pa)'
                        cla(app.UIAxes)
                        hold(app.UIAxes, 'on')
                        sigma_r_node = zeros(numnode,1);
                        node_count = zeros(numnode,1);
                        
                        for iel = 1:numElem
                            nodes_iel = element(iel,:);
                            sigma_r_node(nodes_iel) = sigma_r_node(nodes_iel) + sigma_r(iel);
                            node_count(nodes_iel) = node_count(nodes_iel) + 1;
                        end
                        
                        sigma_r_node = sigma_r_node ./ node_count;
                        switch elemType
                            case 'Q4'
                                patch(app.UIAxes,'Faces', element, ...
                                              'Vertices', node, ...
                                              'FaceVertexCData',  sigma_r_node, ...
                                              'FaceColor', 'int', ...
                                              'EdgeColor', [0.3 0.3 0.3], ...
                                              'LineWidth', 0.5);
                            case 'T3'
                                patch(app.UIAxes,'Faces', element, ...
                                              'Vertices', node, ...
                                              'FaceVertexCData',  sigma_r_node, ...
                                              'FaceColor', 'flat', ...
                                              'EdgeColor', [0.3 0.3 0.3], ...
                                              'LineWidth', 0.5);
                        end
                                  
                        colorbar(app.UIAxes)
                        title(app.UIAxes, '')
                        xlabel(app.UIAxes, 'x')
                        ylabel(app.UIAxes, 'y')
                        axis(app.UIAxes, 'equal')
                        colormap(app.UIAxes, 'jet')
                        view(app.UIAxes, 2)
                        hold(app.UIAxes, 'off')
                    case 'σ_theta (Pa)'
                        cla(app.UIAxes)
                        hold(app.UIAxes, 'on')
                        sigma_theta_node = zeros(numnode,1);
                        node_count = zeros(numnode,1);
                        
                        for iel = 1:numElem
                            nodes_iel = element(iel,:);
                            sigma_theta_node(nodes_iel) = sigma_theta_node(nodes_iel) + sigma_theta(iel);
                            node_count(nodes_iel) = node_count(nodes_iel) + 1;
                        end
                        
                        sigma_theta_node = sigma_theta_node ./ node_count;
                        switch elemType
                            case 'Q4'
                                patch(app.UIAxes,'Faces', element, ...
                                              'Vertices', node, ...
                                              'FaceVertexCData',  sigma_theta_node, ...
                                              'FaceColor', 'int', ...
                                              'EdgeColor', [0.3 0.3 0.3], ...
                                              'LineWidth', 0.5);
                            case 'T3'
                                patch(app.UIAxes,'Faces', element, ...
                                              'Vertices', node, ...
                                              'FaceVertexCData',  sigma_theta_node, ...
                                              'FaceColor', 'flat', ...
                                              'EdgeColor', [0.3 0.3 0.3], ...
                                              'LineWidth', 0.5);
                        end
                                  
                            
                        colorbar(app.UIAxes)
                        title(app.UIAxes, '')
                        xlabel(app.UIAxes, 'x')
                        ylabel(app.UIAxes, 'y')
                        axis(app.UIAxes, 'equal')
                        colormap(app.UIAxes, 'jet')
                        view(app.UIAxes, 2)
                        hold(app.UIAxes, 'off')
                    case 'τ_xy (Pa)'
                        cla(app.UIAxes)
                        hold(app.UIAxes, 'on')
                        tau_xy_node = zeros(numnode,1);
                        node_count = zeros(numnode,1);
                        
                        for iel = 1:numElem
                            nodes_iel = element(iel,:);
                            tau_xy_node(nodes_iel) = tau_xy_node(nodes_iel) + tau_xy_elem(iel);
                            node_count(nodes_iel) = node_count(nodes_iel) + 1;
                        end
                        
                        tau_xy_node = tau_xy_node ./ node_count;
                        switch elemType
                            case 'Q4'
                                patch(app.UIAxes,'Faces', element, ...
                                              'Vertices', node, ...
                                              'FaceVertexCData',  tau_xy_node, ...
                                              'FaceColor', 'int', ...
                                              'EdgeColor', [0.3 0.3 0.3], ...
                                              'LineWidth', 0.5);
                            case 'T3'
                                patch(app.UIAxes,'Faces', element, ...
                                              'Vertices', node, ...
                                              'FaceVertexCData',  tau_xy_node, ...
                                              'FaceColor', 'flat', ...
                                              'EdgeColor', [0.3 0.3 0.3], ...
                                              'LineWidth', 0.5);
                        end
                        colorbar(app.UIAxes)
                        title(app.UIAxes, '')
                        xlabel(app.UIAxes, 'x')
                        ylabel(app.UIAxes, 'y')
                        axis(app.UIAxes, 'equal')
                        colormap(app.UIAxes, 'jet')
                        view(app.UIAxes, 2)
                        hold(app.UIAxes, 'off')
                    case 'Von Mises (Pa)'
                        cla(app.UIAxes)
                        hold(app.UIAxes, 'on')
                        % NODAL AVERAGING FOR VON MISES
                        vm_node = zeros(numnode,1);
                        node_count = zeros(numnode,1);
                        
                        for iel = 1:numElem
                            nodes_iel = element(iel,:);
                            vm_node(nodes_iel) = vm_node(nodes_iel) + vm(iel);
                            node_count(nodes_iel) = node_count(nodes_iel) + 1;
                        end
                        
                        vm_node = vm_node ./ node_count;
                        
                        switch elemType
                            case 'Q4'
                                patch(app.UIAxes,'Faces', element, ...
                                              'Vertices', node, ...
                                              'FaceVertexCData',  vm_node, ...
                                              'FaceColor', 'int', ...
                                              'EdgeColor', [0.3 0.3 0.3], ...
                                              'LineWidth', 0.5);
                            case 'T3'
                                patch(app.UIAxes,'Faces', element, ...
                                              'Vertices', node, ...
                                              'FaceVertexCData',  vm_node, ...
                                              'FaceColor', 'flat', ...
                                              'EdgeColor', [0.3 0.3 0.3], ...
                                              'LineWidth', 0.5);
                        end
                            
                        colorbar(app.UIAxes)
                        title(app.UIAxes, '')
                        xlabel(app.UIAxes, 'x')
                        ylabel(app.UIAxes, 'y')
                        axis(app.UIAxes, 'equal')
                        colormap(app.UIAxes, 'jet')
                        view(app.UIAxes, 2)
                        hold(app.UIAxes, 'off')
                    case 'σ_z (Pa)'
                        cla(app.UIAxes)
                        hold(app.UIAxes, 'on')
                        %NODAL AVERAGING FOR VON MISES
                        sigma_z_node = zeros(numnode,1);
                        node_count = zeros(numnode,1);
                        
                        for iel = 1:numElem
                            nodes_iel = element(iel,:);
                            sigma_z_node(nodes_iel) = sigma_z_node(nodes_iel) + sigma_z(iel);
                            node_count(nodes_iel) = node_count(nodes_iel) + 1;
                        end
                        
                        sigma_z_node = sigma_z_node ./ node_count;
                        round_step = 1e6;
                        sigma_z_round = round(sigma_z_node / round_step)* round_step;
                        switch elemType
                            case 'Q4'
                                patch(app.UIAxes,'Faces', element, ...
                                              'Vertices', node, ...
                                              'FaceVertexCData',  sigma_z_round, ...
                                              'FaceColor', 'int', ...
                                              'EdgeColor', [0.3 0.3 0.3], ...
                                              'LineWidth', 0.5);
                            case 'T3'
                                patch(app.UIAxes,'Faces', element, ...
                                              'Vertices', node, ...
                                              'FaceVertexCData',  sigma_z_round, ...
                                              'FaceColor', 'flat', ...
                                              'EdgeColor', [0.3 0.3 0.3], ...
                                              'LineWidth', 0.5);
                        end
                            
                            
                        colorbar(app.UIAxes)
                        title(app.UIAxes, '')
                        xlabel(app.UIAxes, 'x')
                        ylabel(app.UIAxes, 'y')
                        axis(app.UIAxes, 'equal')
                        colormap(app.UIAxes, 'jet')
                        view(app.UIAxes, 2)
                        hold(app.UIAxes, 'off')
                    case 'u_x (m)'
                        cla(app.UIAxes)
                        hold(app.UIAxes, 'on')
                        % Sử dụng patch để vẽ và blend màu
                        hp = patch(app.UIAxes,'Faces', element, ...
                                   'Vertices', node, ...
                                   'FaceVertexCData', ux, ...    % Gán giá trị ux cho từng nút
                                   'FaceColor', 'interp', ...    % Màu sẽ blend giữa các nút
                                   'EdgeColor', [0.3 0.3 0.3], ... % Hiển thị lưới mờ
                                   'LineWidth', 0.5);
                            
                        colorbar(app.UIAxes)
                        title(app.UIAxes, '')
                        xlabel(app.UIAxes, 'x')
                        ylabel(app.UIAxes, 'y')
                        axis(app.UIAxes, 'equal')
                        colormap(app.UIAxes, 'jet')
                        view(app.UIAxes, 2)
                        hold(app.UIAxes, 'off')

                    case 'u_y (m)'
                        cla(app.UIAxes)
                        hold(app.UIAxes, 'on')
                        % Sử dụng patch để vẽ và blend màu
                        hp = patch(app.UIAxes,'Faces', element, ...
                                   'Vertices', node, ...
                                   'FaceVertexCData', uy, ...    % Gán giá trị ux cho từng nút
                                   'FaceColor', 'interp', ...    % Màu sẽ blend giữa các nút
                                   'EdgeColor', [0.3 0.3 0.3], ... % Hiển thị lưới mờ
                                   'LineWidth', 0.5);
                            
                        colorbar(app.UIAxes)
                        title(app.UIAxes, '')
                        xlabel(app.UIAxes, 'x')
                        ylabel(app.UIAxes, 'y')
                        axis(app.UIAxes, 'equal')
                        colormap(app.UIAxes, 'jet')
                        view(app.UIAxes, 2)
                        hold(app.UIAxes, 'off')
                    case 'u_r (m)'
                        cla(app.UIAxes)
                        hold(app.UIAxes, 'on')
                        % Sử dụng patch để vẽ và blend màu
                        hp = patch(app.UIAxes,'Faces', element, ...
                                   'Vertices', node, ...
                                   'FaceVertexCData', ur, ...    % Gán giá trị ux cho từng nút
                                   'FaceColor', 'interp', ...    % Màu sẽ blend giữa các nút
                                   'EdgeColor', [0.3 0.3 0.3], ... % Hiển thị lưới mờ
                                   'LineWidth', 0.5);
                            
                        colorbar(app.UIAxes)
                        title(app.UIAxes, '')
                        xlabel(app.UIAxes, 'x')
                        ylabel(app.UIAxes, 'y')
                        axis(app.UIAxes, 'equal')
                        colormap(app.UIAxes, 'jet')
                        view(app.UIAxes, 2)
                        hold(app.UIAxes, 'off')

                    case 'u_theta (m)'
                        cla(app.UIAxes)
                        hold(app.UIAxes, 'on')
                        % Sử dụng patch để vẽ và blend màu
                        hp = patch(app.UIAxes,'Faces', element, ...
                                   'Vertices', node, ...
                                   'FaceVertexCData', utheta, ...    % Gán giá trị ux cho từng nút
                                   'FaceColor', 'interp', ...    % Màu sẽ blend giữa các nút
                                   'EdgeColor', [0.3 0.3 0.3], ... % Hiển thị lưới mờ
                                   'LineWidth', 0.5);
                            
                        colorbar(app.UIAxes)
                        title(app.UIAxes, '')
                        xlabel(app.UIAxes, 'x')
                        ylabel(app.UIAxes, 'y')
                        axis(app.UIAxes, 'equal')
                        colormap(app.UIAxes, 'jet')
                        view(app.UIAxes, 2)
                        hold(app.UIAxes, 'off')

 
                    otherwise
                        patch(app.UIAxes, ...
                            'Faces', element, ...
                            'Vertices', node, ...
                            'FaceColor', 'interp', ...
                            'EdgeColor', [0.3 0.3 0.3], ...
                            'LineWidth', 0.5);
                    
                        colorbar(app.UIAxes)
                        title(app.UIAxes, '')
                        xlabel(app.UIAxes, 'x')
                        ylabel(app.UIAxes, 'y')
                        axis(app.UIAxes, 'equal')
                        colormap(app.UIAxes, 'jet')
                        view(app.UIAxes, 2)
                        hold(app.UIAxes, 'off')
                end
            

            
        end

        % Value changed function: TYPEDropDown
        function TYPEDropDownValueChanged(app, event)
            TYPEDropDownClicked(app, event);
        end

        % Button pushed function: RESETButton
        function RESETButtonPushed(app, event)
                    
            % ===== 1. Xóa plot trên UIAxes =====
            cla(app.UIAxes)                 % clear axes
            colorbar(app.UIAxes, 'off')     % tắt colorbar nếu có
            title(app.UIAxes, '')           % xóa title
            xlabel(app.UIAxes, '')
            ylabel(app.UIAxes, '')
        
            % ===== 2. Xóa Label =====
            app.StatusLabel.Text = '';
            % nếu có nhiều label:
            % app.Label2.Text = '';
            % app.Label3.Text = '';
        
            % ===== 3. Xóa TextArea =====
            app.TextArea.Value = {};
        
            % ===== 4. Reset DropDown (nếu muốn) =====
            app.TYPEDropDown.Value = 'None';
        
            % ===== 5. Reset dữ liệu lưu trong app =====
            app.node = [];
            app.element = [];
            app.sigma_x = [];

        end

        % Callback function: MODALButton, MODESpinner
        function MODALButtonPushed(app, event)
            %% input
            E = app.YoungModulusEditField.Value;
            nu = app.PoisonEditField.Value;
            thickness = app.ThicknessEditField.Value;
            r = app.InnerRadiusrEditField.Value;
            R = app.OuterRadiusREditField.Value;
            p_in = app.InnerPressureEditField.Value;
            p_out = app.OuterPressureEditField.Value;
        
            Db=((1-nu)*E/((1+nu)*(1-2*nu)))*[1, nu/(1-nu), 0; 
                                             nu/(1-nu), 1, 0; 
                                             0, 0, (1-2*nu)/(2*(1-nu))];
        
            [node, element, numnode, numElem] = mesh2D('Q4', r, R, 20, 20);
            [K, M] = stiffnessmat2D('Q4', 3, numElem, numnode, node, element, Db, thickness);
            
            numModes = 6;
            [V, D] = eigs(K, M, numModes, 'smallestabs');
            f_hz = sqrt(diag(D)) / (2*pi);
            
            % Giả sử ta đang xem Mode 1
            modeID = round(app.MODESpinner.Value);
            mode_shape = V(:, modeID);
% Tách vector riêng thành chuyển vị x và y cho từng node
        ux_mode = mode_shape(1:2:end);
        uy_mode = mode_shape(2:2:end);
        
        % Tính biên độ biến dạng tổng hợp của từng node
        u_mode_mag = sqrt(ux_mode.^2 + uy_mode.^2);
        
        % =====================================================================
        % ÁP DỤNG THUẬT TOÁN TỈ SỐ BIÊN ĐỘ (TƯƠNG ĐƯƠNG VECTO RIÊNG r)
        % =====================================================================
        % Tìm chuyển vị lớn nhất trong toàn hệ
        max_disp = max(u_mode_mag);
        
        % Chuẩn hóa (Normalize) vector riêng. Lúc này giá trị lớn nhất sẽ là 1 
        % (Tương đương việc đặt X1 = 1 trong hình của bạn)
        ux_normalized = ux_mode / max_disp;
        uy_normalized = uy_mode / max_disp;
        
        % Mảng này chính là tập hợp các tỉ số r của mọi node trên hình so với node max
        r_ratio_array = u_mode_mag / max_disp; 
        
        % Hệ số phóng đại để dễ quan sát (10% bán kính ngoài)
        scaleFactor = 0.1 * R; 

        % =====================================================================
        % BƯỚC 2: KHỞI TẠO ĐỒ HỌA TRÊN APP.UIAXES (Giữ nguyên của bạn)
        % =====================================================================
        cla(app.UIAxes); 
        hold(app.UIAxes, 'on'); 
        axis(app.UIAxes, 'equal'); 
        grid(app.UIAxes, 'on');
        colormap(app.UIAxes, 'jet'); 
    
        hPatch = patch(app.UIAxes, 'Faces', element, ...
                       'Vertices', node, ...
                       'FaceVertexCData', zeros(size(node,1), 1), ... 
                       'FaceColor', 'interp', ...         
                       'EdgeColor', 'k', ... 
                       'EdgeAlpha', 0.3); 
    
        colorbar(app.UIAxes);
        title(app.UIAxes, ['Mode ', num2str(modeID), ' (f = ', num2str(f_hz(modeID), '%.2f'), ' Hz)']);
    
        margin = 0.2 * R;
        xlim(app.UIAxes, [min(node(:,1))-margin, max(node(:,1))+margin]);
        ylim(app.UIAxes, [min(node(:,2))-margin, max(node(:,2))+margin]);

        % =====================================================================
        % BƯỚC 3: THUẬT TOÁN HÌNH ĐỘNG DỰA TRÊN R_RATIO
        % =====================================================================
        nFrames = 60; 
        t_steps = linspace(0, 2*pi, nFrames);
        
        % Vì r_ratio_array đã được chuẩn hóa từ 0 đến 1, ta set clim cố định
        cmin = 0; 
        cmax = 1;
        clim(app.UIAxes, [cmin, cmax]); 
    
        for cycle = 1:10 
            for t = t_steps              
                % Hàm điều hòa
                amp = sin(t);
                
                % Tính tọa độ mới dựa trên vector đã chuẩn hóa r_ratio
                curr_ux = ux_normalized * scaleFactor * amp;
                curr_uy = uy_normalized * scaleFactor * amp;
                new_vertices = [node(:,1) + curr_ux, node(:,2) + curr_uy];
                
                % Màu sắc dựa trên độ lớn của r_ratio tại thời điểm t
                % Dùng hàm abs() vì màu sắc chỉ thể hiện độ lớn ứng suất/biến dạng
                cdata = abs(r_ratio_array * amp);
                
                % Tối ưu hóa: Để MATLAB tự map màu thông qua FaceVertexCData thay vì tự tính RGB
                set(hPatch, 'Vertices', new_vertices, 'FaceVertexCData', cdata);
                
                pause(0.02); 
                drawnow limitrate;
            end
        end
        end

        % Button down function: UIAxes_2
        function UIAxes_2ButtonDown(app, event)
            node = app.node;
            element = app.element;
            numElem = app.numElem;
            R = app.R;
    
            cla(app.UIAxes)
            hold(app.UIAxes, 'on')
        
            for e = 1:numElem
                nodes_e = element(e,:);
                coord = node(nodes_e,:);
    
                coord = [coord; coord(1,:)];
    
                patch(app.UIAxes, coord(:,1), coord(:,2), ...
                    'w', 'EdgeColor', 'k');
            end
    
            plot(app.UIAxes, node(:,1), node(:,2), 'r.')
    
            axis(app.UIAxes, 'equal')
            xlim(app.UIAxes, [0 R])
            ylim(app.UIAxes, [0 R])
    
            hold(app.UIAxes, 'off')

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [1 1 1];
            colormap(app.UIFigure, 'jet');
            app.UIFigure.Position = [100 100 1629 681];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.Pointer = 'hand';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.ButtonDownFcn = createCallbackFcn(app, @UIAxesButtonDown, true);
            app.UIAxes.Position = [296 180 449 358];

            % Create UIAxes_2
            app.UIAxes_2 = uiaxes(app.UIFigure);
            title(app.UIAxes_2, 'Title')
            xlabel(app.UIAxes_2, 'X')
            ylabel(app.UIAxes_2, 'Y')
            zlabel(app.UIAxes_2, 'Z')
            app.UIAxes_2.ButtonDownFcn = createCallbackFcn(app, @UIAxes_2ButtonDown, true);
            app.UIAxes_2.Position = [978 196 449 358];

            % Create YoungModulusEditFieldLabel
            app.YoungModulusEditFieldLabel = uilabel(app.UIFigure);
            app.YoungModulusEditFieldLabel.HorizontalAlignment = 'right';
            app.YoungModulusEditFieldLabel.Position = [50 598 87 22];
            app.YoungModulusEditFieldLabel.Text = 'Young Modulus';

            % Create YoungModulusEditField
            app.YoungModulusEditField = uieditfield(app.UIFigure, 'numeric');
            app.YoungModulusEditField.ValueChangedFcn = createCallbackFcn(app, @YoungModulusEditFieldValueChanged, true);
            app.YoungModulusEditField.Position = [158 598 105 22];
            app.YoungModulusEditField.Value = 200000000000;

            % Create PoisonEditFieldLabel
            app.PoisonEditFieldLabel = uilabel(app.UIFigure);
            app.PoisonEditFieldLabel.HorizontalAlignment = 'right';
            app.PoisonEditFieldLabel.Position = [50 557 42 22];
            app.PoisonEditFieldLabel.Text = 'Poison';

            % Create PoisonEditField
            app.PoisonEditField = uieditfield(app.UIFigure, 'numeric');
            app.PoisonEditField.ValueChangedFcn = createCallbackFcn(app, @PoisonEditFieldValueChanged, true);
            app.PoisonEditField.Position = [158 557 105 22];
            app.PoisonEditField.Value = 0.3;

            % Create ThicknessEditFieldLabel
            app.ThicknessEditFieldLabel = uilabel(app.UIFigure);
            app.ThicknessEditFieldLabel.HorizontalAlignment = 'right';
            app.ThicknessEditFieldLabel.Position = [50 516 59 22];
            app.ThicknessEditFieldLabel.Text = 'Thickness';

            % Create ThicknessEditField
            app.ThicknessEditField = uieditfield(app.UIFigure, 'numeric');
            app.ThicknessEditField.ValueChangedFcn = createCallbackFcn(app, @ThicknessEditFieldValueChanged, true);
            app.ThicknessEditField.Position = [158 516 105 22];
            app.ThicknessEditField.Value = 1;

            % Create InnerRadiusrEditFieldLabel
            app.InnerRadiusrEditFieldLabel = uilabel(app.UIFigure);
            app.InnerRadiusrEditFieldLabel.HorizontalAlignment = 'right';
            app.InnerRadiusrEditFieldLabel.Position = [50 475 80 22];
            app.InnerRadiusrEditFieldLabel.Text = 'Inner Radius r';

            % Create InnerRadiusrEditField
            app.InnerRadiusrEditField = uieditfield(app.UIFigure, 'numeric');
            app.InnerRadiusrEditField.ValueChangedFcn = createCallbackFcn(app, @InnerRadiusrEditFieldValueChanged, true);
            app.InnerRadiusrEditField.Position = [158 475 105 22];
            app.InnerRadiusrEditField.Value = 0.5;

            % Create OuterRadiusREditFieldLabel
            app.OuterRadiusREditFieldLabel = uilabel(app.UIFigure);
            app.OuterRadiusREditFieldLabel.HorizontalAlignment = 'right';
            app.OuterRadiusREditFieldLabel.Position = [50 434 88 22];
            app.OuterRadiusREditFieldLabel.Text = 'Outer Radius R';

            % Create OuterRadiusREditField
            app.OuterRadiusREditField = uieditfield(app.UIFigure, 'numeric');
            app.OuterRadiusREditField.ValueChangedFcn = createCallbackFcn(app, @OuterRadiusREditFieldValueChanged, true);
            app.OuterRadiusREditField.Position = [158 434 105 22];
            app.OuterRadiusREditField.Value = 1;

            % Create InnerPressureEditFieldLabel
            app.InnerPressureEditFieldLabel = uilabel(app.UIFigure);
            app.InnerPressureEditFieldLabel.HorizontalAlignment = 'right';
            app.InnerPressureEditFieldLabel.Position = [50 393 84 22];
            app.InnerPressureEditFieldLabel.Text = 'Inner Pressure';

            % Create InnerPressureEditField
            app.InnerPressureEditField = uieditfield(app.UIFigure, 'numeric');
            app.InnerPressureEditField.ValueChangedFcn = createCallbackFcn(app, @InnerPressureEditFieldValueChanged, true);
            app.InnerPressureEditField.Position = [158 393 105 22];
            app.InnerPressureEditField.Value = 20000000;

            % Create OuterPressureEditFieldLabel
            app.OuterPressureEditFieldLabel = uilabel(app.UIFigure);
            app.OuterPressureEditFieldLabel.HorizontalAlignment = 'right';
            app.OuterPressureEditFieldLabel.Position = [50 352 86 22];
            app.OuterPressureEditFieldLabel.Text = 'Outer Pressure';

            % Create OuterPressureEditField
            app.OuterPressureEditField = uieditfield(app.UIFigure, 'numeric');
            app.OuterPressureEditField.ValueChangedFcn = createCallbackFcn(app, @OuterPressureEditFieldValueChanged, true);
            app.OuterPressureEditField.Position = [158 352 105 22];
            app.OuterPressureEditField.Value = -10000000;

            % Create NodesRadialEditFieldLabel
            app.NodesRadialEditFieldLabel = uilabel(app.UIFigure);
            app.NodesRadialEditFieldLabel.HorizontalAlignment = 'right';
            app.NodesRadialEditFieldLabel.Position = [50 311 77 22];
            app.NodesRadialEditFieldLabel.Text = 'Nodes Radial';

            % Create NodesRadialEditField
            app.NodesRadialEditField = uieditfield(app.UIFigure, 'numeric');
            app.NodesRadialEditField.ValueChangedFcn = createCallbackFcn(app, @NodesRadialEditFieldValueChanged, true);
            app.NodesRadialEditField.Position = [158 311 105 22];
            app.NodesRadialEditField.Value = 11;

            % Create NodesAngularEditFieldLabel
            app.NodesAngularEditFieldLabel = uilabel(app.UIFigure);
            app.NodesAngularEditFieldLabel.HorizontalAlignment = 'right';
            app.NodesAngularEditFieldLabel.Position = [50 270 84 22];
            app.NodesAngularEditFieldLabel.Text = 'Nodes Angular';

            % Create NodesAngularEditField
            app.NodesAngularEditField = uieditfield(app.UIFigure, 'numeric');
            app.NodesAngularEditField.ValueChangedFcn = createCallbackFcn(app, @NodesAngularEditFieldValueChanged, true);
            app.NodesAngularEditField.Position = [158 270 105 22];
            app.NodesAngularEditField.Value = 21;

            % Create ElementTypeSwitchLabel
            app.ElementTypeSwitchLabel = uilabel(app.UIFigure);
            app.ElementTypeSwitchLabel.HorizontalAlignment = 'center';
            app.ElementTypeSwitchLabel.Position = [118 206 78 22];
            app.ElementTypeSwitchLabel.Text = 'Element Type';

            % Create ElementTypeSwitch
            app.ElementTypeSwitch = uiswitch(app.UIFigure, 'slider');
            app.ElementTypeSwitch.Items = {'T3', 'Q4'};
            app.ElementTypeSwitch.ValueChangedFcn = createCallbackFcn(app, @ElementTypeSwitchValueChanged, true);
            app.ElementTypeSwitch.Position = [134 233 45 20];
            app.ElementTypeSwitch.Value = 'T3';

            % Create RUNANALYSISButton
            app.RUNANALYSISButton = uibutton(app.UIFigure, 'push');
            app.RUNANALYSISButton.ButtonPushedFcn = createCallbackFcn(app, @RUNANALYSISButtonPushed, true);
            app.RUNANALYSISButton.BackgroundColor = [0.1608 0.5804 0.1961];
            app.RUNANALYSISButton.FontSize = 20;
            app.RUNANALYSISButton.FontWeight = 'bold';
            app.RUNANALYSISButton.FontColor = [1 1 1];
            app.RUNANALYSISButton.Position = [77 158 162 39];
            app.RUNANALYSISButton.Text = 'RUN ANALYSIS';

            % Create StatusLabel
            app.StatusLabel = uilabel(app.UIFigure);
            app.StatusLabel.HorizontalAlignment = 'center';
            app.StatusLabel.Position = [30 129 254 30];

            % Create TYPEDropDownLabel
            app.TYPEDropDownLabel = uilabel(app.UIFigure);
            app.TYPEDropDownLabel.HorizontalAlignment = 'right';
            app.TYPEDropDownLabel.Position = [340 574 36 22];
            app.TYPEDropDownLabel.Text = 'TYPE';

            % Create TYPEDropDown
            app.TYPEDropDown = uidropdown(app.UIFigure);
            app.TYPEDropDown.Items = {'NONE', 'σ_x (Pa)', 'σ_y (Pa)', 'σ_r (Pa)', 'σ_theta (Pa)', 'τ_xy (Pa)', 'σ_z (Pa)', 'Von Mises (Pa)', 'Total Sigma (Pa)', 'u_x (m)', 'u_y (m)', 'u_r (m)', 'u_theta (m)'};
            app.TYPEDropDown.ValueChangedFcn = createCallbackFcn(app, @TYPEDropDownValueChanged, true);
            app.TYPEDropDown.Position = [391 574 100 22];
            app.TYPEDropDown.Value = 'NONE';

            % Create MESHINGButton
            app.MESHINGButton = uibutton(app.UIFigure, 'push');
            app.MESHINGButton.ButtonPushedFcn = createCallbackFcn(app, @MESHINGButtonPushed, true);
            app.MESHINGButton.BackgroundColor = [0.1608 0.349 0.5804];
            app.MESHINGButton.FontSize = 20;
            app.MESHINGButton.FontWeight = 'bold';
            app.MESHINGButton.FontColor = [1 1 1];
            app.MESHINGButton.Position = [764 557 178 39];
            app.MESHINGButton.Text = 'MESHING';

            % Create APPLIEDFORCEButton
            app.APPLIEDFORCEButton = uibutton(app.UIFigure, 'push');
            app.APPLIEDFORCEButton.ButtonPushedFcn = createCallbackFcn(app, @APPLIEDFORCEButtonPushed, true);
            app.APPLIEDFORCEButton.BackgroundColor = [0.6392 0.4745 0.0902];
            app.APPLIEDFORCEButton.FontSize = 20;
            app.APPLIEDFORCEButton.FontWeight = 'bold';
            app.APPLIEDFORCEButton.FontColor = [1 1 1];
            app.APPLIEDFORCEButton.Position = [764 499 177 39];
            app.APPLIEDFORCEButton.Text = 'APPLIED FORCE';

            % Create FEMPROGRAMINTERFACELabel
            app.FEMPROGRAMINTERFACELabel = uilabel(app.UIFigure);
            app.FEMPROGRAMINTERFACELabel.FontSize = 24;
            app.FEMPROGRAMINTERFACELabel.FontWeight = 'bold';
            app.FEMPROGRAMINTERFACELabel.FontColor = [0.7294 0.2196 0.2196];
            app.FEMPROGRAMINTERFACELabel.Position = [686 631 335 31];
            app.FEMPROGRAMINTERFACELabel.Text = 'FEM PROGRAM INTERFACE';

            % Create TextAreaLabel
            app.TextAreaLabel = uilabel(app.UIFigure);
            app.TextAreaLabel.HorizontalAlignment = 'right';
            app.TextAreaLabel.Position = [38 97 55 22];
            app.TextAreaLabel.Text = 'Text Area';

            % Create TextArea
            app.TextArea = uitextarea(app.UIFigure);
            app.TextArea.ValueChangedFcn = createCallbackFcn(app, @TextAreaValueChanged, true);
            app.TextArea.Position = [108 20 826 101];

            % Create RESETButton
            app.RESETButton = uibutton(app.UIFigure, 'push');
            app.RESETButton.ButtonPushedFcn = createCallbackFcn(app, @RESETButtonPushed, true);
            app.RESETButton.Position = [821 129 100 23];
            app.RESETButton.Text = 'RESET';

            % Create MODALButton
            app.MODALButton = uibutton(app.UIFigure, 'push');
            app.MODALButton.ButtonPushedFcn = createCallbackFcn(app, @MODALButtonPushed, true);
            app.MODALButton.BackgroundColor = [0.6118 0.1333 0.5765];
            app.MODALButton.FontSize = 20;
            app.MODALButton.FontWeight = 'bold';
            app.MODALButton.FontColor = [1 1 1];
            app.MODALButton.Position = [764 437 177 39];
            app.MODALButton.Text = 'MODAL';

            % Create MODESpinnerLabel
            app.MODESpinnerLabel = uilabel(app.UIFigure);
            app.MODESpinnerLabel.HorizontalAlignment = 'right';
            app.MODESpinnerLabel.Position = [777 393 41 22];
            app.MODESpinnerLabel.Text = 'MODE';

            % Create MODESpinner
            app.MODESpinner = uispinner(app.UIFigure);
            app.MODESpinner.Limits = [1 6];
            app.MODESpinner.ValueChangedFcn = createCallbackFcn(app, @MODALButtonPushed, true);
            app.MODESpinner.Position = [833 393 100 22];
            app.MODESpinner.Value = 1;

            % Create TYPEDropDown_2Label
            app.TYPEDropDown_2Label = uilabel(app.UIFigure);
            app.TYPEDropDown_2Label.HorizontalAlignment = 'right';
            app.TYPEDropDown_2Label.Position = [1160 594 36 22];
            app.TYPEDropDown_2Label.Text = 'TYPE';

            % Create TYPEDropDown_2
            app.TYPEDropDown_2 = uidropdown(app.UIFigure);
            app.TYPEDropDown_2.Items = {'NONE', 'σ_x (Pa)', 'σ_y (Pa)', 'σ_r (Pa)', 'σ_theta (Pa)', 'τ_xy (Pa)', 'σ_z (Pa)', 'Von Mises (Pa)', 'Total Sigma (Pa)', 'u_x (m)', 'u_y (m)', 'u_r (m)', 'u_theta (m)'};
            app.TYPEDropDown_2.Position = [1211 594 100 22];
            app.TYPEDropDown_2.Value = 'NONE';

            % Create CONVERGCEButton
            app.CONVERGCEButton = uibutton(app.UIFigure, 'push');
            app.CONVERGCEButton.BackgroundColor = [0.1608 0.349 0.5804];
            app.CONVERGCEButton.FontSize = 20;
            app.CONVERGCEButton.FontWeight = 'bold';
            app.CONVERGCEButton.FontColor = [1 1 1];
            app.CONVERGCEButton.Position = [1431 541 178 39];
            app.CONVERGCEButton.Text = 'CONVERGCE';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = FEM_Interface

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
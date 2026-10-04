%trc function ko co clc
function [W,Q] = quadrature2D(type, quadorder) %đăt tên function trùng với tên file .m

switch type
    case 'Q4'
        rlpt = zeros(quadorder,1);
        rlwt = zeros(quadorder,1);

        switch quadorder
            case 1
                rlpt(1)=0;
                rlwt(1)=2;

            case 2
                rlpt = [ 0.577350269189626;
                        -0.577350269189626];
                rlwt = [1;1];

            case 3
                rlpt = [ 0.774596669241483;
                        -0.774596669241483;
                         0];
                rlwt = [0.555555555555556;
                        0.555555555555556;
                        0.888888888888889];

            case 4
                rlpt = [ 0.861134311594053;
                        -0.861134311594053;
                         0.339981043584856;
                        -0.339981043584856];
                rlwt = [0.347854845137454;
                        0.347854845137454;
                        0.652145154862546;
                        0.652145154862546];
        end

        %gauss 2 chieu = tich 2 gauss 1 chieu 
        n = 1;
        Q = zeros(quadorder^2,2);
        W = zeros(quadorder^2,1);

        for i = 1:quadorder
            for j = 1:quadorder
                Q(n,:) = [rlpt(i), rlpt(j)];
                W(n)   = rlwt(i) * rlwt(j); %weight
                n = n + 1;
            end
        end

    case 'T3' %T3 chi co quad tu 1->3
        switch quadorder
            case 1
                Q = [1/3, 1/3];
                W = 1/2;

            case 3
                Q = [1/6 1/6;
                     2/3 1/6;
                     1/6 2/3];
                W = [1/6;1/6;1/6];
        end

    otherwise
        error('Chỉ nhận T3 hoạc Q4');

end
end
function test_shape_vem2
    % Undeformed vertices
    V0= [ 
         1 3;
         2 3;
         3 3;
         3 2;
         3 1;
         2 1;
         1 1;
         1 2;
        ];
    
    % Deformed vertices
%     V1= [ 
%          1 3;
%          2 3.5;
%          3 3;
%          3 2;
%          3 1;
%          2 0.5;
%          1 1;
%          1 2;
%         ];
    V1= [ 
         1 5;
         2 2.9;
         3 3;
         3 2;
         2.8 1.2;
         2 1;
         1 1;
         1.2 2;
        ];
    
    
    % Edges
    E=[ 1 2;
        2 3;
        3 4;
        4 5;
        5 6;
        6 7;
        7 8;
        8 1 ];
    % Plotting surface
    fig=figure(1)
    clf;
%     plot([V0(E(:,1),1)'; V0(E(:,2),1)'], [V0(E(:,1),2)'; V0(E(:,2),2)'],'Color',[0.5 0.5 0.5]);
%     hold on;
    hold on;
    
    x = 1:0.2:3;
    y = 1:0.2:3;
    [X,Y] = meshgrid(x,y);
    X=[X(:) Y(:)];
    k=6;

        
    p1 = plot(X(:,1),X(:,2),'.','Color','m','MarkerSize',15);
    hold on
    p2 = plot([V0(E(:,1),1)'; V0(E(:,2),1)'], [V0(E(:,1),2)'; V0(E(:,2),2)'],'Color','b');
    xlim([0 4]);
    ylim([0 5]);

    V=V0;
    
    n=50;
    dV = (V1-V0)/n;
    Vt = V0 + dV;
    for ii=1:n

        for i = 1:numel(p2)
            p2(i).XData = [Vt(E(i,1),1) Vt(E(i,2),1)];
            p2(i).YData = [Vt(E(i,1),2) Vt(E(i,2),2)];
        end
        
        X_com = mean(V0,1)';    % undeformed center of mass
        x_com = mean(Vt,1)';    % deformed center of mass


        % Compute monomial basis for each query point
        M=zeros(size(X,1), k);
        for i = 1:size(X,1)
            q=(X(i,:)'-X_com)';
            Mi = [1 q(1) q(2) q(1).^2 q(2).^2 q(1)*q(2)];
            M(i,:) = Mi;
        end

        P = zeros(size(E,1),2,2);
        % Computing shape matching matrix for each shape (edge in this case)
        for i =1:size(E,1)
            q = V0(E(i,:),:)' - X_com;
            p = Vt(E(i,:),:)' - x_com;
            P(i,:,:) = p * q' / (q * q');
        end

        % Compute weighting coefficients
        a = zeros(size(X,1), size(E,1));
        options = optimoptions('lsqlin','Display', 'off');

        for i = 1:size(X,1)
            C = V';
            d = X(i,:)';

            % All weights non-negative
            A = -eye(size(V,1));
            b = zeros(size(V,1),1);

            % Weights sum to one.
            Aeq = ones(1,size(E,1));
            beq = 1;
            
            a(i,:) = lsqlin(C,d,A,b,Aeq,beq, [],[],[], options);
        end

        % Compute deformed positions from 
        for i=1:size(X,1)
            x = [0 0];
            x2 = [0 0];

            for j=1:size(E,1)
                x = x + a(i,j) * Vt(j,:);
                g = squeeze(P(j,:,:)) * (X(i,:)'-X_com) + x_com;
                x2 = x2 + a(i,j) * g';
            end
            p1.XData(i) = x2(1);
            p1.YData(i) = x2(2);
        end
        drawnow;
%         fn=sprintf('output_png\\out2_%03d.png',ii)
%         saveas(fig,fn);
%         
        Vt = Vt + dV;
        
    end
    
    

end
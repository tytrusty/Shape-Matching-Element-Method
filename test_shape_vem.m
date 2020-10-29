function test_shape_vem
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
    figure(1)
    clf;
%     plot([V0(E(:,1),1)'; V0(E(:,2),1)'], [V0(E(:,1),2)'; V0(E(:,2),2)'],'Color',[0.5 0.5 0.5]);
%     hold on;
    plot([V1(E(:,1),1)'; V1(E(:,2),1)'], [V1(E(:,1),2)'; V1(E(:,2),2)'],'Color','b');
    hold on;
    
    V=V0;
    
    X_com = mean(V0,1)';    % undeformed center of mass
    x_com = mean(V1,1)';    % deformed center of mass
    
    x = 1:0.1:3;
    y = 1:0.1:3;
    [X,Y] = meshgrid(x,y);
%     plot(X(:),Y(:),'o','Color',[0.5 0.5 0.5],'MarkerSize',10);
    X=[X(:) Y(:)];
    k=6;
    
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
        p = V1(E(i,:),:)' - x_com;
        P(i,:,:) = p * q' / (q * q');
    end
   
    % Compute weighting coefficients
    a = zeros(size(X,1), size(E,1));
    for i = 1:size(X,1)
        C = V';
        d = X(i,:)';
        
        % All weights non-negative
        A = -eye(size(V,1));
        b = zeros(size(V,1),1);
        
        % Weights sum to one.
        Aeq = ones(1,size(E,1));
        beq = 1;
        lsqlin(C,d,A,b,Aeq,beq)
        a(i,:) = lsqlin(C,d,A,b,Aeq,beq);
    end
    
    % Compute deformed positions from 
    for i=1:size(X,1)
        x = [0 0];
        x2 = [0 0];
        
        for j=1:size(E,1)
            x = x + a(i,j) * V1(j,:);
            p = squeeze(P(j,:,:))
            
            g = squeeze(P(j,:,:)) * (X(i,:)'-X_com) + x_com;
            x2 = x2 + a(i,j) * g;
        end
        plot(x(1),x(2),'.','Color','r','MarkerSize',15);
        plot(x2(1),x2(2),'.','Color','m','MarkerSize',15);
        hold on
    end
        % All undeformed positions
%     n=10;
%     Xx=linspace(1,3,n);
%     Xy=repelem(3,n);
%         plot(X(1),X(2),'.','Color',[0.5 0.5 0.5]);

    plot(X_com(1), X_com(2),'x','MarkerSize',10,'Color','r');
    xlim([0 4]);
    ylim([0 4]);

end
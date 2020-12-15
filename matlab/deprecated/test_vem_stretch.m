function test_vem_stretch
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
        
    V1= [1 3; 4 3; 7 3; 7 2; 7 1; 4 1; 1 1;1 2;]; % stretch x-axis
%     V1= [1 3; 2 3; 7 3; 7 2; 7 1; 2 1; 1 1;1 2;]; % stretch x-axis
    
    % Edges
    E=[ 1 2;
        2 3;
        3 4;
        4 5;
        5 6;
        6 7;
        7 8;
        8 1 ];
    E1=[1 2 3 4 5 6 7 8];
    
%     E1=[ 8 1 2;
%     1 2 3;
%     2 3 4;
%     3 4 5;
%     4 5 6;
%     5 6 7;
%     6 7 8;
%     7 8 1 ];
%     E1=E;

    % Plotting surface
    figure(1)
    clf;

    plot([V0(E(:,1),1)'; V0(E(:,2),1)'], [V0(E(:,1),2)'; V0(E(:,2),2)'],'Color',[0.5 0.5 0.5]);
    hold on;
    plot([V1(E(:,1),1)'; V1(E(:,2),1)'], [V1(E(:,1),2)'; V1(E(:,2),2)'],'Color','b');
    hold on;
    
    V=V0;
    
    X_com = mean(V0,1)';    % undeformed center of mass
    x_com = mean(V1,1)';    % deformed center of mass
    
    x = 1:0.1:3;
    y = 1:0.2:3;
    [X,Y] = meshgrid(x,y);
    X=[X(:) Y(:)];

    k=5;
%     k=2;
    
    P = zeros(size(E1,1),2,k);
    % Computing shape matching matrix for each shape (edge in this case)
    for i =1:size(E1,1)
        q = V0(E1(i,:),:)' - X_com;
        
        q_ = zeros(5, size(q,2));
        q_(1:2,:) = q;
        q_(3,:) = q(1,:).^2;
        q_(4,:) = q(2,:).^2;
        q_(5,:) = q(1,:).*q(2,:);
        q=q_;
        
        [SU, S, SV] = svd(q*q');
        S = diag(S);
        S = 1 ./ max(S, 1e-4);
        Aqq = SV * diag(S) * SU';
        
        
        p = V1(E1(i,:),:)' - x_com;
        P(i,:,:) = (p*q') * Aqq;
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
        a(i,:) = lsqlin(C,d,A,b,Aeq,beq,[],[],[],options);
    end
    
    % Compute deformed positions from 
    for i=1:size(X,1)
        x = [0 0];
        x2 = [0 0];
        
        for j=1:size(E1,1)
            x = x + a(i,j) * V1(j,:);
            
            q = X(i,:)' - X_com;
            q_ = zeros(5, size(q,2));
            q_(1:2,:) = q;
            q_(3,:) = q(1,:).^2;
            q_(4,:) = q(2,:).^2;
            q_(5,:) = q(1,:).*q(2,:);
            q=q_;
        
            g = squeeze(P(j,:,:)) * q + x_com;
            x2 = x2 + 1 * g';
        end
%         plot(x(1),x(2),'.','Color','r','MarkerSize',15);
        plot(x2(1),x2(2),'.','Color','m','MarkerSize',15);
        hold on
    end

    plot(X_com(1), X_com(2),'x','MarkerSize',10,'Color',[0.5 0.5 0.5]);
    plot(x_com(1), x_com(2),'x','MarkerSize',10,'Color','b');
    xlim([0 8]);
%     ylim([0 4]);
    axis equal
    plot([V1(E(:,1),1)'; V1(E(:,2),1)'], [V1(E(:,1),2)'; V1(E(:,2),2)'],'.','MarkerSize', 25, 'Color','b');
end
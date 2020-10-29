function test_vem_v3
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
         1 3.5;
         2 3.5;
         3 3.5;
         3 2;
         3 1;
         2 0.5;
         1 1;
         1 2;
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
    plot([V0(E(:,1),1)'; V0(E(:,2),1)'], [V0(E(:,1),2)'; V0(E(:,2),2)'],'Color',[0.5 0.5 0.5]);
    hold on;
    plot([V1(E(:,1),1)'; V1(E(:,2),1)'], [V1(E(:,1),2)'; V1(E(:,2),2)'],'Color','b');
    hold on;
    
    V=V0;
    
    X_com = mean(V0,1)';    % undeformed center of mass
    
    x = 1:0.1:3;
    y = 1:0.1:3;
    [X,Y] = meshgrid(x,y);
    plot(X(:),Y(:),'o','Color',[0.5 0.5 0.5],'MarkerSize',10);
    X=[X(:) Y(:)];
    k=6;
    
    % Compute monomial basis for each query point
    M=zeros(size(X,1), k);
    for i = 1:size(X,1)
        q=(X(i,:)'-X_com)';
        Mi = [1 q(1) q(2) q(1).^2 q(2).^2 q(1)*q(2)];
        M(i,:) = Mi;
    end
    
    % pseudo-inverse for M
    Pm = pinv(M);
        
    ProjX = zeros(size(E,1),size(X,1),2);
    
    P = zeros(size(E,1), k, 1);
    
    % Compute per-edge projection matrix
    for i = 1:size(E,1)
%         u = [V(E(i,1),:); V(E(i,2),:)]';
        u = V(i,:)';
        
        Pu = pinv(u)';
        
        p = Pm*X*Pu;
        x=M*p*u';
        ProjX(i,:,:) = x;
        P(i,:,:)=p;
    end
    
    % Compute weighting terms
    a = zeros(size(X,1), size(E,1));
    for i = 1:size(X,1)
        C = squeeze(ProjX(:,i,:))';
        d = X(i,:)';
        
        % All weights non-negative
        A = -eye(size(E,1));
        b = zeros(size(E,1));
        
        % Weights sum to one.
        Aeq = ones(1,size(E,1));
        beq = 1;
        
        a(i,:) = lsqlin(C,d,A,b,Aeq,beq);
    end
    
    % Compute deformed positions from 
    for i=1:size(X,1)
        x = [0 0];
        
        for j=1:size(E,1)
%             u = [V1(E(j,1),:); V1(E(j,2),:)]';
            u = V1(j,:); % Deformed DOF values
            size(squeeze(P(j,:,:)))
            x = x + a(i,j) * M(i,:) * squeeze(P(j,:,:))' * u;
            x
        end
        plot(x(1),x(2),'.','Color','r','MarkerSize',15);
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
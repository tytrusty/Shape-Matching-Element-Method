function test_vem_v4
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
         1 8;
         2 3.0;
         3 3.0;
         3 2;
         3 1;
         2 1;
         1 1;
         1 2;
        ];
    
    
    % Edges
    E=[ 1 2 3 4 5 6 7 8];
    % E=[ 1 2 3];

    E1=[ 1 2;
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
    plot([V0(E1(:,1),1)'; V0(E1(:,2),1)'], [V0(E1(:,1),2)'; V0(E1(:,2),2)'],'Color',[0.5 0.5 0.5]);
    hold on;
    plot([V1(E1(:,1),1)'; V1(E1(:,2),1)'], [V1(E1(:,1),2)'; V1(E1(:,2),2)'],'Color','b');
    hold on;
    
    V=V0;
    
    X_com = mean(V0,1)';    % undeformed center of mass
    x_com = mean(V1,1)';    % undeformed center of mass
    
    x = 1:0.2:3;
    y = 1:0.2:3;
    [X,Y] = meshgrid(x,y);
    plot(X(:),Y(:),'o','Color',[0.5 0.5 0.5],'MarkerSize',10);
    X=[X(:) Y(:)];
    k=6;
    
    % Compute monomial basis for each query point
    M=zeros(2*size(X,1), k);
    for i = 1:size(X,1)
        q=(X(i,:)'-X_com)';
        Mi = [1 0; 0 1; q; q(1) 0; 0 q(2); -q(1) q(2)]';
        M(2*i-1:2*i,:) = Mi;
    end
    
    % pseudo-inverse for M
    Pm = pinv(M);
        
    ProjX = zeros(2*size(X,1),size(V,1));
    X_flat = X';
    X_flat = X_flat(:);
    
    P = zeros(size(E,1), k, 2);
    
    % Compute per-edge projection matrix
    for i = 1:size(E,1)
        u = V(E(i,:),:)' ;
        Pu = pinv(u);

        X_ = repelem(X_flat, 1, numel(E(i,:)));
        p = Pm*X_*Pu;
        x=M*p*u;
        norm(X_-x)
        ProjX(:,E(i,:)) = x;
        P(i,:,:)=p;
    end
    
    % Compute weighting terms
    a = zeros(size(X,1), size(E,2));
    for i = 1:size(X,1)
        C = squeeze(ProjX(2*i-1:2*i,:));
        d = X(i,:)';
        
        % All weights non-negative
        A = -eye(size(E,2));
        b = zeros(size(E,2),1);
        
        % Weights sum to one.
        Aeq = ones(1,size(E,2));
        beq = 1;
        a(i,:) = lsqlin(C,d,A,b,Aeq,beq);
    end
    
    % Compute deformed positions from 
    for i=1:size(X,1)
        x = [0 0]';
        
        u = V1(E(1,:),:)';
        
        for j=1:size(E,2)
%             size(squeeze(P(1,:,:)))
%             size(u(:,j))
%             size(M(2*i-1:2*i,:))
            x = x + a(i,j) * M(2*i-1:2*i,:) * squeeze(P(1,:,:)) * u(:,j);
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
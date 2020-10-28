function test_vem_v2
    % Undeformed vertices
    V0= [ 
         1 3;
         2 3;
         3 3;
         3 1;
         2 1 
         1 1;
        ];
    
    % Deformed vertices
    V1= [ 
         1 3;
         2 3.5;
         3 3;
         3 1;
         2 0.5 
         1 1;
        ];
    
    % Edges
    E=[ 1 2;
        2 3;
        3 4;
        4 5;
        5 6;
        6 1 ];
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
    
    % All undeformed positions
%     n=10;
%     Xx=linspace(1,3,n);
%     Xy=repelem(3,n);
    n=1;
    Xx=[2.5];
    Xy=[2.75];
    
    for i = 1:n
        X=[Xx(i) Xy(i)]';
        
        % Compute monomial basis
        q=(X-X_com)';
        M = [1 1; q; q.^2; q.*circshift(q,1)]';
        
        % Computing projection onto polynomial coefficients
        Pm = M'/(M*M'); % pseudo-inverse for M
        
        PU = zeros(size(M,2),size(V,1));
        % Stack DOFs into one column
        for j=1:size(V,1)
            u = V(j,:)';
            ue = [V(E(j,1),:) V(E(j,2),:)];
            Pu = (u'*u)\u'; % pseudo-inverse for U
            P = Pm*X*Pu
            PU(:,j)=P*u;
        end
        
        A = M*PU;
        
        plot(X(1),X(2),'.','Color',[0.5 0.5 0.5]);
        hold on;
    end

        
    % Plotting surface
    %     plot(V(:,1), V(:,2),'.','MarkerSize',15,'Color','b');
    %     hold on;
    plot(X_com(1), X_com(2),'x','MarkerSize',10,'Color','r');
    xlim([0 4]);
    ylim([0 4]);

end
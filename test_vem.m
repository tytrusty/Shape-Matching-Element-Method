function test_vem
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
    n=10;
    Xx=linspace(1,3,n);
    Xy=repelem(2,n);
    
    for i = 1:n
        X=[Xx(i) Xy(i)]';

        q=(X-X_com)';
        M = [1 1; q; q.^2; q.*circshift(q,1)]';
        
        % Stack DOFs into one column
        U=V';
        U=U(:);

        % Computing projection onto polynomial coefficients
        Pu = (U'*U)\U'; % pseudo-inverse for U
        Pm = M'/(M*M'); % pseudo-inverse for M
        P = Pm*X*Pu;

        U1=V1';
        U1=U1(:);
        x = M*P*U1
        
        plot(X(1),X(2),'.','Color',[0.5 0.5 0.5]);
        plot(x(1),x(2),'o','Color','r');
    end

        
    % Plotting surface
    %     plot(V(:,1), V(:,2),'.','MarkerSize',15,'Color','b');
    %     hold on;
    plot(X_com(1), X_com(2),'x','MarkerSize',10,'Color','r');
    xlim([0 4]);
    ylim([0 4]);

end
function nurbs = nurbs_plot(nurbs)
    cm=jet(numel(nurbs));
    for ii = 1:numel(nurbs)
        %         plt = patch('Faces',nurbs{ii}.hires_T,'Vertices',nurbs{ii}.hires_x0', ...

        plt = patch('Faces',nurbs{ii}.T,'Vertices',nurbs{ii}.x0', ...
            'FaceAlpha',0.5,'EdgeColor','none','FaceColor',cm(ii,:));
        hold on;
        nurbs{ii}.plt=plt;
    end
    
    set(gcf,'color','w');
    axis equal
    lighting gouraud;
    
    lightangle(gca,0, 20)
%     lightangle(gca,0, 40)
    %     zlim([-30 50]);
    %     ylim([-20 20]);
    xlabel('x');
    zlabel('z');
    view(-15,30);
%     view(0,0);
%     view(28,13);
end
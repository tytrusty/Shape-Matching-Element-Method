function nurbs = nurbs_plot(nurbs)
    cm=jet(numel(nurbs));
    for ii = 1:numel(nurbs)
        plt = patch('Faces',nurbs{ii}.T,'Vertices',nurbs{ii}.x0', ...
            'FaceAlpha',.9,'EdgeColor','none','FaceColor',cm(ii,:));
        hold on;
        nurbs{ii}.plt=plt;
    end
    
    set(gcf,'color','w');
    axis equal
    lighting gouraud;
%         shading interp
    
%     lightangle(gca,0, 20)
%     zlim([-30 50]);
%     ylim([-20 20]);
%     axis auto;
    xlabel('x');
end
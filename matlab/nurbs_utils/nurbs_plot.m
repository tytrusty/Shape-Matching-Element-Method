function nurbs = nurbs_plot(nurbs)
    cm=jet(numel(nurbs));
    for ii = 1:numel(nurbs)
        xi = reshape(nurbs{ii}.x0, 3, nurbs{ii}.subd(1), nurbs{ii}.subd(2));
        plt = surf(squeeze(xi(1,:,:)),squeeze(xi(2,:,:)),squeeze(xi(3,:,:)),'FaceAlpha',.9,'EdgeColor','black','FaceColor',cm(ii,:));
        hold on;
        nurbs{ii}.plt=plt;
    end
    
    set(gcf,'color','w');
    axis equal
    lighting gouraud;
    %     shading interp
    
    lightangle(gca,0, 20)
%         zlim([-70 150]);
    xlabel('x');
end
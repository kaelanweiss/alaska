function pos = cbarpos(ax,offset,width,varargin)
% Function to calculate position of colorbar to right of existing axes.
%
% pos = cbarpos(ax,offset,width)


% get axes position
if numel(ax) == 1 && isa(ax,'matlab.graphics.axis.Axes')
    ax_pos = get(ax,'position');
elseif numel(ax) > 1 && isa(ax(1),'matlab.graphics.axis.Axes')
    % create new ax_pos that encompasses all tiled axes
    ax_pos = nan(numel(ax),4);
    for i = 1:numel(ax)
        ax_pos(i,:) = ax(i).Position;
    end
    x_right = ax_pos(:,1)+ax_pos(:,3);
    y_top = ax_pos(:,2)+ax_pos(:,4);
    ax_pos = [min(ax_pos(:,1)) min(ax_pos(:,2)) max(x_right)-min(ax_pos(:,1)) max(y_top)-min(ax_pos(:,2))];

else
    ax_pos = ax;
end

% change x_cbar to x+dx+offset
pos = ax_pos;
pos(1) = pos(1)+pos(3)+offset;

% change dx_cbar to width
pos(3) = width;

end


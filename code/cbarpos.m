function pos = cbarpos(ax,offset,width)
% Function to calculate position of colorbar to right of existing axes.
%
% pos = cbarpos(ax,offset,width)


% get axes position
if isa(ax,'matlab.graphics.axis.Axes')
    pos = get(ax,'position');
else
    pos = ax;
end

% change x_cbar to x+dx+offset
pos(1) = pos(1)+pos(3)+offset;

% change dx_cbar to width
pos(3) = width;

end


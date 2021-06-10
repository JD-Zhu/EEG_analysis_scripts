% find_centroid.m 
%
% Given a list of points (specified by their coordinates), find the 
% centre of mass and output its coordinates.
% Works in any number of dimensions. Works for both convex and concave
% (if the points form a concave surface, the convex hull will be used, 
% and the "inside" points will be ignored).
%
% @param P: the list of points (specified by their xyz-coordinates)
%
% Requires: centroid.m
% https://au.mathworks.com/matlabcentral/fileexchange/8514-centroid-of-a-convex-n-dimensional-polyhedron
%
function C = find_centroid (P)
    %P = gallery('uniformdata',30,3,5); % for testing purpose, generate some random points

    % find the convex hull for P (specified as triangle planes formed by indices of points in P)
    % Any of these 3 methods below should work:
    %hull = boundary(P, 0); 
    %hull = convhull(P);
    hull = convhulln(P);

    % remove repeated points in the hull, then find the centroid
    hull_points = unique(hull);
    C = centroid(P(hull_points, :));

    % plot for quality check
    %{
    figure; hold on;
    plot3(P(:,1),P(:,2),P(:,3),'.','MarkerSize',10); grid on; % plot all points
    trisurf(hull, P(:,1),P(:,2),P(:,3),'Facecolor','red','FaceAlpha',0.1); % plot convex hull
    plot3(C(1), C(2), C(3), '*', 'MarkerSize', 10); % plot centroid
    hold off;
    %}
end

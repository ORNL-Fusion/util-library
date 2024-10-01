function area = calc_polygon_area(x, y)
% area = compute_polygon_area(x, y)
% calculate the area of a polygon using Gauss's formula

    n = length(x); % Number of sides
    
    % Check if the loop is already closed
    if x(1) == x(n) && y(1) == y(n)
        n = n - 1; % Reduce the count as the loop is already closed
    else
        % Add the first vertex at the end to complete the loop
        x(end+1) = x(1);
        y(end+1) = y(1);
    end
    
    area = 0;
    for i = 1:n
        area = area + (x(i) * y(i+1) - x(i+1) * y(i));
    end   
    area = abs(area) / 2;
end

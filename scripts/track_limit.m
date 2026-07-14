function [value, isterminal, direction] = track_limit(t, z, limit)
    value = limit - abs(z(1));
    isterminal = 1;
    direction = 0;
end
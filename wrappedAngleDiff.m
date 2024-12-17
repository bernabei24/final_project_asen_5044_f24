function angleDiff = wrappedAngleDiff(angle1, angle2)
    angleDiff = mod(angle1 - angle2 + pi, 2*pi) - pi;
end
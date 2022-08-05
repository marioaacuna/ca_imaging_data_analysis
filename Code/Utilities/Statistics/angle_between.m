function theta = angle_between(v1, v2)
% Returns the angle in radians between vectors 'v1' and 'v2'.
%
%   angle_between([1, 0, 0], [0, 1, 0])
%       1.5707963267948966
%   angle_between([1, 0, 0], [1, 0, 0])
%   	0.0
%   angle_between([1, 0, 0], [-1, 0, 0])
%       3.141592653589793
%
% Reference: https://stackoverflow.com/a/13849249
%

v1_u = unit_vector(v1);
v2_u = unit_vector(v2);
theta = acosd(clip(dot(v1_u, v2_u), -1.0, 1.0));


function uv = unit_vector(vector)
    % Returns the unit vector of the vector
    uv = vector / norm(vector);


function a = clip(a, a_min, a_max)
    a(a < a_min) = a_min;
    a(a > a_max) = a_max;

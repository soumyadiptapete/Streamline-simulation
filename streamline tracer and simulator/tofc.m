function [ tof_d ] = tofc( t_d )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
global e;
b_d=sort(t_d);
for i=1:4
    if b_d(i)>e;
        tof_d=b_d(i);
       break;
    end

end


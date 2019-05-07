function [m_new,n_new,a01,b01 ] = Exit_face(tof1,t111,t121,t211,t221,alpha_d,beta_d,m_d1,n_d1)
global e;
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if abs(tof1-t111)<=e
    % exit on left face
    m_new=m_d1;
    n_new=n_d1-1;
    a01=1;
    b01=beta_d;
end

if abs(tof1-t121)<=e
    %exit on right face
    m_new=m_d1;
    n_new=n_d1+1;
    a01=0;
    b01=beta_d;
end

if abs(tof1-t211)<=e
    % exit on bottom face
    m_new=m_d1-1;
    n_new=n_d1;
    b01=1;
    a01=alpha_d;
end

if abs(tof1-t221)<=e
    %exit on top face
    m_new=m_d1+1;
    n_new=n_d1;
    b01=0;
    a01=alpha_d;
    
end

%    %exit on top right corner
 if (abs(tof1-t121)<=e && abs(tof1-t221)<=e)
    m_new=m_d1+1;
    n_new=n_d1+1;
    b01=0;
    a01=0;
 end
 
 % exit on bottom left corner 
 if abs(tof1-t111)<=e && abs(tof1-t211)<=e
    m_new=m_d1-1;
    n_new=n_d1-1;
    a01=1;
    b01=1;
 end

    
end


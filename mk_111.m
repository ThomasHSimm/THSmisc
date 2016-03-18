%% Plots 111 fcc planes 
%% modified by TH Simm 16
%% input is inv(orientation.Euler) for MTEX or orientation.Euler for
%% Channel5
%% plane111 are 111, -111, 1-11, -1-11
function mk_111(phi1,phi,phi2,plane111)
if nargin==1
    phi2=phi1(3);phi=phi1(2);
    phi1=phi1(1);plane111=1;
elseif nargin==3
    plane111=1;
elseif nargin==0
    phi2=0;phi=0;
    phi1=0;plane111=1;
end 
if plane111==1;figure(101);close;f_handle=figure(101);
elseif plane111==2;figure(102);close;f_handle=figure(102);
elseif plane111==3; figure(103);close; f_handle=figure(103);
elseif plane111==4; figure(104);close; f_handle=figure(104);  
elseif plane111==5 figure(105);close; f_handle=figure(105);
end
%%
         %bottom cube -Top cube -x-side   -x0 side  -yside   -y0 side    %%111 plane  %
xvert=  [ 0 1 1 0      0 1 1 0   1 1 1 1  0 0 0 0   0 1 1 0   0 1 1 0   ...
        1 0 0      1 0 1        1 0 0       1 0 1  ]'; 
        %%111       -111        1-11        -(11-1)
yvert=  [ 0 0 1 1      0 0 1 1   0 1 1 0  0 1 1 0   1 1 1 1   0 0 0 0    ...
        0 1 0      1 0 0       1 0 1       0 1 1 ]';
        %%111       -111       1-11         -(11-1)                             %%24
zvert=( [ 0 0 0 0      1 1 1 1   0 0 1 1  0 0 1 1   0 0 1 1   0 0 1 1    ...
        0 0 1      0 0 1        0 0 1       0 0 1 ]'); 
        %%111       -111        1-11        -(11-1)
v=[xvert yvert zvert];
 
    %bot        %top        %x+             %x-
f=[1 2 3 4;     5 6 7 8;    9 10 11 12;     13 14 15 16; ...
    %y+             %y-             %111
    17 18 19 20;    21 22 23 24      ;25 26 27 27; ];
fc=[1 1 1; 1 1 1;1 1 1;1 1 1; 1 1 1; 1 1 1; 1 0 0; ];%colours

switch plane111
    case 1
        
    case 2
        f(end,:)=[28 29 30 30];
        fc(end,:)=[0 1 1];
    case 3
        f(end,:)=[31 32 33 33];
        fc(end,:)=[1 0 1];
    case 4
        f(end,:)=[34 35 36 36];
        fc(end,:)=[0 0 1];
    case 5
        f(end+1,:)=[28 29 30 30];
        fc(end+1,:)=[0 1 1];  
        f(end+1,:)=[31 32 33 33];
        fc(end+1,:)=[1 0 1];
        f(end+1,:)=[34 35 36 36];
        fc(end+1,:)=[0 0 1];
end    

h=patch('Vertices',v,'Faces',f, 'FaceColor', 'flat','FaceVertexCData',fc, 'FaceAlpha',0.4, 'LineWidth', 3);  
%%
grid,box
axis equal;
% set(gca,'xticklabels','')
% set(gca,'yticklabels','')
% set(gca,'zticklabels','')
xlabel('x','fontsize',16)
ylabel('y','fontsize',16)
%%
rotate0(h,phi1, phi, phi2);
%%
view([0 0 1])
view([.5 .5 0.5])

%%  
    function rotate0(h0, phi1, PHI, phi2)
        rotate(h0,[0 0 1],(phi1*180/pi));%%%rotation about z-axis x->y by phi1
        vec_a=[cos(phi1) sin(phi1) 0];%%%%%%
        rotate(h0,vec_a, phi*180/pi);%%%%rotation about new x axis by phi
        vec_b=[sin(phi1)*sin(phi) -cos(phi1)*sin(phi) cos(phi)]; 
        rotate(h0,vec_b, (phi2)*180/pi);%%%%rotation about the new z-axis given by the z-component of the rotation matrix
 
    end
end

% In this script, we calculated the coordinates of the electrodes used in this project, 
%following the international 10-20 standard as closely as possible.
coord=zeros(30,2);

r=100;
t=linspace(0,2*pi,1000);
head_x      = r.*cos(t);        % Cartesian X of the head
head_y      = r.*sin(t);        % Cartesian Y of the head

figure
plot(head_x,head_y)
axis equal;
hold on;

ri=4/5*r;
head_xi      = ri.*cos(t);        % Cartesian X of the head
head_yi      = ri.*sin(t);        % Cartesian Y of the head

plot(head_xi,head_yi,'--k')
hold on;
plot([-100,100],[0,0],'--b'); hold on;
plot([0,0],[-100,100],'--b'); hold on;

coord(13,1:2)=[0,-ri]; %Oz

ti=linspace(0,2*pi,21);
coord(19,1:2)=[ri*cos(ti(3)),ri*sin(ti(3))]; % F8
coord(16,1:2)=[ri*cos(ti(5)),ri*sin(ti(5))]; % FP2
coord(1,1:2)=[ri*cos(ti(7)),ri*sin(ti(7))]; % FP1
coord(3,1:2)=[ri*cos(ti(9)),ri*sin(ti(9))]; % F7
coord(8,1:2)=[ri*cos(ti(13)),ri*sin(ti(13))]; % P7
coord(10,1:2)=[ri*cos(ti(14)),ri*sin(ti(14))]; % PO7
coord(12,1:2)=[ri*cos(ti(15)),ri*sin(ti(15))]; % O1
coord(30,1:2)=[ri*cos(ti(17)),ri*sin(ti(17))]; % O2
coord(28,1:2)=[ri*cos(ti(18)),ri*sin(ti(18))]; % PO8
coord(26,1:2)=[ri*cos(ti(19)),ri*sin(ti(19))]; % P8
coord(9,1:2)=[r*cos(ti(13)),r*sin(ti(13))]; % P9
coord(27,1:2)=[r*cos(ti(19)),r*sin(ti(19))]; % P10
coord(22,1:2)=[0,0]; % Cz
coord(15,1:2)=[0,-r/5]; % CPz
coord(14,1:2)=[0,-2*r/5]; % Pz
coord(21,1:2)=[0,r/5]; % FCz
coord(17,1:2)=[0,2*r/5]; % Fz
coord(23,1:2)=[2*r/5,0]; % C4
coord(24,1:2)=[3*r/5,0]; % C6
coord(5,1:2)=[-2*r/5,0]; % C3
coord(6,1:2)=[-3*r/5,0]; % C5
coord(2,1:2)=(coord(3,1:2)+coord(17,1:2))/2; % F3
coord(18,1:2)=(coord(19,1:2)+coord(17,1:2))/2; % F4
coord(7,1:2)=(coord(14,1:2)+coord(8,1:2))/2; % P3
coord(25,1:2)=(coord(14,1:2)+coord(26,1:2))/2; % P4
coord(4,1:2)=(coord(2,1:2)+coord(5,1:2))/2; % FC3
coord(20,1:2)=(coord(18,1:2)+coord(23,1:2))/2; % FC4
coord(11,1:2)=(coord(10,1:2)+coord(14,1:2))/2; % PO3
coord(29,1:2)=(coord(14,1:2)+coord(28,1:2))/2; % PO4


% Dummy electrodes
ind=[1:12,14:18,20];
for k=ind
    coord(30+k,1:2)=[r*cos(ti(k)),r*sin(ti(k))]; 
end

scatter(coord(:,1),coord(:,2),'k'); hold on;



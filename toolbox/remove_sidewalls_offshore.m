function [cam,lidar] = remove_sidewalls_offshore(cam,lidar)
% remove side walls & offshore

% LiDAR
lidar.Hs(lidar.y<-12.8,:)=NaN;
lidar.Hs(lidar.y>12.8,:)=NaN;
lidar.Hs(:,lidar.x<24) = NaN;
lidar.See(lidar.y<-12.8,:,:)=NaN;
lidar.See(lidar.y>12.8,:)= NaN;
lidar.See(:,lidar.x<24,:) = NaN;

% Stereo
cam.Hs(cam.y(:,1)<-12.9,:)=NaN;
cam.Hs(cam.y(:,1)>13,:)=NaN;
cam.Hs(:,cam.x(1,:)<28.35) = NaN;
cam.See(cam.y(:,1)<-12.9,:,:)=NaN;
cam.See(cam.y(:,1)>13,:,:)=NaN;
cam.See(:,cam.x(1,:)<28.35,:) = NaN;
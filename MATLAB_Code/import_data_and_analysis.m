%%% DATA PREPARATION %%%
%  To Prepare the Data we must remove all unauthorized tags to ensure we
%  have equal detections from each antenna on the station. 
%%%%%%%%%%%%%%%%%%%%%%%%

%%% WORKSPACE PREPATION AND DATA IMPORT%%%
% We start by clearing all variables and figures as well as our workspace.
% We then continue by importing our two datasets. They were collected on
% two seperate days while utilizing the same calibration technique (pole).
% We use a series of dummy variables to collate the data and then split the
% data up by antenna. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
clf; 
close all;

T = readtable('apr7_tag.csv');
%T.Time.Format = 'MM/dd/yyyy HH:mm:ss';
T1 = table2cell(T);
T3 = T1(:,1:4);

% TA = readtable('CTT_344C4B66_Test_Data_Day_2.csv');
% TA1 = table2cell(TA);
% T4 = TA1(:,1:4);

T2 = [T3];

firstColumn = cell2mat(T1(:,2));
Antenna_1  = T2(firstColumn==1, :);
Antenna_2  = T2(firstColumn==2, :);
Antenna_3  = T2(firstColumn==3, :);
Antenna_4  = T2(firstColumn==4, :);

%Once the data has been inputed we need to convert the times given to us by
%the station to POSIX time. We do this four times once for each antenna. I did this so it is easier to compare apples and apples when comparing Strava time to Station Time.   

holding_matrix = [];

for x=1:numel(Antenna_1)/4
    
    Ant_1_time = posixtime(Antenna_1{x,1});
    holding_matrix=[holding_matrix;Ant_1_time];
end

converted_antenna_1_matrix = [holding_matrix,cell2mat(Antenna_1(:,4))];
clear x; 


 holding_matrix_2 = [];

for x=1:numel(Antenna_2)/4
    
    Ant_2_time = posixtime(Antenna_2{x,1});
    holding_matrix_2=[holding_matrix_2;Ant_2_time];
end

converted_antenna_2_matrix = [holding_matrix_2,cell2mat(Antenna_2(:,4))];
clear x; 

converted_antenna_2_matrix = [holding_matrix_2,cell2mat(Antenna_2(:,4))];
clear x; 
holding_matrix_3 = [];

for x=1:numel(Antenna_3)/4
    
    Ant_3_time = posixtime(Antenna_3{x,1});
    holding_matrix_3=[holding_matrix_3;Ant_3_time];
end

converted_antenna_3_matrix = [holding_matrix_3,cell2mat(Antenna_3(:,4))];
clear x; 
holding_matrix_4 = [];

for x=1:numel(Antenna_4)/4
    
    Ant_4_time = posixtime(Antenna_4{x,1});
    holding_matrix_4=[holding_matrix_4;Ant_4_time];
end

converted_antenna_4_matrix = [holding_matrix_4,cell2mat(Antenna_4(:,4))];

%converted_antenna_1_matrix = [converted_antenna_1_matrix;converted_antenna_2_matrix;converted_antenna_3_matrix;converted_antenna_4_matrix];

%Once we have done that, we hasve to read in our Strava Data. I suppose I
%could have done this before inporting the antenna data. It's
%interchangeable I guess. 

%%%% THE FOLLOWING IS ANTENNA INDEPENDENT %%%%
% P = readtable('Morning_Ride.csv');
% P1 = table2cell(P);
% P2 = [P1(:,2),P1(:,1),P1(:,6),P1(:,7)];
% P4 = posixtime(datetime(P2(:,4),'InputFormat','yyyy/MM/dd HH:mm:ss+SS'));
% A1 = cell2mat(P2(:,1:3));

P = readtable('apr9_gps.csv');
P1 = table2cell(P);
P2 = [P1(:,2),P1(:,1),P1(:,6),P1(:,7)];
P4 = posixtime(datetime(P2(:,4),'InputFormat','yyyy/MM/dd HH:mm:ss+SS'));
A1 = cell2mat(P2(:,1:3));

% P = readtable('Morning_Ride_1.csv');
% P1 = table2cell(P);
% P2 = [P1(:,2),P1(:,1),P1(:,6),P1(:,7)];
% P5 = posixtime(datetime(P2(:,4),'InputFormat','yyyy/MM/dd HH:mm:ss+SS'));
% A2 = cell2mat(P2(:,1:3));

% P = readtable('Morning_Ride_2.csv');
% P1 = table2cell(P);
% P2 = [P1(:,2),P1(:,1),P1(:,6),P1(:,7)];
% P6 = posixtime(datetime(P2(:,4),'InputFormat','yyyy/MM/dd HH:mm:ss+SS'));
% A3 = cell2mat(P2(:,1:3));
% 
% 
% P = readtable('Afternoon_Ride.csv');
% P1 = table2cell(P);
% P2 = [P1(:,2),P1(:,1),P1(:,6),P1(:,7)];
% P7 = posixtime(datetime(P2(:,4),'InputFormat','yyyy/MM/dd HH:mm:ss+SS'));
% A4 = cell2mat(P2(:,1:3));
% 
% P = readtable('Lunch_Ride_1.csv');
% P1 = table2cell(P);
% P2 = [P1(:,2),P1(:,1),P1(:,6),P1(:,7)];
% P8 = posixtime(datetime(P2(:,4),'InputFormat','yyyy/MM/dd HH:mm:ss+SS'));
% A5 = cell2mat(P2(:,1:3));
% 
% P = readtable('Morning_Ride_3.csv');
% P1 = table2cell(P);
% P2 = [P1(:,2),P1(:,1),P1(:,6),P1(:,7)];
% P9 = posixtime(datetime(P2(:,4),'InputFormat','yyyy/MM/dd HH:mm:ss+SS'));
% A6 = cell2mat(P2(:,1:3));
% 
% P = readtable('20210311110346-19125-data.csv');
% P1 = table2cell(P);
% P2 = [P1(:,2),P1(:,1),P1(:,6),P1(:,7)];
% P10 = posixtime(datetime(P2(:,4),'InputFormat','yyyy/MM/dd HH:mm:ss+SS'));
% A7 = cell2mat(P2(:,1:3));
% 
% 
% P = readtable('1100_B3DA0657.csv');
% P1 = table2cell(P);
% P2 = [P1(:,2),P1(:,1),P1(:,6),P1(:,7)];
% P11 = posixtime(datetime(P2(:,4),'InputFormat','yyyy/MM/dd HH:mm:ss+SS'));
% A8 = cell2mat(P2(:,1:3));


%P3 = [P4;P5;P6;P7;P8;P9;P10;P11];
%A = [A1;A2;A3;A4;A5;A6;A7;A8];
P3 = [P4];
A = [A1];
converted_matrix_strava = [A,P3];

antenna_location = [41.520752667,-71.599254667];

%%%% THE ABOVE IS ANTENNA INDEPENDENT %%%%

latlong_of_detections = [];
latlong_of_detections_2 = [];
latlong_of_detections_3  = [];
latlong_of_detections_4 = [];

for z=1:numel(holding_matrix)
    
    for zx=1:numel(converted_matrix_strava(:,4))
        
        if holding_matrix(z) ==  converted_matrix_strava(zx,4)
            latlong_of_detections = [latlong_of_detections;A(zx,:),converted_antenna_1_matrix(z,2)];
            
        end
        
        
    end
        
end

for z=1:numel(holding_matrix_2)
    
    for zx=1:numel(converted_matrix_strava(:,4))
        
        if holding_matrix_2(z) ==  converted_matrix_strava(zx,4)
            latlong_of_detections_2 = [latlong_of_detections_2;A(zx,:),converted_antenna_2_matrix(z,2)];
            
        end
        
        
    end
        
end

for z=1:numel(holding_matrix_3)
    
    for zx=1:numel(converted_matrix_strava(:,4))
        
        if holding_matrix_3(z) ==  converted_matrix_strava(zx,4)
            latlong_of_detections_3 = [latlong_of_detections_3;A(zx,:),converted_antenna_3_matrix(z,2)];
            
        end
        
        
    end
        
end

for z=1:numel(holding_matrix_4)
    
    for zx=1:numel(converted_matrix_strava(:,4))
        
        if holding_matrix_4(z) ==  converted_matrix_strava(zx,4)
            latlong_of_detections_4 = [latlong_of_detections_4;A(zx,:),converted_antenna_4_matrix(z,2)];
            
        end
        
        
    end
        
end

xyz = lla2ecef([latlong_of_detections(:,1),latlong_of_detections(:,2),latlong_of_detections(:,3)]);
xyz_1 =lla2ecef([latlong_of_detections_2(:,1),latlong_of_detections_2(:,2),latlong_of_detections_2(:,3)]);
xyz_2 =lla2ecef([latlong_of_detections_3(:,1),latlong_of_detections_3(:,2),latlong_of_detections_3(:,3)]);
xyz_3 =lla2ecef([latlong_of_detections_4(:,1),latlong_of_detections_4(:,2),latlong_of_detections_4(:,3)]);

base = lla2ecef([41.520752667,-71.599254667,40.0]);
normal_xyz = xyz - base;
normal_xyz_1 = xyz_1 - base; 
normal_xyz_2 = xyz_2 - base; 
normal_xyz_3 = xyz_3 - base;

figure(1); hold on;
scatter(normal_xyz(:,1),normal_xyz(:,2),30,latlong_of_detections(:,4));
scatter(normal_xyz_1(:,1),normal_xyz_1(:,2),30,latlong_of_detections_2(:,4));
scatter(normal_xyz_2(:,1),normal_xyz_2(:,2),30,latlong_of_detections_3(:,4));
scatter(normal_xyz_3(:,1),normal_xyz_3(:,2),30,latlong_of_detections_4(:,4));

x =  normal_xyz(:,1);
y =  normal_xyz(:,2);
z = latlong_of_detections(:,3) - 36;

[azimuth,elevation,r] = cart2sph(x,y,z); %note all angles are in radians

x =  normal_xyz_1(:,1);
y =  normal_xyz_1(:,2);
z = latlong_of_detections_2(:,3) - 36;

[azimuth1,elevation1,r1] = cart2sph(x,y,z); %note all angles are in radians

x =  normal_xyz_2(:,1);
y =  normal_xyz_2(:,2);
z = latlong_of_detections_3(:,3) - 36;

[azimuth2,elevation2,r2] = cart2sph(x,y,z); %note all angles are in radians


x =  normal_xyz_3(:,1);
y =  normal_xyz_3(:,2);
z = latlong_of_detections_4(:,3) - 36;

[azimuth3,elevation3,r3] = cart2sph(x,y,z); %note all angles are in radians

azimuth = azimuth;
azimuth1 = azimuth1 - deg2rad(90);
azimuth2 = azimuth2 + deg2rad(180);
azimuth3 = azimuth3 + deg2rad(90);

total_detections_rssi = [latlong_of_detections(:,4);latlong_of_detections_2(:,4);latlong_of_detections_3(:,4);latlong_of_detections_4(:,4)];
combined_azimuth = [azimuth;azimuth1;azimuth2;azimuth3];
combined_radius = [r;r1;r2;r3];
combined_elevation = [elevation;elevation1;elevation2;elevation3];
hold off;
figure(2)
% test_matrix = [rad2deg(azimuth),normalize(latlong_of_detections_3(:,4))]
%test_matrix = [rad2deg(azimuth),normalize(r)]
k = find(total_detections_rssi(:,1) <= -105);
%norm_rssi = normalize(latlong_of_detections(:,4));
%rssi = latlong_of_detections_3(:,4)
test = [combined_azimuth(k),combined_radius(k)];
test = sortrows(test)
%norm_radius = normalize(r);
polarscatter(test(:,1),test(:,2),30);
%f = fit(test(:,1),test(:,2),'sin3');
%a = linspace(-pi,pi,100);
%hold on;
%polarplot(a,f(a))

% color = strcat('#', dec2hex(floor(rand * 2^16), 6));
% dt = delaunayTriangulation(x,z,y);
% tetramesh(dt, 'FaceColor', color);



%  for alph = min(finalSolution(:,4)):max(finalSolution(:,4))
%      holding_matrix = [];
%      for beta = 1:numel(finalSolution(:,4))
%          
%          if finalSolution(beta,4) == -53
%              holding_matrix = [holding_matrix;x(beta),y(beta),z(beta)];
%          end
%      end
%      try
%          color = strcat('#', dec2hex(floor(rand * 2^16), 6));
%          dt = delaunayTriangulation(holding_matrix(:,1),holding_matrix(:,2),holding_matrix(:,3));
%          tetramesh(dt, 'FaceColor', color);
%          hold on;
%      catch
%          warning('Triagnualtion Error')
%      end
%  end


%finalSolution = [finalSolution(mod(finalSolution(:,1)+160,360)<40,:)]
% 
% clf;
%  f2 = fit(finalSolution(:,1),finalSolution(:,4),'exp2')
%  f1 = fit(finalSolution(:,1),finalSolution(:,4),'exp1')
%  plot(f2,'green');
%  hold on;
%  scatter(finalSolution(:,1),finalSolution(:,4));
%  title('Signal Strength vs. Distance Plot')
%  ylabel('dbm')
%  xlabel('Distance (m)')
% % legend('First Order Exponetial', 'Raw Data','Second Order Exponetial')
% % 
% 
% % C = abs(latlong_of_detections(:,4));

figure(3); 
 geoscatter(latlong_of_detections(:,1),latlong_of_detections(:,2), 30, latlong_of_detections(:,4),'filled');
 hold on;
 geoscatter(latlong_of_detections_2(:,1),latlong_of_detections_2(:,2), 30, latlong_of_detections_2(:,4),'filled');
 geoscatter(latlong_of_detections_3(:,1),latlong_of_detections_3(:,2), 30, latlong_of_detections_3(:,4),'filled');
 geoscatter(latlong_of_detections_4(:,1),latlong_of_detections_4(:,2), 30, latlong_of_detections_4(:,4),'filled');

 geobasemap satellite

%plot(1:numel(z),z) plot for alt
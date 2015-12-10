function [] = tunnelPlot()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TUNNELPLOT Summary:                                                      %
%   Function for plotting the grids based on the tunnel geometry.          %
%                                                                          %
% Geometry:                                                                %
%         ______________________________________________________________   %
%      ))                                                               |  %
% Shock)) <---------------------- Tunnel Length ----------------------->|  %
%      )) ______________________________________________________________|  %
%                                                                          %
%                                                                          %
% EECS587-TunnelShock          Dec. 05, 2015                   Matt Graham %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Import data
fID = fopen('X.txt');
H = textscan(fID,'%s %f %s %f',1);
cols = H{2};
rows = H{4};
X = zeros(rows,cols);
for j=1:rows    
    A = textscan(fID,'%f',cols);
    X(j,:) = A{:}';
end
fclose(fID);
fID = fopen('Y.txt');
Y = zeros(rows,cols);
for j=1:rows
    A = textscan(fID,'%f',cols);
    Y(j,:) = A{:}';
end
fclose(fID);
fID = fopen('RHO.txt');
rho = zeros(rows,cols);
for j=1:rows
    A = textscan(fID,'%f',cols);
    rho(j,:) = A{:}';
end
fclose(fID);
% Refined Cells
fID = fopen('refineX.txt');
H = textscan(fID,'%s %f %s %f',1);
cols = H{2};
rows = H{4};
refineX = zeros(rows,cols);
for j=1:rows    
    A = textscan(fID,'%f',cols);
    refineX(j,:) = A{:}';
end
fclose(fID);
fID = fopen('refineY.txt');
refineY = zeros(rows,cols);
for j=1:rows
    A = textscan(fID,'%f',cols);
    refineY(j,:) = A{:}';
end
fclose(fID);
fID = fopen('refineRHO.txt');
refinerho = zeros(rows,cols);
for j=1:rows
    A = textscan(fID,'%f',cols);
    refinerho(j,:) = A{:}';
end
fclose(fID);

% Plotting
figure(1)
subplot(2,1,1)
contourf(X,Y,rho)
xlabel('X')
ylabel('Y')
colorbar
subplot(2,1,2)
meshc(X,Y,rho)
xlabel('X')
ylabel('Y')
zlabel('Pressure')
% suptitle(['t = ' num2str(dt*n)])
% print('-depsc',filename)

figure(2)
subplot(2,1,1)
contourf(refineX,refineY,refinerho)
xlabel('X')
ylabel('Y')
colorbar
subplot(2,1,2)
meshc(refineX,refineY,refinerho)
xlabel('X')
ylabel('Y')
zlabel('Pressure')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Abaqus simulation  data post-processing

%       FFT analysis for the displacement field
%                                 Sifan Yin
%                             02/01/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all;close all;
data_name = 'SY_TR035RR075MR7_OR54_Hex_implicit';Radius = 5.4;
rpt =[data_name,'.rpt'];
data = [data_name,'_data.txt']
filename_original_png = [data_name,'_RealSpace_0202.png']; 
filename_Fspace_png = [data_name,'_FourierSpace_0202.png']; 
finp = fopen(rpt,'r+');
fop = fopen(data,'w+');
k = 0;

while 1
    k = k+1;
    line{k} = fgetl(finp);
    if strncmp(line{k},'    PART-1-1',12)
        line{k};
        fprintf(fop,'%s\r\n',line{k});
    end
    if strcmp(line{k},'  Part Instance     Node ID  Attached elements  U, Magnitude')
        wow = 1
        break
    end
end
fclose('all');
 


filename = data;
formatSpec = '%*24s%14f%14f%14f%14f%14f%f%[^\n\r]';
fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
DATA = table(dataArray{1:end-1}, 'VariableNames', {'VarName3','VarName4','VarName5','VarName6','VarName7','VarName8'});
M = table2array(DATA);
clearvars filename formatSpec fileID dataArray ans DATA;
close all;

node_num = length(M)        
node_X = M(:,1);  node_Y = M(:,2);  node_Z = M(:,3);
node_x = M(:,4);  node_y = M(:,5);   node_z = M(:,6);
node_R1= sqrt(node_X(1)^2+node_Y(1)^2)
node_R2= sqrt(node_X(1000)^2+node_Y(1000)^2)


for ii = 1:length(M)
    X = node_X(ii);
    Y = node_Y(ii);
    Z = node_Z(ii);
    x = node_x(ii); 
    y = node_y(ii);
    z = node_z(ii);
    u_x = x-X;
    u_y = y-Y;
    u_z(ii) = z-Z;
    complex = X+1i*Y;
    theta = angle(complex);
    Rot_matrix = [cos(theta) sin(theta);-sin(theta) cos(theta)];
    
    node_Theta(ii) = Radius*theta;
    u_rtheta = Rot_matrix*[u_x;u_y];
    u_r(ii) = u_rtheta(1);
    u_theta(ii) = u_rtheta(2);
end

node_Theta = node_Theta';
u_r = u_r';
u_theta = u_theta';
u_z = u_z';
node_cylinder = [node_Z node_Theta u_r u_theta u_z];
node_cylinder_sort = sortrows(node_cylinder);
Z_node_sort = node_cylinder_sort(:,1);
Theta_node_sort = node_cylinder_sort(:,2);
Ur_sort = node_cylinder_sort(:,3);
Utheta_sort = node_cylinder_sort(:,4);
Uz_sort = node_cylinder_sort(:,5);
[Z_value, start, floor_row] = unique(node_cylinder_sort(:,1),'rows'); 
 Z_floor = length(Z_value)
floor_num = node_num/Z_floor 
mesh_row = Z_floor;
mesh_colume = floor_num;
Length_Z = max(Z_value)
Z_grid =kron(Z_value,ones(1,mesh_colume));   
Theta_grid = reshape(Theta_node_sort,mesh_colume,mesh_row);
Theta_grid = Theta_grid';
Length_th = max(Theta_node_sort)
Ur_grid = reshape(Ur_sort,mesh_colume,mesh_row);
Ur_grid = Ur_grid';
Utheta_grid = reshape(Utheta_sort,mesh_colume,mesh_row);
Utheta_grid = Utheta_grid';
Uz_grid = reshape(Uz_sort,mesh_colume,mesh_row);
Uz_grid = Uz_grid';
U_ampl = sqrt(Ur_grid.^2+Utheta_grid.^2+Uz_grid.^2);

fig1 = figure(1)
set(fig1,'color',[1 1 1]); 
set(gcf, 'position',[50 50 800 600]);
Theta_grid_double = [Theta_grid(1:end-1,:);Theta_grid];      size(Theta_grid_double)
Z_grid_double = [Z_grid(1:end-1,:)-Length_Z;Z_grid];          size(Z_grid_double)
Ur_grid_double = [Ur_grid(1:end-1,:);Ur_grid];                     size(Ur_grid_double)
U_ampl_double = [U_ampl(1:end-1,:);U_ampl];    
surf(Theta_grid_double,Z_grid_double,Ur_grid_double,'edgecolor','none')
xlabel('\theta\cdot R_i' )
ylabel('z')
colormap;
c = colorbar
c.Label.String = 'Displacement';
shading interp
axis equal
title('U_r in real space')
view(90,90)
set(gca,'Fontsize',22,'Fontname','Times New Roman')
print(filename_original_png,'-dpng')

N_th = 128; L_th = Length_th*2;
N_z = 128; L_z = Length_Z*2;
Theta_interp = L_th/N_th*[-N_th/2:N_th/2-1];
Z_interp = L_z/N_z*[-N_z/2:N_z/2-1];
k_th = 2*pi/L_th*[0:N_th/2-1 -N_th/2:-1];
k_z =  2*pi/L_z*[ 0:N_z/2-1  -N_z/2:-1];
[Kth,Kz] = meshgrid(k_th, k_z);
[Theta_grid_std, Z_grid_std] = meshgrid(Theta_interp, Z_interp);
Ur_grid_std = griddata(Theta_grid_double, Z_grid_double, Ur_grid_double, Theta_grid_std, Z_grid_std,'cubic');
 

fig3 = figure(3)
set(fig3,'color',[1 1 1]); 
set(gcf, 'position',[50 50 700 600]);
Urh_real= abs(real(fft2(Ur_grid_std)/(N_th*N_z)));
Urh_real(1) = 0;
surf(Kth,Kz,Urh_real,'edgecolor','none')
xlabel('k_\theta')
ylabel('k_z')
xlim([-10 10])
ylim([-10 10])
colormap;
c = colorbar
c.Label.String = 'fft(U_r)';
shading interp
title('U_r in Fourier space')
axis equal
view(0,90)
set(gca,'Fontsize',22,'Fontname','Times New Roman')
print(filename_Fspace_png,'-dpng')









%% gut tensile test analysis
% hasreet gill, 2018
% adapted from nanstron_gut_compress by nandan nerurkar

% set scale (mm per pixel)
% choice = questdlg('set scale?');
% switch choice
%     case 'Yes'
%         msgbox('choose scale photo');
%         [filename, pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif;*.lsm','All Image Files';...
%           '*.*','All Files' },'mytitle');
%         cd(pathname);
%         im = imread(filename);
%         figure, imshow(im);
%         h = imdistline;
%         pause
%         scale = 1/(getDistance(h)); 
%     case 'No'
%         scale = inputdlg('enter number of pixels per mm');
%         scale = 1/(str2num(scale{1}));      
% end
clear
npix = 207; % 10.0=526 8.0=426, 6.3=326
scale = 1/npix;

% set constants
% determine lever velocity from calibration video
choice = questdlg('set step size?');
switch choice
    case 'Yes'
        msgbox('choose calibration video');
        disp = levercal(scale);
        for x = 1:length(disp)-1; 
           diff(x) = disp(x+1) - disp(x);       
        end
        stepsize = mean(nonzeros(diff));
    case 'No'
        stepsize = 0.10;
end


% stepsize = 0.0133;
E = 4.1e11; % elastic modulus of tungsten (Pa)
r = 0.10/2; % radius of lever (mm)
% e14.5 mgfp hg - endmes = 0.05/2; mesmus, whole = 0.1/2; 
% e14.5 mgfp mg = 0.05/2 

I = (pi*(r.^4))/4; % formula for area moment of inertia of of a circle (mm^4)
L = 43; % length of lever (mm)
% e14.5 mgfp hg - endmes = 45; mesmus, whole = ?; 
% e14.5 mgfp mg = 45

kb = (3*E*I)/(L.^3); % Bending stiffness of the rod = 3EI/L^3 (Pa*mm)
ri = 0.068; % inner radius (whole tube) - e8 hind = 0.041, e8 fore = 0.068
ro = 0.128; % outer radius (whole tube) - e8 hind = 0.146, e8 fore = 0.128
% l = 0.22;
% msgbox('measure tissue width')
% [filename, pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif;*.lsm','All Image Files';...
%   '*.*','All Files' },'mytitle');
% im2 = imread(filename);
% figure, imshow(im2);
% w1 = imdistline;
% w2 = imdistline;
% pause
% w1 = round(getDistance(w1));
% w2 = round(getDistance(w2));
% w = mean([w1,w2]);
% scale2 = 0.13495/512;
% w = (w*scale2); % Width of one arm of tissue ring (mm)

% w = 0.02; 
% e17 - 0.16; 
% e14.5 mgfp hg - endmes = 0.15; mesmus = 0.25; whole = 0.29 
% e14.5 mgfp mg - end = 0.015; endmes = 0.08

% perform tensile test
% nvids = inputdlg('how many videos to analyze?');
% nvids = str2num(nvids{1});

%[name, datastruct] = tensiletest(scale, stepsize, l, kb, w, nvids);
% [name, datastruct] = tensiletest(scale, stepsize, l, kb, w);
[name, datastruct] = tensiletest(scale, stepsize, kb, ri, ro);
% save(['datastruct_',name,'.mat'],'datastruct')

%% calculate elastic moduli
% 
% %%%%%% if no datastruct
% nvids = inputdlg('how many samples?');
% nvids = str2num(nvids{1});
% 
% [filename, pathname] = uigetfile('*.xls');
% name = strsplit(pathname,filesep);
% name = name{end-1};
% [~,sheet_name]=xlsfinfo(filename);
% for i = 1:nvids;
%     alldata = xlsread(filename,sheet_name{i});
%     datastruct.avgstrain{i} = alldata(4,:);
%     datastruct.stress{i} = alldata(6,:);
% end
% 
% %%%%%%

% stravg = 0.2;

% plot and save results
close all
% figure;
% for k = 1:nvids
%     plot((datastruct.avgnormdist{k}*1000),datastruct.force{k})
%     hold on
% end
% xlabel('distance (um)')
% ylabel('force (uN)')
% legend('1','2','3','4','5','6','7','8')

figure,

%     [value, idx] = min(abs(datastruct.avgstrain-stravg));
%     a = datastruct.stress(idx) - datastruct.stress(1);
%     b = datastruct.avgstrain(idx) - datastruct.avgstrain(1);
%     emod = a/b;
%     smod = emod./(2*(1+0.5));
%     
%     plot(datastruct.avgstrain,datastruct.stress)
    plot(datastruct.time,datastruct.stress)
% % for k = 1:nvids
%     
%     [value, idx] = min(abs(datastruct.avgstrain{k}-stravg));
%     a = datastruct.stress{k}(idx) - datastruct.stress{k}(1);
%     b = datastruct.avgstrain{k}(idx) - datastruct.avgstrain{k}(1);
%     emod(k) = a/b;
%     smod(k) = emod(k)./(2*(1+0.5));
%     
%     plot(datastruct.avgstrain{k},datastruct.stress{k})
%     hold on
% end

% xlabel('average strain')
xlabel('time')
ylabel('stress (Pa)')
% legend('1','2','3','4','5','6','7','8')
saveas(gcf,[name, '_stressvtime.png'])

% titles = {'young modulus';'shear modulus'};
% xlswrite([name,'_moduli','.xls'],titles,1,'A1')
% xlswrite([name,'_moduli','.xls'],emod,1,'B1')
% xlswrite([name,'_moduli','.xls'],smod,1,'B2')

figure,
plot(datastruct.time,datastruct.avgstrain)
xlabel('time')
ylabel('average strain')
saveas(gcf,[name, '_strainvtime.png'])




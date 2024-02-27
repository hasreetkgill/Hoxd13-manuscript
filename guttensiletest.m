%% gut tensile test analysis
% hasreet gill, 2018
% adapted from nanstron_gut_compress by nandan nerurkar

% set scale (mm per pixel)
choice = questdlg('set scale?');
switch choice
     case 'Yes'
         msgbox('choose scale photo'); % use photo of ruler taken at mag corresponding to test videos
        [filename, pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif;*.lsm','All Image Files';...
           '*.*','All Files' },'mytitle');
        cd(pathname);
        im = imread(filename);
        figure, imshow(im);
        h = imdistline; % draw line that measures 1 mm
        pause
        scale = 1/(getDistance(h)); 
    case 'No'
        scale = inputdlg('enter number of pixels per mm');
        scale = 1/(str2num(scale{1}));      
end

% use the following if not measuring scale from an image
npix = 207; % known pixels/mm for relevant magnifications: 10.0=526 8.0=426, 6.3=326
scale = 1/npix;

% set constants
% determine lever velocity from calibration video (lever moving alone without sample)
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


% if not calculating stepsize, enter manually
stepsize = 0.0133;

% other constants
E = 4.1e11; % elastic modulus of tungsten (Pa)
r = 0.10/2; % radius of lever (mm)
I = (pi*(r.^4))/4; % formula for area moment of inertia of of a circle (mm^4)
L = 43; % length of lever (mm)
kb = (3*E*I)/(L.^3); % Bending stiffness of the rod = 3EI/L^3 (Pa*mm)
ri = 0.068; % inner radius (whole tube) - e8 hind = 0.041, e8 fore = 0.068
ro = 0.128; % outer radius (whole tube) - e8 hind = 0.146, e8 fore = 0.128
l = 0.22; % tissue ring thickness (long axis)

% measure tissue width from transverse section images - take two measurements
msgbox('measure tissue width')
[filename, pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif;*.lsm','All Image Files';...
  '*.*','All Files' },'mytitle');
im2 = imread(filename);
figure, imshow(im2);
w1 = imdistline;
w2 = imdistline;
pause
w1 = round(getDistance(w1));
w2 = round(getDistance(w2));
w = mean([w1,w2]);
scale2 = 0.13495/512; % scale of transverse section image pixels/mm
w = (w*scale2); % Width of one arm of tissue ring (mm)

% if not measuring tissue width, enter manually
w = 0.02;

% perform tensile test

% if analyzing multiple videos
% nvids = inputdlg('how many videos to analyze?');
% nvids = str2num(nvids{1});
% [name, datastruct] = tensiletest(scale, stepsize, l, kb, w, nvids); % multiple videos

% one video
[name, datastruct] = tensiletest(scale, stepsize, l, kb, w); % one video
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

stravg = 0.2; % average strain used to calculate elastic modulus

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
% legend('1','2','3','4','5','6','7','8') % multiple videos
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




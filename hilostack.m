%% Initial setup
clear;

bar = waitbar(0,'Define folders and parameters');

% Define image folder and names
[field_corr,dir_corr] = uigetfile(...
    '*.tif',...
    'Please select flat field illumination correction file');
[fnames_uni,indir_uni] = uigetfile(...
    '*.tif',...
    'Please select ALL UNIFORM files to process',...
    'MultiSelect','on');
[fnames_spe,indir_spe] = uigetfile(...
    '*.tif',...
    'Please select ALL SPECKLE files to process',...
    'MultiSelect','on');
outdir = uigetdir('Please select output directory for HiLo images');
opts.Interpreter = 'tex';
opts.Resize = 'on';
opts.WindowStyle = 'normal';
input_parameters = ...                  % Set optical sectioning thickness
    str2double(...                      % parameter and low frequency
    inputdlg(...                        % scaling factor in dialog box
    {'Optical sectioning parameter',...
    'Low frequency scaling factor',...
    'Lateral resolution in px/\mum'},...
    'Define process parameters',...
    1,{'1','1','4.456'},opts));
%% Define parameters 

tic

time = {'Create filter',...
    'Pre-allocation',...
    'Reading images and median filtering',...
    'Bandpass filtering',...
    'Contrast evaluation',...
    'Filtering and scaling',...
    'Reconstruction and Saving',...
    'Total';...
    0,0,0,0,0,0,0,0};

stck_sz = numel(fnames_uni);            % Number of slices
info = imfinfo([indir_spe,'\',...
    char(fnames_spe(1))]);              % Getting image properties
res = input_parameters(3)*1e6*0.0254;   % Lateral resolution
w = info.Width;                         % get image width
h = info.Height;                        % get image height
sigmaBP = w/(10*input_parameters(1));	% band pass filter frequency
kc = nearest(sigmaBP*0.18);             % cut-off frequency between hp and
                                        % lp filter
sigmaLP = kc*2/2.355;                   % Finding sigma value for low pass
                                        % and high pass filter
lambda = nearest(w/(2*kc));             % Side length of contrast
                                        % evalutation mask by rounding to
                                        % nearest odd integer
if mod(lambda,2) == 0                   % lambda must be odd, this is to
                                        % check if lambda is even and add
                                        % 1 to make it odd
    lambda = lambda+1;
else
end
if field_corr == 0                      % Check if correction file is selected
    gauss_corr = 1;
else
    gauss_corr = single(imread([dir_corr,field_corr]));	% Create gaussian field illumination correction matrix if file was selected
    gauss_corr = gpuArray(gauss_corr);
    gauss_corr = gauss_corr+0.2*max(gauss_corr(:));
    gauss_corr = gauss_corr/max(gauss_corr(:));
end
nh = gpuArray(ones(lambda));            % create kernel for local contrast
h = h+2*lambda;                         % increase image size by lambda for
w = w+2*lambda;                         % padding
Nk = sum(nh(:));                        % evalutation
%% Creating filters, pre-allocating and reading image stacks

% Create band pass, high and low pass filters
waitbar(1/10,bar,'Creating filters')

lp = lpgauss(h,w,sigmaLP);
hp = hpgauss(h,w,sigmaLP);
bp = bpgauss(h,w,sigmaBP);
bp = bp/max(bp(:));

time{2,1} = round(toc/60,2);
tic

% Preallocate stacks
waitbar(2/10,bar,'Pre-allocating stacks')

uni = zeros(h,w,stck_sz,'single');
dif = zeros(h,w,stck_sz,'single');
weight = zeros(h,w,stck_sz,'single');

time{2,2} = round(toc/60,2);
tic

for    i=1:stck_sz
    
    % Read images into stack
    waitbar(3/10,bar,['Reading image pair ' int2str(i) ' of ' int2str(stck_sz)])
    
    u = medfilt2(...
         padarray(...
          gpuArray(single(imread([indir_uni,'\',char(fnames_uni(i))])))./gauss_corr,...
         [lambda lambda],'symmetric'),...
        [3 3]);
    s = medfilt2(...
         padarray(...
          gpuArray(single(imread([indir_spe,'\',char(fnames_spe(i))])))./gauss_corr,...
         [lambda lambda],'symmetric'),...
        [3 3]);

    uni(:,:,i) = gather(u);
    dif(:,:,i) = gather(medfilt2(s-u,[5 5]));
end
clear s u gauss_corr

time{2,3} = round(toc/60,2);
tic
%% Bandpass filtering difference image stack

waitbar(4/10,bar,'Bandpass filtering difference image stack')

dif_bp = real(ifft2(fft2(dif).*bp));
dif_min = min(dif_bp(:));
dif_bp = dif_bp-dif_min;

time{2,4} = round(toc/60,2);
tic
%% Contrast evaluation

waitbar(5/10,bar,'Contrast evaluation')

for k = 1:stck_sz
    
    waitbar(5/10,bar,['Contrast evaluation ',int2str(k),' of ',int2str(stck_sz)])

    d = gpuArray(dif_bp(:,:,k));
    mu = imfilter(d,nh,'same')/Nk;
    sq = imfilter(d.^2,nh,'same');
    sd = real(sqrt(complex((sq-Nk*mu.^2)/(Nk-1))));
    weight(:,:,k) = gather(sd./mu);


end

weight = weight-min(weight(:));
intermediate_weight = weight(lambda+1:end-lambda,lambda+1:end-lambda,:);
weight = weight/max(intermediate_weight(:));

waitbar(6/10,bar,'Weighting')

Lo = fft2(weight.*uni);
time{2,5} = round(toc/60,2);
tic
%% Filtering

waitbar(7/10,bar,'Filtering and low frequency scaling')

Hi = real(ifft2(fft2(uni).*hp));
Hi = Hi(lambda+1:end-lambda,lambda+1:end-lambda,:);
Lo = real(ifft2(Lo.*lp));
Lo = Lo(lambda+1:end-lambda,lambda+1:end-lambda,:);
%% Scaling

nabla = abs(fft2(Hi(1,kc+1,:)))./(abs(fft2(Lo(1,kc+1,:))));
nabla_mu = mean(nabla);
nabla_med = median(nabla);
nabla_min = min(nabla);
LoWeight = num2cell(num2str(input_parameters(2)*nabla_min));
LoWeight = str2double(inputdlg({['Minimum ratio: ' num2str(input_parameters(2)*nabla_min,3)...
    ', Median ratio: ' num2str(input_parameters(2)*nabla_med,3)...
    ', Mean ratio: ' num2str(input_parameters(2)*nabla_mu,3)]},...
    'Choose final low scaling factor',1,LoWeight));
Lo = Lo*LoWeight;
time{2,6} = round(toc/60,2);
tic
%% Reconstruction and 16bit conversion

waitbar(8/10,bar,'Reconstructing and converting to 16bit')

HiLo = Hi+Lo;
HiLo = HiLo/max(HiLo(:));
HiLo = im2uint16(HiLo);
mkdir(outdir,strcat('HiLo_',int2str(input_parameters(1)),'_',num2str(LoWeight,2)));
%% Saving

for j=1:stck_sz
    
    waitbar(9/10,bar,['Saving ' int2str(j) ' of ' int2str(stck_sz)])
    
    imwrite(HiLo(:,:,j),[outdir,...
        '\HiLo_',...
        int2str(input_parameters(1)),'_',num2str(LoWeight,2),'\HiLo',int2str(j),'.tif'],...
        'Resolution',res,'Compression','none')
    
end
time{2,7} = round(toc/60,2);
time{2,8} = round(sum(cell2mat(time(2,1:7))),2);
waitbar(1,bar,['DONE!!! ','Total time taken: ','~',num2str(round(cell2mat(time(2,8)))),' Minutes'])
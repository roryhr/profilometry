function im_analyze5() 
% im_analyze4.m
% clear
% -------------- Program Variables --------------------------------%
im_dir = 'Trial 32_8\';   % Folder where images or video are stored
% im_dir = 'Nov 6\';
% im_list = 'All';        % Images or frames to read. 1st image sets reference plane
im_list = [10:30];  
plots_on = 2;           	% Generate diagnostic plots

% ----------- Experimental Parameters -----------------------------%
L_0 = (28+4)*25.4;      % Distance from object to camera & projector: mm
D = 4.75*25.4;             % Distance from camera to projector: mm
% f_0 = 13.5/25.4;        % Fundamental frequency of observed grating, gratings/mm
% pixel_pitch = 190/25.4;	% pixels / mm
pixel_pitch = 6.59;         % Trial 30 -- pixels / mm
% win_length = 40;        % = 1/standard deviation: pixels^-1


% -------- Data analysis variables -----------------------------%
ROI = [200, 90, 1050, 990];	% X Start (1st column), Y start (1st row),
%                          	% X End (last column), Y end (last row)
% ROI = 'off';
rotate_on = 0;      % 0 if gratings are horizontal (default), 1 rotates images
high_pass_wid = 80;         % Width of high pass gaussian filter
gauss_filter_ratio = 1/3;   % Ratio of standard deviation to wavenumber

% ------------ Read images ------------------------------%
files = dir(im_dir);
if strcmpi(files(3).name(end-2:end),'avi') 
    if ischar(im_list)       
        im_list = 1:length(files);      % Create a list of all frames in video to import
    end
      
    num_images = length(im_list);
    readerobj = mmreader(strcat(im_dir, files(3).name));
    im.name = readerobj.Name; 
    for i = 1:num_images
        temp = read(readerobj,im_list(i));  % Read frame  
        im(i).image = temp(:,:,1);      % for B&W images only use R of RGB
    end
else 
    if ischar(im_list)       
        im_list = 1:(length(files)-2);	% Create a list of all images in directory
    end
    num_images = length(im_list);
    for i = 1:num_images
        im(i).image = imread(strcat(im_dir, files(im_list(i)+2).name));  
        im(i).name = files(im_list(i)+2).name;
    end
    
end

% ------------- Region of Interest ------------------------%
if plots_on >= 1
    subplot(1,2,1)
    imshow(im(1).image, [0 mean(max(im(1).image))])
    title('Reference image')
end
if ~ischar(ROI(1))
    if plots_on
        hold on
        plot([ROI(1),ROI(3),ROI(3),ROI(1),ROI(1)],[ROI(2),ROI(2),...
            ROI(4),ROI(4),ROI(2)], 'Color','red', 'LineWidth',1)  % Rectangle over ROI
        hold off
    end

    % Discard subset of image not in ROI    
    for i = 1:num_images
        im(i).image = im(i).image(ROI(2):ROI(4),ROI(1):ROI(3)); 
    end
end

if rotate_on
    for i = 1:num_images
        im(i).image = im(i).image';
    end
end
    

% ----------- Analysis & processing ------------------------%
[N c] = size(im(1).image);      % N is the length of the column vectors = fft length
% NFFT = 2^nextpow2(N);           % Next power of 2
mid_x = floor(c/2);
mid_y = floor(N/2);
x_p = 0:c-1;                    % X Pixels
y_p = 0:N-1;                    % Y Pixels

y_range = floor(.45*N:.55*N);
mich_contrast = double(max(im(1).image(y_range,mid_x)) - min(im(1).image(y_range,mid_x)))...
    / double(max(im(1).image(y_range,mid_x)) + min(im(1).image(y_range,mid_x)));


if plots_on >= 1
    subplot(1,2,1)
    hold on
    plot([ROI(1)+mid_x,ROI(1)+mid_x],[ROI(2)+y_range(1),ROI(2)+y_range(end)],...
        'Color','blue', 'LineWidth',1); 
    hold off
    
    subplot(1,2,2);
    plot(y_range,im(1).image(y_range,mid_x))
   	title(sprintf('Michelson Contrast = %1.2f', mich_contrast))
    xlabel('Y (pixels)')
    ylabel('Intensity')
%     xlim([1 2])
end

% % Remove average
% im_1 = im_1 - mean(mean(im_1));
% im_2 = im_2 - mean(mean(im_2));

% ---------- Take Fourier Transform -----------
filt = fspecial('disk',3);      % Initialize averaging filter

for i = 1:num_images
%  	im(i).image = imfilter(im(i).image, filt);  % Filter image first 
    im(i).fft = fft(imfilter(double(im(i).image),filt));
%     im(i).image = fft(im(i).image);     % Fourier transform columns
end


% --------- Plot power spectrum ---------------
f0 = (-N/2:N/2-1);

% % Plot single-sided amplitude spectrum.
% f = 1/2*linspace(0,1,NFFT/2+1);
% plot(f,2*abs(im(1).fft(1:NFFT/2+1,mid_x))) 
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')

if plots_on >= 3
    % 2D Plot
    mesh(x_p,f0,fftshift(abs(im(1).fft)))
end

% ---------- Window (filter) data ----------
% High-pass filter (remove lowest frequencies)
win = 1-gausswin(N,high_pass_wid);
if plots_on >= 2
    plot(fftshift(abs(im(1).fft(:,mid_x))))
    hold on
    plot(win*max(abs(im(1).fft(:,mid_x))))
    hold off
end

for i = 1:num_images
    im(i).fft = ifftshift(fftshift(im(i).fft).*repmat(win,1,c));
end
% im_1 = ifftshift(fftshift(im_1).*repmat(win,1,c)); 
% im_2 = ifftshift(fftshift(im_2).*repmat(win,1,c));

% Identify where to filter using reference image (this happens once)
[C I] = max(abs(im(1).fft(:,mid_x)));
% shift_length = I;
shift_length = I-round(N/2);
win_length = (I-1)*gauss_filter_ratio;       % Standard deviation of window
win = circshift(gausswin(N,(N/2)/(win_length)),shift_length); % Generate and shift gaussian window
win = win/max(win);

if plots_on >= 1
    plot(abs(im(1).fft(:,mid_x))/N)
    hold on
    plot(win*max(abs(im(1).fft(:,mid_x)))/N)
    title(sprintf('Gaussian Filter with Sigma = %1.2f', win_length))
    xlim([0 3*I])
    hold off
end

% Apply window to fft data
for i = 1:num_images
    im(i).fft = im(i).fft.*repmat(win,1,c);
end

% ----- Transform back to real space ------- 
filt = fspecial('disk',5); % Initialize averaging filter

for i = 1:num_images
    im(i).phase = ifft(im(i).fft);
end

f_0 = (I/N)*pixel_pitch;        % Gratings / pixel * pixel/ mm
% ---- Calculate delta_phi (change of phase with respect to reference)---
for i = 2:num_images
    im(i).phase = angle(im(i).phase.*conj(im(1).phase));

    im(i).phase = imfilter(im(i).phase, filt); % Image filter 

    im(i).phase = unwrap(im(i).phase);          % Remove jumps along columns
    im(i).phase = unwrap(im(i).phase')';        % Check for jumps along rows too

    
    im(i).h = L_0 * im(i).phase./(im(i).phase-2*pi*f_0*D); % Compute height from phase

    if rotate_on==1
        im(i).h = im(i).h';  % flip back to original orientation
    end
    
    figure(1)
    mesh(x_p/pixel_pitch, y_p/pixel_pitch, im(i).h)
    title(sprintf('Image %g of %g', i-1, num_images-1))
    xlabel('X (mm)')
    ylabel('Y (mm)')
    zlabel('Z (mm)')
    zlim([-4 4])

    figure(2)
    
    subplot(2,2,1)
    plot(f0,fftshift(abs(im(i).fft(:,mid_x))))
    title('Filtered Frequency')
    subplot(2,2,2)
    plot(y_p,im(i).phase(:,mid_x))
    title('Unwrapped Phase')
    ylabel('Radians')
    subplot(2,2,3)
    plot(y_p/pixel_pitch,im(i).h(:,mid_x))
   	title('Middle X Slice')
    xlabel('Y (mm)')
    subplot(2,2,4)
   	plot(x_p/pixel_pitch,im(i).h(mid_y,:))
    title('Middle Y')
    ylabel('Height (mm)')
    xlabel('X (mm)')
    
end
% zlim([0 .3])
% figure
% plot(y_p/pixel_pitch,im.h(mid_x,:))
end


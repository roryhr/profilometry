function grating_gen()
% grating_gen.m
% Generates a grating image given an image resolution and spatial frequency
full_im_name = '.\Grating Images\test72.png';
N_x = 1080; %Resolution of image in pixels
aspect_ratio = 16/9;
N_y = aspect_ratio*N_x;
f_0 = N_x/12;
%---------------------------------------%


% ----------- Build Matrix -- 1D variation in x -----------------------%
col = 0:(N_x-1);
I_1D = cos(2*pi*col*f_0/N_x)/2.1+.5;
I_2D = ones(N_y,N_x)*diag(I_1D);
% I_2D = I_2D'; 

% imshow(I_2D', [0 2])
% figure
% scatter(1:length(I_1D),I_1D)
% hold on 
% plot(I_1D)
% hold off
% xlim([0,2*N_x/f_0])

% --------- Save Image --------------------------%
imwrite(I_2D',full_im_name,'png','bitdepth',16)       % Transpose for vertical grating

% -------- Open Image and plot profile ---------------------------%
im = imread(full_im_name);
y_range = [1,2*N_y/f_0];
x_range = floor(N_y/2);

imshow(im);
hold on
plot([x_range(1),x_range(1)],[y_range(1),y_range(2)], 'Color','red', 'LineWidth',1)  % Rectangle over R
hold off

figure
scatter(y_range(1):y_range(2),im(y_range(1):y_range(2),x_range))
hold on
plot(y_range(1):y_range(2),im(y_range(1):y_range(2),x_range))
hold off

xlabel('Pixels in Y')
ylabel('Intensity')



end



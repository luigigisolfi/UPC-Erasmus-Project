% Specify the directory where the images are saved
imageDir = '/Users/luigigisolfi/Documents/mcodes/mcodes/jupiter_gif'; % Change this to the directory where your images are saved
outputGif = '/Users/luigigisolfi/Documents/mcodes/mcodes/jupiter_gif/lol.gif'; % Name of the output GIF file

% List all the image files in the directory
imageFiles = dir(fullfile(imageDir, '*.png')); % Adjust the extension if your images are not PNGs

% Sort the files by name to ensure correct order
[~, idx] = sort({imageFiles.date});
imageFiles = imageFiles(idx);

% Create the GIF
for s = 1:length(imageFiles)
    % Read the image
    imageFiles(s)
    img = imread(fullfile(imageDir, imageFiles(s).name));
    % Convert the image to indexed format
    [imind, cm] = rgb2ind(img, 256);
    % Write to the GIF file
    if s == 1
        imwrite(imind, cm, outputGif, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
    else
        imwrite(imind, cm, outputGif, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
end
%% Find true spacing
img = rgb2gray(IMG_0064);
img = double(img);
imagesm=img(2100-160*3:2100+160*3,3300-221*3:3300+221*3);
imagesm = imagesm / max(imagesm(:));
[lx, ly] = size(imagesm);
score = 0;
spacing_x = 0;
for ii = 1:100
    newimg = imtranslate(imagesm, [ii, 0]);
    product = imagesm .* newimg;
    product = sum(product(:)) / (lx * (ly - ii));
    if product > score
        score = product;
        spacing_x = ii;
    end
end
score = 0;
spacing_y = 0;
for jj = 1:100
    newimg = imtranslate(imagesm, [0, jj]);
    product = imagesm .* newimg;
    product = sum(product(:)) / (ly * (lx - jj));
    if product > score
        score = product;
        spacing_y = jj;
    end
end

%% Find real position
diff_score = 10000;
pos = [0, 0];
for xx = 1:spacing_x - 1
    for yy = 1:spacing_y - 1
        mask = ones(lx, ly);
        mask(xx:spacing_x:lx, yy:spacing_y:ly) = 0;
        diff = mean(mean(abs(mask - imagesm)));
        if diff < diff_score
            diff_score = diff;
            pos = [xx, yy];
            real_mask = mask;
        end
    end
end
        

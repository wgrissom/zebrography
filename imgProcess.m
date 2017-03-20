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
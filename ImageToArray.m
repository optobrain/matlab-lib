%% to use imshow but keep X and Y axis

function I = ImageToArray(img)

    [ny,nx,nch] = size(img);  
    I = zeros(nx,ny,nch,class(img));
    for ich=1:nch
        I(:,:,ich) = img(:,:,ich)';
    end
    I = I(:,ny:-1:1,:);

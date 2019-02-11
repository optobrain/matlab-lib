%% to use imshow but keep X and Y axis

function img = ArrayToImage(I)

    [nx,ny,nch] = size(I);  
    img = zeros(ny,nx,nch,class(I));
    for ich=1:nch
        img(:,:,ich) = I(:,:,ich)';
    end
    img = img(ny:-1:1,:,:);

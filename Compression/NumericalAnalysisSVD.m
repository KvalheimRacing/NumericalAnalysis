% Function for using SVD to factorize and then compress images
function [] = NumericalAnalysisSVD(r)


    format shortG % Pretty format
    close all     % Close all previous figures

    fprintf('\nChecking for input validity')
    if ~(floor(r)==r)
        error('r = %d is not an integer. Please provide an integer.\n', r)
        exit
    end


    fprintf('\nConverting images to grey scale')
    im1_grey = rgb2gray(imread('chessboard.png'));
    im2_grey = rgb2gray(imread('jellyfish.jpg'));
    im3_grey = rgb2gray(imread('new_york.jpg'));

    fprintf('\nConverting images to double between 0 and 1')
    im1 = im2double(im1_grey);
    im2 = im2double(im2_grey);
    im3 = im2double(im3_grey);

    fprintf('\nStoring grey scale images in current foler')
    imwrite (im1, 'im1_grey.png');
    imwrite (im2, 'im2_grey.jpg');
    imwrite (im3, 'im3_grey.jpg');

    fprintf('\nGetting info about images before compression')
    chessinfo = imfinfo('im1_grey.png');
    jellyinfo = imfinfo('im2_grey.jpg');
    nyorkinfo = imfinfo('im3_grey.jpg');

    fprintf('\nDecomposing grayscale images with SVD')
    [U1,S1,V1] = svd(im1);
    [U2,S2,V2] = svd(im2);
    [U3,S3,V3] = svd(im3);

    fprintf('\nGetting matrix dimensions')
    [mU1,nU1] = size(U1); % 1280x1280
    [mU2,nU2] = size(U2); % 1280x1280
    [mU3,nU3] = size(U3); % 1280x1280
    [mS1,nS1] = size(S1); % 1280x1236
    [mS2,nS2] = size(S2); % 1280x1920
    [mS3,nS3] = size(S3); % 1280x1920
    [mV1,nV1] = size(V1); % 1236x1236
    [mV2,nV2] = size(V2); % 1920x1920
    [mV3,nV3] = size(V3); % 1920x1920

    % Don't run program if r is too big
    if r > min([mS1,nS1,mS2,nS2,mS3,nS3])
        warning(['With r = %d, you are choosing more singular values ',...
                 'than you have matrix dimensions for. Impossible\n'], r)
        return
    end


    % Let now m=n=r be the constants determening the compression
    fprintf('\nCompressing images as a funcnction of n=m=r')
    im1_compressed = U1(1:mU1,1:r)*S1(1:r,1:r)*V1(1:nV1,1:r)';
    im2_compressed = U2(1:mU2,1:r)*S2(1:r,1:r)*V2(1:nV2,1:r)';
    im3_compressed = U3(1:mU3,1:r)*S3(1:r,1:r)*V3(1:nV3,1:r)';

    fprintf('\nStoring compressed grey scale images in current foler')
    imwrite (im1_compressed, 'im1_grey_compressed.png');
    imwrite (im2_compressed, 'im2_grey_compressed.jpg');
    imwrite (im3_compressed, 'im3_grey_compressed.jpg');

    fprintf('\nGetting info about picures after compression')
    chessinfo_compressed = imfinfo('im1_grey_compressed.png');
    jellyinfo_compressed = imfinfo('im2_grey_compressed.jpg');
    nyorkinfo_compressed = imfinfo('im3_grey_compressed.jpg');

    fprintf('\nThen compression ratios are given as\n')
    fprintf('   Compression ratio for chessboard image is %f, so we have compressed by %f%%\n',...
                                         (chessinfo.FileSize/chessinfo_compressed.FileSize),...
                      (100.0-((1.0/(chessinfo.FileSize/chessinfo_compressed.FileSize))*100.0)))
    fprintf('   Compression ratio for jelly fish image is %f, so we have compressed by %f%%\n',...
                                         (jellyinfo.FileSize/jellyinfo_compressed.FileSize),...
                            (100-((1/(jellyinfo.FileSize/jellyinfo_compressed.FileSize))*100)))
    fprintf('   Compression ratio for New York__ image is %f, so we have compressed by %f%%\n',...
                                         (nyorkinfo.FileSize/nyorkinfo_compressed.FileSize),...
                            (100-((1/(nyorkinfo.FileSize/nyorkinfo_compressed.FileSize))*100)))


    fprintf('Plotting original images\n')
    dataset4window = figure('Name','Awesome Chessboard','NumberTitle','off');
    movegui(dataset4window,'southeast')
    imshow(im1, 'InitialMagnification', 30)
    title(['Chessboard in Grayscale, before compression'])

    dataset5window = figure('Name','Awesome Jellyfish','NumberTitle','off');
    movegui(dataset5window,'southwest')
    imshow(im2, 'InitialMagnification', 25)
    title(['Jellyfish in Grayscale, before compression'])

    dataset6window = figure('Name','New York City','NumberTitle','off');
    movegui(dataset6window,'northwest')
    imshow(im3, 'InitialMagnification', 25)
    title(['New York in Grayscale, before compression'])


    fprintf('Plotting compressed images\n')
    dataset7window = figure('Name','Awesome Compressed Chessboard','NumberTitle','off');
    movegui(dataset7window,'southeast')
    imshow(im1_compressed, 'InitialMagnification', 30)
    title(['Chessboard in Grayscale, after compression'])

    dataset8window = figure('Name','Awesome Compressed Jellyfish','NumberTitle','off');
    movegui(dataset8window,'southwest')
    imshow(im2_compressed, 'InitialMagnification', 25)
    title(['Jellyfish in Grayscale, after compression'])

    dataset9window = figure('Name','Compressed New York City','NumberTitle','off');
    movegui(dataset9window,'northwest')
    imshow(im3_compressed, 'InitialMagnification', 25)
    title(['New York in Grayscale, after compression'])


    fprintf('Plotting singular values\n\n')
    dataset1window = figure('Name','Singular Values of Chessboard image','NumberTitle','off');
    movegui(dataset1window,'south')
    plot(log(svd(im1)))
    title(['Log of ' num2str(nnz(svd(im1))) ' Non Zero Singluar Values in decending order' ])

    dataset2window = figure('Name','Singular Values of Jellyfish image','NumberTitle','off');
    movegui(dataset2window,'north')
    plot(log(svd(im2)))
    title(['Log of ' num2str(nnz(svd(im2))) ' Non Zero Singluar Values in decending order' ])

    dataset3window = figure('Name','Singular Values of New York image','NumberTitle','off');
    movegui(dataset3window,'northeast')
    plot(log(svd(im3)))
    title(['Log of ' num2str(nnz(svd(im3))) ' Non Zero Singluar Values in decending order' ])


    % We fint that the highest reasonable compression for the chessboard image
    % is given by r=2. For the jellyfish image we find that r=70 gives a reasonable image.
    % lastly for the picture of new york, r=250 is a good compression.
    % The high compression rate on the chessboard image is due to the difference in the
    % singular values, beeing a lot that is very close to 0. The high singular
    % values are right at the start, so for the chessboard, two is all we need.


end

function [ newimage ] = ComputeHomography_2015CSB1032( img1, img2 )
    
    imshow(img1);
    title('Please select 4 points in first image for homography');
    img1_coordinates = zeros(4,2);
    for i=1:4
        [x_coord, y_coord] = ginput(1);
        img1_coordinates(i, :) = [x_coord, y_coord];
        hold on;
        plot(img1_coordinates(i, 1), img1_coordinates(i, 2), 'g+', 'MarkerSize', 9);
    end
    hold off;

    imshow(img2);
    title('Please select 4 points in second image for homography');
    img2_coordinates = zeros(4,2);
    for i=1:4
        [x_coord, y_coord] = ginput(1);
        img2_coordinates(i, :) = [x_coord, y_coord];
        hold on;
        plot(img2_coordinates(i, 1), img2_coordinates(i, 2), 'g+', 'MarkerSize', 9);
    end
    hold off;
    
    A = zeros(8,9);
    for i=1:4
        row1 = [img1_coordinates(i, 1) img1_coordinates(i, 2) 1 0 0 0 -(img1_coordinates(i, 1)*img2_coordinates(i, 1)) -(img1_coordinates(i, 2)*img2_coordinates(i, 1)) -(img2_coordinates(i, 1))];
        row2 = [0 0 0 img1_coordinates(i, 1) img1_coordinates(i, 2) 1 -(img1_coordinates(i, 1)*img2_coordinates(i, 2)) -(img1_coordinates(i, 2)*img2_coordinates(i, 2)) -(img2_coordinates(i, 2))];
        A((2*i)-1, :) = row1;
        A(2*i, :) = row2;
    end
    
    % Finding eigenvector of transpose(A)*A with smallest eigenvalue
    A_transpose = A';
    X = A_transpose*A;
    
    [eig_vectors, eig_values] = eig(X);
    [~, min_index] = min(min(eig_values));
    
    H_vertical = eig_vectors(:, min_index);
    H_matrix = reshape(H_vertical, 3, 3);
    H_matrix = inv(H_matrix')';
    projective_transform = projective2d(H_matrix);

    [warped_img, ref_img2] = imwarp(img2, projective_transform);
    imwrite(warped_img, 'Warped_image.jpg');
    
    warped_coordinates = zeros(4,2);
    for i=1:4
        [homography_x, homography_y] = transformPointsForward(projective_transform, img2_coordinates(i,1), img2_coordinates(i,2));
        finalx = homography_x - ref_img2.XWorldLimits(1);
        finaly = homography_y - ref_img2.YWorldLimits(1);
        warped_coordinates(i,:) = [finalx, finaly];
    end
    
    left_image1 = min(img1_coordinates(:,1));
    left_image2 = min(warped_coordinates(:,1));
    left_newimage = max(left_image1, left_image2);

    right_image1 = size(img1, 2) - max(img1_coordinates(:,1));
    right_image2 = size(warped_img, 2) - max(warped_coordinates(:,1));
    right_newimage = max(right_image1, right_image2);
    
    top_image1 = min(img1_coordinates(:,2));
    top_image2 = min(warped_coordinates(:,2));
    top_newimage = max(top_image1, top_image2);

    bottom_image1 = size(img1, 1) - max(img1_coordinates(:,2));
    bottom_image2 = size(warped_img, 1) - max(warped_coordinates(:,2));
    bottom_newimage = max(bottom_image1, bottom_image2);
    
    center_lr = max(warped_coordinates(:,1)) - min(warped_coordinates(:,1));
    center_tb = max(warped_coordinates(:,2)) - min(warped_coordinates(:,2));
    
    newimage_size = [int64(top_newimage+center_tb+bottom_newimage), int64(left_newimage+center_lr+right_newimage)];
    newimage = uint8( zeros(newimage_size(1), newimage_size(2), 3) );
    
    x1 = int64(top_newimage-top_image1+1);
    y1 = int64(left_newimage-left_image1+1);
    x2 = x1+size(img1, 1)-1;
    y2 = y1+size(img1, 2)-1;
            
    xx1 = int64(top_newimage-top_image2+1);
    yy1 = int64(left_newimage-left_image2+1);
    xx2 = xx1+size(warped_img, 1)-1;
    yy2 = yy1+size(warped_img, 2)-1;

    newimage(xx1:xx2, yy1:yy2, :) = warped_img;

    for i=x1:min(x2, size(newimage,2))
        for j=y1: min(y2, size(newimage,1))
            if (newimage(i,j,1)==0 && newimage(i,j,2)==0 && newimage(i,j,3)==0)
                newimage(i,j,:)=img1(i-x1+1,j-y1+1,:);
            else
                newimage(i,j,:) = ( newimage(i,j,:)/2 )+ ( img1(i-x1+1,j-y1+1,:)/2 );
            end
        end
    end
    imwrite(newimage, 'stitched_image.jpg');
    imshow(newimage);
    title('Stitched Image');
    
end
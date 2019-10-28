function [Processed_image]=Pre_processing(Concat_Image,perform_adjustment,perform_histeq,perform_background_removal,perform_decorrelation,perform_smoothing,adjust_param,sigma_param)
%Pre_processing : clean and pre-process the loaded images
%   Several cleaning steps necessary for the analysis : adjustment of
%   brightness and contrast, histogramm equalisation across chanels and
%   top-hat filter to remove unequal illunination

if nargin < 6 
    perform_adjustment = true;
    perform_histeq = false;
    perform_background_removal = true;
    perform_decorrelation = true;
    perform_smoothing = true;

end


if nargin < 8
    adjust_param=0.99999;
    sigma_param= 30;
end

numstack = size(Concat_Image,3);
l = size(Concat_Image,4);


%%Remove non specific signal
if perform_decorrelation
 disp("Perfoming decorrelation stretch");

    for k=1:numstack
        corrected_stacks = decorrstretch(squeeze(Concat_Image(:,:,k,:)));
        Concat_Image(:,:,k,:) = corrected_stacks;
    end
end


%%Performing background removal using Gaussian smoothing
if perform_background_removal
    disp("Perfoming background removal");
    for k = 1:l 
        for i = 1:numstack
        bg_image = imgaussfilt(Concat_Image(:,:,i,k),sigma_param);
        Concat_Image(:,:,i,k) = Concat_Image(:,:,i,k) - bg_image;
        end
    end
end

Concat_Image(Concat_Image<0) = 0;


%%Performing intensity adjustment 

if perform_adjustment
    disp("Perfoming contrast/intensity adjustment");
    for k = 1:l 
        Concat_Image(:,:,:,k)=imadjustn(Concat_Image(:,:,:,k),double([ quantile(reshape(Concat_Image(:,:,:,k),[],1),1-adjust_param) quantile(reshape(Concat_Image(:,:,:,k),[],1),adjust_param)]));
    end
end

%%Performing 3D smoothing using a small sigma value
if perform_smoothing 
   disp("Perfoming 3D image smoothing");
   for k=1:l
   Concat_Image(:,:,:,k) = imgaussfilt3(Concat_Image(:,:,:,k),1);
   end
end

if perform_histeq
    disp("Perfoming histogram equalisation across chanels");
    for k=1:l 
        Concat_Image(:,:,:,k)=imhistmatchn(Concat_Image(:,:,:,k),Concat_Image(:,:,:,1));
    end
end


Processed_image=single(Concat_Image);
disp("Pre-processing finished !");

end


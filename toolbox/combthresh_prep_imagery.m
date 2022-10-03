function [IMavg]  = combthresh_prep_imagery(odir,imageno,Y,imthreshmethod)
% find average value for imagery  to compute threshold from

if imthreshmethod == 1
    I = [];
    for i = 1:length(imageno)
        
        imagefile = fullfile(odir,['c2_',sprintf('%05d',imageno(i)),'.tiff']);
        IM = double(rgb2gray(imread(imagefile)))/255;
        if i==1
            I = IM;
        else
            I = I+IM;
        end
    end
    
    imreg = I(:,70:120,1)/length(imageno);
    imreg(imreg>0.45)=NaN;
    IMavg = squeeze(nanmean(imreg,2));
    IMavg(IMavg<0.22)=0.22;
    IMavg = repmat(IMavg,1,length(xtemp));
elseif imthreshmethod ==  2
    for i = 1:length(imageno)
        imagefile = fullfile(odir,['c2_',sprintf('%05d',imageno(i)),'.tiff']);
        I(:,:,i) = double(rgb2gray(imread(imagefile)))/255;
    end
    
    imagetemp = movmean(I,5*8*2,3,'omitnan');% 10 waves
    [iy(1)] = find(Y(:,1)==-13.05);
    [iy(2)] = find(Y(:,1)==13.05);
    
    for i = 1:size(I,3)
        imagetemp(iy(1):iy(2),:,i) = smoothdata(imagetemp(iy(1):iy(2),:,i),1,'movmean',[6 6],'omitnan');%[5 5]%  just this one!
        imagetemp(iy(1):iy(2),:,i) = smoothdata(imagetemp(iy(1):iy(2),:,i),2,'movmean',[6 6],'omitnan');% [3 3]
    end
    for i = 1:20
        imagetemp(:,:,i) = imagetemp(:,:,20);
        imagetemp(:,:,end-i+1) = imagetemp(:,:,end-20);
    end
    IMavg = imagetemp;
end
end
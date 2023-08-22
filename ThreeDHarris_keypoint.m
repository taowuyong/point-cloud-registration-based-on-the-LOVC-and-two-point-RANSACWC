function [keypoint] = ThreeDHarris_keypoint(PC,Rk)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[n m]=size(PC);
[idx dist]=rangesearch(PC,PC,Rk);
MnormalV=[];
for i=1:n
    KNN=PC(idx{i},:);
    [h1 l1]=size(KNN);
    Cov=(1/h1)*(KNN-ones(h1,1)*mean(KNN))'*(KNN-ones(h1,1)*mean(KNN));
    [U S V]=svd(Cov);
    normalV=V(:,3);
    MnormalV=[MnormalV normalV];
end
keypoint=[];
for i=1:n
    nKNN=MnormalV(:,idx{i});
    [h1 l1]=size(nKNN);
    COV=(1/l1)*(nKNN)*(nKNN)';
    RH=det(COV);
    if RH>0.01                    %parameter
        keypoint=[keypoint;PC(i,:)];
    end
end
end
function [ DV ] = LOVC( KNN,keypoint,V,RR )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[h l]=size(KNN);
KNNt=(KNN-ones(h,1)*keypoint)*V;
wbin=15;              %%%%²ÎÊý
Lbin=2*RR/wbin;
for i=1:wbin
    for j=1:wbin
        for k=1:wbin
            idd=find(KNNt(:,1)>=-RR+(i-1)*Lbin & KNNt(:,1)<-RR+i*Lbin & KNNt(:,2)>=-RR+(j-1)*Lbin & KNNt(:,2)<-RR+j*Lbin & KNNt(:,3)>=-RR+(k-1)*Lbin & KNNt(:,3)<-RR+k*Lbin);
            if isempty(idd)
                D(i,j,k)=0;
            else
                D(i,j,k)=1;
            end
        end
    end
end
[m n l]=size(D);
DV=[];
for i=1:m
    DD=D(:,:,i);
    VF=[];
    for j=1:m
        for k=j+1:m
            if DD(k,j)==0 && DD(j,k)==0
                vf=0;
            end
            if DD(k,j)==0 && DD(j,k)==1
                vf=1;
            end
            if DD(k,j)==1 && DD(j,k)==0
                vf=1;
            end  
            if DD(k,j)==1 && DD(j,k)==1
                vf=0;
            end
            VF=[VF vf];
        end
    end
    DV=[DV VF];
end
end


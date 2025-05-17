%point cloud registration based on the LOVC and 2-point RANSACWC
str='E:\compile document\matlab\data\Indoor and outdoor dataset\scan000-scan001.txt';
fid = fopen(str,'r');
D = textscan(fid, '%f%f%f%f');
fclose(fid);
Ttrue=[D{1} D{2} D{3} D{4}];
pcloud1=pcread('E:\compile document\matlab\data\Indoor and outdoor dataset\scan000.ply');
PC2=pcloud1.Location;
pcloud2=pcread('E:\compile document\matlab\data\Indoor and outdoor dataset\scan001.ply');
PC1=pcloud2.Location;
%     plot3(PC1(:,1),PC1(:,2),PC1(:,3),'.b','MarkerSize',1);
%     hold on;
%     plot3(PC2(:,1),PC2(:,2),PC2(:,3),'.r','MarkerSize',1);
%     set(gca,'DataAspectRatio',[1 1 1]);
%     axis off
%     pr=0.022431135;          %apartment4  apartment6
% pr=0.025664071;          %boardroom0  boardroom3
%     pr=0.12308180;          %City1  City2
    pr=0.25640265;             %scan000  scan001
[n m]=size(PC2);
Rk=5*pr;
[keypoint1] = ThreeDHarris_keypoint(PC1,Rk);
[keypoint2] = ThreeDHarris_keypoint(PC2,Rk);
%     plot3(keypoint1(:,1),keypoint1(:,2),keypoint1(:,3),'.b','MarkerSize',1);
%     hold on;
%     plot3(keypoint2(:,1),keypoint2(:,2),keypoint2(:,3),'.r','MarkerSize',1);
%     set(gca,'DataAspectRatio',[1 1 1]);
%     axis off
RR=15*pr;
[n1 m1]=size(keypoint1);
[idx,dist]=rangesearch(PC1,keypoint1,RR);
MV1=[];
MDV1=[];
for i=1:n1
    KNN=PC1(idx{i},:);
    d=dist{i};
    [V] = LRF_TriLCI(KNN,RR,d,keypoint1(i,:));
    MV1=[MV1;V];
    [DV] = LOVC(KNN,keypoint1(i,:),V,RR);
    MDV1=[MDV1;DV];
end
[n2 m2]=size(keypoint2);
[idx,dist]=rangesearch(PC2,keypoint2,RR);
MV2=[];
MDV2=[];
for i=1:n2
    KNN=PC2(idx{i},:);
    d=dist{i};
    [V] = LRF_TriLCI(KNN,RR,d,keypoint2(i,:));
    MV2=[MV2;V];
    [DV] = LOVC(KNN,keypoint2(i,:),V,RR);
    MDV2=[MDV2;DV];
end
[idxx distt]=knnsearch(MDV1,MDV2,'Dist','hamming','k',2);
Mmatch=[];
for i=1:n2
    if distt(i,1)/distt(i,2)<=0.9
        match=[i idxx(i,1)];
        Mmatch=[Mmatch;match];
    end
end
[n3 m3]=size(Mmatch);
%2-point RANSACWC
tic
ninlier0=0;
for i=1:2000
    Sidx=randperm(n3,2);
    p1=keypoint2(Mmatch(Sidx(1),1),:);
    p2=keypoint2(Mmatch(Sidx(2),1),:);
    q1=keypoint1(Mmatch(Sidx(1),2),:);
    q2=keypoint1(Mmatch(Sidx(2),2),:);
    Vp1=MV2(3*Mmatch(Sidx(1),1)-2:3*Mmatch(Sidx(1),1),:);
    Vp2=MV2(3*Mmatch(Sidx(2),1)-2:3*Mmatch(Sidx(2),1),:);
    Vq1=MV1(3*Mmatch(Sidx(1),2)-2:3*Mmatch(Sidx(1),2),:);
    Vq2=MV1(3*Mmatch(Sidx(2),2)-2:3*Mmatch(Sidx(2),2),:);
    d1=norm(p1-p2);
    d2=norm(q1-q2);
    d12=abs(Vp1(:,3)'*(p2-p1)');
    d13=abs(Vp2(:,3)'*(p1-p2)');
    d22=abs(Vq1(:,3)'*(q2-q1)');
    d23=abs(Vq2(:,3)'*(q1-q2)');
    datatheta1=real(acos((trace(Vp1*inv(Vp2))-1)/2)*(180/pi));
    datatheta2=real(acos((trace(Vq1*inv(Vq2))-1)/2)*(180/pi));
    datatheta11=real(acos(Vp1(:,3)'*Vp2(:,3))*(180/pi));
    datatheta22=real(acos(Vq1(:,3)'*Vq2(:,3))*(180/pi));
    if abs(d1-d2)<2*pr && abs(datatheta1-datatheta2)<10 && abs(datatheta11-datatheta22)<10 && abs(d12-d22)<2*pr && abs(d13-d23)<2*pr
        p11=p1+Vp1(:,3)';
        p22=p2+Vp2(:,3)';
        q11=q1+Vq1(:,3)';
        q22=q2+Vq2(:,3)';
        A=[p1;p2;p11;p22];
        Y=[q1;q2;q11;q22];
        [h,l]=size(A);
        uA=[mean(A(:,1)) mean(A(:,2)) mean(A(:,3))];
        uY=[mean(Y(:,1)) mean(Y(:,2)) mean(Y(:,3))];
        H=zeros(3);
        for j=1:h
            H=H+(A(j,:)-uA)'*(Y(j,:)-uY);
        end
        [U S V]=svd(H);
        D=diag([1 1 det(U*V')]);
        R1=V*D*U';
        t1=mean([q1;q2])-mean([p1;p2])*R1';
        PC2t=PC2*R1'+ones(n,1)*t1;
        [idx1 dist1]=knnsearch(PC1,PC2t,'k',1);
        iddist1=find(dist1<3*pr);
        ninlier=length(iddist1);
        if ninlier>ninlier0
            R=R1;
            t=t1;
            ninlier0=ninlier;
        end
    end
end
PC2t=PC2*R'+ones(n,1)*t;
T=[R t';zeros(1,3) 1];
time=toc;
Rtrue=Ttrue(1:3,1:3);
ttrue=Ttrue(1:3,4);
Rerror=real(acos((trace(Rtrue*inv(R))-1)/2)*(180/pi));
terror=norm(ttrue-t');
plot3(PC1(:,1),PC1(:,2),PC1(:,3),'.b','MarkerSize',1);
hold on;
plot3(PC2t(:,1),PC2t(:,2),PC2t(:,3),'.r','MarkerSize',1);
set(gca,'DataAspectRatio',[1 1 1]);
axis off

v = VideoWriter('RC.mp4');
v.FrameRate = 30;
open(v)
A = VideoReader('Movie S7.mp4');
B = VideoReader('Movie S8v2.avi');
C = VideoReader('Movie S9.mp4');
D = [];
fc=0;
while hasFrame(B)
    fc=fc+1;
    A1 = readFrame(A);
    B1 = readFrame(B);
    C1 = readFrame(C);
    D(:,:,1) = [A1(:,:,1);B1(:,:,1);C1(:,:,1)];
    D(:,:,2) = [A1(:,:,2);B1(:,:,2);C1(:,:,2)];
    D(:,:,3) = [A1(:,:,3);B1(:,:,3);C1(:,:,3)];
    sc = 1.55;
    cn = 56/34;
    cn = cn*sc/1.4;
    fs = 20;
    DD = imresize(uint8(D),[568,640]*sc);
    imshow(DD)
    text(10,(10+170)*sc,'Active','fontsize',fs,'fontname','Arial')
    text(10,(10+162+192)*sc,'Reverse','fontsize',fs,'fontname','Arial')
    text(10,(10+170+192+184)*sc,'Lorenz','fontsize',fs,'fontname','Arial')
    hold on
    plot([170,269]*sc/1.4,[788,788]*sc/1.4,'k','linewidth',3)
    text(180*sc/1.4,(10+170+192+184)*sc-10*sc/1.4,'60 cm','fontsize',fs,'fontname','Arial')
    DDD = getframe;
    DDD = DDD.cdata;
    DDD = DDD(1:end-1,:,:);
    p1 = (1080-size(DDD,1))/2;
    p2 = (1920-size(DDD,2))/2;
    DDD = padarray(DDD,[p1,p2]);
    writeVideo(v,DDD)
    clf;
%     pause(.01)
end
close(v)
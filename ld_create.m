clear;clc
close all;
precision=0.5;%激光雷达精度 
T=1;%扫描间隔
maxTime=100;%总时间
s=3;%每次逆时针旋转随机0~s度
clutterNum=50 %杂波最大值
Scenario.Z{1}=cell(1,maxTime)
F = [1, T, 0, 0
    0, 1, 0, 0
    0, 0, 1, T
    0, 0, 0, 1]; %目标转移矩阵
W=8;
H=5;
V=8;%初始x方向速度
a=[-50,V,-125+200,0]';%分别为[x轴起始坐标，x方向速度，y起始坐标，y方向速度] 矩形左下角坐标及参数
b=[a(1)+W,a(2),a(3),a(4)]';
c=[a(1)+W,a(2),a(3)+H,a(4)]';
d=[a(1),a(2),a(3)+H,a(4)]';
for i=1:T:maxTime
    %绘图
    figure(1);
    line([a(1) b(1)],[a(3)-200 b(3)-200]);
    line([b(1) c(1)],[b(3)-200 c(3)-200]);
    line([c(1) d(1)],[c(3)-200 d(3)-200]);
    line([d(1) a(1)],[d(3)-200 a(3)-200]);
    hold on;
    for j=1:precision:179 %激光雷达扫描扇形的范围
        k1=tand(j);
        xs=[a(1),b(1),c(1),d(1),a(1)];
        ys=[a(3),b(3),c(3),d(3),a(3)];
        temporary=[];
        %y1=linspace(0,10);
        %x1=y1/k1;
        %plot(x1,y1);
        hold on;
        for k=2:5
            if((xs(k)-xs(k-1))~=0)
                k2=(ys(k)-ys(k-1))/(xs(k)-xs(k-1));
                u=[(-xs(k)*k2+ys(k))/(k1-k2),k1*(((-xs(k)*k2+ys(k))/(k1-k2)))]';
                u(3)=sqrt(u(1)^2+u(2)^2);
                if (u(1)<=max(xs(k),xs(k-1))+0.01 && u(1)>=min(xs(k),xs(k-1))-0.01 && u(2)<=max(ys(k),ys(k-1))+0.01 && u(2)>=min(ys(k),ys(k-1))-0.01)
                    temporary=[temporary,u];
                end
            end
            if((xs(k)-xs(k-1))==0)
                u=[xs(k),k1*xs(k)]';
                u(3)=sqrt(u(1)^2+u(2)^2);
                if (u(1)<=max(xs(k),xs(k-1))+0.01 && u(1)>=min(xs(k),xs(k-1))-0.01 && u(2)<=max(ys(k),ys(k-1))+0.01 && u(2)>=min(ys(k),ys(k-1))-0.01)
                    temporary=[temporary,u];
                end
            end
        end
        if ~isempty(temporary)
            temporarys=temporary';
            temporarys=sortrows(temporarys,3);
            hold on;
            temporarys(1,2)=temporarys(1,2)-200;
            plot(temporarys(1,1),temporarys(1,2),'*');
            Scenario.Z{1}{i}=[Scenario.Z{1}{i},temporarys(1,1:2)'];
        end
    end
    %状态改变（主运动）
    a=F*a;
    b=F*b;
    c=F*c;
    d=F*d;
    %转动
    rx0=(a(1)+c(1))/2;
    ry0=(a(3)+c(3))/2;
%    ra=rand();
%     if ra>=0.66
%         rs=5;
%     elseif ra>=0.33 && ra<0.66
%         rs=-5;
%     else
%         rs=0;
%    end
    rs=s;%*rand();
    M=[cosd(rs),-sind(rs)
            sind(rs),cosd(rs)];%转动矩阵
    va=[a(2),a(4)]';
    va=M*va;
    a(2)=va(1);
    a(4)=va(2);
    b(2)=va(1);
    b(4)=va(2);
    c(2)=va(1);
    c(4)=va(2);
    d(2)=va(1);
    d(4)=va(2);
    ax=(a(1)-rx0)*cosd(rs)-(a(3)-ry0)*sind(rs)+rx0;
    a(3)=(a(1)-rx0)*sind(rs)+(a(3)-ry0)*cosd(rs)+ry0;
    bx=(b(1)-rx0)*cosd(rs)-(b(3)-ry0)*sind(rs)+rx0;
    b(3)=(b(1)-rx0)*sind(rs)+(b(3)-ry0)*cosd(rs)+ry0;
    cx=(c(1)-rx0)*cosd(rs)-(c(3)-ry0)*sind(rs)+rx0;
    c(3)=(c(1)-rx0)*sind(rs)+(c(3)-ry0)*cosd(rs)+ry0;
    dx=(d(1)-rx0)*cosd(rs)-(d(3)-ry0)*sind(rs)+rx0;
    d(3)=(d(1)-rx0)*sind(rs)+(d(3)-ry0)*cosd(rs)+ry0;
    a(1)=ax;
    b(1)=bx;
    c(1)=cx;
    d(1)=dx;
    ClutterNum=rand()*clutterNum;
%     for j=1:ClutterNum %杂波生成
%         cx=-200+400*rand();
%         cy=-200+400*rand();
%         Scenario.Z{1}{i}=[Scenario.Z{1}{i},[cx,cy]'];
%     end
end
% figure(2)
% for i = 1:size(Scenario.Z{1},2)
%     pause(0.1)
%     if ~isempty(Scenario.Z{1}{i})
%         plot(Scenario.Z{1}{i}(1,:),Scenario.Z{1}{i}(2,:),'b.');
%         hold on;
%     end
% end
Scenario.false_alarm_rate=60;
Scenario.detection_prob=0.9;
save ('64673','Scenario') ;
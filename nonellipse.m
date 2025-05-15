
function nonellipse(numZ,detection_prob,false_alarm_rate,plotline)
isellipse=1
% numZ=1
%椭圆生成
% detection_prob=0.99;%目标漏检概率
precision=0.08;%激光雷达精度
Scenario.false_alarm_rate=false_alarm_rate; %泊松杂波强度
detection_prob_ld=0.99; %激光雷达每次漏检概率
clutter=1;%杂波
u_gauss=0.05; %高斯噪声强度

if isellipse
    for z=1:numZ
        
        T=1;%扫描间隔
        maxTime=100;%总时间
        a=6;
        b=4;%椭圆的长短半轴
        rs=3;%每次逆时针旋转随机0~s度
        v=[8 0];%初始x方向速度
        theta=0;%椭圆一初始角度
        theta2=90;%椭圆二初始角度
        x0=300;
        y0=250;
        x02=300;
        y02=250;
        syms x k;
        Scenario.Z{z}=cell(1,maxTime);
        for i=1:T:maxTime
            [3 z i]
            
            %%%%%%%%%%%%% 绘制椭圆
            theta_vals = linspace(0, 2*pi, 100);
            ellipse_x = a*cos(theta_vals)*cosd(theta) - b*sin(theta_vals)*sind(theta)+x0;
            ellipse_y = a*cos(theta_vals)*sind(theta) + b*sin(theta_vals)*cosd(theta)+y0;
            
            ellipse_x2 = a*cos(theta_vals)*cosd(theta2) - b*sin(theta_vals)*sind(theta2)+x02;
            ellipse_y2 = a*cos(theta_vals)*sind(theta2) + b*sin(theta_vals)*cosd(theta2)+y02;
            
            plot(ellipse_x, ellipse_y, 'b-', 'LineWidth', 2); 
            hold on
            plot(ellipse_x2, ellipse_y2, 'b-', 'LineWidth', 2); 
            hold on
            axis  equal;
            %%%%%%%%%%%%%%%%%%%%%%存储真实数据
            if z==1
                targetTracks(1).x(1:4,i)=[(x0+x02)/2 (y0+y02)/2 v]';
                theta_ = linspace(0, 2*pi, 10000);
                ellipse_x = a*cos(theta_)*cosd(theta) - b*sin(theta_)*sind(theta)+x0;
                ellipse_y = a*cos(theta_)*sind(theta) + b*sin(theta_)*cosd(theta)+y0;
                
                ellipse_x2 = a*cos(theta_)*cosd(theta2) - b*sin(theta_)*sind(theta2)+x0;
                ellipse_y2 = a*cos(theta_)*sind(theta2) + b*sin(theta_)*cosd(theta2)+y02;
                
                X=cov([ellipse_x ellipse_x2],[ellipse_y ellipse_y2]);
                targetTracks(1).X(:,:,i)=X;
            end
            
            ellipse_eq = ((x-x0)*cosd(theta) + (k*x-y0)*sind(theta))^2/a^2 + ((x-x0)*sind(theta) - (k*x-y0)*cosd(theta))^2/b^2 - 1;
            x_sol = solve(ellipse_eq,x);
            
            ellipse_eq2 = ((x-x02)*cosd(theta2) + (k*x-y02)*sind(theta2))^2/a^2 + ((x-x02)*sind(theta2) - (k*x-y02)*cosd(theta2))^2/b^2 - 1;
            x_sol2 = solve(ellipse_eq2,x);
            
            k_=atand((y0+y02)/(x0+x02));
            if rand()<detection_prob
                for j=k_-5:precision:k_+5 %激光雷达扫描扇形的范围
                    %[3 z i j]
                    tk=tand(j);
%                     y1=linspace(200,600);
%                     x1=y1/tk;
%                     plot(x1,y1);
                    hold on;

                    
                    x_solv=subs(x_sol,k,tk);
                    y_solv=tk*x_solv;
                    
                    x_solv2=subs(x_sol2,k,tk);
                    y_solv2=tk*x_solv2;
                    
                    points1 = [double(x_solv), double(y_solv)];
                    points2 = [double(x_solv2), double(y_solv2)];
                    pointsu = [points1' points2'];
                    points=[];
                    hpoints=[];
                    for l=1:length(pointsu(1,:))
                        if isreal(pointsu(1,l))
                            points(1:2,end+1)=pointsu(1:2,l);
                            hpoints(end+1)= hypot(points(1,end),points(1,end));
                        end
                    end

                    if ~isempty(points)
                        [~,index]=sort(hpoints);
                        in=index(1);
%                         if hypot(points(1,1),points(1,2))<hypot(points(2,1),points(2,2))
%                             temporary=points(1,:)';
%                         else
%                             temporary=points(2,:)';
%                         end
                        temporary=points(:,in);
                        if rand()<detection_prob_ld
%                             figure(1);
%                             plot(temporary(1,1), temporary(2,1),'or'); % 绘制交点
%                             hold on
                            temporary(1)=temporary(1)+randn()*u_gauss;
                            temporary(2)=temporary(2)+randn()*u_gauss;
                            plot(temporary(1),temporary(2),'or')
                            Scenario.Z{z}{i}=[Scenario.Z{z}{i} temporary];
                        end
                     end
                end
            end
            %转态改变
            x0=x0+v(1);
            y0=y0+v(2);
            x02=x02+v(1);
            y02=y02+v(2);
            theta=theta+rs;
            theta2=theta2+rs;
            %转动

            M=[cosd(rs),-sind(rs)
                sind(rs),cosd(rs)];%转动矩阵
            v=(M*v')';%速度转动
            
            if clutter
                clutterNum=poissrnd(Scenario.false_alarm_rate);%杂波生成
                for j=1:clutterNum
                    cx=100+600*rand();
                    cy=100+600*rand();
                    Scenario.Z{z}{i}=[Scenario.Z{z}{i},[cx,cy]'];
                end
            end
        end
    targetTracks(1).birthTime=1;
    targetTracks(1).deathTime=100;
    end
end
if plotline
    figure(2)
    for i = 1:size(Scenario.Z{1},2)
        %pause(0.1)
        if ~isempty(Scenario.Z{1}{i})
            plot(Scenario.Z{1}{i}(1,:),Scenario.Z{1}{i}(2,:),'b.');
            hold on;
        end
    end
end

Scenario.detection_prob=detection_prob;
save ('64673e','Scenario') ;
save ('64673targetTrackse','targetTracks')
end
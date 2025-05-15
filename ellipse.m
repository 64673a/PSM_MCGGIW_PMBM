
function ellipse(numZ,detection_prob,false_alarm_rate,plotline)
isellipse=1;

precision=0.04;%Lidar accuracy
Scenario.false_alarm_rate=false_alarm_rate; %Clutter rate 
detection_prob_ld=0.99; %Single measurement detection probability
clutter=1;%杂波
u_gauss=0.02; %Gaussian noise intensity

if isellipse
    for z=1:numZ
        
        T=1;%time step
        maxTime=100;%总时间
        a=4;
        b=2.5;%The major and minor axes of an ellipse
        rs=4;
        v=[8 0];%Initial velocity
        theta=0;%Initial ellipse angle
        x0=300;
        y0=250;
        syms x k;
        Scenario.Z{z}=cell(1,maxTime);
        for i=1:T:maxTime
            [3 z i]
            if plotline
                figure(11)
                axis([100,600,100,600]);
                theta_vals = linspace(0, 2*pi, 100);
                ellipse_x = a*cos(theta_vals)*cosd(theta) - b*sin(theta_vals)*sind(theta)+x0;
                ellipse_y = a*cos(theta_vals)*sind(theta) + b*sin(theta_vals)*cosd(theta)+y0;

                plot(ellipse_x, ellipse_y, 'c-', 'LineWidth', 1.5); 
                hold on

                axis  equal;
            end
            if z==1
                targetTracks(1).x(1:4,i)=[x0 y0 v]';
                theta_ = linspace(0, 2*pi, 10000);
                ellipse_x = a*cos(theta_)*cosd(theta) - b*sin(theta_)*sind(theta)+x0;
                ellipse_y = a*cos(theta_)*sind(theta) + b*sin(theta_)*cosd(theta)+y0;

                
                X=cov(ellipse_x,ellipse_y);
                targetTracks(1).X(:,:,i)=X;
            end
            
            ellipse_eq = ((x-x0)*cosd(theta) + (k*x-y0)*sind(theta))^2/a^2 + ((x-x0)*sind(theta) - (k*x-y0)*cosd(theta))^2/b^2 - 1;
            x_sol = solve(ellipse_eq,x);

            k_=atand(y0/x0);
            if rand()<detection_prob || i<10
                for j=k_-1:precision:k_+1
                    tk=tand(j);
                    
                    x_solv=subs(x_sol,k,tk);
                    y_solv=tk*x_solv;
                    
                    points = [double(x_solv), double(y_solv)];
                    if isreal(points)
                        if hypot(points(1,1),points(1,2))<hypot(points(2,1),points(2,2))
                            temporary=points(1,:)';
                        else
                            temporary=points(2,:)';
                        end
                        if rand()<detection_prob_ld
                            temporary(1)=temporary(1)+randn()*u_gauss;
                            temporary(2)=temporary(2)+randn()*u_gauss;
                            Scenario.Z{z}{i}=[Scenario.Z{z}{i} temporary];
                        end
                    end
                end
            end
            %转态改变
            x0=x0+v(1);
            y0=y0+v(2);
            if i==50
                rs=rs+1;
            end
            theta=theta+rs;
            %转动
    
            M=[cosd(rs),-sind(rs)
                sind(rs),cosd(rs)];%转动矩阵
            v=(M*v')';%速度转动
            
            if clutter
                clutterNum=poissrnd(Scenario.false_alarm_rate);%杂波生成
                for j=1:clutterNum
                    cx=100+500*rand();
                    cy=100+500*rand();
                    Scenario.Z{z}{i}=[Scenario.Z{z}{i},[cx,cy]'];
                end
            end
        end
    targetTracks(1).birthTime=1;
    targetTracks(1).deathTime=100;
    end
end

Scenario.detection_prob=detection_prob;
save ('64673e','Scenario') ;
save ('64673targetTrackse','targetTracks')
end
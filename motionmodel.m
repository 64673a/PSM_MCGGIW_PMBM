classdef motionmodel
    %MOTIONMODEL is a class containing different motion models
    
    methods (Static)
        function obj = cvmodel(T,sigma)
            %CVMODEL creates a 2D nearly constant velocity model 创建二维近似等速模型
            %INPUT:     T: sampling time --- scalar 采样间隔
            %           sigma: standard deviation of motion noise --- scalar 运动噪声标准差
            %OUTPUT:    obj.d: object state dimension --- scalar 状态维度 
            %           obj.F: function handle return a motion transition matrix --- 2 x 2 matrix 运动过渡矩阵
            %           obj.Q: motion noise covariance --- 4 x 4 matrix 运动噪声协方差
            %           obj.f: function handle return state prediction --- 4 x 1 vector 状态预测
            % NOTE: the motion model assumes that the state vector x consist of the following states: 运动模型假设状态向量x由以下状态组成:
            %           px          X-position x坐标
            %           py          Y-position y坐标
            %           vx          X-velocity x方向v
            %           vy          Y-velocity y方向v
            obj.d = 4;
            obj.F{1} = @(x) [
                1 0 T 0;
                0 1 0 T;
                0 0 1 0;
                0 0 0 1];
            obj.PPPF = @(x) [
                1 0 T 0;
                0 1 0 T;
                0 0 1 0;
                0 0 0 1]
            obj.Q = sigma^2*[
                T^4/4   0       T^3/2   0;
                0       T^4/4   0       T^3/2;
                T^3/2   0       T^2     0;
                0       T^3/2   0       T^2;
                ];
            obj.f = @(x) obj.F(x)*x;
        end
        
        function obj = immmodel(T,sigma)
            %CVMODEL creates a 2D nearly constant velocity model 创建二维近似等速模型
            %INPUT:     T: sampling time --- scalar 采样间隔
            %           sigma: standard deviation of motion noise --- scalar 运动噪声标准差
            %OUTPUT:    obj.d: object state dimension --- scalar 状态维度 
            %           obj.F: function handle return a motion transition matrix --- 2 x 2 matrix 运动过渡矩阵
            %           obj.Q: motion noise covariance --- 4 x 4 matrix 运动噪声协方差
            %           obj.f: function handle return state prediction --- 4 x 1 vector 状态预测
            % NOTE: the motion model assumes that the state vector x consist of the following states: 运动模型假设状态向量x由以下状态组成:
            %           px          X-position x坐标
            %           py          Y-position y坐标
            %           vx          X-velocity x方向v
            %           vy          Y-velocity y方向v
            w=6*pi/180;
            obj.w=w;
            obj.d = 4;

            obj.F{1}= @(x) [
                                  1 0 T 0;
                                  0 1 0 T;
                                  0 0 1 0;
                                  0 0 0 1];
            obj.F{2} = @(x) [
                                  1 0 sin(w*T)/w (cos(w*T)-1)/w;
                                  0 1 (1-cos(w*T))/w sin(w*T)/w;
                                  0 0 cos(w*T) -sin(w*T);
                                  0 0 sin(w*T) cos(w*T)];
            obj.F{3}= @(x) [
                                  1 0 sin(-w*T)/-w (cos(-w*T)-1)/-w;
                                  0 1 (1-cos(-w*T))/-w sin(-w*T)/-w;
                                  0 0 cos(-w*T) -sin(-w*T);
                                  0 0 sin(-w*T) cos(-w*T)];
            obj.Q = sigma^2*[
                T^4/4   0       T^3/2   0;
                0       T^4/4   0       T^3/2;
                T^3/2   0       T^2     0;
                0       T^3/2   0       T^2;
                ];
           obj.PPPF = @(x) [
                                  1 0 T 0;
                                  0 1 0 T;
                                  0 0 1 0;
                                  0 0 0 1];
        end
        
        function obj = ctmodel(T,sigmaV,sigmaOmega)
            %CTMODEL creates a 2D coordinate turn model with nearly constant polar velocity and turn rate
            %INPUT:     T: sampling time --- scalar
            %           sigmaV: standard deviation of motion noise added to polar velocity --- scalar 极速运动噪声的标准偏差
            %           sigmaOmega: standard deviation of motion noise added to turn rate --- scalar 添加到转动率的运动噪声的标准偏差
            %OUTPUT:    obj.d: object state dimension --- scalar
            %           obj.F: function handle return a motion Jacobian matrix --- 5 x 5 matrix
            %           obj.f: function handle return state prediction --- 5 x 1 vector
            %           obj.Q: motion noise covariance --- 5 x 5 matrix
            % NOTE: the motion model assumes that the state vector x consist of the following states:
            %           px          X-position
            %           py          Y-position
            %           v           velocity
            %           phi         heading
            %           omega       turn-rate
            obj.d = 5;
            obj.f = @(x) x + [
                T*x(3)*cos(x(4));
                T*x(3)*sin(x(4));
                0;
                T*x(5);
                0];
            obj.F = @(x) [
                1 0 T*cos(x(4)) -T*x(3)*sin(x(4)) 0;
                0 1 T*sin(x(4)) T*x(3)*cos(x(4))  0;
                0 0 1           0                 0;
                0 0 0           1                 T;
                0 0 0           0                 1
                ];
            G = [zeros(2,2); 1 0; 0 0; 0 1];
            obj.Q = G*diag([sigmaV^2 sigmaOmega^2])*G';
        end

    end
end
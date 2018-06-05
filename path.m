function [x_vector,y_vector,theta_vector] = path()
    global t0;
    global tf;
    
    % multiple assignment
    % dx0 is velocity in x direction at t0
    [t0,tc,tf] = deal(0,2.5,5);
    [x0,y0,dx0,dy0] = deal(0,0,0,0);
    [xc,yc,dxc,dyc] = deal(5,5,1,sqrt(3));
    [xf,yf,dxf,dyf] = deal(10,10,2,2);
    
    % Bernstein coefficients and derivatives 
    % B0(u(tc)) -> B(0,mu(tc))
    range = 0:5;
    Bt0 = arrayfun(@B,range,repmat(mu(t0),1,6));
    Btc = arrayfun(@B,range,repmat(mu(tc),1,6));
    Btf = arrayfun(@B,range,repmat(mu(tf),1,6));
    dBt0 = arrayfun(@dB,range,repmat(mu(t0),1,6));
    dBtc = arrayfun(@dB,range,repmat(mu(tc),1,6));
    dBtf = arrayfun(@dB,range,repmat(mu(tf),1,6));
    
    % Coefficient(Wx) to find out x(t) -> Ax = Bx*Wx 
    Wx = zeros(1,6);
    [Wx(1),Wx(6)] = deal(x0,xf);
    Bx = zeros(4,4);
    Bx(1,1:4) = Btc(2:5);
    Bx(2,1:4) = dBt0(2:5);
    Bx(3,1:4) = dBtf(2:5);
    Bx(4,1:4) = dBtc(2:5);
    
    Ax = zeros(4,1);
    Ax(1) = xc - Wx(1)*Btc(1) - Wx(6)*Btc(6);
    Ax(2) = dx0 - Wx(1)*dBt0(1) - Wx(6)*dBt0(6);
    Ax(3) = dxf - Wx(1)*dBtf(1) - Wx(6)*dBtf(6);
    Ax(4) = dxc - Wx(1)*dBtc(1) - Wx(6)*dBtc(6);
    
    Wx(2:5) = pinv(Bx)*Ax;
    disp(Wx);
    
    % Coefficient(Wk) to find out y(t) -> Ak = Bk*Wk
    Wk = zeros(1,6);
    [Wk(1),Wk(6)] = deal(0,dyf/dxf);
    
    % Bk
    Bk = zeros(4,4);
    Bk(1,1:4) = dBt0(2:5);
    Bk(2,1:4) = dBtf(2:5);
    uptc = mu(tc);
    Bk(3,1:4) = [F(uptc,1,Wx),F(uptc,2,Wx),F(uptc,3,Wx),F(uptc,4,Wx)];
    uptf = mu(tf);
    Bk(4,1:4) = [F(uptf,1,Wx),F(uptf,2,Wx),F(uptf,3,Wx),F(uptf,4,Wx)];
    
    % Ak
    Ak = zeros(4,1);
    Ak(1) = sec(0)^2 - Wk(1)*dBt0(1) - Wk(6)*dBt0(6);
    Ak(2) = sec(atan(dyf/dxf))^2 - Wk(1)*dBtf(1) - Wk(6)*dBtf(6);
    Ak(3) = yc - Wk(1)*F(uptc,0,Wx) - Wk(6)*F(uptc,5,Wx);
    Ak(4) = yf - Wk(1)*F(uptf,0,Wx) - Wk(6)*F(uptf,5,Wx);
    
    Wk(2:5) = pinv(Bk)*Ak;
    disp(Wk);
    
%     Final Wx for original case = [0 0 10 1 8 10];
%     Final Wk for original case = [0 1 1.2788 4.3196 -1 1];
    
    % x(t),y(t),theta(t)
    syms t;
    time = 0:0.1:5;
    
    syms x(t);
    sum = sym('0');
    for i = 0:5
        sum = sum + Wx(i+1) * B(i,mu(t)); 
    end
    x(t) = sum;
    x_vector = double(x(time));
    
    syms y(t);
    sum = sym('0');
    for i = 0:5
        sum = sum + Wk(i+1) * F(mu(t),i,Wx);
    end
    y(t) = y0 + sum;
    y_vector = double(y(time));
    
    syms theta(t);
    sum = sym('0');
    for i = 0:5
       sum = sum + Wk(i+1)*B(i,mu(t)); 
    end
    theta(t) = atan(sum);
    theta_vector = rad2deg(double(theta(time)));
    
    % data saving
    save('data','x_vector','y_vector','theta_vector');
end

function m = mu(t)
    global t0;
    global tf;
    m = (t - t0)/(tf - t0); 
end

function Bi = B(i,m)
    % wrt variable mu(t)
    Bi = nchoosek(5,i) * ((1-m).^(5-i)) * (m.^i);
end

function dBi = dB(i,mu)
    global t0;
    global tf;
    syms f(m)
    f(m) = ((1-m).^(5-i)) * (m.^i);
    % mu in terms function of t
    df = diff(f,m) / (tf - t0);
    dBi = double(nchoosek(5,i) * df(mu));
end

function dBi = dB_sym(i,mu)
    % symbolic dB for integration
    global t0;
    global tf;
    syms f(m)
    f(m) = ((1-m).^(5-i)) * (m.^i);
    % mu in terms function of t
    df = diff(f,m) / (tf - t0);
    dBi = (nchoosek(5,i) * df(mu));
end

function integration = F(Upper,j,Wx)
    % definite integral: t0 -> Upper = mu(t)
    global t0;
    global tf;
    syms m;
    sum = sym('0');
    for i = 0:5
       sum = sum + (Wx(i+1) * dB_sym(i,m));
    end
    f = (tf - t0) * sum * B(j,m);
    integration = int(f,mu(t0),Upper);
end



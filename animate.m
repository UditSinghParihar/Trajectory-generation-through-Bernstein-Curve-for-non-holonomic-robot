load('data','x_vector','y_vector','theta_vector');
[x,y,theta] = deal(x_vector,y_vector,theta_vector);

len = 1;
for i = 1:1:51
    X = [x(i) - len*cos(deg2rad(theta(i))), x(i) + len*cos(deg2rad(theta(i)))];
    Y = [y(i) - len*sin(deg2rad(theta(i))), y(i) + len*sin(deg2rad(theta(i)))];
    plot(x,y,X,Y,'--');
    set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 3);
    axis([-1 11 -1 11]);
    drawnow;
    pause(0.1);
end

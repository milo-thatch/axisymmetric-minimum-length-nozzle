function results=wall_minus(theta0,point0,point1,point3)
    global th
    %theta0 is a scalar
    %point0 is an array
    %point1 and point3 are arrays
    m13=(point3(4)-point1(4))/(point3(3)-point1(3));
    A=@(INC) [1,-m13;1,-tan(0.5*(theta0+INC))]; %2x2 matrix
    B=@(INC) [point3(4)-m13*point3(3);point0(2)-tan(0.5*(theta0+INC))*point0(1)];
    
    tol=1e3;
    theta_new=theta0;
    point_new=[point0(1);point0(2)];
    while(tol>th)
        old_point=point_new;
        old_theta=theta_new;
        point_new=A(theta_new)\B(theta_new); %array
        theta_new=point3(1)+sqrt(((point_new(1)-point1(3))^2+(point_new(2)-point1(4))^2)/...
            ((point3(3)-point1(3))^2+(point3(4)-point1(4))^2))*(point1(1)-point3(1));
        tol=max(abs(point_new(1)-old_point(1)),abs(point_new(2)-old_point(2)));
        tol=max(tol,abs(theta_new-old_theta));
    end
    results=[theta_new,point_new'];
end
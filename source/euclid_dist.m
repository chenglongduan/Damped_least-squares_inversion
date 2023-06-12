function dist = euclid_dist( A, B, mode )
%Calculate the Euclid distance from point A to point B (two-dimension)
%Inputs:   A=[x1 y1] or [x1 x2] or [x1 x2...xn]; 
%          B=[x2 y2] or [y1 y2] or [y1 y2...yn].
%          mode = 1 : Coordinate mode,i.e.,A=[x1 y1]; B=[x2 y2].
%               = 2 : Component mode,i.e.,A=[x1 x2]; B=[y1 y2].
%               = 3 : Component vector,i.e.,A=[x1 x2...xn];B=[y1 y2...yn].
%Outputs:  Distance between A and B.

%-----By Chenglong Duan,Nanjing University,2015.-----

if nargin<3
    error('Lack inputs!');
end

if mode==1
    dist=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
end

if mode==2
    dist=sqrt((A(1)-A(2))^2+(B(1)-B(2))^2);
end

if mode==3
    dist=zeros(1,length(A)-1);
    for i=1:length(A)-1
       dist(i)=sqrt((A(i)-A(i+1))^2+(B(i)-B(i+1))^2);
    end
end

end


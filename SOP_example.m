%Maintaining symmetry, I want to try a bigger intercept for the gamma
%function PLUS probabilistic one-sided breaking (same gamma as lobby4).

%gamma: 1.25 + cn^.2
%cn = .00166
% tn = (8*(1 + .00156^.1)-5)/(68-8*(1 + .00156^.1)) = 0.12592

tn = (8*(1.25 + .00166^.2)-5)/(68-8*(1.25 + .00166^.2));

a=.2;
tns=tn;
ta = [0:.001:.1259];
tad = size(ta,2);
c=[0:.00005:.005];
tac = size(c,2);
[ta,c]=meshgrid(ta,c);

lhs = (8*(1.25 + (c).^a)-5).*(tn-ta) + (4*(1.25 + (c).^a)-34).*(tn^2-ta.^2)-3.*(tns-ta)+9.*(tns^2-ta.^2);

%add e to gamma
e1=.25; e2=.25;
up = (8*(1.25 + (c).^a +e1)-5).*(tn-ta) + (4*(1.25 + (c).^a +e1)-34).*(tn^2-ta.^2)-3.*(tns-ta)+9.*(tns^2-ta.^2);
low = (8*(1.25 + (c).^a -e2)-5).*(tn-ta) + (4*(1.25 + (c).^a -e2)-34).*(tn^2-ta.^2)-3.*(tns-ta)+9.*(tns^2-ta.^2);
spread = up - low;
b = up./spread;
for j=1:tad
    if lhs(1,j)< 0
        b(1,j) = 0;
        else
        b(1,j) = 1;
     end; 
     for i=2:101
        if up(i,j) < 0
        b(i,j) = 0;
        elseif low(i,j) > 0
        b(i,j) = 1;
      end
    end
end

l = b;
l2 = 4/49*((1+tn).^2);
l3 = l.*l2;
l4 = 4/49*((1+ta).^2);
l5 = l3 + (1-l).*l4 - c;
[M,I] = max(l5,[],1);
opt_C = (I -1)/20000;
%plot(TA,opt_C);


TA = [0:.001:.125];

for j=1:tad
  B(j) = b(I(j),j);
end

E=B;
E2 = 21/49*(2-2*tn^2);
E3 = E*E2;
E4 = 21/49.*(2-2*TA.^2);
E5 = (1-E).*E4;
E6 = E3 + (1-E).*E4 ;
%plot(TA,E6);

[M2,I2] = max(E6); %I2 gives E's choice of tau^a (63 ==> .062)
ta_star = (I2 - 1)/1000
c_star = I(I2) %c_star indicates what contributions the lobby chooses
e_welfare = E6(I2)
c_star_ = (c_star -1)/20000
b_star = b(c_star,I2) %break probability
gamma = 1.25 + (c_star_).^a

%e =[.00:.01:.25];
%t =[.113 .111 .109 .107 .105 .104 .104 .104 .105 .106 .107 .108 .109 .111 .096 .096 .091 .09 .086 .087 .083 .081 .08 .08 .078 .078];
%plot(e,t);
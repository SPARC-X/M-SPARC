function Ylm = sphericalHarmonics(X,Y,Z,l,m,kind)
% @brief    sphericalHarmonics(X,Y,Z,l,m,kind) calculates the sperical harmonics
%           Ylm (real or complex) at the given positions (X,Y,Z)
%
% @param X      The x coordinates of the positions 
% @param Y      The y coordinates of the positions 
% @param Z      The z coordinates of the positions 
% @param l      The azimuthal quantum number
% @param m      The magnetic quantum number
% @param kind   Optional arguments: 'real'    - Compute real spherical harmonics
%                                   'complex' - Compute complex spherical harmonics
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%============================================================================================

szX = size(X); szY = size(Y); szZ = size(Z);
if ~isequal(szX,szY,szZ)
	error('sphericalHarmonics(X,Y,Z,l,m,kind): dimensions of <X> <Y> <Z> must match');
end

isInt = l - round(l);
if (isInt ~= 0)||(l < 0)
	error('sphericalHarmonics(X,Y,Z,l,m,kind): <l> must be a nonnegative integer');
end

isInt = m - round(m);
if (isInt ~= 0)||(m < -l)||(m > l)
	error('sphericalHarmonics(X,Y,Z,l,m,kind): <m> must be a integer between [-l,l]');
end

if strcmp(kind,'real')
	isRealSH = true;
elseif strcmp(kind,'complex')
	isRealSH = false;
else
	error('sphericalHarmonics(X,Y,Z,l,m,kind): <kind> must be either ''real'' or ''complex''');
end

if (isRealSH)
	Ylm = RealSphericalHarmonics(X,Y,Z,l,m);
else
	Ylm = ComplexSphericalHarmonics(X,Y,Z,l,m);
end

end




% ------------- Calculator for Real Spherical Harmonics
function Ylm = RealSphericalHarmonics(X,Y,Z,l,m)
% RealSphericalHarmonics calculates real spherical harmonics

if (l > 6)
	error('Only <l> less than or equal to 6 supported.');
end

r = sqrt(X.^2 + Y.^2 + Z.^2) ;

if (l == 0)
	% l=0
	C00 = 0.282094791773878; 	% 0.5*sqrt(1/pi)
	Ylm = C00 * ones(size(X));
elseif (l == 1)
	% l=1
	C1m1 = 0.488602511902920; 	% sqrt(3/(4*pi))
	C10 = 0.488602511902920; 	% sqrt(3/(4*pi))
	C1p1 = 0.488602511902920; 	% sqrt(3/(4*pi))
	if(m==-1)
		Ylm = C1m1 * (Y ./ r);
	elseif(m==0)
		Ylm = C10 * (Z ./ r);
	elseif(m==1)
		Ylm = C1p1 * (X ./ r);
	end
elseif(l == 2)
	% l=2
	C2m2 = 1.092548430592079; 	% 0.5*sqrt(15/pi)
	C2m1 = 1.092548430592079; 	% 0.5*sqrt(15/pi) 
	C20 =  0.315391565252520; 	% 0.25*sqrt(5/pi)
	C2p1 = 1.092548430592079; 	% 0.5*sqrt(15/pi) 
	C2p2 =  0.546274215296040;	% 0.25*sqrt(15/pi)
	if(m==-2)
		Ylm = C2m2*(X.*Y)./(r.*r);
	elseif(m==-1)
		Ylm = C2m1*(Y.*Z)./(r.*r);
	elseif(m==0)
		Ylm = C20*(-X.*X - Y.*Y + 2*Z.*Z)./(r.*r);
	elseif(m==1)
		Ylm = C2p1*(Z.*X)./(r.*r);
	elseif(m==2)
		Ylm = C2p2*(X.*X - Y.*Y)./(r.*r);
	end		
elseif(l == 3)
	% l=3
	C3m3 =  0.590043589926644;	% 0.25*sqrt(35/(2*pi))  
	C3m2 = 2.890611442640554;	% 0.5*sqrt(105/(pi))
	C3m1 = 0.457045799464466;	% 0.25*sqrt(21/(2*pi))
	C30 =  0.373176332590115; 	% 0.25*sqrt(7/pi)
	C3p1 =  0.457045799464466; 	% 0.25*sqrt(21/(2*pi))
	C3p2 = 1.445305721320277; 	% 0.25*sqrt(105/(pi))
	C3p3 = 0.590043589926644; 	% 0.25*sqrt(35/(2*pi))
	if(m==-3)
		Ylm = C3m3*(3*X.*X - Y.*Y).*Y./(r.*r.*r);
	elseif(m==-2)
		Ylm = C3m2*(X.*Y.*Z)./(r.*r.*r);
	elseif(m==-1)
		Ylm = C3m1*Y.*(4*Z.*Z - X.*X - Y.*Y)./(r.*r.*r);
	elseif(m==0)
		Ylm = C30*Z.*(2*Z.*Z-3*X.*X-3*Y.*Y)./(r.*r.*r);
	elseif(m==1)
		Ylm = C3p1*X.*(4*Z.*Z - X.*X - Y.*Y)./(r.*r.*r);
	elseif(m==2)
		Ylm = C3p2*Z.*(X.*X - Y.*Y)./(r.*r.*r);
	elseif(m==3)
		Ylm = C3p3*X.*(X.*X-3*Y.*Y)./(r.*r.*r);
	end	
elseif(l == 4)
	if(m==-4)
		Ylm=(3.0/4.0)*sqrt(35.0/pi)*(X.*Y.*(X.*X-Y.*Y))./(r.^4);
	elseif(m==-3)
		Ylm=(3.0/4.0)*sqrt(35.0/(2.0*pi))*(3.0*X.*X-Y.*Y).*Y.*Z./(r.^4);
	elseif(m==-2)
		Ylm=(3.0/4.0)*sqrt(5.0/pi)*X.*Y.*(7.0*Z.*Z-r.*r)./(r.^4);
	elseif(m==-1)
		Ylm=(3.0/4.0)*sqrt(5.0/(2.0*pi))*Y.*Z.*(7.0*Z.*Z-3.0*r.*r)./(r.^4);
	elseif(m==0)
		Ylm=(3.0/16.0)*sqrt(1.0/pi)*(35.0*Z.^4-30.0*Z.*Z.*r.*r+3.0*r.^4)./(r.^4);
	elseif(m==1)
		Ylm=(3.0/4.0)*sqrt(5.0/(2.0*pi))*X.*Z.*(7.0*Z.*Z-3.0*r.*r)./(r.^4);
	elseif(m==2)
		Ylm=(3.0/8.0)*sqrt(5.0/(pi))*(X.*X-Y.*Y).*(7.0*Z.*Z-r.*r)./(r.^4);
	elseif(m==3)
		Ylm=(3.0/4.0)*sqrt(35.0/(2.0*pi))*(X.*X-3.0*Y.*Y).*X.*Z./(r.^4);
	elseif(m==4)
		Ylm=(3.0/16.0)*sqrt(35.0/pi)*(X.*X.*(X.*X-3.0*Y.*Y) - Y.*Y.*(3.0*X.*X-Y.*Y))./(r.^4);
	end
elseif(l == 5)
	p = sqrt(X.*X+Y.*Y);                                     
	if(m==-5)
		Ylm = (3.0*sqrt(2*77/pi)/32.0)*(8.0*X.^4.*Y-4.0*X.*X.*Y.^3 + 4.0*Y.^5-3.0*Y.*p.^4)./(r.^5);
	elseif(m==-4)
		Ylm = (3.0/16.0)*sqrt(385.0/pi)*(4.0*X.^3.*Y - 4.0*X.*Y.^3).*Z./(r.^5);
	elseif(m==-3)
		Ylm = (sqrt(2.0*385.0/pi)/32.0)*(3.0*Y.*p.*p - 4.0*Y.^3).*(9*Z.*Z-r.*r)./(r.^5);
	elseif(m==-2)
		Ylm = (1.0/8.0)*sqrt(1155.0/pi)*2.0*X.*Y.*(3.0*Z.^3-Z.*r.*r)./(r.^5);
	elseif(m==-1)
		Ylm = (1.0/16.0)*sqrt(165.0/pi)*Y.*(21.0*Z.^4 - 14.0*r.*r.*Z.*Z+r.^4)./(r.^5);
	elseif(m==0)
		Ylm = (1.0/16.0)*sqrt(11.0/pi)*(63.0*Z.^5 -70.0*Z.^3.*r.*r + 15.0*Z.*r.^4)./(r.^5);
	elseif(m==1)
		Ylm = (1.0/16.0)*sqrt(165.0/pi)*X.*(21.0*Z.^4 - 14.0*r.*r.*Z.*Z+r.^4)./(r.^5);
	elseif(m==2)
		Ylm = (1.0/8.0)*sqrt(1155.0/pi)*(X.*X-Y.*Y).*(3.0*Z.^3 - r.*r.*Z)./(r.^5);
	elseif(m==3)
		Ylm = (sqrt(2.0*385.0/pi)/32.0)*(4.0*X.^3-3.0*p.*p.*X).*(9.0*Z.*Z-r.*r)./(r.^5);
	elseif(m==4)
		Ylm = (3.0/16.0)*sqrt(385.0/pi)*(4.0*(X.^4+Y.^4)-3.0*p.^4).*Z./(r.^5);
		%Ylm = (3.0/16.0)*sqrt(385.0/pi)*(X.^4+Y.^4-6*X.^2.*Y.^2).*Z./(r.^5);
	elseif(m==5)
		Ylm = (3.0*sqrt(2.0)/32.0)*sqrt(77.0/pi)*(4.0*X.^5 + 8.0*X.*Y.^4 -4.0*X.^3.*Y.*Y -3.0*X.*p.^4)./(r.^5);	
	end
elseif(l == 6)
	p = sqrt(X.*X+Y.*Y);                           
	if(m==-6)
		Ylm = (sqrt(2.0*3003.0/pi)/64.0)*(12.0*X.^5.*Y+12.0*X.*Y.^5 - 8.0*X.^3.*Y.^3-6.0*X.*Y.*p.^4)./(r.^6);
	elseif(m==-5)
		Ylm = (3.0/32.0)*sqrt(2.0*1001.0/pi)*(8.0*X.^4.*Y - 4.0*X.*X.*Y.^3 + 4.0*Y.^5 -3.0*Y.*p.^4).*Z./(r.^6);
	elseif(m==-4)
		Ylm = (3.0/32.0)*sqrt(91.0/pi)*(4.0*X.^3.*Y -4.0*X.*Y.^3).*(11.0*Z.*Z-r.*r)./(r.^6);
	elseif(m==-3)
		Ylm = (sqrt(2.0*1365.0/pi)/32.0)*(-4.0*Y.^3 + 3.0*Y.*p.*p).*(11.0*Z.^3 - 3.0*Z.*r.*r)./(r.^6);
	elseif(m==-2)
		Ylm = (sqrt(2.0*1365/pi)/64.0)*(2.0*X.*Y).*(33.0*Z.^4-18.0*Z.*Z.*r.*r + r.^4)./(r.^6);
	elseif(m==-1)
		Ylm = (sqrt(273.0/pi)/16.0)*Y.*(33.0*Z.^5-30.0*Z.^3.*r.*r +5.0*Z.*r.^4)./(r.^6);
	elseif(m==0)
		Ylm = (sqrt(13.0/pi)/32.0)*(231.0*Z.^6-315*Z.^4.*r.*r + 105.0*Z.*Z.*r.^4 -5.0*r.^6)./(r.^6);
	elseif(m==1)
		Ylm = (sqrt(273.0/pi)/16.0)*X.*(33.0*Z.^5-30.0*Z.^3.*r.*r +5*Z.*r.^4)./(r.^6);
	elseif(m==2)
		Ylm = (sqrt(2.0*1365/pi)/64.0)*(X.*X-Y.*Y).*(33.0*Z.^4 - 18.0*Z.*Z.*r.*r + r.^4)./(r.^6);
	elseif(m==3)
		Ylm = (sqrt(2.0*1365.0/pi)/32.0)*(4.0*X.^3 -3.0*X.*p.*p).*(11.0*Z.^3 - 3.0*Z.*r.*r)./(r.^6);
	elseif(m==4)
		Ylm = (3.0/32.0)*sqrt(91.0/pi)*(4.0*X.^4+4.0*Y.^4 -3.0*p.^4).*(11.0*Z.*Z -r.*r)./(r.^6);
	elseif(m==5)
		Ylm = (3.0/32.0)*sqrt(2.0*1001.0/pi)*(4.0*X.^5 + 8.0*X.*Y.^4-4.0*X.^3.*Y.*Y-3.0*X.*p.^4).*Z./(r.^6);
	elseif(m==6)
		Ylm = (sqrt(2.0*3003.0/pi)/64.0)*(4.0*X.^6-4.0*Y.^6 +12.0*X.*X.*Y.^4-12.0*X.^4.*Y.*Y + 3.0*Y.*Y.*p.^4-3.0*X.*X.*p.^4)./(r.^6);	
	end
else
	error('<l> must be an integer between [0,6]');
end


if (l > 0)
	Ylm(abs(r) < 1e-10) = 0;
end

end



% ------------- Calculator for Complex Spherical Harmonics
function Ylm = ComplexSphericalHarmonics(X,Y,Z,l,m)
% ComplexSphericalHarmonics calculates complex spherical harmonics


%fprintf('In Complex spherical harmonics\n');

r = sqrt(X.*X + Y.*Y + Z.*Z);
iota = sqrt(-1);

if l == 0
	
	C00 = 0.282094791773878; % 0.5*sqrt(1/pi)
	Ylm = C00*ones(size(X));	
elseif l == 1	
%fprintf('In Complex spherical harmonics with l = %d \n',l);
	C11 = 0.345494149471335; % 0.5 * sqrt(3/(2*pi))
	C10 = 0.488602511902920; % sqrt(3/(4*pi))
	
	if m==-1
		Ylm = C11 * ((X-iota*Y)./r);
	elseif m==0
		Ylm = C10*(Z./r);
	elseif m==1
		Ylm = -C11 * ((X+iota*Y)./r);
	else
		error('ERROR:incorrect l: %d,m: %d \n',l,m);
	end
elseif l == 2
%fprintf('In Complex spherical harmonics with l = %d \n',l);
	C22 = 0.386274202023190;
	C21 = 0.772548404046379;
	C20 = 0.315391565252520;
	
	r2 = r.^2;

	if m==-2
		Ylm = C22 * ((X-iota*Y).^2)./r2;
	elseif m==-1
		Ylm = C21 * ((X-iota*Y).*Z)./r2;
	elseif m==0
		Ylm = C20 * (2*Z.*Z - X.*X - Y.*Y)./r2;
	elseif m==1
		Ylm = -C21 * ((X+iota*Y).*Z)./r2;
	elseif m==2
		Ylm = C22 * ((X+iota*Y).^2)./r2;
	else
		error('ERROR: incorrect l: %d, m: %d \n',l,m);
	end
elseif l == 3
	C33 = 0.417223823632784;
	C32 = 1.021985476433282;
	C31 = 0.323180184114151;
	C30 = 0.373176332590115;

	r3 = r.^3;

	if m==-3
		Ylm = C33 * ((X-iota*Y).^3)./r3;
	elseif m==-2
		Ylm = C32 * (((X-iota*Y).^2).*Z)./r3;
	elseif m==-1
		Ylm = C31 * (X-iota*Y).*(4*Z.*Z - X.*X - Y.*Y)./r3;
	elseif m==0
		Ylm = C30 * Z.*(2*Z.*Z - 3*X.*X - 3*Y.*Y)./r3;
	elseif m==1
		Ylm = -C31 * (X+iota*Y).*(4*Z.*Z - X.*X - Y.*Y)./r3;
	elseif m==2
		Ylm = C32 * (((X+iota*Y).^2).*Z)./r3;
	elseif m==3
		Ylm = -C33 * ((X+iota*Y).^3)./r3;
	else
		error('ERROR: incorrect l: %d, m: %d \n',l,m);
	end    
elseif l == 4
	C44 = 0.442532692444983;
	C43 = 1.251671470898352;
	C42 = 0.334523271778645;
	C41 = 0.473087347878780;
	C40 = 0.105785546915204; 

	r4 = r.^4;

	if m==-4
		Ylm = C44 * ((X-iota*Y).^4)./r4;
	elseif m==-3
		Ylm = C43 * (((X-iota*Y).^3).*Z)./r4;
	elseif m==-2
		Ylm = C42 * (((X-iota*Y).^2).*(7*Z.^2 - r.^2))./r4;
	elseif m==-1
		Ylm = C41 * ((X-iota*Y).*Z.*(7*Z.^2 - 3*r.^2))./r4;
	elseif m==0
		Ylm = C40 * (35*Z.^4 - 30*(Z.^2).*(r.^2) + 3*r.^4)/r4;
	elseif m==1
		Ylm = -C41 * ((X+iota*Y).*Z.*(7*Z.^2 - 3*r.^2))./r4;
	elseif m==2
		Ylm = C42 * (((X+iota*Y).^2).*(7*Z.^2 - r.^2))./r4;
	elseif m==3
		Ylm = -C43 * (((X+iota*Y).^3).*Z)./r4;
	elseif m==4
		Ylm = C44 * ((X+iota*Y).^4)./r4;
	else
		error('ERROR: incorrect l: %d, m: %d \n',l,m);
	end 
elseif l == 5
	C55 = 0.464132203440858;
	C54 = 1.467714898305751; 
	C53 = 0.345943719146840;
	C52 = 1.694771183260899;
	C51 = 0.320281648576215;
	C50 = 0.116950322453424;

	r5 = r.^5;

	if m==-5
		Ylm = C55 * ((X-iota*Y).^5)./r5;
	elseif m==-4
		Ylm = C54 * (((X-iota*Y).^4).*Z)./r5;
	elseif m==-3
		Ylm = C53 * (((X-iota*Y).^3).*(9*Z.^2 - r.^2))./r5;
	elseif m==-2
		Ylm = C52 * (((X-iota*Y).^2).*(3*Z.^3 - Z.*(r.^2)))./r5;
	elseif m==-1
		Ylm = C51 * ((X-iota*Y).*(21*Z.^4 - 14*Z.*Z.*r.*r + r.^4))./r5;
	elseif m==0
		Ylm = C50 * (63*Z.^5 - 70*(Z.^3).*(r.^2) + 15*Z.*(r.^4))./r5;
	elseif m==1
		Ylm = -C51 * ((X+iota*Y).*(21*Z.^4 - 14*Z.*Z.*r.*r + r.^4))./r5;
	elseif m==2
		Ylm = C52 * (((X+iota*Y).^2).*(3*Z.^3 - Z.*(r.^2)))./r5;
	elseif m==3
		Ylm = -C53 * (((X+iota*Y).^3).*(9*Z.^2 - r.^2))./r5;
	elseif m==4
		Ylm = C54 * (((X+iota*Y).^4).*Z)./r5;
	elseif m==5
		Ylm = -C55 * ((X+iota*Y).^5)./r5;
	else
		error('ERROR: incorrect l: %d, m: %d \n',l,m);
	end
elseif l == 6
	C66 = 0.483084113580066; 
	C65 = 1.673452458100098;
	C64 = 0.356781262853998; 
	C63 = 0.651390485867716; 
	C62 = 0.325695242933858; 
	C61 = 0.411975516301141;
	C60 = 0.063569202267628;

	r6 = r.^6;

	if m==-6
		Ylm = C66 * ((X-iota*Y).^6)./r6;
	elseif m==-5
		Ylm = C65 * (((X-iota*Y).^5).*Z)./r6;
	elseif m==-4
		Ylm = C64 * (((X-iota*Y).^4).*(11*Z.^2 - r.^2))./r6;
	elseif m==-3
		Ylm = C63 * (((X-iota*Y).^3).*(11*Z.^3 - 3*Z.*(r.^2)))./r6;
	elseif m==-2
		Ylm = C62 * (((X-iota*Y).^2).*(33*Z.^4 - 18*(Z.^2).*(r.^2) + r.^4))./r6;
	elseif m==-1	
		Ylm = C61 * ((X-iota*Y).*(33*Z.^5 - 30*(Z.^3).*(r.^2) + 5*Z.*(r.^4)))./r6;
	elseif m==0
		Ylm = C60 * (231*Z.^6 - 315*(Z.^4).*(r.^2) + 105*(Z.^2).*(r.^4) - 5*r6)./r6;
	elseif m==1
		Ylm = -C61 * ((X+iota*Y).*(33*Z.^5 - 30*(Z.^3).*(r.^2) + 5*Z.*(r.^4)))./r6;
	elseif m==2
		Ylm = C62 * (((X+iota*Y).^2).*(33*Z.^4 - 18*(Z.^2).*(r.^2) + r.^4))./r6;
	elseif m==3
		Ylm = -C63 * (((X+iota*Y).^3).*(11*Z.^3 - 3*Z.*(r.^2)))./r6;
	elseif m==4
		Ylm = C64 * (((X+iota*Y).^4).*(11*Z.^2 - r.^2))./r6;
	elseif m==5
		Ylm = -C65 * (((X+iota*Y).^5).*Z)./r6;
	elseif m==6
		Ylm = C66 * ((X+iota*Y).^6)./r6;
	else
		error('Error: incorrect l: %d, m: %d \n',l,m);
	end
else
	error('Only l=0,1,supported. Input l:%d\n',l);
end 

if l > 0
	Ylm(abs(r)<(1e-10))=0;
end

end

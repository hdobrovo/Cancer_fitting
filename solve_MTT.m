#!/usr/bin/octave
# Exemple de solution d'equations differentielles

global odata;
global par
#odata = load("DOX.txt");

function eps = give_eps( drg, ic, emax, m )
    if( nargin == 3 )
        m = 1;
    end;
    eps = emax .* drg.^m ./ (drg.^m+ic.^m);
endfunction

function xdot = log_growth(x,t)
	global par;
	l = par(1); K = par(2);
	V = 1; 

	xdot = zeros(1,1);
	xdot(V) = l*x(V)*(1-x(V)/K);
endfunction

# Params
p = [0.651726 245211 0.00539614 0.812666]; #MCF-7 lambda
#p = [0.727047 492117 0.013937 1]; #HeLa lambda
#p = [0.727047 492117 0.00349155 1]; #Hela K
#p = [0.651726 245211 0.0018641 0.989111]; #MCF-7 lambda

x0 = [5000];

# Now plot the results of the fit
fd = fopen( "data/MCF7_aMTT.dat", "w" );

# Model solution
mtt=[];
drug_conc = p(3)*logspace(-3,4,100);
tvec = [0:0.01:2]';
epse = give_eps( drug_conc, p(3), p(4) );
par = [ p(1), p(2)];
y_base = lsode("log_growth",x0,tvec);
for i = 1:length(drug_conc)
	par = [ p(1)*(1-epse(i)), p(2)];
	y = lsode("log_growth",x0,tvec);
	mtt = [mtt; drug_conc(i) y(end)/y_base(end)];
end;

fprintf( fd, "%g %g\n", mtt' );
fclose(fd);


% The base of this file is from the page of CONICOPF at the address https://github.com/cbingane/conicopf
% It has been modified to create a DAT file in the Pyomo format from the
% MATPOWER problems. These DAT files can be used for the DiCARP software
% package. 


% INPUTS
%   casedata: MATPOWER case
%   model: either 0 for loss minimization or 1 for cost minimization
%   filename: The name of the file for the output DAT file, for example
%   case9.dat. 

function [n, slack, angslack, pL, qL, gs, bs, vl, vu,...
    nGen, pGl, pGu, qGl, qGu, c2, c1, c0, busgen,...
    nBranch, from, to, y, bsh, tap, shift, su, dl, du,...
    incidentF, incidentT, edges] = opf_data_file(casedata, model, filename)
casedata = loadcase(casedata);
mpc = ext2int(casedata);
fileID = fopen(filename, 'w');
%% baseMVA
base = mpc.baseMVA;

%% bus data
bus = mpc.bus;
n = size(bus,1);

% slack
slack = find(bus(:,2) == 3); angslack = bus(slack,9)*pi/180;
% loads
pL = sparse(bus(:,3)/base); qL = sparse(bus(:,4)/base);
% bounds Voltages
vl = bus(:,13); vu = bus(:,12);
% shunt compensators
gs = sparse(bus(:,5)/base); bs = sparse(bus(:,6)/base);
%
%% branch data
branch = mpc.branch; nBranch = size(branch,1);
e1=[26 29 45 89 90 91 92 93];
branch(e1,:)
branch(89,:)
% from to
from = branch(:,1); to = branch(:,2);
% series admittance
y = 1./(branch(:,3) + 1j*branch(:,4));
% shunt conductance
bsh = branch(:,5);
% tap
tap = sparse(branch(:,9));
% shift
shift = sparse(branch(:,10)*pi/180);
% flow limit
su = sparse(branch(:,6)/base);
% angle difference
dl = sparse(branch(:,12)*pi/180); du = sparse(branch(:,13)*pi/180);
% incidence
incidentF = sparse(1:nBranch,from,ones(nBranch,1),nBranch,n);
incidentT = sparse(1:nBranch,to,ones(nBranch,1),nBranch,n);
% adjacence
edges = unique(sort([from to],2), 'rows');
%
%% generator bounds
gen = mpc.gen; nGen = size(gen,1);
pGl = gen(:,10)/base; pGu = gen(:,9)/base;
qGl = gen(:,5)/base; qGu = gen(:,4)/base;
busgen = sparse(1:nGen,gen(:,1),ones(nGen,1),nGen,n);
%% generator cost coefficients
gencost = mpc.gencost;
if size(gencost,2) == 7
    c2 = sparse(gencost(:,5))*base^2; c1 = sparse(gencost(:,6))*base; c0 = sparse(gencost(:,7));
end
if size(gencost,2) == 6
    c2 = sparse(nGen,1)*base^2; c1 = sparse(gencost(:,5))*base; c0 = sparse(gencost(:,6));
end
if model == 0
    c2 = sparse(nGen,1)*base^2; c1 = ones(nGen,1)*base; c0 = sparse(nGen,1);
end


Yft = makeYft_c(nBranch,y,bsh,tap,shift);


fprintf(fileID, 'param baseMVA := %d; \n', base);
fprintf(fileID, 'param nv := %d; \n', n);
fprintf(fileID, 'param ne := %d; \n', nBranch);
fprintf(fileID, 'param ng := %d; \n\n', nGen);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, 'set gi :=');
gen_n = [];
cur = gen(1,1);
n_cur = 0;
for k=1:nGen
    if gen(k,1)==cur
        n_cur = n_cur+1;
    else
        cur = gen(k,1);
        n_cur = 1;
    end
    gen_n=[gen_n;n_cur];
end

for k=1:nGen
    fprintf(fileID, ' %d  %d ', gen(k,1), gen_n(k));
end
fprintf(fileID, '; \n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, 'param slack := %d; \n', slack);
fprintf(fileID, 'param angslack := %.8f; \n\n', angslack);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, 'param pL :=');
for k=1:n
    fprintf(fileID, '%d %.8f ', k, full(pL(k)));
end
fprintf(fileID, '; \n');
fprintf(fileID, 'param qL :=');
for k=1:n
    fprintf(fileID, '%d %.8f ', k, full(qL(k)));
end
fprintf(fileID, '; \n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, 'param vL :=');
for k=1:n
    fprintf(fileID, '%d %.8f ', k, full(vl(k)));
end
fprintf(fileID, '; \n');
fprintf(fileID, 'param vU :=');
for k=1:n
    fprintf(fileID, '%d %.8f ', k, full(vu(k)));
end
fprintf(fileID, '; \n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, 'param gs :=');
for k=1:n
    fprintf(fileID, '%d %.8f ', k, full(gs(k)));
end
fprintf(fileID, '; \n');
fprintf(fileID, 'param bs :=');
for k=1:n
    fprintf(fileID, '%d %.8f ', k, full(bs(k)));
end
fprintf(fileID, '; \n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, 'param fr :=');
for k=1:nBranch
    fprintf(fileID, '%d %d ', k, full(from(k)));
end
fprintf(fileID, '; \n');
fprintf(fileID, 'param to :=');
for k=1:nBranch
    fprintf(fileID, '%d %d ', k, full(to(k)));
end
fprintf(fileID, '; \n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, 'param Y : Y11r Y11i Y12r Y12i Y21r Y21i Y22r Y22i := \n');
for k=1:nBranch
    fprintf(fileID, '%d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f \n', k, real(Yft{k}(1,1)), imag(Yft{k}(1,1)), real(Yft{k}(1,2)), imag(Yft{k}(1,2)), real(Yft{k}(2,1)), imag(Yft{k}(2,1)), real(Yft{k}(2,2)), imag(Yft{k}(2,2)));
end
fprintf(fileID, '; \n');



fprintf(fileID, 'param e_con :=');
for k=1:nBranch
    fprintf(fileID, '%d %.8f ', k, real(y(k)));
end
fprintf(fileID, '; \n');
fprintf(fileID, 'param e_sus :=');
for k=1:nBranch
    fprintf(fileID, '%d %.8f ', k, imag(y(k)));
end
fprintf(fileID, '; \n');
fprintf(fileID, 'param e_bsh :=');
for k=1:nBranch
    fprintf(fileID, '%d %.8f ', k, (bsh(k)));
end
fprintf(fileID, '; \n');
fprintf(fileID, 'param e_su :=');
for k=1:nBranch
    fprintf(fileID, '%d %.8f ', k, full(su(k)));
end
fprintf(fileID, '; \n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, 'param e_dl :=');
for k=1:nBranch
    fprintf(fileID, '%d %.8f ', k, full(dl(k)));
end
fprintf(fileID, '; \n');
fprintf(fileID, 'param e_du :=');
for k=1:nBranch
    fprintf(fileID, '%d %.8f ', k, full(du(k)));
end
fprintf(fileID, '; \n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, 'param pGl :=');
for k=1:nGen
    fprintf(fileID, '%d %d %.8f ', gen(k,1), gen_n(k), pGl(k));
end
fprintf(fileID, '; \n');
fprintf(fileID, 'param pGu :=');
for k=1:nGen
    fprintf(fileID, '%d %d %.8f ', gen(k,1), gen_n(k), pGu(k));
end
fprintf(fileID, '; \n');
fprintf(fileID, 'param qGl :=');
for k=1:nGen
    fprintf(fileID, '%d %d %.8f ', gen(k,1), gen_n(k), qGl(k));
end
fprintf(fileID, '; \n');
fprintf(fileID, 'param qGu :=');
for k=1:nGen
    fprintf(fileID, '%d %d %.8f ', gen(k,1), gen_n(k), qGu(k));
end
fprintf(fileID, '; \n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, 'param c2 :=');
for k=1:nGen
    fprintf(fileID, '%d %d %.8f ', gen(k,1), gen_n(k), full(c2(k)));
end
fprintf(fileID, '; \n');
fprintf(fileID, 'param c1 :=');
for k=1:nGen
    fprintf(fileID, '%d %d %.8f ', gen(k,1), gen_n(k), full(c1(k)));
end
fprintf(fileID, '; \n');
fprintf(fileID, 'param c0 :=');
for k=1:nGen
    fprintf(fileID, '%d %d %.8f ', gen(k,1), gen_n(k), full(c0(k)));
end
fprintf(fileID, '; \n');


fclose(fileID);
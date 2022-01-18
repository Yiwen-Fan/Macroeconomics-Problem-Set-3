clear;
clc;

global CHI_C CHI_N VARPHI THETA BETA cPHIpi ELB OMEGA TAU PIstar;
CHI_C = 1; CHI_N = 0.5; VARPHI = 100;
THETA = 6; BETA = 0.99; cPHIpi = 2; ELB = 1; OMEGA= 0.5;
TAU = 1/(THETA-1);

%Set the range of PIstar
PIstarVals = 0.9975:0.0005:1.0025;
NVals   = size(PIstarVals,2);

%Order of the output variables of fsolve function:
fsolveVar = ["C","Y","N","w","Pi","R"];
Nvars = size(fsolveVar,2);

fun1 = @TAYLOR;
fun2 = @TAYLOR_ZLB;
x0 = ones(Nvars,1);
result = zeros(Nvars, NVals);
result_ZLB = zeros(Nvars, NVals);

for i = 1:NVals
    PIstar = PIstarVals(i);
    result(1:Nvars, i) = fsolve(fun1,x0);
    result_ZLB(1:Nvars, i) = fsolve(fun2,x0);
end

% Welfare
if CHI_C == 1
    V_STD = (log(result(find(fsolveVar=="C"), 1:NVals)) - (result(find(fsolveVar=="N"), 1:NVals).^(1+CHI_N))/(1+CHI_N))./(1-BETA);
    V_ZLB = (log(result_ZLB(find(fsolveVar=="C"), 1:NVals)) - (result_ZLB(find(fsolveVar=="N"), 1:NVals).^(1+CHI_N))/(1+CHI_N))./(1-BETA);
else
    V_STD = (((result(find(fsolveVar=="C"), 1:NVals)).^(1-CHI_C))/(1-CHI_C) - (result(find(fsolveVar=="N"), 1:NVals).^(1+CHI_N))/(1+CHI_N))./(1-BETA);
    V_ZLB = (((result_ZLB(find(fsolveVar=="C"), 1:NVals)).^(1-CHI_C))/(1-CHI_C) - (result_ZLB(find(fsolveVar=="N"), 1:NVals).^(1+CHI_N))/(1+CHI_N))./(1-BETA);
end

%Price Adjustment Cost (Y-C)
AC_STD = result(find(fsolveVar=="Y"), 1:NVals) - result(find(fsolveVar=="C"), 1:NVals);
AC_ZLB = result_ZLB(find(fsolveVar=="Y"), 1:NVals) - result_ZLB(find(fsolveVar=="C"), 1:NVals);

figvar = ["Pi","R","Y","C","w","AC","V"];
figtitle = ["Inflation Rates (\Pi)","Policy Rate (R)","Output (Y)","Consumption (C)","Real Wage (w)",...
    "Price Adjustment Cost (Y-C)","Welfare (V)"];

set(0,'DefaultAxesTitleFontWeight','normal');
fig(1) = figure(1);


for ifig=1:1:size(figvar,2)
    subplot(3,3,ifig);
    if ifig < 6
        plot(0:NVals-1, result(find(fsolveVar==figvar(ifig)), 1:NVals), 'black','LineWidth',2); hold on
        plot(0:NVals-1, result_ZLB(find(fsolveVar==figvar(ifig)), 1:NVals), 'red','LineWidth',2);
        if figvar(ifig) == "R"
            vecRstar = zeros(1,NVals); vecRstar(1,1:NVals) = 1/BETA;
            plot(0:NVals-1, vecRstar(1,1:NVals), 'k--'); hold off
        elseif figvar(ifig) == "Pi" | figvar(ifig) == "w"
            plot(0:NVals-1, ones(1,NVals), 'k--'); hold off
        elseif figvar(ifig) == "Y"
            ylim([0.9998 1.0004])
        end
    elseif figvar(ifig) == "AC"
        plot(0:NVals-1, AC_STD(1, 1:NVals), 'black','LineWidth',2); hold on
        plot(0:NVals-1, AC_ZLB(1, 1:NVals), 'red','LineWidth',2); hold off
        ylim([-0.01, 0.01])
    elseif figvar(ifig) == "V"
        L(1)=plot(0:NVals-1, V_STD, 'black', 'LineWidth',2); hold on
        L(2)=plot(0:NVals-1, V_ZLB, 'red', 'LineWidth',2); hold off
    end
    title(figtitle(ifig))
    xlabel("Inflation Target (\Pi^{targ})")
    xticklabels([min(PIstarVals) (min(PIstarVals)+max(PIstarVals))/2 max(PIstarVals)])
end

sh=subplot(3,3,9);
p=get(sh,'position');
lh=legend(sh,L,"Standard Equilibrium","Deflationary Equilibrium");
set(lh,'position',p);
axis(sh,'off');
sgtitle('Steady State for the Inflation Target \Pi^{targ} \in [0.9975, 1.0025]')
saveas(gcf,'Figure32.jpg');


function F = TAYLOR(x)

global CHI_C CHI_N VARPHI THETA BETA cPHIpi ELB OMEGA TAU PIstar;

Css=x(1);
Yss=x(2);
Nss=x(3);
wss=x(4);
Piss=x(5);
Rss=x(6);

% FONCs w.r.t. Lagrange multipliers (i.e. private sector equilibrium conditions)
F(1) = Rss - 1/BETA*Piss;
F(2) = wss - (Nss^CHI_N * Css^CHI_C);
F(3) = (Yss*(Css^(-CHI_C)))*(VARPHI*(Piss^(1-OMEGA) - 1)*Piss^(1-OMEGA) - (1+TAU)*(1-THETA) - THETA*wss) - (BETA*(Yss*(Css^(-CHI_C)))*VARPHI*(Piss^(1-OMEGA)-1)*Piss^(1-OMEGA));
F(4) = Yss - (Css + (VARPHI / 2)*((Piss^(1-OMEGA)-1)^2)*Yss);
F(5) = Yss - Nss;
% FONCs w.r.t. private sector variables
F(6) = Rss - PIstar/BETA*(Piss/PIstar)^cPHIpi;
%F(6) = Rss - ELB;
end

function F = TAYLOR_ZLB(x)

global CHI_C CHI_N VARPHI THETA BETA cPHIpi ELB OMEGA TAU PIstar;

Css=x(1);
Yss=x(2);
Nss=x(3);
wss=x(4);
Piss=x(5);
Rss=x(6);

% FONCs w.r.t. Lagrange multipliers (i.e. private sector equilibrium conditions)
F(1) = Rss - 1/BETA*Piss;
F(2) = wss - (Nss^CHI_N * Css^CHI_C);
F(3) = (Yss*(Css^(-CHI_C)))*(VARPHI*(Piss^(1-OMEGA) - 1)*Piss^(1-OMEGA) - (1+TAU)*(1-THETA) - THETA*wss) - (BETA*(Yss*(Css^(-CHI_C)))*VARPHI*(Piss^(1-OMEGA)-1)*Piss^(1-OMEGA));
F(4) = Yss - (Css + (VARPHI / 2)*((Piss^(1-OMEGA)-1)^2)*Yss);
F(5) = Yss - Nss;
% FONCs w.r.t. private sector variables
%F(6) = Rss - PIstar/BETA*(Piss/PIstar)^cPHIpi;
F(6) = Rss - ELB;
end


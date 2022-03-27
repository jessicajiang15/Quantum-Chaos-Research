(* ::Package:: *)

(* ::Input::Initialization:: *)
 


(* ::Input::Initialization:: *)
currTime=0;
totalSteps=100;
tf=30;
ti=0.1;
\[HBar]=1;
\[Omega]=0.5;
m=1;

xmin=-5;xmax=5;pmin=-5;pmax=5;
nbox=40;
npoints=2;
cutoff=0.0001;
nmax=10;
L=1;


(* ::Input::Initialization:: *)
V[x_]=1/2 m \[Omega]^2 x^2;


(* ::Input::Initialization:: *)
hamonicOscEnergy[n_]:=\[HBar] \[Omega](1/2+n);
calculateTimeFactor[n_, t_]:=E^(-I hamonicOscEnergy[n] t/\[HBar])
generateNthHarmonicOscillatorEigenstateWithT[n_]:=1/Sqrt[2^n n!] ((m  \[Omega])/(\[Pi] \[HBar]))^(1/4) E^(-(m \[Omega] x^2)/(2 \[HBar])) HermiteH[n, Sqrt[m \[Omega]/\[HBar]]x]calculateTimeFactor[n, t];


(* ::Input::Initialization:: *)
generateNthInfiniteSquareWell[n_]:=Sqrt[2/L]Sin[n x/L]


(* ::Input::Initialization:: *)
currTime=0;
tStep=0.1;
\[HBar]=1;
\[Omega]=0.5;
m=1;
nMin=1;
nMax=10;
num=10;

(*when := you have a different definition*)
generateNthHarmonicOscillatorEigenstate[n_]:=1/Sqrt[2^n n!] ((m  \[Omega])/(\[Pi] \[HBar]))^(1/4) E^(-(m \[Omega] x^2)/(2 \[HBar])) HermiteH[n, Sqrt[m \[Omega]/\[HBar]]x];



(*Function to take in wavefunction and convert into a Wigner Function

to ask: is this valid syntax*)
WignerFunction[f_, t_]:=1/(\[Pi] \[HBar]) Integrate[Conjugate[f[(x+y), t]]f[(x-y), t]E^(2 I p y/\[HBar]), {y, -\[Infinity], \[Infinity]},Assumptions->{x\[Element]Reals, p\[Element]Reals, t\[Element] Reals}]

hamonicOscEnergy[n_]:=\[HBar] \[Omega](1/2+n);
calculateTimeFactor[n_, t_]:=E^(-I hamonicOscEnergy[n] t/\[HBar])

generateNthHarmonicOscillatorEigenstateWithT[n_]:=1/Sqrt[2^n n!] ((m  \[Omega])/(\[Pi] \[HBar]))^(1/4) E^(-(m \[Omega] x^2)/(2 \[HBar])) HermiteH[n, Sqrt[m \[Omega]/\[HBar]]x]calculateTimeFactor[n, t];


A=1/Sqrt[2]
B=1/Sqrt[2]
C1=1/Sqrt[3]
D1=Sqrt[2/3]
\[CapitalPsi]1[x_, t_]=A*generateNthHarmonicOscillatorEigenstateWithT[1]+B*generateNthHarmonicOscillatorEigenstateWithT[2];

\[CapitalPsi]2[x_, t_]=C1*generateNthHarmonicOscillatorEigenstateWithT[1]+D1*generateNthHarmonicOscillatorEigenstateWithT[2];
WignerFunction[f_, t_]:=1/(\[Pi] \[HBar]) Integrate[Conjugate[f[(x+y), t]]f[(x-y), t]E^(2 I p y/\[HBar]), {y, -\[Infinity], \[Infinity]},Assumptions->{x\[Element]Reals, p\[Element]Reals, t\[Element] Reals}]
WignerFunctionWDomain[f_, t_]:=1/(\[Pi] \[HBar]) Integrate[Conjugate[f[(x+y), t]]f[(x-y), t]E^(2 I p y/\[HBar]), {y, -\[Infinity], \[Infinity]},Assumptions->{x\[Element]Reals, p\[Element]Reals, t\[Element] Reals}]


(* ::Input::Initialization:: *)
A


(* ::Input::Initialization:: *)
ti=0.5
W1[x_, p_,t_]=WignerFunction[\[CapitalPsi]1, ti]
W2[x_, p_,t_]=WignerFunction[\[CapitalPsi]2, ti]


(* ::Input::Initialization:: *)
Clear[p,x,W]
timeEvolve[f_,V_, ti_,tf_, nmax_]:=NDSolve[{D[W[x, p, t],t]==-p/m D[W[x, p, t],x]+Sum[(1/2)^(2n) 1/(2n+1)! (D[V[x],{x,2n+1} ])(D[W[x,p,t],{p,2n+1}]),{n,0,nmax}], W[xmin,p,t]==0,W[xmax,p,t]==0, W[x,pmin,t]==0,W[x,pmax,t]==0,W[x,p,0]==f[x,p, 0]},W[x,p,t],{x,xmin,xmax},{p,pmin,pmax},{t,ti, tf}];


(* ::Input::Initialization:: *)
snap=timeEvolve[W1, V, ti, tf, nmax]
snap1=timeEvolve[W2, V,ti,  tf, nmax]


(* ::Input::Initialization:: *)
f[x_,p_,t_]=W[x,p,t]/.snap[[1]]
g[x_,p_,t_]=W[x,p,t]/.snap1[[1]]


(* ::Input::Initialization:: *)
Plot3D[f[x,p,0],{x,xmin,xmax},{p,pmin,pmax},PlotRange->All]


(* ::Input::Initialization:: *)



(* ::Input::Initialization:: *)
temp=Table[Plot3D[f[x,p,t],{x,xmin,xmax},{p,pmin,pmax},PlotRange->All],{t,0,tf}]


(* ::Input::Initialization:: *)
Export["thetestt.gif",temp]

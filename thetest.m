(* ::Package:: *)

(* ::Input::Initialization:: *)
 


(* ::Input::Initialization:: *)
currTime=0;
totalSteps=100;
tf=10.1;
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
calculateEMD[finitnorm_, ffinalnorm_, xmin_, xmax_,ymin_, ymax_, nbox_, gridboxcutoff_]:=Module[{dx, dy, diffarray, outboxes, inboxes, nout, nin, nvars, mat0, supplyamount, demandamount, mat, distance, k, l, m, n, i, j},dx=(xmax-xmin)/nbox;dy=(ymax-ymin)/nbox;
f1[x_,p_]=finitnorm;
f2[x_,p_]=ffinalnorm;
Print["hello"];
Monitor[diffarray=Table[NIntegrate[Re[f[x,p,0.1]-g[x,p,0.1]],{x,xmin+(i-1)*dx,xmin+i*dx},{p,pmin+(j-1)*dy,pmin+j*dy},Method->{"AdaptiveMonteCarlo",MaxPoints->10^7}],{i,1,nbox},{j,1,nbox}],i];
Print["done"];
outboxes={};inboxes={}; 
For[i=1,i<=nbox,i++,
For[j=1,j<=nbox,j++,
If[diffarray[[i,j]]>gridboxcutoff,AppendTo[outboxes,{xmin+(i-.5)*dx,ymin+(j-.5)*dy}]];
If[diffarray[[i,j]]<-gridboxcutoff,AppendTo[inboxes,{xmin+(i-.5)*dx,ymin+(j-.5)*dy}]]]];
nout=Length[outboxes];nin=Length[inboxes];
nvars=nout*nin;
supplyamount ={}; demandamount = {};
For[i=1,i<=nbox,i++,
For[j=1,j<=nbox,j++,
If[diffarray[[i,j]]>gridboxcutoff,AppendTo[supplyamount,diffarray[[i,j]]]];
If[diffarray[[i,j]]<-gridboxcutoff,AppendTo[demandamount,diffarray[[i,j]]]];
]];
mat = Table[
EuclideanDistance[outboxes[[o]], inboxes[[q]]],{o,Length@outboxes},{q, Length@inboxes}];
mat0 = Table[0,{Length@outboxes+Length@inboxes},{Length@outboxes+Length@inboxes}];
For[i=1,i<=Length@outboxes,i++,
For[j=Length@outboxes+1,j<=(Length@inboxes+Length@outboxes),j++,
mat0[[i,j]]=mat[[i,j-Length@outboxes]]
]];

mat0//MatrixForm;

(*Then put in the numbers for the upper left section*)
For[k=1,k<=Length@outboxes,k++,
For[l=k+1,l<=Length@outboxes,l++,
mat0[[k,l]]=EuclideanDistance[outboxes[[k]], outboxes[[l]]]
]];

mat0//MatrixForm;
(*For the lower right..*)
For[m=Length@outboxes+1,m<=(Length@inboxes+Length@outboxes-1),m++,
For[n=m+1,n<=(Length@inboxes+Length@outboxes),n++,
mat0[[m,n]]=EuclideanDistance[inboxes[[m-Length@outboxes]], inboxes[[n-Length@outboxes]]]
]];

mat0//MatrixForm;
distance=FindMinimumCostFlow[ mat0, Join[supplyamount,demandamount]];
Print["Done one"];
Return[distance];
]




(* ::Input::Initialization:: *)
nbox=40;
gridboxcutoff=.0001;
distances=Table[calculateEMD[f[x,p,t], g[x, p, t], xmin, xmax,pmin, pmax, nbox, gridboxcutoff],{t,ti, tf,1}]
plt=ListPlot[distances]
Export["plotnbox40.gif",plt]
Export["numbersnbox40.csv",distances,"CSV"]



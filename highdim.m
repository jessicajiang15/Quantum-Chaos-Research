(* ::Package:: *)

(* ::Input::Initialization:: *)
begint=AbsoluteTime[]


(* ::Title::Initialization:: *)
(*Numerical Calculation: FindMininumCostFlow*)


(* ::Text::Initialization:: *)
(*This notebook is designed to take two distributions of sand in a specified box, Subscript[f, init](x,y) and Subscript[f, final](x,y) and calculate the cost of moving the sand from the initial to the final distribution. Cost is defined as the amount of sand moved times the distance.*)


(* ::Text::Initialization:: *)
(*The FindMinimumCostFlow function will output the minimum total cost of flowing moving from supply to demand. *)


(* ::Section::Initialization:: *)
(*Enter Parameters*)


(* ::Text::Initialization:: *)
(*Define the initial and final sand configurations. (These will be automatically normalized below.)*)


(* ::Input::Initialization::Plain:: *)
x1=0;y1=0;
x2=0.1;y2=0;
\[Rho]=1; p=1; \[Theta]=0;


(* ::Input::Initialization:: *)
aa=((Cos[\[Theta]])^2/p^2+(p^2 (Sin[\[Theta]])^2)/\[Rho]^2);
bb=2(Sin[\[Theta]])*(Cos[\[Theta]])*(1/p^2-p^2/\[Rho]^2);
cc=((Sin[\[Theta]])^2/p^2+(p^2 (Cos[\[Theta]])^2)/\[Rho]^2);
finit[x_,y_,z_,w_]=1/\[Pi] E^(-(x)^2-(y)^2/\[Rho]^2-z^2-w^2)
(*ffinal[x_,y_]=1/\[Pi]E^((-(x-x2)^2/(p^2))-((p^2*(y-y2)^2)/(\[Rho]^2)) )*)
ffinal[x_,y_,z_,w_]=1/\[Pi] E^-(aa*(x-x2)^2+bb*(x-x2)*(y-y2)+cc*(y-y2)^2+z^2+w^2)


(* ::Input::Initialization:: *)
Plot3D[{ffinal[x,y,0,0],finit[x,y,0,0]},{x,-8,8},{y,-5,5},PlotRange->All]


(* ::Text::Initialization:: *)
(*Set parameters for the sandbox*)
(*  xmin, xmax, ymin, and ymax define the borders of the sandbox*)
(*  nbox is the number of gridpoints on a side. The sandbox will be divided into an nbox x nbox grid and sand will be moved between those gridpoints.*)


(* ::Input::Initialization:: *)
xmin=-4;xmax=4;ymin=-4;ymax=4;zmin=-4;zmax=4;wmin=-4;wmax=4;
(*running 40 nboxes causes findminimumcostflow to stall. 50 nboxes finishes running.*)
nbox=5;


(* ::Input::Initialization:: *)
zmin


(* ::Text::Initialization:: *)
(*Set parameters for the optimization*)
(*  gridboxcutoff is the minimum amount of sand that a box must contain to be included. Specifically, if the difference between initial and final amounts of sand in that box is less than gridboxcutoff then the box is excluded from the calculation. Set this to 0 to include all boxes.*)


(* ::Input::Initialization:: *)
gridboxcutoff=.0001/(nbox^2/20)^2
(*scale to nbox
divide by nbox square*)


(* ::Input::Initialization:: *)
finitnorm[x,y,z,y]


(* ::Input::Initialization:: *)
finit[x,y,z,w]
ffinal[x,y,z,w]


(* ::Input::Initialization:: *)
Method->{"AdaptiveQuasiMonteCarlo",MaxPoints->10^10}


(* ::Section::Initialization:: *)
(*Calculate Derived Quantities*)


(* ::Text::Initialization:: *)
(*Normalize the initial and final sand distributions within the total box.*)


(* ::Input::Initialization:: *)
finitnorm[x_,y_,z_,w_]=finit[x,y,z,w]/NIntegrate[finit[x,y,z,w],{x,xmin,xmax},{y,ymin,ymax},{z,zmin,zmax},{w,wmin,wmax},Method->{"GlobalAdaptive"}];
ffinalnorm[x_,y_,z_,w_]=ffinal[x,y,z,w]/NIntegrate[ffinal[x,y,z,w],{x,xmin,xmax},{y,ymin,ymax}, {z,zmin,zmax},{w,wmin,wmax},Method->{"GlobalAdaptive"}];


(* ::Input::Initialization:: *)
finitnorm[x,y,z,w]


(* ::Input::Initialization:: *)
ffinalnorm[x,y,z,w]


(* ::Text::Initialization:: *)
(*Calculate the amount of sand that needs to be moved from each box of the grid and store the results in an nbox x nbox grid called diffarray.*)
(*This is the difference between the initial and final amounts of sand, so a positive amount means sand needs to be moved away from this box (sources) and a negative amount means sand needs to be moved into this box (targets)*)


(* ::Input::Initialization:: *)
dx=(xmax-xmin)/nbox;dy=(ymax-ymin)/nbox;dz=(zmax-zmin)/nbox;dw=(wmax-wmin)/nbox;
diffarray=ParallelTable[NIntegrate[finitnorm[x,y,z,w]-ffinalnorm[x,y,z,w],{x,xmin+(i-1)*dx,xmin+i*dx},{y,ymin+(j-1)*dy,ymin+j*dy},{z,zmin+(k-1)*dz,zmin+k*dz},{w,wmin+(l-1)*dw,wmin+l*dw},Method->{"GlobalAdaptive"}],{i,1,nbox},{j,1,nbox},{k,1,nbox},{l,1,nbox}];


(* ::Input::Initialization:: *)
diffarray


(* ::Input::Initialization:: *)
diffarray


(* ::Section::Initialization:: *)
(*Create Necessary Inputs to run FindMinimumCostFlow*)


(* ::Text::Initialization:: *)
(**)
(*Create lists of all the boxes from which sand needs to be moved out (outboxes) and to which sand needs to be moved in (inboxes).*)
(*Each entry in one of those lists will be a list of two numbers: the x and y coordinates of the center of the grid box.*)


(* ::Input::Initialization:: *)
Dimensions[diffarray]


(* ::Input::Initialization:: *)
gridboxcutoff=0.0001


(* ::Input::Initialization:: *)
outboxes={};inboxes={}; supplyamount ={}; demandamount = {};
For[i=1,i<=nbox,i++,
For[j=1,j<=nbox,j++,
For[k=1,k<=nbox,k++,
For[l=1,l<=nbox,l++,
If[Abs[diffarray[[i,j,k,l]]]>gridboxcutoff,AppendTo[outboxes,{xmin+(i-.5)*dx,ymin+(j-.5)*dy,zmin+(k-.5)*dz,wmin+(l-.5)*dw}];
AppendTo[supplyamount,diffarray[[i,j,k,l]]]
];
If[diffarray[[i,j,k,l]]<-gridboxcutoff,AppendTo[inboxes,{xmin+(i-.5)*dx,ymin+(j-.5)*dy,zmin+(i-.5)*dz,wmin+(i-.5)*dw}];
AppendTo[demandamount,diffarray[[i,j,k,l]]]
]

]
]
]];





(* ::Input::Initialization:: *)
diffarray[[1,1,1,1]]


(* ::Input::Initialization:: *)
inboxes


(* ::Text::Initialization:: *)
(*The variables nout and nin are the numbers of boxes from which sand will be moved out or in.*)
(* Their product, nvars, represents the total number of possible sand movements.*)


(* ::Input::Initialization:: *)
nout=Length[outboxes];nin=Length[inboxes];
nvars=nout*nin;


(* ::Text::Initialization:: *)
(* *)
(*Creating a graph of vertices and edges contributing to the flow (bipartite graph of two sets of vertices: source set and target sets) *)
(*EdgeCost: Assigning cost to each edge which is calculated by the Euclidean Distance between the two vertices. *)


(* ::Input::Initialization:: *)



(* ::Input::Initialization:: *)
nout+nin


(* ::Input::Initialization:: *)



(* ::Input::Initialization:: *)
diffarray


(* ::Input::Initialization:: *)
outboxes


(* ::Section::Initialization:: *)
(*New codes begin here*)


(* ::Input::Initialization:: *)
Clear[mat0,i,j,mat,k,l,m,n]


(* ::Input::Initialization:: *)
(*Create the cost matrix with rows being all supplying boxes and columns being all demanding boxes.*)
mat = Table[
EuclideanDistance[outboxes[[o]], inboxes[[q]]],{o,Length@outboxes},{q, Length@inboxes}];

(*FindMinimumCostFlow method requires the columns and the rows to contain all boxes,no matter supplying or demanding, though I think the method will not actually use the diagonally lower half of the matrix, so we can keep elements in that part quatiling to 0. 
Now create a zero matrix with # of columns=# of rows=# of all boxes, and put the cost matrix from previous section onto the upper right part of the zero matrix,where the rows are supplying boxes and columns are demanding boxes.*)
mat0 = Table[0,{Length@outboxes+Length@inboxes},{Length@outboxes+Length@inboxes}];

For[i=1,i<=Length@outboxes,i++,
For[j=Length@outboxes+1,j<=(Length@inboxes+Length@outboxes),j++,
mat0[[i,j]]=mat[[i,j-Length@outboxes]]
]]

For[m=Length@outboxes+1,m<=(Length@inboxes+Length@outboxes-1),m++,
For[i=1,i<=Length@outboxes,i++,
mat0[[i,m]]=mat[[i,m-Length@outboxes]]
]]

mat0//MatrixForm;


(* ::Input::Initialization:: *)
Dimensions[mat0]


(* ::Input::Initialization:: *)
(*Then put in the numbers for the upper left section*)
For[k=1,k<=Length@outboxes,k++,
For[l=k+1,l<=Length@outboxes,l++,
mat0[[k,l]]=EuclideanDistance[outboxes[[k]], outboxes[[l]]]
]]

mat0//MatrixForm;



(* ::Input::Initialization:: *)
(*For the lower right..*)
For[m=Length@outboxes+1,m<=(Length@inboxes+Length@outboxes-1),m++,
For[n=m+1,n<=(Length@inboxes+Length@outboxes),n++,
mat0[[m,n]]=EuclideanDistance[inboxes[[m-Length@outboxes]], inboxes[[n-Length@outboxes]]]
]]

mat0//MatrixForm;


(* ::Input::Initialization:: *)
Dimensions[mat0]


(* ::Input::Initialization:: *)
Export["supplyamount.csv",supplyamount,"CSV"]
Export["demandamount.csv",demandamount,"CSV"]
Export["mat0.csv",mat0,"CSV"]


(* ::Input::Initialization:: *)
(*Find the minimum cost flow*)

distance=FindMinimumCostFlow[ mat0, Join[supplyamount,demandamount]]


(* ::Input::Initialization:: *)
totaltime=AbsoluteTime[]-begint;


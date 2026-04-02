(* ::Package:: *)

BeginPackage["IdealGas`"]
CreateSimulationBox::usage="Input the linear dimension of the simulation box (# cells, length 2.22). Identifies voxels & neighbourhoods."
TorusDistance::usage="calculates distance between particles on a torus"
PopulateSimulationBox::usage="populates simulation box with randomly placed particles"
TorusExtractTrajectories::usage="time evolves the state for a specified number of steps"
ConstructCollisionArray::usage="constructs collision array from trajectory and state"
IdentifyVoxel::usage="test"
ExtractVoxelDirectory::usage="test"
NeighborhoodVoxels::usage="test"
Begin["Private`"]


\[CapitalDelta]l=20/9;
PrototypeNeighborhood=Tuples[{0,-\[CapitalDelta]l,\[CapitalDelta]l},3];

LoopBack3[{x_,y_,z_},L_]:=Module[{\[DoubleStruckP],correction},
\[DoubleStruckP]=2Abs[#]/L&/@{x,y,z}//Floor;
correction=-L \[DoubleStruckP] Sign[{x,y,z}];
{x,y,z}+correction
];

TorusDistance[{{v1_,v2_,v3_},{w1_,w2_,w3_}},L_]:=Block[{\[CapitalDelta]x,\[CapitalDelta]y,\[CapitalDelta]z},
\[CapitalDelta]x=Min[Abs[v1-w1],L-Abs[v1-w1]];
\[CapitalDelta]y=Min[Abs[v2-w2],L-Abs[v2-w2]];
\[CapitalDelta]z=Min[Abs[v3-w3],L-Abs[v3-w3]];
Sqrt[# . #]&@{\[CapitalDelta]x,\[CapitalDelta]y,\[CapitalDelta]z}
];

CreateSimulationBox[cells_]:=Block[{\[ScriptCapitalL],DisplacementVectors,LinearBoundaries,VoxelBoundaries,VoxelCentres,Neighborhoods},
\[ScriptCapitalL]=\[CapitalDelta]l*cells;
DisplacementVectors=\[ScriptCapitalL] Range[-1/2,1/2,1/cells]//Tuples[#,3]&;
Neighborhoods=Map[LoopBack3[#,\[ScriptCapitalL]+\[ScriptCapitalL]/cells]&,((PrototypeNeighborhood+ConstantArray[#,27])&/@DisplacementVectors),{2}];
LinearBoundaries=Range[-\[ScriptCapitalL]/2-\[ScriptCapitalL]/(2*cells),\[ScriptCapitalL]/2+\[ScriptCapitalL]/(2*cells),\[CapitalDelta]l]//Partition[#,2,1]&;
VoxelBoundaries=Tuples[LinearBoundaries,3];
VoxelCentres=Flatten[Neighborhoods,1]//DeleteDuplicates;
<|"VoxelBoundaries"->VoxelBoundaries,"VoxelCentres"->VoxelCentres,"Neighborhoods"->Neighborhoods,"Length"->\[ScriptCapitalL]|>
];

PopulateSimulationBox[simulationbox_,\[ScriptCapitalN]_]:=Module[{\[ScriptCapitalL],\[DoubleStruckX],\[DoubleStruckS]tate,\[DoubleStruckV]},
\[ScriptCapitalL]=simulationbox["Length"];
\[DoubleStruckX]=RandomReal[{-\[ScriptCapitalL]/2,\[ScriptCapitalL]/2},3\[ScriptCapitalN]]; 
\[DoubleStruckS]tate=Join[{0},\[DoubleStruckX]];
\[DoubleStruckV]=RandomPoint[Sphere[{0,0,0}],\[ScriptCapitalN]]//Flatten;
Append[simulationbox,{"Positions"->\[DoubleStruckX],"State"->\[DoubleStruckS]tate,"Momenta"->\[DoubleStruckV]}]
];

LoopBack32[{x_,y_,z_},L_]:=Module[{\[DoubleStruckP],correction},
\[DoubleStruckP]=2Abs[#]/L&/@{x,y,z}//Floor;
correction=-L \[DoubleStruckP] Sign[{x,y,z}]
]; (*implements periodic boundary conditions on a single particle*)

LoopBack[state_,L_]:=Module[{working,parts,corrections},
working=Delete[state,1];
parts=Partition[working,3];
corrections=LoopBack32[#,L]&/@parts;
Join[{0},Flatten[corrections,1]]+state
];(*implements periodic boundary conditions on all particles*)

TorusEulerUpdate[state_,\[DoubleStruckV]_,L_]:=Module[{\[CapitalDelta]t=.25},
LoopBack[state+\[CapitalDelta]t Join[{1},\[DoubleStruckV]],L]
];(*Euler steps the state subject to periodic boundary conditions*)

TorusExtractTrajectories[\[DoubleStruckS]tate_,steps_]:=Module[{history,\[CapitalDelta]t,state,\[DoubleStruckV],\[ScriptCapitalL],\[ScriptCapitalN]},
\[ScriptCapitalL]=\[DoubleStruckS]tate["Length"];
\[DoubleStruckV]=\[DoubleStruckS]tate["Momenta"];
\[ScriptCapitalN]=Partition[\[DoubleStruckV],3]//Length;
state=Join[{0.},\[DoubleStruckS]tate["Positions"]];
history=NestList[TorusEulerUpdate[#,\[DoubleStruckV],\[ScriptCapitalL]]&,state,steps];
\[CapitalDelta]t=history[[{1,2}]]//Differences[#][[1,1]]&;
history=Delete[#,1]&/@history;
history=Partition[#,3]&/@history;
Table[
<|#[[1]]->#[[2]]&/@(Transpose@{Range[\[ScriptCapitalN]],history[[t]]})|>
,{t,1,steps,1}]
];

IdentifyVoxel[{x_,y_,z_},VoxelCentres_,VoxelBoundaries_]:=Block[{cTb,X,VoxelAssoc},
cTb=Transpose@{VoxelCentres,VoxelBoundaries};
X=Cases[cTb,ctb_/;ctb[[2,1,1]]<=x<=ctb[[2,1,2]]&&ctb[[2,2,1]]<=y<=ctb[[2,2,2]]&&ctb[[2,3,1]]<=z<=ctb[[2,3,2]]][[1,1]];
VoxelAssoc=<|{#[[2]]->#[[1]],#[[1]]->#[[2]]}&/@(Transpose@{Range[Length[VoxelCentres]],VoxelCentres})|>;
VoxelAssoc[X]
];

ExtractVoxelDirectory[trajectories_,VoxelCentres_,VoxelBoundaries_]:=Map[IdentifyVoxel[#,VoxelCentres,VoxelCentres]&,trajectories,{2}];

NeighborhoodVoxels[Neighborhoods_,VoxelCentres_,VoxelBoundaries_]:=<|#[[1]]->#[[2]]&/@(Transpose@{Range[Length[Neighborhoods]],Map[IdentifyVoxel[#,VoxelCentres,VoxelBoundaries]&,Neighborhoods,{2}]})|>;

InvertForKeys[x_,y_]:=Select[x,MemberQ[y,#]&]//Keys;

IdentifyOccupiedNeighborhoods[voxeldirectory_,NeighborhoodVoxels_]:=Block[{voxels,min},
min=Table[
voxels=NeighborhoodVoxels[[i]];
Cases[voxeldirectory,n_/;MemberQ[voxels,n]]
,{i,729}];
InvertForKeys[voxeldirectory,#]&/@min
];

IdentifyPairs[OccupiedNeighbourhoods_]:=Subsets[#,{2}]&/@OccupiedNeighbourhoods//Flatten[#,1]&//Map[Sort,#,{1}]&//DeleteDuplicates;

ComputeDistances[pairs_,trajectories_,\[ScriptCapitalL]_]:=TorusDistance[#,\[ScriptCapitalL]]&/@((trajectories/@#)&/@pairs);
ConstructSparseArray[pairs_,distances_,\[ScriptCapitalN]_]:=#[[1]]->#[[2]]&/@Transpose@{pairs,HeavisideTheta[-distances+2]}//SparseArray[#,{\[ScriptCapitalN],\[ScriptCapitalN]},\[Infinity]]&//Normal[#]/.{0->\[Infinity]}&//SparseArray;

SparseMin[ListOfArrays_]:=SparseArray[MapThread[Min,Normal/@ListOfArrays,2]]

ConstructCollisionArray[trajectory_,\[DoubleStruckS]tate_]:=Module[{voxeldirectory,ONtimeline,PairTimeline,DistanceTimeline,SparseArrays,tWeights,Neighborhoods,NVoxels,\[ScriptCapitalL],\[ScriptCapitalN],Vc,Vb},
\[ScriptCapitalL]=\[DoubleStruckS]tate["Length"];\[ScriptCapitalN]=\[DoubleStruckS]tate["Positions"]//Partition[#,3]&//Length;
Print["extracted length and particle number..."];
Vc=\[DoubleStruckS]tate["VoxelCentres"];Vb=\[DoubleStruckS]tate["VoxelBoundaries"];
Print["extracted voxel centres and boundaries..."];
Neighborhoods=\[DoubleStruckS]tate["Neighborhoods"];
Print["extracted all neighborhoods..."];
voxeldirectory=ExtractVoxelDirectory[trajectory,Vb,Vc];
Print["extracted voxel directory..."];
NVoxels=NeighborhoodVoxels[Neighborhoods,Vc,Vb];
Print["extracted all neighborhood\[LeftArrow]voxel relations..."];
ONtimeline=IdentifyOccupiedNeighborhoods[#,NVoxels]&/@voxeldirectory;
Print["identified occupied neighborhoods for all time..."];
PairTimeline=IdentifyPairs[#]&/@ONtimeline;
Print["identified relevant pairs for all time..."];
DistanceTimeline=ComputeDistances[#[[1]],#[[2]],\[ScriptCapitalL]]&/@(Transpose@{PairTimeline,trajectory});
Print["computed pairwise distances..."];
SparseArrays=ConstructSparseArray[#[[1]],#[[2]],\[ScriptCapitalN]]&/@(Transpose@{PairTimeline,DistanceTimeline});
Print["constructed sparse arrays..."];
tWeights=.25Range[Length[SparseArrays]];
SparseMin[tWeights SparseArrays]//Normal[#]/.{\[Infinity]->0.}&
];


End[]
EndPackage[]

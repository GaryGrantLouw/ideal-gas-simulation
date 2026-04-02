(* ::Package:: *)

BeginPackage["AlphaGlue`"]
CreateSimulationBox::usage="Input the linear dimension of the simulation box (# cells, length 2.22). Identifies voxels & neighbourhoods."
FindNeighbors::usage="identifies neighbors of a voxel"
TorusDistance::usage="calculates distance between particles on a torus"
PopulateSimulationBox::usage="populates simulation box with randomly placed particles"
Begin["Private`"]


\[CapitalDelta]l=20/9;
PrototypeNeighborhood=Tuples[{0.,-\[CapitalDelta]l,\[CapitalDelta]l},3];

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

FindNeighbors[voxel_,voxelcentres_,\[CapitalDelta]l_,\[ScriptCapitalL]_]:=Module[{pairs,neighborlist},
pairs={voxel,#}&/@voxelcentres;
neighborlist=Cases[pairs,n_/;TorusDistance[n,\[ScriptCapitalL]]<=Sqrt[3]\[CapitalDelta]l];
neighborlist[[All,2]]
];

CreateSimulationBox[cells_]:=Block[{\[ScriptCapitalL],DisplacementVectors,LinearBoundaries,VoxelBoundaries,VoxelCentres,Neighborhoods},
\[ScriptCapitalL]=\[CapitalDelta]l*cells;
DisplacementVectors=\[ScriptCapitalL] Range[-.5,.5,1/cells]//Tuples[#,3]&;
Neighborhoods=Map[LoopBack3[#,\[ScriptCapitalL]+\[ScriptCapitalL]/cells]&,((PrototypeNeighborhood+ConstantArray[#,27])&/@DisplacementVectors),{2}];
LinearBoundaries=Range[-\[ScriptCapitalL]/2,\[ScriptCapitalL]/2,\[CapitalDelta]l]//Partition[#,2,1]&;
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


End[]
EndPackage[]

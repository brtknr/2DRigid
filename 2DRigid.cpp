/*
 Program Name        2D Rigid Frame Analysis Tool
 
 Program Objective   Analyse a 2D Rigid Frame from the data fed through the input files
 
 Units        This program makes use of the following standard units
 Delta        m
 Force        N
 However, this is not a concrete rule. If the units are consistent, its fine.
 
 Input Files         Members.txt      ~ Contains the information about various members
 - The first line contains the ID of the last member.
 - Data Fields: MemID, NodeID[0], NodeID[1], E, A, I
 - The member of MemID lies between NodeID[0] and NodeID[1] with E, A and I
 
 Nodes.txt        ~ Contains the information about various nodes
 - The first line contains the ID of the last node.
 - Data Fields: NodeID, X co-ordinate, Y co-ordinate
 - The node of NodeID has the coordinate of (X co-ordinate, Y co-ordinate)
 
 Loads.txt        ~ Contains the information about the loads experienced by various nodes
 - Data Fields: NodeID, Dir, Magnitude
 - The last line contains the terminating value of -1
 - The node of NodeID experiences Magnitude of force/moment in the following Dir
 0: X axis      (Rightwards +ve)
 1: Y axis      (Upwards    +ve
 2: Rotational  (CCW        +ve)
 
 Supports.txt     ~ Contains the information about the supports available to various nodes
 - Data Fields: NodeID, Dir
 - The last line contains the terminating value of -1
 - The node of NodeID has linear/rotational support available in the following Dir
 0: Horizontal fixity
 1: Vertical fixity
 2: Rotational fixity
 
 Output Files        Drawing.dxf      ~ Drawing of the structure and its deformed shape as well as bending moment diagram
 Data.txt         ~ Tension experienced by the members in the structure
 
 Date Completed      29 April 2010
 
 Total Time Spent    ~ 70 hours 
 */

// Include System Libraries
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
using namespace std;

// Load constants
#include "Constants.h"

// Define Procedures
void Invert(int, double[][MaxDoF + 1], double[][MaxDoF + 1], double[][MaxDoF + 1]);
int Opposite(int);
void OutputStiffnessMatrix(int, double[][MaxDoF + 1], string);
void DrawStructure(void);
void DXFSetUp(void);
void DXFMember(int);
void DXFBMD(int);
void DXFLoad(int, int);
void DXFMemberNo(int);
void DXFNodeNo(int);
void DXFNodeCircle(int);
void DXFFinishOff(void);

// Define Variables
int Ends[2][MaxMemID + 1],	              // Refers to the 2 Node IDs of a member
Fixed[MaxDoF + 1],	                  // Refers to the boolean for fixity in the DoF
MemID,				                  // Temporary value that refers to the ID of the member
LastMemID,			                  // Value of the last member ID specified in the data file
NodeID,			                      // Temporary value that refers to the ID of the node
LastNodeID,			                  // Value of the last node ID specified in the data file
LastDoF,			                  // The last DoF in the matrices
Dir;			                      // The direction of the action where [0]: X axis, [1]: Y axis, [2]: Rotational

double LengthSQR,		              // Temporary Delta Squared
LL,                                   // LL value for the matrix
LargestLength,                        // Largest dimension in the drawing
LargestDelta,                         // Largest change in the drawing
LargestLoad,                          // Largest load in the drawing
LargestMoment,                        // Largest moment in the drawing

MemberScaleFactor,                    // Scale factor for the members
MomentScaleFactor,                       // Scale factor for the BMD
LoadScaleFactor,                      // Scale factor for the indication of load
NodeRadius,                           // Size of the node radius as % of the size of the largest member
TextSize,                             // Size of the text as a % of the size of the largest member

ThisStiffness,                        // The elementary stiffness
Tension,                              // The member tension

Elongation,                           // The member elongation
Displacement,                         // The relative displacement
YoungsModulus,                        // The young's modulus
Area,                                 // The area
SMoA,                                 // The second moment of area
Coord[2][MaxNodeID + 1],	          // The co-ordinate matrix for the nodes where [0]: X axis, [1]: Y axis
Load[MaxDoF + 1],	                  // The load matrix
Delta[MaxDoF + 1],	                  // The displacement matrix
K[MaxDoF + 1][MaxDoF + 1],	          // The global stiffness matrix
T[MaxDoF + 1][MaxDoF + 1],            // The global transformed matrix
C[MaxDoF + 1][MaxDoF + 1],            // The global check matrix
M[6][6],                              // The elemental stiffness matrix
Forces[6][MaxMemID + 1],              // The forces in the structure
EA[MaxMemID + 1],		              // The EA matrix
EI[MaxMemID + 1],		              // The EI matrix
DeltaLen[2][MaxMemID + 1],	          // The coordinate difference matrix
MemLen[MaxMemID + 1];	              // The member length matrix

stringstream MemIDStr;                    // MemID in text for output

// Initialize Inputs
ifstream MembersInput;
ifstream NodesInput;
ifstream LoadsInput;
ifstream SupportsInput;

// Initialize Outputs
ofstream DrawingOutput;
ofstream DataOutput;
ofstream StiffOutput;

int main(void)
{
    // MAIN PROCEDURE
    // ---------------------------------------------------------------------------------------------------------------------------------
    cout << "2D RIGID JOINTED FRAME ANALYSIS\n\n";
    cout << "Initiating program...\n";
    
    // Check compatiability with the maximum assigned DoF
    if (3 * MaxNodeID + 2 > MaxDoF) {
    	cout << "Problem with DoF sizes...\n";
    	return 0;
    }
    
    char buffer[1000];
    std::cout << "Current directory is: " << getcwd(buffer, 1000) << "\n";    
    
    // Retrieve name of the structure
    string Folder;
    cout << "Name the folder containing structural data: ";   
    cin >> Folder;
    
    string MemberFile   = Folder + "/Members.txt";
    string NodesFile    = Folder + "/Nodes.txt";
    string LoadsFile    = Folder + "/Loads.txt";
    string SupportsFile = Folder + "/Supports.txt";
    
    string DrawingFile  = Folder + "/Drawing.dxf";
    string DataFile     = Folder + "/Data.csv";
    string StiffFile    = Folder + "/Stiffness.txt";
    
    // Open files
    MembersInput.open(MemberFile.data());
    NodesInput.open(NodesFile.data());
    LoadsInput.open(LoadsFile.data());
    SupportsInput.open(SupportsFile.data());
    DrawingOutput.open(DrawingFile.data());
    DataOutput.open(DataFile.data());
    StiffOutput.open(StiffFile.data());
    
    // Check if the files exist
    if (!MembersInput.is_open() || !NodesInput.is_open() || !LoadsInput.is_open() || !SupportsInput.is_open()){
    	cout << "All or part of the structure data is missing.\n";
        sleep(3);
    	return 0;     
    }
    
    // MEMBERS
    // -----------------------------------------------------------------------------------------------------------------------------
    cout << "Loading Members.txt file...\n";
    
    // Read the last member ID on the file from the first line and check if it is within allowed limit
    MembersInput >> LastMemID;
    if (LastMemID > MaxMemID) {
    	cout << "There are too many members.\n";
        sleep(3);
    	return 0;
    }
    // Read the member attributes: Member ID; Node 0, Node 1, Member Stiffness
    for (int i = 0; i <= LastMemID; i++) {
    	MembersInput >> MemID;
    	MembersInput >> Ends[0][MemID] >> Ends[1][MemID];
        MembersInput >> YoungsModulus >> Area >> SMoA;
        EA[MemID] = YoungsModulus * Area;
        EI[MemID] = YoungsModulus * SMoA;
    }
    // Close the member data file
    MembersInput.close();
    
    // NODES
    // -----------------------------------------------------------------------------------------------------------------------------
    cout << "Loading Node.txt file...\n";
    
    // Read the last node ID on the file from the first line and check if it is within allowed limit  
    NodesInput >> LastNodeID;
    if (LastNodeID > MaxNodeID) {
    	cout << "There are too many nodes.\n";
  	    system("pause");
    	return 0;
    }
    // Read the node attributes: Node ID; X co-ordinate, Y co-ordinate
    for (int i = 0; i <= LastNodeID; i++) {
    	NodesInput >> NodeID;
    	NodesInput >> Coord[0][NodeID] >> Coord[1][NodeID];
    }
    // Close the node data file
    NodesInput.close();
    
    // DEGREE OF FREEDOM
    // -----------------------------------------------------------------------------------------------------------------------------    
    cout << "Determining degree of freedom...\n";
    
    // Calculate the sum total degree of freedom
    LastDoF = 3 * LastNodeID + 2;
    // Zero the DoF dependent matrices
    for (int ThisDoF = 0; ThisDoF <= LastDoF; ThisDoF++) {
    	Load[ThisDoF] = 0.0;
    	Delta[ThisDoF] = 0.0;
    	Fixed[ThisDoF] = 0;
    	for (int ThatDoF = 0; ThatDoF <= LastDoF; ThatDoF++) {
    	    K[ThisDoF][ThatDoF] = 0.0;
    	    C[ThisDoF][ThatDoF] = 0.0;
    	}
    }
    
    // LOADS
    // -----------------------------------------------------------------------------------------------------------------------------
    cout << "Loading Loads.txt file...\n";
    
    // Read the load attributes: Node ID; Dir; Magnitude
    for (;;) {
    	LoadsInput >> NodeID;  
        // Break the loop if the node ID is negative
    	if (NodeID < 0) break;
    	LoadsInput >> Dir;
    	LoadsInput >> Load[3 * NodeID + Dir];
    }
    
    // SUPPORTS
    // -----------------------------------------------------------------------------------------------------------------------------
    cout << "Loading Supports.txt file...\n";
    
    // Read the support attributes: Node ID; Dir
    for (;;) {
    	SupportsInput >> NodeID;
        // Break the loop if the node ID is negative
    	if (NodeID < 0) break;  
    	SupportsInput >> Dir;
        // Mark the particular DoF
    	Fixed[3 * NodeID + Dir] = 1;
        // Zero the load if the node is fixed
    	Load[3 * NodeID + Dir] = 0.0;	
    }
        
    // GLOBAL STIFFNESS MATRIX FOR A 2D RIGID FRAME
    // -----------------------------------------------------------------------------------------------------------------------------
    cout << "Building the global stiffness matrix...\n";
    
    // Load a member at a time
    for (int MemID = 0; MemID <= LastMemID; MemID++) {
        // Zero the length squared variable
    	LengthSQR = 0.0;
    	// Load the same axis of either end of the member
    	for (int ThisAxis = 0; ThisAxis <= 1; ThisAxis++) {
            // Work out the difference in length in that axis
    	    DeltaLen[ThisAxis][MemID] = Coord[ThisAxis][Ends[1][MemID]] - Coord[ThisAxis][Ends[0][MemID]];
            // Add the square of the difference to the total sum
    	    LengthSQR += pow(DeltaLen[ThisAxis][MemID], 2);
        }
    	// Determine the member length using Pythagoras
    	MemLen[MemID] = sqrt(LengthSQR);
    	
    	// Determine the member with the largest length
        if (MemLen[MemID] > LargestLength) {
            LargestLength = MemLen[MemID]; 
        }
        
    	// EA/L
    	for (int ThisAxis = 0; ThisAxis <= 1; ThisAxis++) {
    	    for (int ThatAxis = 0; ThatAxis <= 1; ThatAxis++) {
                LL = (DeltaLen[ThisAxis][MemID] * DeltaLen[ThatAxis][MemID]) / pow(MemLen[MemID], 2);
        		ThisStiffness = (EA[MemID] * LL) / MemLen[MemID];
        		for (int ThisEnd = 0; ThisEnd <= 1; ThisEnd++) {
        		    for (int ThatEnd = 0; ThatEnd <= 1; ThatEnd++) {
            			int ThisDoF = 3 * Ends[ThisEnd][MemID] + ThisAxis;
            			int ThatDoF = 3 * Ends[ThatEnd][MemID] + ThatAxis;
            			if (ThisEnd == ThatEnd)
            			    K[ThisDoF][ThatDoF] += ThisStiffness;
            			else
            			    K[ThisDoF][ThatDoF] -= ThisStiffness;
        		    }
        		}
    	    }
    	}
        
        // 12EI/L^3
    	for (int ThisAxis = 0; ThisAxis <= 1; ThisAxis++) {
    	    for (int ThatAxis = 0; ThatAxis <= 1; ThatAxis++) {
                LL = (DeltaLen[Opposite(ThisAxis)][MemID] * DeltaLen[Opposite(ThatAxis)][MemID]) / pow(MemLen[MemID], 2);
        		ThisStiffness = (12.0 * EI[MemID] * LL) / pow(MemLen[MemID], 3);
        		for (int ThisEnd = 0; ThisEnd <= 1; ThisEnd++) {
        		    for (int ThatEnd = 0; ThatEnd <= 1; ThatEnd++) {
            			int ThisDoF = 3 * Ends[ThisEnd][MemID] + ThisAxis;
            			int ThatDoF = 3 * Ends[ThatEnd][MemID] + ThatAxis;
            			if (ThisEnd == ThatEnd) {
            			    K[ThisDoF][ThatDoF] += ThisStiffness * pow(-1.0, ThisAxis + ThatAxis);
            			} else {
            			    K[ThisDoF][ThatDoF] -= ThisStiffness * pow(-1.0, ThisAxis + ThatAxis);
                        }
        		    }
        		}
            }
    	}
        
        // 6EI/L^2
    	for (int ThisAxis = 0; ThisAxis <= 1; ThisAxis++) {
            LL = DeltaLen[Opposite(ThisAxis)][MemID] / MemLen[MemID];
    		ThisStiffness = 6.0 * EI[MemID] * LL / pow(MemLen[MemID], 2);
    		for (int ThisEnd = 0; ThisEnd <= 1; ThisEnd++) {
    		    for (int ThatEnd = 0; ThatEnd <= 1; ThatEnd++) {
        			int ThisDoF = 3 * Ends[ThisEnd][MemID];
        			int ThatDoF = 3 * Ends[ThatEnd][MemID];
                    if (ThisEnd == ThatEnd){
        			    K[ThisDoF + 2][ThatDoF + ThisAxis] -= ThisStiffness * pow(-1.0, ThisEnd + ThisAxis) ;
        			    K[ThisDoF + ThisAxis][ThatDoF + 2] -= ThisStiffness * pow(-1.0, ThisEnd + ThisAxis);
        			} else {
        			    K[ThisDoF + 2][ThatDoF + ThisAxis] += ThisStiffness * pow(-1.0, ThisEnd + ThisAxis);
        			    K[ThisDoF + ThisAxis][ThatDoF + 2] -= ThisStiffness * pow(-1.0, ThisEnd + ThisAxis);
                    }
    		    }
    		}
    	}
        
        // 2EI/L
	    ThisStiffness = 2.0 * EI[MemID] / MemLen[MemID];
		for (int ThisEnd = 0; ThisEnd <= 1; ThisEnd++) {
	        for (int ThatEnd = 0; ThatEnd <= 1; ThatEnd++) {
    			int ThisDoF = 3 * Ends[ThisEnd][MemID] + 2;
    			int ThatDoF = 3 * Ends[ThatEnd][MemID] + 2;
    			if (ThisEnd == ThatEnd)
    			    K[ThisDoF][ThatDoF] += ThisStiffness * 2.0;
    			else
    			    K[ThisDoF][ThatDoF] += ThisStiffness;
		    }
		}
		MemIDStr.str("");
        MemIDStr << MemID;
	    OutputStiffnessMatrix(LastDoF,K,"GSM after accounting Member " + MemIDStr.str());
        
    }
    
    // ACCOUNT FOR SUPPORTS
    // -----------------------------------------------------------------------------------------------------------------------------
    cout << "Accounting for supports...\n";
    OutputStiffnessMatrix(LastDoF, K,"GSM before accounting for Supports");
    
    for (int ThisDoF = 0; ThisDoF <= LastDoF; ThisDoF++) {
    	if (Fixed[ThisDoF] == 1) {
    	    for (int ThatDoF = 0; ThatDoF <= LastDoF; ThatDoF++) {
    		    K[ThatDoF][ThisDoF] = 0.0;
    		    K[ThisDoF][ThatDoF] = 0.0;
    	    }
    	    K[ThisDoF][ThisDoF] = 1.0;
    	}
    }
    
    OutputStiffnessMatrix(LastDoF, K, "GSM after accounting for Supports");
    
    // INVERT MATRIX
    // -----------------------------------------------------------------------------------------------------------------------------
    cout << "Inverting matrix...\n";
    
    Invert(LastDoF, K, T, C);
    
    OutputStiffnessMatrix(LastDoF, K, "Inverted Matrix");
    
    OutputStiffnessMatrix(LastDoF, C, "Check Matrix");
    
    // DISPLACEMENTS AND ROTATIONS
    // -----------------------------------------------------------------------------------------------------------------------------
    cout << "Working out displacements and rotations...\n";
    
    for (int ThisDoF = 0; ThisDoF <= LastDoF; ThisDoF++) {
    	for (int ThatDoF = 0; ThatDoF <= LastDoF; ThatDoF++) {
    	    Delta[ThisDoF] += T[ThisDoF][ThatDoF] * Load[ThatDoF];
    	}
        // Determine if it is the largest delta
        if (abs(Delta[ThisDoF]) > LargestDelta) {
            LargestDelta = abs(Delta[ThisDoF]); 
        }
        // Determine the largest load
        if (abs(Load[ThisDoF]) > LargestLoad) {
            LargestLoad = abs(Load[ThisDoF]);
        }
    }    
    
    // FORCES AND MOMENTS
    // -----------------------------------------------------------------------------------------------------------------------------
    cout << "Working out shear force and rotation...\n";
    
    // Load a member at a time
    for (int MemID = 0; MemID <= LastMemID; MemID++) {
        
        for (int ThisDoF = 0; ThisDoF <= 5; ThisDoF++) for (int ThatDoF = 0; ThatDoF <= 5; ThatDoF++) M[ThisDoF][ThatDoF] = 0.0;
        
    	// EA/L
    	for (int ThisAxis = 0; ThisAxis <= 1; ThisAxis++) {
    	    for (int ThatAxis = 0; ThatAxis <= 1; ThatAxis++) {
                LL = (DeltaLen[ThisAxis][MemID] * DeltaLen[ThatAxis][MemID]) / pow(MemLen[MemID], 2);
        		ThisStiffness = (EA[MemID] * LL) / MemLen[MemID];
        		for (int ThisEnd = 0; ThisEnd <= 1; ThisEnd++) {
        		    for (int ThatEnd = 0; ThatEnd <= 1; ThatEnd++) {
            			int ThisDoF = 3 * ThisEnd + ThisAxis;
            			int ThatDoF = 3 * ThatEnd + ThatAxis;
            			if (ThisEnd == ThatEnd)
            			    M[ThisDoF][ThatDoF] += ThisStiffness;
            			else
            			    M[ThisDoF][ThatDoF] -= ThisStiffness;
        		    }
        		}
    	    }
    	}
        
        // 12EI/L^3
    	for (int ThisAxis = 0; ThisAxis <= 1; ThisAxis++) {
    	    for (int ThatAxis = 0; ThatAxis <= 1; ThatAxis++) {
                LL = (DeltaLen[Opposite(ThisAxis)][MemID] * DeltaLen[Opposite(ThatAxis)][MemID]) / pow(MemLen[MemID], 2);
        		ThisStiffness = (12.0 * EI[MemID] * LL) / pow(MemLen[MemID], 3);
        		for (int ThisEnd = 0; ThisEnd <= 1; ThisEnd++) {
        		    for (int ThatEnd = 0; ThatEnd <= 1; ThatEnd++) {
            			int ThisDoF = 3 * ThisEnd + ThisAxis;
            			int ThatDoF = 3 * ThatEnd + ThatAxis;
            			if (ThisEnd == ThatEnd) {
            			    M[ThisDoF][ThatDoF] += ThisStiffness * pow(-1.0, ThisAxis + ThatAxis);
            			} else {
            			    M[ThisDoF][ThatDoF] -= ThisStiffness * pow(-1.0, ThisAxis + ThatAxis);
                        }
        		    }
        		}
            }
    	}
        
        // 6EI/L^2
    	for (int ThisAxis = 0; ThisAxis <= 1; ThisAxis++) {
            LL = DeltaLen[Opposite(ThisAxis)][MemID] / MemLen[MemID];
    		ThisStiffness = 6.0 * EI[MemID] * LL / pow(MemLen[MemID], 2);
    		for (int ThisEnd = 0; ThisEnd <= 1; ThisEnd++) {
    		    for (int ThatEnd = 0; ThatEnd <= 1; ThatEnd++) {
        			int ThisDoF = 3 * ThisEnd;
        			int ThatDoF = 3 * ThatEnd;
                    if (ThisEnd == ThatEnd){
        			    M[ThisDoF + 2][ThatDoF + ThisAxis] -= ThisStiffness * pow(-1.0, ThisEnd + ThisAxis) ;
        			    M[ThisDoF + ThisAxis][ThatDoF + 2] -= ThisStiffness * pow(-1.0, ThisEnd + ThisAxis);
        			} else {
        			    M[ThisDoF + 2][ThatDoF + ThisAxis] += ThisStiffness * pow(-1.0, ThisEnd + ThisAxis);
        			    M[ThisDoF + ThisAxis][ThatDoF + 2] -= ThisStiffness * pow(-1.0, ThisEnd + ThisAxis);
                    }
    		    }
    		}
    	}
        
        // 2EI/L
	    ThisStiffness = 2.0 * EI[MemID] / MemLen[MemID];
		for (int ThisEnd = 0; ThisEnd <= 1; ThisEnd++) {
	        for (int ThatEnd = 0; ThatEnd <= 1; ThatEnd++) {
    			int ThisDoF = 3 * ThisEnd + 2;
    			int ThatDoF = 3 * ThatEnd + 2;
    			if (ThisEnd == ThatEnd)
    			    M[ThisDoF][ThatDoF] += ThisStiffness * 2.0;
    			else
    			    M[ThisDoF][ThatDoF] += ThisStiffness;
		    }
		}
        
        for (int ThisEnd = 0; ThisEnd <= 1; ThisEnd++) {
            for (int ThisAxis = 0; ThisAxis <= 2; ThisAxis++) {
//              int ThisDoF = 3 * Ends[ThisEnd][MemID] + ThisAxis;
                for (int ThatEnd = 0; ThatEnd <= 1; ThatEnd++) {                             
           	        for (int ThatAxis = 0; ThatAxis <= 2; ThatAxis++) {
                        int ThatDoF = 3 * Ends[ThatEnd][MemID] + ThatAxis; 
                        Forces[ThisEnd * 3 + ThisAxis][MemID] += M[3 * ThisEnd + ThisAxis][3 * ThatEnd + ThatAxis] * Delta[ThatDoF];
                    }
                }
            }
   	    }		
    }
    
    
    // DATA OUTPUT
    // -----------------------------------------------------------------------------------------------------------------------------
    cout << "Writing Data.txt file...\n";
    
	DataOutput.precision(3);
    DataOutput.setf(ios::fixed, ios::floatfield);
   	DataOutput << "Member" 
    << ",Tension" 
    << ",Shear Force" 
    << ",Displacement" 
    << ",Elongation"
    << ",Moment"
    << ",Rotation"
    << ",@"
    << ",Moment"
    << ",Rotation"
    << ",@\n";    
   	DataOutput << "ID" 
    << ",N" 
    << ",N" 
    << ",mm" 
    << ",mm"
    << ",Nmm"
    << ",rad"
    << ",Node"
    << ",Nmm"
    << ",rad"
    << ",Node\n";                             
    
    for (int MemID = 0; MemID <= LastMemID; MemID++) {
        int ThisDoF      = 3 * Ends[0][MemID];
        int ThatDoF      = 3 * Ends[1][MemID];
        Displacement     = (DeltaLen[1][MemID] * (Delta[ThatDoF + 0] - Delta[ThisDoF + 0]) 
                            +  DeltaLen[0][MemID] * (Delta[ThatDoF + 1] - Delta[ThisDoF + 1]))
        /  MemLen[MemID];
        Elongation       = (DeltaLen[0][MemID] * (Delta[ThatDoF + 0] - Delta[ThisDoF + 0])
                            +  DeltaLen[1][MemID] * (Delta[ThatDoF + 1] - Delta[ThisDoF + 1]))
        /  MemLen[MemID];
        Tension          = Elongation * (EA[MemID] / MemLen[MemID]);    
    	DataOutput << MemID
        << "," << Tension 
        << "," << (Forces[2][MemID] + Forces[5][MemID]) / MemLen[MemID] 
        << "," << Displacement 
        << "," << Elongation 
        << "," << Forces[2][MemID] 
        << "," << Delta[ThisDoF + 2]
        << "," << Ends[0][MemID]
        << "," << Forces[5][MemID] 
        << "," << Delta[ThatDoF + 2] 
        << "," << Ends[1][MemID] <<"\n";
        
        // Determine the largest load
        if (abs(Forces[2][MemID]) > LargestMoment) {
            LargestMoment = abs(Forces[2][MemID]);
        }        
        if (abs(Forces[5][MemID]) > LargestMoment) {
            LargestMoment = abs(Forces[5][MemID]);
        }        
        
    }
    // Close the output files
    DataOutput.close();
    StiffOutput.close();

    // DRAWING OUTPUT
    // -----------------------------------------------------------------------------------------------------------------------------
    
    MemberScaleFactor = MemberScale * LargestLength/LargestDelta;
    LoadScaleFactor = LoadScale * LargestLength/LargestLoad;
    MomentScaleFactor = MomentScale * LargestLength/LargestMoment;    
    TextSize = TextScale * LargestLength;
    NodeRadius = NodeScale * LargestLength;

    cout << "LargestLoad: "<< LargestDelta << "\n";
    cout << "LargestLoad: "<< LargestLoad << "\n";
    
    cout << "Writing Drawing.dxf file...\n";      
    
    DrawStructure();
    
    // END OF MAIN PROCEDURE
    // -----------------------------------------------------------------------------------------------------------------------------
    cout << "End of program.\n";
    return 0;
}


// SUB PROCEDURES
// -----------------------------------------------------------------------------------------------------------------------------

int Opposite(int Value)
{
    if (Value == 0) return 1;
    else            return 0;
}

void OutputStiffnessMatrix(int Max, double MDArray[][MaxDoF + 1], string Title)
{
    if (LastDoF < MaxDoFOutput){
        StiffOutput.precision(2);
        StiffOutput.setf(ios::fixed, ios::floatfield);
        StiffOutput << Title << "\n";
        for (int i = 0; i <= Max; i++) {
            StiffOutput << i;
            for (int j = 0; j <= Max; j++) {
                StiffOutput << "\t";
                StiffOutput << MDArray[i][j];
            }
            StiffOutput << "\n";
        }
      	StiffOutput << "\n";
    }
}

void DXFSetUp(void)
{
    DrawingOutput << "0\n" << "SECTION\n" << "2\n" << "ENTITIES\n";
}

void DrawStructure(void)
{
    DXFSetUp();
    for (int MemID = 0; MemID <= LastMemID; MemID += 1) {
        DXFBMD(MemID);   
    }
    for (int MemID = 0; MemID <= LastMemID; MemID += 1) {
        DXFMember(MemID);
    }
    for (int NodeID = 0; NodeID <= LastNodeID; NodeID += 1) {
        DXFNodeCircle(NodeID);
    }
    for (int MemID = 0; MemID <= LastMemID; MemID += 1) {
        DXFMemberNo(MemID);
    }
    for (int NodeID = 0; NodeID <= LastNodeID; NodeID += 1) {
        DXFNodeNo(NodeID);
    }
    for (int NodeID = 0; NodeID <= LastNodeID; NodeID += 1) {
        for (int Dir = 0; Dir <= 2; Dir++) {
//          cout << NodeID << ": " << Dir << ": " << Load[3 * NodeID + Dir] << "\n";          
            if (Load[3 * NodeID + Dir] != 0.0) {
                DXFLoad(NodeID, Dir);
            }
        }
    }
    DXFFinishOff();    
}

void DXFMember(int ThisMem)
{
    DrawingOutput << "0\nLINE\n8\nUndeflected\n" 
    << "10\n" << Coord[0][Ends[0][ThisMem]] << "\n" 
    << "20\n" << Coord[1][Ends[0][ThisMem]] << "\n" 
    << "11\n" << Coord[0][Ends[1][ThisMem]] << "\n" 
    << "21\n" << Coord[1][Ends[1][ThisMem]] << "\n" 
    << "62\n0\n";
    DrawingOutput << "0\nLINE\n8\nDeflected\n" 
    << "10\n" << Coord[0][Ends[0][ThisMem]] + MemberScaleFactor * Delta[3 * Ends[0][ThisMem] + 0] << "\n" 
    << "20\n" << Coord[1][Ends[0][ThisMem]] + MemberScaleFactor * Delta[3 * Ends[0][ThisMem] + 1] << "\n" 
    << "11\n" << Coord[0][Ends[1][ThisMem]] + MemberScaleFactor * Delta[3 * Ends[1][ThisMem] + 0] << "\n" 
    << "21\n" << Coord[1][Ends[1][ThisMem]] + MemberScaleFactor * Delta[3 * Ends[1][ThisMem] + 1] << "\n" 
    << "62\n1\n";
}

void DXFNodeCircle(int ThisNode)
{
    DrawingOutput << "0\nCIRCLE\n8\nUndeflected\n" 
    << "10\n" << Coord[0][ThisNode] << "\n" 
    << "20\n" << Coord[1][ThisNode] << "\n" 
    << "40\n" << NodeRadius << "\n" 
    << "62\n0\n";
    DrawingOutput << "0\nCIRCLE\n8\nDeflected\n" 
    << "10\n" << Coord[0][ThisNode] + MemberScaleFactor * Delta[3 * ThisNode + 0]<< "\n" 
    << "20\n" << Coord[1][ThisNode] + MemberScaleFactor * Delta[3 * ThisNode + 1]<< "\n" 
    << "40\n" << NodeRadius  << "\n" 
    << "62\n1\n";
}

void DXFBMD(int ThisMem)
{
    double Angle     = atan(-DeltaLen[0][ThisMem]/DeltaLen[1][ThisMem]);
    int Reverse      = 1;
    if (Coord[0][Ends[1][ThisMem]] < Coord[0][Ends[0][ThisMem]]) Reverse = -1;
    if (Coord[1][Ends[1][ThisMem]] < Coord[1][Ends[0][ThisMem]]) Reverse = -1;
    DrawingOutput << "0\nLINE\n8\nBMD\n" 
    << "10\n" << Coord[0][Ends[0][ThisMem]] << "\n" 
    << "20\n" << Coord[1][Ends[0][ThisMem]] << "\n" 
    << "11\n" << Coord[0][Ends[0][ThisMem]] - cos(Angle) * Forces[2][ThisMem] * Reverse * MomentScaleFactor << "\n" 
    << "21\n" << Coord[1][Ends[0][ThisMem]] - sin(Angle) * Forces[2][ThisMem] * Reverse * MomentScaleFactor << "\n" 
    << "62\n4\n";
    DrawingOutput << "0\nLINE\n8\nBMD\n" 
    << "10\n" << Coord[0][Ends[0][ThisMem]] - cos(Angle) * Forces[2][ThisMem] * Reverse * MomentScaleFactor << "\n" 
    << "20\n" << Coord[1][Ends[0][ThisMem]] - sin(Angle) * Forces[2][ThisMem] * Reverse * MomentScaleFactor << "\n" 
    << "11\n" << Coord[0][Ends[1][ThisMem]] + cos(Angle) * Forces[5][ThisMem] * Reverse * MomentScaleFactor << "\n" 
    << "21\n" << Coord[1][Ends[1][ThisMem]] + sin(Angle) * Forces[5][ThisMem] * Reverse * MomentScaleFactor << "\n" 
    << "62\n4\n";
    DrawingOutput << "0\nLINE\n8\nBMD\n" 
    << "10\n" << Coord[0][Ends[1][ThisMem]] + cos(Angle) * Forces[5][ThisMem] * Reverse * MomentScaleFactor << "\n" 
    << "20\n" << Coord[1][Ends[1][ThisMem]] + sin(Angle) * Forces[5][ThisMem] * Reverse * MomentScaleFactor << "\n" 
    << "11\n" << Coord[0][Ends[1][ThisMem]] << "\n" 
    << "21\n" << Coord[1][Ends[1][ThisMem]] << "\n" 
    << "62\n4\n";
}

void DXFLoad(int ThisNode, int ThisDirection)
{
    if (ThisDirection == 2){
        int StartAngle = 340;
        int EndAngle = 120;
        double Radius = Load[3 * ThisNode + ThisDirection] * LoadScaleFactor;
        if (Radius < 0) {
            Radius *= -1;        
            StartAngle = 240;
            EndAngle = 20;
        }
        DrawingOutput << "0\nARC\n8\nInformation\n" 
        << "10\n" << Coord[0][ThisNode] << "\n" 
        << "20\n" << Coord[1][ThisNode] << "\n"
        << "40\n" << Radius << "\n" 
        << "50\n" << StartAngle << "\n" 
        << "51\n" << EndAngle << "\n"                    
        << "62\n3\n";
    } else { 
        DrawingOutput << "0\nLINE\n8\nInformation\n" 
        << "10\n" << Coord[0][ThisNode] << "\n" 
        << "20\n" << Coord[1][ThisNode] << "\n";
        if (ThisDirection == 0)
            DrawingOutput << "11\n"
            << Coord[0][ThisNode] - Load[3 * ThisNode + ThisDirection] * LoadScaleFactor << "\n"    
            << "21\n" << Coord[1][ThisNode] << "\n";
        else if (ThisDirection == 1)
            DrawingOutput << "11\n"
            << Coord[0][ThisNode] << "\n"
            << "21\n" << Coord[1][ThisNode] - Load[3 * ThisNode + ThisDirection] * LoadScaleFactor << "\n";      
        DrawingOutput << "62\n1\n";
    }
}

void DXFMemberNo(int ThisMem)
{
    DrawingOutput << "0\nTEXT\n8\nInformation\n" 
    << "1\n" << ThisMem << "\n"
    << "10\n" << (Coord[0][Ends[0][ThisMem]] + Coord[0][Ends[1][ThisMem]]) / 2.0 << "\n" 
    << "20\n" << (Coord[1][Ends[0][ThisMem]] + Coord[1][Ends[1][ThisMem]]) / 2.0 << "\n" 
    << "40\n" << TextSize << "\n" 
    << "62\n0\n";
}

void DXFNodeNo(int ThisNode)
{
    DrawingOutput << "0\nTEXT\n8\nInformation\n" 
    << "1\n" << ThisNode << "\n"
    << "10\n" << Coord[0][ThisNode] + TextSize << "\n"
    << "20\n" << Coord[1][ThisNode] + TextSize << "\n"
    << "40\n" << TextSize << "\n" 
    << "62\n5\n";
}

void DXFFinishOff(void)
{
    DrawingOutput << "0\n" << "ENDSEC\n" << "0\n" << "EOF\n";
    DrawingOutput.close();
}

// END OF SUB PROCEDURES
// -----------------------------------------------------------------------------------------------------------------------------
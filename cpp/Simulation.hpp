#ifndef SIMULATIONDEF
#define SIMULATIONDEF

#define NUMTYPE_IS_FLOAT true
#define PERIODIC_BOUNDARY_CONDITIONS true

#include <iostream> 		// Input and output
#include <iomanip>			// Input and output formatting
#include <fstream>			// File streams

#include <armadillo> 		// Matrix library

#include <random> 			// Random numbers
#include <curand.h>
#include <curand_kernel.h>
#include <math.h>		 	// Mathmatical constants
#include <algorithm>

#include <vector>			// C++ standard vector
#include <string.h> 		// Strings

#include <cassert> 			// Assertions

#include <sys/types.h> 		// Packages for the directory
#include <sys/stat.h>		//    information,
#include <dirent.h>		 	//    handling and etc.
#include <ctime>			// Time functions

#if NUMTYPE_IS_FLOAT
typedef float numtype;
typedef arma::fmat mattype;
typedef arma::fcolvec colvectype;
#else
typedef double numtype;
typedef arma::mat mattype;
typedef arma::colvec colvectype;
#endif

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>

/* Class to contain simulation parameters, and which drives the simulation forward */
class Chains
{
private:
	numtype P_0; 				// [1/µm^3] The density of invading phages in the simulation initially

	int N_max; 					//          The maximum number of cells in the simulation
	int M_max; 					//          The maximum number of phages in the simulation

	numtype dT;					// [hour]   Size of the time step

	int nSamp;	 				//          Number of samples to save per simulation hour

	numtype L; 			    	// [µm]     Length of boundary condition box (x direction)

	numtype margin; 			// 			Number of target typical lengthscales to simulate
	numtype T_i; 				// [hour]   Time when the phage infections begins (-1 disables phage infection)

	numtype g; 					// [1/hour] Growth rate for the cells

	numtype Ld; 				// [µm]     The length scale for division (Typical volume (0.6 - 0.7) µm^3)
	numtype R;					// [µm]     The "radius" of the cells

	numtype k_rep; 				// [N*m]    Parameter for repulsive potential
	numtype k_att; 				// [N*m]    Parameter for attraction potential
	numtype k_int; 				// [N*m]    Parameter for internal spring potential
	numtype k_pull; 			// [N]      Parameter for colony formation (colony gravity)

	numtype gamma; 				//          Probability to infect cell (also contains adsorption)
	int beta;			 		//          Multiplication factor phage
	numtype delta; 				// [1/hour] Rate of phage decay

	numtype r; 					// [1/hour] Rate of lysis

	numtype eta; 				//          Amount of division noise (0.5+eta) * V and (0.5-eta) * V
	numtype nu;					//          Amount of displacement noise of the new poles (width of gaussian)

	numtype bendingAngle; 		// [rad]		The allowed angle between two cells

	numtype D_B; 				// [µm^2/hour] Diffusion constant for the cells
	numtype D_P; 				// [µm^2/hour] Diffusion constant for the phage

	numtype Time; 				// Integer to keep track of time
	numtype RunTime; 			// Variable tracking total run time

	bool lockCells; 			// Boolean to lock the cells (stop updating)
	bool allCaptured; 			// Boolean to stop phage updating (all are captured)

	int debug;		 			// Set the debug level
	bool exit;		 			// Exit bool; terminates execution
	bool firstRun; 				// Bool to indicate if this run is the first

	bool wellMixed; 			// Bool to indicate if bacteria should be well mixed

	bool exportAny;				//
	bool exportCellData;	 	// Booleans to control the export output
	bool exportColonySize; 		//
	bool exportPhageData;		//

	bool ready; 				// Boolean to indicate whether the data is ready on the GPU

	numtype rngSeed;			// The seed for the random number generator
	std::mt19937 rng; 			// Mersenne twister, random number generator
	std::uniform_real_distribution<numtype> rand;
	std::normal_distribution<numtype> randn;

	std::ofstream f_cells;		// Filestream to save configuration of cells
	std::ofstream f_colonySize; // Filestream to save size of colony (Volume, N, N_Lys, N_Lyt)
	std::ofstream f_phages;		// Filestream to save configuration of phages
	std::ofstream f_log;		// Filestream to save log.txt

	std::string path; 			// Sets the path to store in

	// CUDA related variables
	thrust::host_vector<numtype>   h_cells;
	thrust::device_vector<numtype> d_cells;
	thrust::device_vector<numtype> d_cells_new;

	thrust::host_vector<numtype>   h_phages;
	thrust::device_vector<numtype> d_phages;

	thrust::host_vector<int>   h_active;
	thrust::device_vector<int> d_active;

	thrust::host_vector<curandState>   h_rng_state;
	thrust::device_vector<curandState> d_rng_state;

	int phagesBlockSize;
	int phagesGridSize;

	int cellsBlockSize;
	int cellsGridSize;

public:
	// Constructers
	explicit Chains(int N_max);										// Direct constructer
	explicit Chains(Chains &other); 							    // Copy constructor

	// Drivers
	int Run(numtype T); 											// Controls the evaluation of the simulation

private:
	void Initialize(); 												// Initialize the simulation

public:
	void SpawnPhages(); 						            		// Spawns phages according to spawning rules
	void LyseCell(int I);											// Replace cell I with beta phages

	void AutoScale();												// Auto scale the simulation space
	void Relax();													// Auto relaxes the bacteria
	void Equilibrate(numtype T);									// Equilibrates the phage distribution

public:
	// Settings
	void TimeStep(numtype dT); 										// Set the size of the time-step
	void SetLength(numtype L); 										// Set the length of the simulation space
	void SetMargin(numtype margin); 								// Set the margin (number of lengthscales to simulate)

	void PhageInvasionStartTime(numtype T_i);						// Sets the time when the phages should start infecting
	void PhageInitialDensity(numtype P_0);							// Sets initial density of the phages (1/µm^3)
	void PhageDiffusionConstant(numtype D_P);						// Sets the diffusion constant of the phages
	void PhageInfectionRate(numtype r);								// Sets rate of the infection increaasing in stage
	void PhageDecayRate(numtype delta);								// Set the decay rate of the phages
	void PhageBurstSize(int beta);									// Set the size of the bursts
	void PhageAdsorptionParameter(numtype gamma); 					// Changes the adsorption parameter gamma

	void CellDiffusionConstant(numtype D_B);		 				// Sets the diffusion constant of the bacteria
	void CellRepulsiveParameter(numtype k_rep);						// Sets the strength of the repulsive potential
	void CellAttractiveParameter(numtype k_att); 					// Sets the strength of the attractive potential
	void CellInternalParameter(numtype k_int);	 					// Sets the strength of the internal potential
	void CellLength(numtype Ld); 									// Sets the (division)length of the bacteria
	void CellRadius(numtype R);										// Sets the radius of the bactera
	void CellBendingAngle(numtype bendingAngle); 					// Sets the bending angle between cells
	void CellLock();												// Locks the cells in their current configuration
	void WellMixed();												// Sets the bacteria to be well mixed
	void SphericalColony(numtype k_pull);							// Sets the bacteria to form a spherical colony (with force k_pull)

	void MaxPhageCount(int M_max); 									// Sets the maximum number of phages in the simulation

	void SetRngSeed(int n); 										// Sets the seed of the random number generator

private:
	void deb(const std::string &input, int n); 						// Debug function (prints "input")
	void error(const std::string &input);			 				// Error function (prints "input" and set exit = true)
	void warning(const std::string &input);			 				// Warning function (prints "input")
	void WriteLog();											 	// Write a log.txt file

public:
	void Quiet();			 										// Set debug level to 0
	void Debug(int n); 												// Set the debug level to n

	void SetSamples(int nSamp); 									// Set the number of output samples

	// File outputs
	void ExportCellData();	 										//
	void ExportColonySize(); 										// Sets booleans for export functions
	void ExportPhageData();											//

	void ExportCellDataNow();
	void ExportPhageDataNow();

private:
	void ExportData(numtype t);										// Master function to export the data
	void f_ExportCellData(numtype t);								// Export the position and size of the cells
	void f_ExportColonySize(numtype t); 							// Export the volume of colony and number of cells.
	void f_ExportPhageData(numtype t);								// Export the position and size of the phages

	// Data handling
	void OpenFileStream(std::ofstream &stream, 						// Open filstream if not allready opened
											std::string &fileName);
	std::string GeneratePath(); 									// Generates a save path for datafiles

public:
	void SetFolderNumber(int number);			 					// Sets the folder number (useful when running parallel code)
	void SetPath(const std::string &path); 							// Sets the output path (useful when running parallel code)

	// Get properties
	std::string GetPath(); 											// Returns the save path
	int GetTime();				 									// Returns the time
	numtype GetDeltaT();	 										// Returns the time-step dT

	// Operators
	Chains &operator=(const Chains &rhs); 							// Copy assignment

	// Clean up
	void DeleteFolder(); 											// Delete the data folder
private:
	void DeleteFolderTree(const char *directory_name); 				// Delete folders recursively

public:
	~Chains(); 														// Destructor
};

#endif

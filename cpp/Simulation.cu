#include "Simulation.hpp"

using namespace std;
using namespace arma;

#include "Simulation_kernels.cu.h"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
// Put these after kernel calls
//gpuErrchk( cudaPeekAtLastError() );
//gpuErrchk( cudaDeviceSynchronize() );


inline numtype cpu_sqrt(numtype x){
	#if NUMTYPE_IS_FLOAT
	return sqrtf(x);
	#else
	return sqrt(x);
	#endif
}

inline numtype cpu_sqr(numtype x){
	return x * x;
}

inline numtype cpu_pow(numtype x, numtype y){
	#if NUMTYPE_IS_FLOAT
	return powf(x, y);
	#else
	return pow(x, y);
	#endif
}

inline numtype cpu_sin(numtype x){
	#if NUMTYPE_IS_FLOAT
	return sinf(x);
	#else
	return sin(x);
	#endif
}

inline numtype cpu_cos(numtype x){
	#if NUMTYPE_IS_FLOAT
	return cosf(x);
	#else
	return cos(x);
	#endif
}

inline numtype cpu_log(numtype x){
	#if NUMTYPE_IS_FLOAT
	return logf(x);
	#else
	return log(x);
	#endif
}

// Constructors /////////////////////////////////////////////////////////////////////////
// Direct constructor
Chains::Chains(int N_max) {

    // Store the maximum number of cells in the simulation
    this->N_max = N_max;

    // Set some default parameters (initialize some default objects)
    P_0                 = 0;        // [1/µm^2] The density of invading phages in the simulation initially

    M_max               = 1e6;      // Maximum number of phages in simulation

    dT                  = 1e-6;     // [hour]   Size of the time step

    nSamp               = 100;      // Number of samples to save per simulation hour
    L                   = -1;       // [µm]     Length of boundary condition box (x direction) (L = -1 sets auto scaling)

    margin              = 1;        // 			Number of target typical length scales to simulate

    Ld                  = 3.00;     // [µm]     The length scale for division (Typical volume 1.33 µm^3)
    R                   = 0.45;     // [µm]     The "radius" of the cells

    k_int               = 500;      // [N * m]    Parameter for internal spring potential
    k_rep               = k_int/2;  // [N * m]    Parameter for repulsive potential
    k_att               = k_int/4;  // [N * m]    Parameter for attraction potential
    k_pull              = 0;        // [N]      Parameter for colony formation (colony gravity)

    gamma               = 1/dT;     //          Probability to infect cell
    beta                = 100;      //          Multiplication factor phage
    delta               = 0.003;    // [1/hour] Rate of phage decay

    r                   = 1/0.5;    // [1/hour] Rate of lysis
    T_i                 = -1;       // [hours]  Time when the phage infections begins (less than 0 disables phage infection)

    eta                 = 0.1;      // Amount of division noise along length axis (width of gaussian)
    nu                  = 0.05;     // Amount of displacement noise of the new poles (width of gaussian)

	bendingAngle        = M_PIl/6;	// [rad]		The allowed angle between two cells

    D_B                 = 0;        // [µm^2/hour] Diffusion constant for the cells
    D_P                 = 13000;    // [µm^2/hour] Diffusion constant for the phage

    Time                = 0.0;      // Counter for how many time steps have passed
    RunTime             = 0.0;  	// Variable tracking total run time

    lockCells           = false; 	// Boolean to lock the cells (stop updating)
    allCaptured         = true;	    // Boolean to stop phage updating (all are captured)

    debug               = 1;        // The amount of information to print to terminal
    exit                = false;    // Boolean to control early exit
    firstRun            = true;     // Bool to indicate if this run is the first (i.e. first time we write data)

    wellMixed           = false;    // Bool to indicate if bacteria should be well mixed

    exportAny           = false;    //
    exportCellData      = false;    // Booleans to control the export output
    exportColonySize    = false;    //
    exportPhageData     = false;    //

    ready               = false;    // Boolean to indicate whether the data is ready on the GPU

    rngSeed = -1;                   // Random number seed  ( set to -1 if unused )

}


// Copy constructor
Chains::Chains(Chains& other) {

    P_0                 = other.P_0;                        // [1/µm^2] The density of invading phages in the simulation initially

    N_max               = other.N_max;                      // Maximum number of cells in simulation
    M_max               = other.M_max;                      // Maximum number of phages in simulation

    dT                  = other.dT;                         // Size of the time step

    nSamp               = other.nSamp;                      // Number of samples to save
    L                   = other.L;                          // [µm]     Length of boundary condition box

    margin              = other.margin;

    Ld                  = other.Ld;                         // [µm]     The length scale for division (Typical volume x.x µm^3)
    R                   = other.R;                          // [µm]     The "radius" of the cells

    k_rep               = other.k_rep;                      // [N * m]    Parameter for repulsive potential
    k_att               = other.k_att;                      // [N * m]    Parameter for attraction potential
    k_int               = other.k_int;                      // [N * m]    Parameter for internal spring potential
    k_pull              = other.k_pull;

    gamma               = other.gamma;                      //          Probability to infect cell
    beta                = other.beta;
    delta               = other.delta;                      // [1/hour] Rate of phage decay

    r                   = other.r;
    T_i                 = other.T_i;                        // [hours]  Time when the phage infections begins (less than 0 disables phage infection)

    eta                 = other.eta;                        // Amount of division noise along length axis (width of gaussian)
    nu                  = other.nu;                         // Amount of displacement noise of the new poles (width of gaussian)

    bendingAngle        = other.bendingAngle;

    D_B                 = other.D_B;                        // [µm^2/hour] Diffusion constant for the cells
    D_P                 = other.D_P;                        // [µm^2/hour] Diffusion constant for the phage

    Time                = other.Time;                       // Counter for how many time steps have passed
    RunTime             = other.RunTime;        			// Variable tracking total run time

    lockCells           = other.lockCells;
    allCaptured         = other.allCaptured;

    debug               = other.debug;                      // The amount of information to print to terminal
    exit                = other.exit;                       // Boolean to control early exit
    firstRun            = other.firstRun;                   // Bool to indicate if this run is the first

    wellMixed           = other.wellMixed;

    exportAny           = other.exportAny;                  //
    exportCellData      = other.exportCellData;             // Booleans to control the export output
    exportColonySize    = other.exportColonySize;           //
    exportPhageData     = other.exportPhageData;            //

    ready               = false;                            // Boolean to indicate whether the data is ready on the GPU
    other.ready         = false;

    rngSeed             = other.rngSeed;                    // The seed for the random number generator

    // Copy random number generator
    rng = other.rng;

    // Copy the configuration arguments
    cellsBlockSize = other.cellsBlockSize;
    cellsGridSize  = other.cellsGridSize;

    phagesBlockSize = other.phagesBlockSize;
    phagesGridSize  = other.phagesGridSize;

    // Copy the device pointers
    d_cells        = other.d_cells;
    d_cells_new    = other.d_cells_new;
    d_phages       = other.d_phages;
    d_active       = other.d_active;

    // Copy GPU data to host
    h_cells        = other.d_cells;
    other.h_cells  = other.d_cells;

    h_phages       = other.d_phages;
    other.h_phages = other.d_phages;

    h_active       = other.d_active;
    other.h_active = other.d_active;

    h_rng_state       = other.d_rng_state;
    other.h_rng_state = other.d_rng_state;

}


// Controls the evaluation of the simulation
int Chains::Run(numtype T) {

    if (T < dT) {
        error("Cannot run simulation for less than dT!");
        return 1;
    }

    if (exit) {
        error("Cannot run: exit flag is set!");
        return 1;
    }

    // Get start time
    time_t  tic;
    time(&tic);

    // Things to run only when simulation is initialized
    if (Time == 0.0) {

        // Initialize the simulation matrices
        Initialize();
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );

        if (exit) return 1;

    } else {

        // Delete Completed.txt
        string path_s = path + "/Completed.txt";
        struct stat info;
        if (stat(path_s.c_str(), &info) == 0 && S_ISREG(info.st_mode)) unlink(path_s.c_str());
    }

    // Check if data is loaded and ready
    if (not ready) {

        // Ensure GPU data is current
        d_cells     = h_cells;
        d_phages    = h_phages;
        d_active    = h_active;
        d_rng_state = h_rng_state;

    }

    // Check if it is time to spawn phages
    SpawnPhages();

    // Check if export has been enabled, and if so, generate a path
    if (exportAny) path = GeneratePath();

    // Check if we have written data before
    if ((firstRun) and ((f_cells.is_open()) or (f_phages.is_open()) or (f_colonySize.is_open()))) { // This can never be true.
        firstRun = false;
    }

    // Write the reproducible command to log.txt
    if (exportAny) WriteLog();

    // Export the start configuration
    if (firstRun) {

        if ( exportCellData          and (debug > 0)) {cout << "\tExporting Cell Position Data" << endl;}
        if ( exportColonySize        and (debug > 0)) {cout << "\tExporting Colony Size" << endl;}
        if ( exportPhageData         and (debug > 0)) {cout << "\tExporting Phage Position Data" << endl;}

        // Export data
        if (exportAny) ExportData(Time);

    }

    // Run the time evolution

    // Determine the number of time steps between samplings
    int nStepsPerSample = (int)round(1/(nSamp * dT));

    // Determine the number of samples to take
    int nSamplings = nSamp * T;

    // Store the current time
    numtype Time_0 = Time;

    // Loop over samplings
    int t = 0;
    for (int n = 1; n <= nSamplings; n++) {

        // Check for exit flag
        if (exit) { break; }

        // Run time inside samples
        while (n * T / nSamplings - t * dT > dT / 2) {

            // Count steps
            t++;

            // Check for exit flag
            if (exit) { break; }

            // Compute the time
            Time = Time_0 + t * dT;

            // Update remaining phages
            if (not allCaptured) {
                PhageUpdateKernel<<<phagesGridSize, phagesBlockSize>>>(thrust::raw_pointer_cast(&d_phages[0]),
                    thrust::raw_pointer_cast(&d_active[0]),
                    thrust::raw_pointer_cast(&d_cells[0]),
                    cpu_sqrt(2 * D_P * dT),
                    L,
                    R,
                    thrust::raw_pointer_cast(&d_rng_state[0]),
                    h_active.size(),
                    N_max,
                    false);
            }

            // Update cells
            if (not lockCells) {

                CellUpdateKernel<<<cellsGridSize, cellsBlockSize>>>(thrust::raw_pointer_cast(&d_cells[0]),
                    thrust::raw_pointer_cast(&d_cells_new[0]),
                    Ld,
                    R,
                    k_rep,
                    k_att,
                    k_int,
                    k_pull,
                    dT,
                    L,
                    N_max);

                thrust::swap(d_cells, d_cells_new);
            }

            // Spawn phages
            SpawnPhages();
        }

        // Export the data
        if (exportAny) ExportData(Time);

        // Show progress bar
        if ((n > 0) and (debug > 0)) {
            cout << "\t[";
            int pos = 60 * static_cast<float>(n) / static_cast<float>(nSamplings);;
            for (int i = 0; i < 60; ++i) {
                if (i <= pos) cout << ".";
                else cout << " ";
            }
            cout << "] " << "\r";
            cout.flush();
        }
    }

    // Get stop time
    time_t  toc;
    time(&toc);

    // Calculate time difference
    RunTime += difftime(toc, tic);
    float seconds = RunTime;
    float hours   = floor(seconds/3600);
    float minutes = floor(seconds/60);
    minutes -= hours * 60;
    seconds -= minutes * 60 + hours * 3600;

    if (debug > 0) {
        cout << endl;
        cout << "Simulation complete after ";
        if (hours > 0.0)   cout << hours   << " hours and ";
        if (minutes > 0.0) cout << minutes << " minutes and ";
        cout  << seconds << " seconds." << endl;
    }

    std::ofstream f_out;
    f_out.open(GetPath() + "/Completed.txt", fstream::trunc);
    f_out << "\tSimulation complete after ";
    if (hours > 0.0)   f_out << hours   << " hours and ";
    if (minutes > 0.0) f_out << minutes << " minutes and ";
    f_out  << seconds << " seconds." << "\n";
    f_out.flush();
    f_out.close();

    // Write success to log
    if (exit) {
        f_log << ">>Simulation completed with exit flag<<" << endl;
    }

    if (exit) {
        return 1;
    } else {
        return 0;
    }
}


// Initialize the simulation
void Chains::Initialize() {
    deb("Initializing", 1);

    // Set the random number generator seed
    if (rngSeed >= 0.0) {
        rng.seed( rngSeed );
    } else {
        static std::random_device rd;
        rng.seed(rd());
    }

    // Initialize the cells
    deb("- Spawning cells", 1);

    // Compute GPU block and grid size
    cellsBlockSize = 256;
    cellsGridSize = (N_max + cellsBlockSize - 1) / cellsBlockSize;

    // Generate new cell vector
    h_cells.reserve(N_max * 8);

    // Initialize chains
    if (not wellMixed) {

        // Keep track of the angles
        numtype theta;
        numtype phi;

        // And keep track of the center
        numtype center[3] = {0.0, 0.0, 0.0};
        for (int n = 0; n < N_max; n++) {

            // Allocate coordinates
            numtype xP, yP, zP;
            numtype xQ, yQ, zQ;
            numtype xR, yR, zR;
            numtype xS, yS, zS;

            // Generate location 1
            if (n == 0) {
                // First iteration, start at (0, 0, 0)
                xP = 0.0;
                yP = 0.0;
                zP = 0.0;

                // Choose random angle
                theta = bendingAngle * rand(rng) + M_PIl / 4;
                phi   =    2 * M_PIl * rand(rng);

                // Generate location 2 coordinates
                xQ = xP + Ld * sin(theta) * cos(phi);
                yQ = yP + Ld * sin(theta) * sin(phi);
                zQ = zP + Ld * cos(theta);

            } else {

                // Choose random angle (giving the vector r)
                theta =     M_PIl * rand(rng);
                phi   = 2 * M_PIl * rand(rng);

                // Compute r dot PQ
                numtype rdotPQ = cpu_sin(theta) * cpu_cos(phi) * (xQ - xP) + cpu_sin(theta) * cpu_sin(phi) * (yQ - yP) + cpu_cos(theta) * (zQ - zP);
                //                              xr                  xPQ                    yr                   yPQ            zr            zPQ

                // Generate U-vector which is perpendicular to PQ
                numtype xU, yU, zU;
                xU = cpu_sin(theta) * cpu_cos(phi) - rdotPQ * (xQ - xP) / cpu_pow(Ld, 2);
                yU = cpu_sin(theta) * cpu_sin(phi) - rdotPQ * (yQ - yP) / cpu_pow(Ld, 2);
                zU = cpu_cos(theta)                - rdotPQ * (zQ - zP) / cpu_pow(Ld, 2);

                // Normalize U
                numtype normU = cpu_sqrt(xU * xU + yU * yU + zU * zU);
                xU /= normU;
                yU /= normU;
                zU /= normU;

                // Choose a random angle to rotate PQ around U
                theta = bendingAngle * rand(rng);

                // Compute rotation (store in temporarily as vector S)
                numtype ct = cpu_cos(theta);
                numtype st = cpu_sin(theta);

                xS = (ct + xU * xU * (1 - ct))           * (xQ - xP) +      (xU * yU * (1 - ct) - zU * st) * (yQ - yP) +      (xU * zU * (1 - ct) + yU * st) * (zQ - zP);
                yS = (     xU * yU * (1 - ct) + zU * st) * (xQ - xP) + (ct + yU * yU * (1 - ct))           * (yQ - yP) +      (yU * zU * (1 - ct) - xU * st) * (zQ - zP);
                zS = (     xU * zU * (1 - ct) - yU * st) * (xQ - xP) +      (yU * zU * (1 - ct) + xU * st) * (yQ - yP) + (ct + zU * zU * (1 - ct))           * (zQ - zP);

                // Move vector R, a distance of two cell radii along the average angle
                xR = xQ + (xQ - xP + xS ) * R / Ld;
                yR = yQ + (yQ - yP + yS ) * R / Ld;
                zR = zQ + (zQ - zP + zS ) * R / Ld;

                // Generate the location of the second pole
                xS += xR;
                yS += yR;
                zS += zR;

                // Store coordinates in xP and xQ
                xP = xR;
                yP = yR;
                zP = zR;

                xQ = xS;
                yQ = yS;
                zQ = zS;
            }

            // Compute center location
            numtype x = 0.5 * (xP + xQ);
            numtype y = 0.5 * (yP + yQ);
            numtype z = 0.5 * (zP + zQ);

            // Store the center
            center[0] += x / N_max;
            center[1] += y / N_max;
            center[2] += z / N_max;

            // Add cell to system
            h_cells.push_back(xP);
            h_cells.push_back(yP);
            h_cells.push_back(zP);
            h_cells.push_back(xQ);
            h_cells.push_back(yQ);
            h_cells.push_back(zQ);

            // Connect the cell to the next cell
            h_cells.push_back(0);
            h_cells.push_back(0);

            // Connect the cells
            if (n > 0) {
                // Connect the cell to the previous cell
                h_cells[8*n+6] = -n;

                // Connect the previous cell to the cell
                h_cells[8*(n-1)+7] = n + 1;
            }

        }

        // Center the chain on (0, 0, 0)
        for (int n = 0; n < N_max; n++) {

            h_cells[8*n+0] -= center[0];
            h_cells[8*n+1] -= center[1];
            h_cells[8*n+2] -= center[2];
            h_cells[8*n+3] -= center[0];
            h_cells[8*n+4] -= center[1];
            h_cells[8*n+5] -= center[2];

        }

        // Update L values (Include safety margin)
        if (L == -1) {

            // Determine radius of smallest sphere that can encapsulate the sphere
            numtype m = 0;

            // Loop over cells
            for (int n = 0; n < N_max; n++) {
                numtype xP = h_cells[8*n+0];
                numtype yP = h_cells[8*n+1];
                numtype zP = h_cells[8*n+2];
                numtype xQ = h_cells[8*n+3];
                numtype yQ = h_cells[8*n+4];
                numtype zQ = h_cells[8*n+5];

                // Compute reach of P and Q coordinates
                numtype rP = xP * xP + yP * yP + zP * zP;
                numtype rQ = xQ * xQ + yQ * yQ + zQ * zQ;

                // Store largest
                if (rP > m) m = rP;
                if (rQ > m) m = rQ;
            }

            // Convert to radius
            m = cpu_sqrt(m)+R;

            // Set new L value
            L = 2 * m * margin;
        }

        // Center the chain on (L, L, L) / 2
        for (int n = 0; n < N_max; n++) {
            for (int j = 0; j < 6; j++) {
                h_cells[8*n+j] += L / 2;
            }
        }

    } else {    // Initialize well mixed

        // Declare variables
        numtype xP, yP, zP;
        numtype xQ, yQ, zQ;
        numtype theta;
        numtype phi;

        // Spawn bacteria
        for (int n = 0; n < N_max; n++) {

            // Draw random location for cell
            xP = -1;
            yP = -1;
            zP = -1;
            xQ = -1;
            yQ = -1;
            zQ = -1;

            // While any are out of bounds
            while ((xP < 0) or (xP > L) or (yP < 0) or (yP > L) or (zP < 0) or (zP > L) or (xQ < 0) or (xQ > L) or (yQ < 0) or (yQ > L) or (zQ < 0) or (zQ > L)) {

                // Generate coordinates
                xP = rand(rng) * L;
                yP = rand(rng) * L;
                zP = rand(rng) * L;

                 // Choose random angle
                theta =     M_PIl * rand(rng);
                phi   = 2 * M_PIl * rand(rng);

                // Generate location 2 coordinates
                xQ = xP + Ld * sin(theta) * cos(phi);
                yQ = yP + Ld * sin(theta) * sin(phi);
                zQ = zP + Ld * cos(theta);

            }

            // Add cell to system
            h_cells.push_back(xP);
            h_cells.push_back(yP);
            h_cells.push_back(zP);
            h_cells.push_back(xQ);
            h_cells.push_back(yQ);
            h_cells.push_back(zQ);

            // Connect the cell to the previous cell
            h_cells.push_back(0);
            h_cells.push_back(0);
        }
    }

    // Copy cell data to GPU
    d_cells     = h_cells;
    d_cells_new = d_cells;

    // Compute GPU block and grid size
    cellsBlockSize = 256;
    cellsGridSize = (N_max + cellsBlockSize - 1) / cellsBlockSize;

    // Set the ready flag
    ready = true;

    // Give warnings
    if (k_rep * dT > 0.05) {
        cout << "Time step might not be small enough! (k_rep * dT = " << k_rep * dT << ")" << endl;
    }
    if (k_att * dT > 0.05) {
        cout << "Time step might not be small enough! (k_att * dT = " << k_att * dT << ")" << endl;
    }
    if (k_int * dT > 0.05) {
        cout << "Time step might not be small enough! (k_int * dT = " << k_int * dT << ")" << endl;
    }
}


// Simulation functions /////////////////////////////////////////////////////////////////

// Spawns phages according to spawning rules
void Chains::SpawnPhages() {

    if (Time < T_i) {return;}

    if (P_0 != 0.0) {
        deb("Spawning Phages", 1);

        // Compute the number of phages and allocate space
        int M = 0;
        if (P_0 > 0) {
            M = (int)(round(P_0 * cpu_pow(L, 3) / 1e12));
        } else {
            M = (int)(round(-P_0));
        }

        // Reset P_0
        P_0 = 0.0;

        // Create thrust array to store phage locations in
        h_phages.reserve(M);
        h_active.reserve(M);
	    phagesBlockSize = 256;
        phagesGridSize = (M + phagesBlockSize - 1) / phagesBlockSize;

        // Spawn phages uniformly within the space
        for (int m = 0; m < M; m++) {

            // Allocate coordinates
            numtype x = -1;
            numtype y = -1;
            numtype z = -1;

            // Use hit and miss method to generate uniformly distributed phages
            while ((x < 0) or (x > L) or (y < 0) or (y > L) or (z < 0) or (z > L)) {
                x = rand(rng) * L;
                y = rand(rng) * L;
                z = rand(rng) * L;
            }

            // Spawn the new phage
            h_phages.push_back(x);
            h_phages.push_back(y);
            h_phages.push_back(z);
            h_active.push_back(1);

        }

        // Copy phage data to GPU
        d_phages = h_phages;
        d_active = h_active;

        // Initialize rng on device
        d_rng_state.resize(M);
        initRNG<<<phagesGridSize, phagesBlockSize>>>(thrust::raw_pointer_cast(&d_rng_state[0]), M);

        // Set the capture boolean
        allCaptured = false;

        // Check for failed phage spawning
        if (M == 0) {
            allCaptured = true;
        }

        // Set the ready flag
        ready = true;

        deb("- Done!", 1);
    }
}

// Replace cell I with beta phages
void Chains::LyseCell(int I) {

    deb("Lysing cell", 1);

    if (I >= N_max) {
        error("Cannot lyse cell, I >= N_max");
    }

    // Copy bacteria to Host
    h_cells  = d_cells;

    // Allocate space for phages
    int M = h_phages.size();
    if (M + 3 * beta > h_phages.capacity()) {
        h_phages.reserve(M + 3 * beta);
        h_active.reserve(M + beta);
    }

    // Create thrust array to store phage locations in
    phagesBlockSize = 256;
    phagesGridSize = (M + 3 * beta + phagesBlockSize - 1) / phagesBlockSize;

    // Extract information of cell I
    numtype xP  = h_cells[8 * I + 0];
    numtype yP  = h_cells[8 * I + 1];
    numtype zP  = h_cells[8 * I + 2];

    numtype xQ  = h_cells[8 * I + 3];
    numtype yQ  = h_cells[8 * I + 4];
    numtype zQ  = h_cells[8 * I + 5];

    // Define end points of the cell
    colvectype P = {xP, yP, zP};
    colvectype Q = {xQ, yQ, zQ};

    // Define coordinate system with PQ as Z axis
    mattype PQ = zeros<mattype>(3, 3);
    PQ.col(0) = Q - P;

    mattype q, r;
    qr( q, r, PQ );

    mattype XYZ = join_horiz( q.cols(1, 2), PQ.col(0) / norm(PQ.col(0)));

    // Define phage location vectors
    colvectype T;

    // Spawn new phages (Uniformly on the interior of the bursted cell)
    for (int b = 0; b < beta; b++) {

        // Generate segment of cell where phage is spawned
        numtype r = rand(rng);

        // Generate a distance
        numtype d = R * rand(rng);

        // Generate location for phage
        numtype theta = M_PIl / 2 * rand(rng);
        numtype phi   = M_PIl * 2 * rand(rng);

        // Determine where on the cell the phage is located
        if (r < (2 * R / 3) / (Ld + 4 * R / 3)) {          // Phage is located in the top half-sphere

            // Determine translation vectors
            T = P;
            T += cos(theta) * cos(phi) * d * XYZ.col(0) + cos(theta) * sin(phi) * d * XYZ.col(1) - sin(theta) * d * XYZ.col(2);


        } else if (r < (4 * R / 3) / (Ld + 4 * R / 3)) {   // Phage is located in the bottom half-sphere

            // Determine translation vector
            T = Q;
            T += cos(theta) * cos(phi) * d * XYZ.col(0) + cos(theta) * sin(phi) * d * XYZ.col(1) + sin(theta) * d * XYZ.col(2);


        } else {                                            // Phage is located along the cylindrical part

            // Determine translation vector
            T = P + rand(rng) * (Q-P);
            T += cos(phi) * d * XYZ.col(0) + sin(phi) * d * XYZ.col(1);

        }

        // Add new phage
        h_phages.push_back(T(0));
        h_phages.push_back(T(1));
        h_phages.push_back(T(2));
        h_active.push_back(1);

    }

    // Set the capture boolean
    if (beta > 0) allCaptured = false;

    // Initialize rng on device
    d_rng_state.reserve(M + 3 * beta);
    initRNG<<<phagesGridSize, phagesBlockSize>>>(thrust::raw_pointer_cast(&d_rng_state[0]), M + 3 * beta);

    // Remove cell I
    d_cells.erase(d_cells.begin() + 8 * I, d_cells.begin() + 8 * (I + 1));
    N_max--;

    // Copy to the GPU
    d_phages = h_phages;
    d_active = h_active;

    // Set the ready flag
    ready = true;

    deb("- Done!", 1);
}

// Autoscale the simulation space
void Chains::AutoScale() {

    // Store P_0 and T_i value
    numtype P_0 = this->P_0;
    numtype T_i = this->T_i;

    // Copy cells to CPU
    h_cells = d_cells;

    // Set L to be smallest posable box
    numtype L_s = 0;

    // If simulating a chain, use small test box
    if (not wellMixed) {

        // Determine most extreme coordinate
        for (int n = 0; n < N_max; n++) {
            for (int j = 0; j < 6; j++) {
                if (abs(h_cells[8 * n + j] - L/2) > L_s) L_s = abs(h_cells[8 * n + j] - L / 2);
            }
        }

        // Add R to extent
        L_s += R;

        // Use the smallest L value
        numtype L_test = max(20.0, 6 * L_s);

        // Center the chain on (L, L, L)
        for (int n = 0; n < N_max; n++) {
            for (int j = 0; j < 6; j++) {
                h_cells[8 * n + j] += (L_test - L) / 2;
            }
        }

        // Update L values
        L = L_test;

    } else {

        // Set L_s value when well-mixed
        L_s = L / 4.0;
    }

    // Copy cells to GPU
    d_cells = h_cells;


    // Auto scale with P_0 = -1e5 first, then -1e6, and then -1e7 if it fails /////////////////////
    numtype M;
    numtype T;
    for (int i = 0; i < 3; i ++) {

        // Overwrite P_0 and T_i value
        if (i == 1)      this->P_0 = P_0 * 10;
        else if (i == 2) this->P_0 = P_0 * 100;
        this->T_i = 0;

        // Spawn phages
        SpawnPhages();

        // Count the number of phages
        int nPhages = h_phages.size() / 3;

        // Take steps until enough phage are adsorbed
        int n = 0;                     // Number of iterations run through
        int t_step = round(1e-2/dT);   // Time between checks

        T = 0.0;        // Elapsed simulation time
        M = 1.0;        // Fraction of free phages (relative to start)

        while ((T < 0.1) and (M > 0.3679)) {
            for (int t = 0; t < t_step; t++) {

                // Increment counter
                n++;

                // Update remaining phages
                PhageUpdateKernel<<<phagesGridSize, phagesBlockSize>>>(thrust::raw_pointer_cast(&d_phages[0]),
                    thrust::raw_pointer_cast(&d_active[0]),
                    thrust::raw_pointer_cast(&d_cells[0]),
                    cpu_sqrt(2 * D_P * dT),
                    L,
                    R,
                    thrust::raw_pointer_cast(&d_rng_state[0]),
                    nPhages,
                    N_max,
                    false);
            }

            // Compute the time
            T = n * dT;

            // Update M
            M = static_cast<numtype>(thrust::reduce(d_active.begin(), d_active.end())) / static_cast<numtype>(h_active.size());

        }

        // Make sure some phages have been adsorbed
        if (M == 1.0) {
            if (i == 0) {
                warning("No hits during autoscaling! Retrying...");

                // Clean up
                h_phages.clear();
                d_phages.clear();

            } else {
                error("No hits during autoscaling!");
            }

        } else {
            break;
        }
    }

    // Estimate adsorption rate eta
    numtype eta = - cpu_pow(L, 3) * cpu_log(M) / T;

    // Set new volume size
    numtype V = - eta / cpu_log(0.9); // Scale so that 90% of phage remain free

    // Scale lengths
    numtype L_new;

    // Set new size
    L_new = cpu_pow(V, 1.0/3.0);
    if (L_new < 2 * 2 * L_s) {

        // If margin is too small, keep the same density of test phages
        P_0 = - P_0 * 1e12 / cpu_pow(L_new, 3.0);

        // Scale the space to be larger
        L_new = 2 * 2 * L_s;

    }

    // If simulating a chain, re-center it
    if (not wellMixed) {

        // Center the chain on (L_new, L_new, L_new) / 2
        for (int n = 0; n < N_max; n++) {
            for (int j = 0; j < 6; j++) {
                h_cells[8 * n + j] += (L_new - L) / 2;
            }
        }

    } else {

        // Redraw well mixed bacteria
        h_cells.clear();

        // Declare variables
        numtype xP, yP, zP;
        numtype xQ, yQ, zQ;
        numtype theta;
        numtype phi;

        // Spawn bacteria
        for (int n = 0; n < N_max; n++) {

            // Allocate coordinates
            xP = -1;
            yP = -1;
            zP = -1;
            xQ = -1;
            yQ = -1;
            zQ = -1;

            // While any are out of bounds
            while ((xP < 0) or (xP > L_new) or (yP < 0) or (yP > L_new) or (zP < 0) or (zP > L_new) or (xQ < 0) or (xQ > L_new) or (yQ < 0) or (yQ > L_new) or (zQ < 0) or (zQ > L_new)) {

                // Generate coordinates
                xP = rand(rng) * L_new;
                yP = rand(rng) * L_new;
                zP = rand(rng) * L_new;

                // Choose random angle
                theta =     M_PIl * rand(rng);
                phi   = 2 * M_PIl * rand(rng);

                // Generate location 2 coordinates
                xQ = xP + Ld * sin(theta) * cos(phi);
                yQ = yP + Ld * sin(theta) * sin(phi);
                zQ = zP + Ld * cos(theta);

            }

            // Add cell to system
            h_cells.push_back(xP);
            h_cells.push_back(yP);
            h_cells.push_back(zP);
            h_cells.push_back(xQ);
            h_cells.push_back(yQ);
            h_cells.push_back(zQ);
            h_cells.push_back(0);
            h_cells.push_back(0);
        }
    }

    // Store new L values
    L = L_new;

    // Update GPU values
    d_cells = h_cells;

    // Clean up
    h_phages.clear();
    d_phages.clear();

    h_active.clear();
    d_active.clear();

    allCaptured = true;

    // Reset P_0 and T_i value
    this->P_0 = P_0;
    this->T_i = T_i;

    // Update log.txt
    if (exportAny) WriteLog();
}


// Auto relaxes the bacteria
void Chains::Relax() {

    // Relax for 0.1 hour
    Run(0.1);

    // If no overlap possible, return
    if (N_max < 2) return;

    // Create vector for overlaps
    thrust::device_vector<numtype> d_overlaps;
    d_overlaps.resize(N_max - 1);

    // Detect current overlap
    CellOverlaps<<<cellsGridSize, cellsBlockSize>>>(
        thrust::raw_pointer_cast(&d_cells[0]),
        thrust::raw_pointer_cast(&d_overlaps[0]),
        R,
        N_max);

    numtype maxOverlap = *(thrust::max_element(d_overlaps.begin(), d_overlaps.end()));

    // While overlap exists, advance time
    while (maxOverlap / R > 0.01) {

        // Relax for 0.1 additional hour
        Run(0.1);

        // Detect current overlap
        CellOverlaps<<<cellsGridSize, cellsBlockSize>>>(
            thrust::raw_pointer_cast(&d_cells[0]),
            thrust::raw_pointer_cast(&d_overlaps[0]),
            R,
            N_max);

        // Determine the overlap
        numtype maxOverlap_new = *(thrust::max_element(d_overlaps.begin(), d_overlaps.end()));

        // Detect convergence
        if (abs(maxOverlap - maxOverlap_new) < 1e-8) {
            break;
        }

        // Save overlap
        maxOverlap = maxOverlap_new;

    }
}


// Equilibrates the phage distribution
void Chains::Equilibrate(numtype T) {

    if (T < dT) {
        return;
    }

    if (Time < T_i) {
        error("Cannot equilibrate before phage invasion begins!");
        return;
    }

    if (exit) {
        error("Cannot equilibrate: exit flag is set!");
        return;
    }

    // Things to run only when simulation is initialized
    if (Time == 0.0) {

        // Initialize the simulation matrices
        Initialize();

        if (exit) return;

    } else {

        // Delete Completed.txt
        string path_s = path + "/Completed.txt";
        struct stat info;
        if (stat(path_s.c_str(), &info) == 0 && S_ISREG(info.st_mode)) unlink(path_s.c_str());

    }

    // Check if it is time to spawn phages
    SpawnPhages();

    // Write the reproducible command to log.txt
    if (exportAny) WriteLog();


    // Run the time evolution

    // Determine the number of time steps to run
    int nSteps = (int)round(T / dT);

    // Store the current time
    numtype Time_0 = Time;

    // Loop over steps
    int t = 0;
    while (T - t * dT > dT / 2) {

        // Check for exit flag
        if (exit) { break; }

        // Count steps
        t++;

        // Compute the time
        Time = Time_0 + t * dT;

        // Update remaining phages
        PhageUpdateKernel<<<phagesGridSize, phagesBlockSize>>>(thrust::raw_pointer_cast(&d_phages[0]),
            thrust::raw_pointer_cast(&d_active[0]),
            thrust::raw_pointer_cast(&d_cells[0]),
            cpu_sqrt(2 * D_P * dT),
            L,
            R,
            thrust::raw_pointer_cast(&d_rng_state[0]),
            h_active.size(),
            N_max,
            true);

        // Spawn phages
        SpawnPhages();
    }

}

// Settings /////////////////////////////////////////////////////////////////////////////
// Set the size of the time-step
void Chains::TimeStep(numtype dT) {this->dT = dT;}


// Set the length of the simulation space
void Chains::SetLength(numtype L) {
    this->L = L;
}


// Set the margin (number of length scales to simulate)
void Chains::SetMargin(numtype margin) {this->margin = margin;}


// Sets the time when the phages should start infecting
void Chains::PhageInvasionStartTime(numtype T_i) {this->T_i = T_i;}


// Sets initial density of the phages (1/µm^3)
void Chains::PhageInitialDensity(numtype P_0) {this->P_0 = P_0;}


// Sets the diffusion constant of the phages
void Chains::PhageDiffusionConstant(numtype D_P) { this->D_P = D_P;}


// Sets rate of the infection increasing in stage
void Chains::PhageInfectionRate(numtype r) {this->r = r;}


// Set the decay rate of the phages
void Chains::PhageDecayRate(numtype delta) {this->delta = delta;}


// Set the size of the bursts
void Chains::PhageBurstSize(int beta) {this->beta = beta;}


// Changes the adsorption parameter gamma
void Chains::PhageAdsorptionParameter(numtype gamma) {this->gamma = gamma;}


// Sets the diffusion constant of the bacteria
void Chains::CellDiffusionConstant(numtype D_B) {this->D_B = D_B;}


// Sets the strength of the repulsive potential
void Chains::CellRepulsiveParameter(numtype k_rep) {this->k_rep = k_rep;}


// Sets the strength of the attractive potential
void Chains::CellAttractiveParameter(numtype k_att) {this->k_att = k_att;}


// Sets the strength of the internal potential
void Chains::CellInternalParameter(numtype k_int) {this->k_int = k_int;}


// Sets the (division)length of the bacteria
void Chains::CellLength(numtype Ld) {this->Ld = Ld;}


// Sets the radius of the bacteria
void Chains::CellRadius(numtype R) {this->R = R;}


// Sets the bending angle between cells
void Chains::CellBendingAngle(numtype bendingAngle){this->bendingAngle = bendingAngle;}


// Locks the cells in their current configuration
void Chains::CellLock() {lockCells = true;}


// Sets the bacteria to be well mixed
void Chains::WellMixed() {wellMixed = true;}


// Sets the bacteria to form a spherical colony (with force k_pull)
void Chains::SphericalColony(numtype k_pull) {this->k_pull = k_pull;}


// Sets the maximum number of phages in the simulation
void Chains::MaxPhageCount(int M_max) {this->M_max = M_max;}

// Helping functions ////////////////////////////////////////////////////////////////////

// Sets the seed of the random number generator
void Chains::SetRngSeed(int n) {
    rngSeed = n;
}

// Debug function (prints "input")
void Chains::deb(const std::string& input, int n) {

    if (debug >= n) {
        std::stringstream stream;
        stream << "<<" << input << ">>" << endl;
        cout << stream.str();
        cout.flush();
    }
}

// Error function (prints "input" and set exit = true)
void Chains::error(const std::string& input) {
    if (not f_log.is_open()) f_log.open(path + "/log.txt", fstream::app);
    cerr << "\t>> " << input << " Exiting...<<" << endl;
    f_log << ">> " << input << " Exiting...<<" << endl;
    exit = true;
}

// Warning function (prints "input")
void Chains::warning(const std::string &input) {
    if (not f_log.is_open()) f_log.open(path + "/log.txt", fstream::app);
    cerr << "\t>> " << input << "<<" << endl;
    f_log << ">> " << input << "<<" << endl;
}

// Write a log.txt file
void Chains::WriteLog() {
    deb("Writing log", 1);
    if ((not f_log.is_open()) and (not exit)) {

        // Open the file stream and write the command
        f_log.open(path + "/log.txt", fstream::out);
        deb("- " + path + "/log.txt", 2);

        // Physical parameters
        f_log << "P_0 = "      << P_0      << "\n";        // Density of invading phages
        f_log << "N_max = "    << N_max    << "\n";        // Maximal allowed cells
        f_log << "M_max = "    << M_max    << "\n";        // Maximum allowed phages

        f_log << "T_i = "      << T_i      << "\n";        // Phage infection time
        f_log << "L = "        << L        << "\n";        // Length of boundary condition box
        f_log << "k_rep = "    << k_rep    << "\n";        // Strength of the repulsion potential
        f_log << "k_att = "    << k_att    << "\n";        // Strength of the attraction potential
        f_log << "k_int = "    << k_int    << "\n";        // Strength of the internal spring potential
        f_log << "k_pull = "   << k_pull   << "\n";        // Parameter for colony formation (colony gravity)

        f_log << "R = "        << R        << "\n";        // Critical radius
        f_log << "Ld = "       << Ld       << "\n";        // The length scale for division
        f_log << "D_B = "      << D_B      << "\n";        // Bacteria diffusion
        f_log << "D_P = "      << D_P      << "\n";        // Phage diffusion
        f_log << "gamma = "    << gamma    << "\n";        // Probability to infect cell
        f_log << "beta = "     << beta     << "\n";        // Multiplication factor phage
        f_log << "delta = "    << delta    << "\n";        // Rate of phage decay
        f_log << "r = "        << r        << "\n";        // Rate of lysis

        f_log << "eta = "      << eta      << "\n";        // Amount of division noise
        f_log << "nu = "       << nu       << "\n";        // Amount of displacement noise

        f_log << "theta = "    << bendingAngle << "\n";    // [rad]		The allowed angle between two cells

        // Non physical parameters
        f_log << "dT = "       << dT       << "\n";        // Time step size
        f_log << "nSamp = "    << nSamp    << "\n";        // Number of samples to save per simulation hour
        f_log << "margin = "   << margin   << "\n";        // The number of length scales to simulate
        f_log << "rngSeed = "  << rngSeed  << "\n";        // Random Number Generator seed
        f_log << "wellMixed = "<< wellMixed<< "\n";        // Output the well mixed boolean

        f_log << "debug = "             << debug            << "\n";
        f_log << "exportAny = "         << exportAny        << "\n";
        f_log << "exportCellData = "    << exportCellData   << "\n";
        f_log << "exportColonySize = "  << exportColonySize << "\n";
        f_log << "exportPhageData = "   << exportPhageData  << endl;

        f_log.close(); // By closing it, we allow it to be overwritten by next .Run()
    }
}


// Set debug level to 0
void Chains::Quiet() { debug=0; };


// Set the debug level to n
void Chains::Debug(int n) {debug = n;}


// Set the number of output samples
void Chains::SetSamples(int nSamp) {this->nSamp = nSamp;};


// File outputs /////////////////////////////////////////////////////////////////////////
// Sets boolean for export function
void Chains::ExportCellData()         { exportCellData          = true; exportAny = true; };
void Chains::ExportColonySize()       { exportColonySize        = true; exportAny = true; };
void Chains::ExportPhageData()        { exportPhageData         = true; exportAny = true; };


void Chains::ExportCellDataNow()      { f_ExportCellData(Time); };
void Chains::ExportPhageDataNow()     { f_ExportPhageData(Time); };

// Master function to export the data
void Chains::ExportData(numtype t) {
    deb("Exporting data", 2);

    // Export the data
    if ( exportCellData ) {
        f_ExportCellData(t);
    }
    if ( exportColonySize ) {
        f_ExportColonySize(t);
    }
    if ( exportPhageData ) {
        f_ExportPhageData(t);
    }

    deb("- Done", 2);
}


// Export the position and size of the cells
void Chains::f_ExportCellData(numtype t) {

    // Verify the file stream is open
    string fileName = "CellData";
    OpenFileStream(f_cells, fileName);

    // Copy from GPU
    h_cells = d_cells;

    // Loop over cells in simulation
    for (int i = 0; i < N_max; i++) {

        // Output format:   (T    x    y    z   r   s)
        f_cells << fixed    << setprecision(5);
        f_cells << setw(8)  << t << "\t";
        f_cells << fixed    << setprecision(8);
        f_cells << setw(12) << h_cells[8 * i + 0] << "\t";
        f_cells << setw(12) << h_cells[8 * i + 1] << "\t";
        f_cells << setw(12) << h_cells[8 * i + 2] << "\t";
        f_cells << setw(12) << h_cells[8 * i + 3] << "\t";
        f_cells << setw(12) << h_cells[8 * i + 4] << "\t";
        f_cells << setw(12) << h_cells[8 * i + 5] << "\t";
        f_cells << fixed    << setprecision(0);
        f_cells << setw(12) << h_cells[8 * i + 6] << "\t";
        f_cells << setw(12) << h_cells[8 * i + 7] << endl;
    }
}


// Export the position and size of the phages
void Chains::f_ExportPhageData(numtype t) {

    // Verify the file stream is open
    string fileName = "PhageData";
    OpenFileStream(f_phages, fileName);

    // Copy from GPU
    h_phages = d_phages;

    // Loop over phages in simulation
    for (int i = 0; i < h_phages.size() / 3; i++) {

        // Output format:   T    x    y    z
        f_phages << fixed    << setprecision(5);
        f_phages << setw(8)  << t << "\t";
        f_phages << setw(12) << h_phages[3 * i + 0] << "\t";
        f_phages << setw(12) << h_phages[3 * i + 1] << "\t";
        f_phages << setw(12) << h_phages[3 * i + 2] << endl;
    }
}


// Export the volume of colony and number of cells.
void Chains::f_ExportColonySize(numtype t) {

    // Verify the file stream is open
    string fileName = "ColonySize";
    OpenFileStream(f_colonySize, fileName);

    // Writes the time, total volume, number of cells,
    // number of lytic stage cells, number of lysogenic cells

    //cout << "d_active.size() = " << d_active.size() << endl;
    //cout << "phagesBlockSkize = " << phagesBlockSize << endl;
    //cout << "phagesGridSize = " << phagesGridSize << endl;
    int M = 0;
    if ((not allCaptured) and (h_active.size() > 0)) {
        M = thrust::reduce(d_active.begin(), d_active.end());
        if (M == 0) allCaptured = true; // If all are captured, stop phage updating
    }

    f_colonySize << fixed    << setprecision(5);
    f_colonySize << setw(8)  << t           << "\t";
    f_colonySize << setw(8)  << M           << endl;
}


// Data Handling ////////////////////////////////////////////////////////////////////////
// Open file-stream if not already opened
void Chains::OpenFileStream(ofstream& stream, string& fileName) {

    // Check that if file stream is open.
    if ((not stream.is_open()) and (not exit)) {

        // Debug info
        if (debug > 0) {cout << "\t\tSaving data to file: " << path << "/" << fileName << ".txt" << endl;}

        // Construct path
        string streamPath;
        streamPath = path + "/" + fileName + ".txt";

        // Open the file stream
        if (firstRun) {
            stream.open(streamPath, fstream::trunc);
        } else {
            stream.open(streamPath, fstream::app);
        }

        // Check stream is open
        if ((not exit) and (not stream.is_open())) {
            cerr << "\t>>Could not open filestream \"" << streamPath << "\"! Exiting..<<" << endl;
            f_log <<  ">>Could not open filestream \"" << streamPath << "\"! Exiting..<<" << endl;
            exit = true;
        };

        // If appending to existing file, do not rewrite the meta data
        if (not firstRun) {
            return;
        }
    }
}


// Generates a save path for data-files
string Chains::GeneratePath() {

    // Generate a directory path
    string prefix = "data";    // Data folder name

    // Create the path variable
    string path_s;

    // Create info object
    struct stat info;

    // Check if user has specified numbered folder
    if (path.empty()) {

        // Add the prefix
        path_s += prefix;

        // Check if path exists
        if (not(stat(path_s.c_str(), &info) == 0 && S_ISDIR(info.st_mode))) {
            // Create path if it does not exist
            mkdir(path_s.c_str(), 0700);
        }

        // Loop over folders in date folder, to find current number
        int currentNumerateFolder = 1;
        DIR * dir;
        if ((dir = opendir (path_s.c_str())) != NULL) {
            struct dirent * ent;
            while ((ent = readdir (dir)) != NULL) {
                if (ent->d_type == DT_DIR) {
                    // Skip . or ..
                    if (ent->d_name[0] == '.') {continue;}
                    currentNumerateFolder++;        // Increment folder number
                }
            }
            closedir (dir);
        }

        // Append numerate folder
        path_s += "/";
        path_s += to_string(currentNumerateFolder);

        // Check if path exists
        if (not(stat(path_s.c_str(), &info) == 0 && S_ISDIR(info.st_mode))) {
            // Create path if it does not exist
            mkdir(path_s.c_str(), 0700);
        }

    } else {    // User has specified a path

        // This path maybe more than one layer deep, so attempt to make it recursively
        int len = path.length();

        // Check if prefix has been added
        if ((len < 4) or (not ((path[0] == 'd') and (path[1] == 'a') and (path[2] == 't') and (path[3] == 'a')))) {
            path_s += prefix;
            path_s += '/';
        }

        // Boolean to see name of first folder
        bool firstFolder = true;

        string folder = "";
        for (int i = 0; i < len; i++) {
            folder += path[i]; // Append char to folder name

            // If separator is found or if end of path is reached, construct folder
            if ((path[i] == '/') or (i == len - 1)) {

                // If separator is found, remove it:
                if (path[i] == '/') folder.pop_back();

                // Check if this is the first sub-folder
                if (firstFolder) {
                    firstFolder = false;

                    // Check if first folder contains date format
                    if (not ((folder.length() == 10) and(folder[4] == '-') and (folder[7] == '-'))) {

                        // Check if path exists
                        if (not(stat(path_s.c_str(), &info) == 0 && S_ISDIR(info.st_mode))) {
                            // Create path if it does not exist
                            mkdir(path_s.c_str(), 0700);
                        }
                    }
                }

                // Append folder to path
                path_s += folder;

                // Make folder
                if (not(stat(path_s.c_str(), &info) == 0 && S_ISDIR(info.st_mode)))
                { // Create path if it does not exist
                    mkdir(path_s.c_str(), 0700);
                }

                // Append a separator
                path_s += "/";

                // Reset folder
                folder = "";
            }
        }

        // Remove last separator
        if (path_s.back() == '/') path_s.pop_back();

    }

    // // Generate state folder  TODO: Change such that it works on the GPU
    // string path_ss = path_s + "/state";
    // if (not(stat(path_ss.c_str(), &info) == 0 && S_ISDIR(info.st_mode)))
    // { // Create path if it does not exist
    //     mkdir(path_ss.c_str(), 0700);
    // }

    // Return the generated path
    return path_s;
}


// Sets the folder number (useful when running parallel code)
void Chains::SetFolderNumber(int number) {path = to_string(number);};


// Sets the folder path (useful when running parallel code)
void Chains::SetPath(const std::string& path) {
    exportAny = true;
    this->path = path;
}


// Get properties ///////////////////////////////////////////////////////////////////////
// Returns the save path
std::string Chains::GetPath() {
    return path;
}


// Returns the time
int Chains::GetTime() {
    return Time;
}


// Returns the time-step dT
numtype Chains::GetDeltaT() {
    return dT;
}

// Clean up /////////////////////////////////////////////////////////////////////////////
// Delete the data folder
void Chains::DeleteFolder() {
    DeleteFolderTree(path.c_str());
}


// Delete folders recursively
void Chains::DeleteFolderTree(const char * directory_name) {

    DIR *            dp;
    struct dirent *  ep;
    char            p_buf[512] = {0};


    dp = opendir(directory_name);

    while ((ep = readdir(dp)) != NULL) {
        // Skip self dir "."
        if (strcmp(ep->d_name, ".") == 0 || strcmp(ep->d_name, "..") == 0) continue;

        sprintf(p_buf, "%s/%s", directory_name, ep->d_name);

        // Is the path a folder?
        struct stat s_buf;
        int IsDirectory = -1;
        if (stat(p_buf, &s_buf)){
            IsDirectory = 0;
        } else {
            IsDirectory = S_ISDIR(s_buf.st_mode);
        }

        // If it is a folder, go recursively into
        if (IsDirectory) {
            DeleteFolderTree(p_buf);
        } else {    // Else delete the file
            unlink(p_buf);
        }
    }

    closedir(dp);
    rmdir(directory_name);
}


// Destructor
Chains::~Chains() {

    // Close file-streams
    if (f_cells.is_open()) {
        f_cells.close();
    }
    if (f_colonySize.is_open()) {
        f_colonySize.close();
    }
    if (f_phages.is_open()) {
        f_phages.close();
    }
    if (f_log.is_open()) {
        f_log.close();
    }
}

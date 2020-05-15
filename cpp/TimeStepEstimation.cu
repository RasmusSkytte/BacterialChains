#include "Simulation.hpp"

using namespace std;

inline numtype cpu_pow(numtype x, numtype y){
	#if NUMTYPE_IS_FLOAT
	return powf(x,y);
	#else
	return pow(x,y);
	#endif
}

int main(int argc, char** argv){

	// Length of the chains
    int N = 10;
    numtype theta = 0.6*M_PIl;

    // Default parameters
    numtype dT_relax = 1e-4;     // Set the time scale during relaxation

    numtype T = 1;               // Set the length of the experiment
    std::vector<numtype> dT = {cpu_pow(10,-4.00), cpu_pow(10,-4.50), cpu_pow(10,-5.00), cpu_pow(10,-5.50), cpu_pow(10,-6.00), cpu_pow(10,-6.50), cpu_pow(10,-7.00), cpu_pow(10,-7.50)};

    int repeats = 3;

    // Streptobacillus
    numtype R  = 0.45;
    numtype Ld = 1.5;

    int k = 0;
    // Loop over time steps
    for (int t = 0; t < dT.size(); t++ ) {
        for (int r = 0; r < repeats; r++) {

            // Increment counter
            k++;

            // Print the process bar
            cout << "\r[" << k << " / " << dT.size() * repeats << "]" << flush;

            // Set the path
            string pathName = "TimeStepEstimation";

            pathName += "/dT_1e";
            char buffer[80];                                  // Create a buffer to store the date
            snprintf(buffer, sizeof(buffer), "%.2f", log10(dT[t]));
            pathName += string(buffer);

            pathName += "/repeat_";
            pathName += to_string(r);

            // Check if data is already made
            string path_s = "data/"; // Data folder name
            path_s += pathName;
            path_s += "/Completed.txt";

            // Check if run exists and is completed
            struct stat info;
            if (stat(path_s.c_str(), &info) == 0 && S_ISREG(info.st_mode)) {
                continue;
            }

            // Load simulation module
            Chains s(N);
            s.SetRngSeed(r);
            s.Debug(0);
            s.SetPath(pathName);

            s.SetSamples(100);
            s.SetLength(350);

            // Set the cell dimensions
            s.CellLength(Ld);
            s.CellRadius(R);

            // Set the angle between cells
            s.CellBendingAngle(theta);

            // Let the cells relax
            s.TimeStep(dT_relax);
            s.Relax();

            // Lock the configuration
            s.CellLock();

            // Configure the phage invasion
            s.PhageInvasionStartTime(0.0);
            s.PhageInitialDensity(-1e5);

            // Set the data to export
            s.ExportColonySize();

            // Set the time step
            s.TimeStep(dT[t]);

            // Run the experiment
            s.Run(T);

        }
    }

    cout << "\rDone!                   " << endl;
	return 0;
}

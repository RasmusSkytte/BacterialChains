#include "Simulation.hpp"

using namespace std;

int main(int argc, char** argv){

    // Length of the chains
    std::vector<int>     N    = {     1,     3,    10,    32,   100,   316,  1000,  3162};
    std::vector<int>     runs = {     9,     9,     9,     9,     9,     9,     9,     9};

    std::vector<numtype> theta   = {1.0, 0.8, 0.6, 0.4, 0.2};

    // Number of chains to generate for persistence length measurements
    int repeats = 100;

    // Default parameters
    numtype dT_relax = 1e-4;     // Set the time scale during relaxation

    numtype T  = 1;              // Set the length of the experiment
    numtype dT = 1e-7;           // Set the time-step of the experiment

    // Streptobacillus
    numtype R  = 0.45;
    numtype Ld = 1.5;

    int k = 0;
    // Loop over chain lengths
    for (int n = 0; n < N.size(); n++) {
        for (int t = 0; t < theta.size(); t++) {
            for (int r = 0; r < repeats; r++) {

                // Increment counter
                k++;

                // Print the process bar
                cout << "\r[" << k << " / " << theta.size() * N.size() * repeats << "]" << flush;

                // Set the path
                string pathName = "Chain";
                pathName += "/N_" + to_string(N[n]);

                pathName += "/theta_";
                char buffer[80]; // Create a buffer to store the date¨
                snprintf(buffer, sizeof(buffer), "%.3f", theta[t]);
                pathName += string(buffer);
                pathName += "_pi";

                pathName += "/repeat_" + to_string(r);

                // Check if data is already made
                string path_s1 = "data/" + pathName + "/Completed.txt";
                string path_s2 = "data/" + pathName + "/ColonySize.txt";

                // Check if runs exists and is completed
                bool exists = false;
                struct stat info;
                if ((stat(path_s1.c_str(), &info) == 0 && S_ISREG(info.st_mode)) and (stat(path_s2.c_str(), &info) == 0 && S_ISREG(info.st_mode))) { // Base run is completed

                    // Count sub runs that are completed
                    int completed = 0;

                    // Check if lysis runs are completed
                    if ((r == 0) and (N[n] > 1)) {
                        for (int i = 0; i < min(N[n], 5); i++) {
                            string test_path = "data/" + pathName + "/lysis_" + to_string(i) + "/Completed.txt";
                            if (stat(test_path.c_str(), &info) == 0 && S_ISREG(info.st_mode)) {
                                completed++;
                            }
                        }
                    }

                    // Check if any runs need to be completed
                    if ((r == 0) and (N[n] > 1)) {  // Lysis runs still need to be run
                        if (completed == min(N[n], 5)) continue;
                    } else {    // No runs are needed
                        continue;
                    }

                    // Mark the base run as existing
                    exists = true;
                }

                // Load simulation module
                Chains s(N[n]);
                s.SetRngSeed(r + 100 * t + 1000 * (2 * n));
                s.Debug(0);

                // Set the path
                if (not exists) s.SetPath(pathName);

                s.SetSamples(100);

                // Set the cell dimensions
                s.CellLength(Ld);
                s.CellRadius(R);

                // Set the angle between cells
                s.CellBendingAngle(theta[t] * M_PIl);

                // Let the cells relax
                s.TimeStep(dT_relax);
                s.Relax();

                // Check if the run is not a phage run
                if (r >= runs[n]) {
                    s.ExportCellDataNow();
                    continue;
                }

                // Set the data to export
                if (not exists) s.ExportColonySize();

                // Configure the phage invasion
                s.PhageInvasionStartTime(0.0);
                s.PhageInitialDensity(-1e5);

                // Lock the configuration
                s.CellLock();

                // Configure the time step
                s.TimeStep(dT);

                // Autoscale simulation
                s.AutoScale();

                // Run lysis experiment (for first repeat only)
                if ((N[n] > 1) && (r == 0)) {

                    // Select lysis site
                    for (int i = 0; i < min(N[n], 5); i++) {

                        // Create off-spring simulation
                        Chains t(s);

                        // Increase the sampling
                        t.SetSamples(10000);

                        // Configure phage attack
                        t.PhageInitialDensity(0);
                        t.PhageBurstSize(1e4);

                        if (i == 0) {
                            t.LyseCell(0);
                        } else {
                            t.LyseCell(i * (N[n] - 1) / (min(N[n], 5) - 1));
                        }

                        // Change the pathname
                        string path_s = pathName + "/lysis_" + to_string(i);
                        t.SetPath(path_s);

                        // Check if data is already made
                        path_s = "data/" + path_s + "/Completed.txt";

                        // Check if run exists and is completed
                        struct stat info;
                        if (stat(path_s.c_str(), &info) == 0 && S_ISREG(info.st_mode)) continue;

                        // Set the data to export
                        t.ExportColonySize();

                        // Run the experiment
                        t.Run(T);

                        // Store the final configuration
                        t.ExportCellDataNow();
                        t.ExportPhageDataNow();

                    }
                }

                // Skip the run if already completed
                if (exists) continue;

                // Run the Experiment
                s.Run(T);

                // Store the final configuration
                s.ExportCellDataNow();
                s.ExportPhageDataNow();
            }
        }
    }

    cout << "\rDone!                   " << endl;
	return 0;
}

#include "Simulation.hpp"

using namespace std;

int main(int argc, char** argv){

    // Length of the chains
    std::vector<int>     N       = {     1,     3,    10,    32,   100,   316,  1000};
    std::vector<int>     repeats = {     3,     3,     3,     3,     3,     3,     3};

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
        for (int r = 0; r < repeats[n]; r++) {

            // Increment counter
            k++;

            // Print the process bar
            cout << "\r[" << k << " / " << std::accumulate(repeats.begin(), repeats.end(), 0) << "]" << flush;

            // Set the path
            string pathName = "WellMixed";
            pathName += "/N_" + to_string(N[n]);
            pathName += "/repeat_" + to_string(r);

            // Check if data is already made
            string path_s = "data/" + pathName + "/Completed.txt";

            // Check if run exists and is completed
            bool exists = false;
            struct stat info;
            if (stat(path_s.c_str(), &info) == 0 && S_ISREG(info.st_mode)) {
                if ((r == 0) and (N[n] > 1)) {
                    int completed = 0;
                    for (int i = 0; i < min(N[n], 5); i++) {
                        string test_path = "data/" + pathName + "/lysis_" + to_string(i) + "/Completed.txt";
                        if (stat(test_path.c_str(), &info) == 0 && S_ISREG(info.st_mode)) {
                            completed++;
                        }
                    }
                    if (completed == min(N[n], 5)) continue;
                } else {
                    continue;
                }
                exists = true;
            }

            // Load simulation module
            Chains s(N[n]);
            s.SetRngSeed(r +  800 * (2 * n));
            s.Debug(0);

            // Set the path
            if (not exists) s.SetPath(pathName);

            s.SetSamples(100);

            // Set the cell dimensions
            s.CellLength(Ld);
            s.CellRadius(R);

            // Set the cells to be dispersed
            s.WellMixed();
            s.SetLength(100);

            // Let the cells relax
            s.TimeStep(dT_relax);
            s.Relax();

            // Configure the phage invasion
            s.PhageInvasionStartTime(0.0);
            s.PhageInitialDensity(-1e5);

            // Set the data to export
            if (not exists) s.ExportColonySize();

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

            // Run the experiment
            s.Run(T);

            // Store the final configuration
            s.ExportCellDataNow();
            s.ExportPhageDataNow();
        }
    }

    cout << "\rDone!                   " << endl;
	return 0;
}

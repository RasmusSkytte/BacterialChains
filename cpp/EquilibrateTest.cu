#include "Simulation.hpp"

using namespace std;

int main(int argc, char** argv){

    // Length of the chains
    std::vector<int> N = {1, 32};
    int repeats = 3;
    numtype theta = 0.6;

    // Default parameters
    numtype dT_relax = 1e-4;     // Set the time scale during relaxation

    numtype T = 1;               // Set the length of the experiment
    numtype dT = 1e-7;           // Set the time-step of the experiment

    std::vector<numtype> T_Eq = { 0, 1, 2, 3, 4, 5, 6 };

    // Streptobacillus
    numtype R  = 0.45;
    numtype Ld = 1.5;

    int k = 0;
    // Loop over repeats
    for (int n = 0; n < N.size(); n++) {
        for (int r = 0; r < repeats; r++) {

            // Increment counter
            k++;

            // Print the process bar
            cout << "\r[" << k << " / " << T_Eq.size() * N.size() * repeats << "]" << flush;

            // Set the path
            string pathName = "EquilibriumTest";

            // Check if run exists and is completed
            struct stat info;
            char buffer[80]; // Create a buffer
            int completed = 0;
            for (int m = 0; m < T_Eq.size(); m++) {
                snprintf(buffer, sizeof(buffer), "%.0f", T_Eq[m]);
                string test_path = "data/" + pathName + "/T_" + string(buffer) + "/N_" + to_string(N[n]) + "/repeat_" + to_string(r) + "/Completed.txt";
                if (stat(test_path.c_str(), &info) == 0 && S_ISREG(info.st_mode)) {
                    completed++;
                }
            }
            if (completed == T_Eq.size()) continue;


            // Load simulation module
            Chains s(N[n]);
            s.SetRngSeed(r);
            s.Debug(0);

            s.SetSamples(100);

            // Set the cell dimensions
            s.CellLength(Ld);
            s.CellRadius(R);

            // Set the angle between cells
            s.CellBendingAngle(theta * M_PIl);

            // Let the cells relax
            s.TimeStep(dT_relax);
            s.Relax();

            // Configure the phage invasion
            s.PhageInvasionStartTime(0.0);
            s.PhageInitialDensity(-1e5);

            // Lock the configuration
            s.CellLock();

            // Configure the time step
            s.TimeStep(dT);

            // Autoscale simulation
            s.AutoScale();

            // Create off-spring simulation
            Chains t(s);

            // Run without equilibration time
            snprintf(buffer, sizeof(buffer), "%.0f", T_Eq[0]);
            string path = pathName + "/T_" + string(buffer) + "/N_" + to_string(N[n]) + "/repeat_" + to_string(r);

            // Check if run is completed
            string path_s = "data/" + path + "/Completed.txt";

            // Check if run exists and is completed
            if (not (stat(path_s.c_str(), &info) == 0 && S_ISREG(info.st_mode))) {

                // Set the path name
                s.SetPath(path);

                // Set the data to export
                s.ExportColonySize();

                // Run the experiment
                s.Run(T);

                // Store the final configuration
                s.ExportCellDataNow();
                s.ExportPhageDataNow();

            }

            // Loop over equilibration times
            for (int m = 1; m < T_Eq.size(); m++) {

                // Increment counter
                k++;

                // Print the process bar
                cout << "\r[" << k << " / " << T_Eq.size() * N.size() * repeats << "]" << flush;

                // Let the master simulation equilibrate for T time
                t.Equilibrate(T_Eq[m] - T_Eq[m-1]);

                // Make copy
                Chains u(t);

                // Set new path
                snprintf(buffer, sizeof(buffer), "%.0f", T_Eq[m]);
                path = pathName + "/T_" + string(buffer) + "/N_" + to_string(N[n]) + "/repeat_" + to_string(r);

                // Check if run is completed
                string path_s = "data/" + path + "/Completed.txt";

                // Check if run exists and is completed
                if (not (stat(path_s.c_str(), &info) == 0 && S_ISREG(info.st_mode))) {

                    // Set the path name
                    u.SetPath(path);

                    // Set the data to export
                    u.ExportColonySize();

                    // Run the experiment
                    u.Run(T);

                    // Store the final configuration
                    u.ExportCellDataNow();
                    u.ExportPhageDataNow();

                }
            }
        }
    }

    cout << "\rDone!                   " << endl;
	return 0;
}

#include "Simulation.hpp"

using namespace std;

int main(int argc, char** argv){

	// Length of the chains
    int N = 100;

    // Default parameters
    numtype dT_relax = 1e-4;     // Set the time scale during relaxation

    std::vector<numtype> theta  = {1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0};

    int repeats = 9;

    // Streptobacillus
    numtype R  = 0.45;
    numtype Ld = 1.5;

    // Loop over time steps
    for (int t = 0; t < theta.size(); t++) {
        for (int r = 0; r < repeats; r++) {

            // Set the path
            string pathName = "BendingAngleValidation/";
            pathName += "theta_";
            char buffer[80]; // Create a buffer to store the date¨
            snprintf(buffer, sizeof(buffer), "%.3f", theta[t]);
            pathName += string(buffer);
            pathName += "_pi";

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

            // Print the process bar
            cout << "\r[" << t * repeats + r + 1 << " / " << theta.size() * repeats << "]" << flush;

            // Load simulation module
            Chains s(N);
            s.SetRngSeed(r);
            s.Debug(0);
            s.SetPath(pathName);

            // Set the cell dimensions
            s.CellLength(Ld);
            s.CellRadius(R);

            // Set the angle between cells
            s.CellBendingAngle(theta[t] * M_PIl);

            // Configure the phage invasion
            s.PhageInitialDensity(0);

            // Let the cells relax
            s.TimeStep(dT_relax);
            s.Relax();

            // Export final configuration
            s.ExportCellDataNow();

        }
    }

    cout << "\rDone!                   " << endl;
	return 0;
}

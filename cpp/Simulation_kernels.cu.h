const int TILE_DIM = 100;

inline __device__ numtype gpu_pow(numtype x, numtype y){
    #if NUMTYPE_IS_FLOAT
    return powf(x,y);
    #else
    return pow(x,y);
    #endif
}

inline __device__ numtype gpu_sqrt(numtype x){
    #if NUMTYPE_IS_FLOAT
    return sqrtf(x);
    #else
    return sqrt(x);
    #endif
}

inline __device__ numtype gpu_max(numtype x, numtype y){
    #if NUMTYPE_IS_FLOAT
    return fmaxf(x, y);
    #else
    return fmax(x, y);
    #endif
}

inline __device__ numtype gpu_min(numtype x, numtype y){
    #if NUMTYPE_IS_FLOAT
    return fminf(x, y);
    #else
    return fmin(x, y);
    #endif
}

inline __device__ numtype gpu_randn(curandState* rng){
    #if NUMTYPE_IS_FLOAT
    return curand_normal(rng);
    #else
    return curand_normal_double(rng);
    #endif
}

inline __device__ numtype gpu_rand(curandState* rng){
    #if NUMTYPE_IS_FLOAT
    return curand_uniform(rng);
    #else
    return curand_uniform_double(rng);
    #endif
}

__global__ void initRNG(curandState* state, int N){

  int i = blockDim.x*blockIdx.x + threadIdx.x;

  if (i < N) {
    curand_init(i, 0, 0, &state[i]);
  }
}

__device__ void Shortest_Distance_Vector(numtype* u, numtype* xU, numtype* yU, numtype* zU, numtype xP, numtype yP, numtype zP, numtype xQ, numtype yQ, numtype zQ, numtype xR, numtype yR, numtype zR, numtype xS, numtype yS, numtype zS ) {

    // Vector between Q and P
    numtype xQP = xQ - xP;
    numtype yQP = yQ - yP;
    numtype zQP = zQ - zP;

    // Vector between S and R
    numtype xSR = xS - xR;
    numtype ySR = yS - yR;
    numtype zSR = zS - zR;

    // Vector between P and R
    numtype xPR = xP - xR;
    numtype yPR = yP - yR;
    numtype zPR = zP - zR;

    // Compute the relevant dot products
    numtype dQP_QP = xQP*xQP + yQP*yQP + zQP*zQP;
    numtype dSR_SR = xSR*xSR + ySR*ySR + zSR*zSR;
    numtype dQP_SR = xQP*xSR + yQP*ySR + zQP*zSR;
    numtype dQP_PR = xQP*xPR + yQP*yPR + zQP*zPR;
    numtype dSR_PR = xSR*xPR + ySR*yPR + zSR*zPR;

    // Compute denominator
    numtype h = dQP_QP*dSR_SR - dQP_SR*dQP_SR;
    numtype hu = h;
    numtype hv = h;

    // Compute u and v
           *u = dQP_SR*dSR_PR - dSR_SR*dQP_PR;
    numtype v = dQP_QP*dSR_PR - dQP_SR*dQP_PR;

    // If u or v lands outside segment, convert back to line segment
    if (*u < 0.0) {
        *u = 0.0;
        v  = dSR_PR;
        hv = dSR_SR;
    } else if (*u > h) {
        *u = h;
        v  = dSR_PR + dQP_SR;
        hv = dSR_SR;
    }

    if (v < 0.0) {
        v = 0.0;

        if (dQP_PR > 0.0) {
            *u = 0.0;
        } else if (-dQP_PR > dQP_QP) {
            *u = h;
        } else {
            *u = -dQP_PR;
            hu = dQP_QP;
        }
    } else if (v > hv) {
        v = hv;
        if ((dQP_SR - dQP_PR) < 0.0) {
            *u = 0;
        } else if ((dQP_SR - dQP_PR) > dQP_QP) {
            *u = h;
        } else {
            *u = dQP_SR - dQP_PR;
            hu = dQP_QP;
        }
    }

    // Compute u and v
    *u = *u / hu;
    v  =  v / hv;

    // Compute the shortest distance vector
    *xU = xPR + *u * xQP - v * xSR;
    *yU = yPR + *u * yQP - v  *ySR;
    *zU = zPR + *u * zQP - v * zSR;
}

__global__ void PhageUpdateKernel(numtype* phages, int* active, numtype* cells, numtype sigma, numtype L, numtype R, curandState* rng, int M, int N, bool redistribute){

    int i = blockDim.x*blockIdx.x + threadIdx.x;

    if ((i < M) && (active[i]==1)) {

        // Store phage coordinates
        numtype xR = phages[3*i];
        numtype yR = phages[3*i+1];
        numtype zR = phages[3*i+2];

        // Add gaussian noise to the coordinates
        numtype xS = xR + sigma*gpu_randn(&rng[i]);
        numtype yS = yR + sigma*gpu_randn(&rng[i]);
        numtype zS = zR + sigma*gpu_randn(&rng[i]);

        // Loop over cells to check for overlap
        #pragma unroll
        for (int n = 0; n < N; n++) {

            // Store cell coordinates (P and Q vector)
            numtype xP = cells[8*n];
            numtype yP = cells[8*n+1];
            numtype zP = cells[8*n+2];
            numtype xQ = cells[8*n+3];
            numtype yQ = cells[8*n+4];
            numtype zQ = cells[8*n+5];

            // Vector between Q and P
            numtype xQP = xQ - xP;
            numtype yQP = yQ - yP;
            numtype zQP = zQ - zP;

            // Vector between S and P
            numtype xSP = xS - xP;
            numtype ySP = yS - yP;
            numtype zSP = zS - zP;

            // Compute the relevant dot products
            numtype dSP_QP = xSP*xQP + ySP*yQP + zSP*zQP;
            numtype dQP_QP = xQP*xQP + yQP*yQP + zQP*zQP;

            // Determine u
            numtype u;
            if ( dSP_QP <= 0 ) { // Closest point is P
                u = 0;
            } else if ( dQP_QP <= dSP_QP ) { // Closest point is Q
                u = 1;
            } else { // Closest point is on PQ
                u = dSP_QP / dQP_QP;
            }

            // Calculate the shortest distance vector
            numtype xU = xSP - u*xQP;
            numtype yU = ySP - u*yQP;
            numtype zU = zSP - u*zQP;

            // Calculate the distance
            numtype d = gpu_sqrt( xU*xU + yU*yU + zU*zU );

            // Check for overlap
            if (d < R) {

                if (not redistribute) {

                    // Set inactive
                    active[i] = 0;

                } else {

                    // Move to new location
                    xS = L * gpu_rand(&rng[i]);
                    yS = L * gpu_rand(&rng[i]);
                    zS = L * gpu_rand(&rng[i]);

                }

            }
        }

        // Apply boundary conditions
        #if PERIODIC_BOUNDARY_CONDITIONS
        // Check X-direction
        numtype s;
        s = (xS > 0) ? 1 : -1;
        if (abs(xS-L/2) >  L/2) {
            xS = (abs(xS) - L) * s;
        }

        // Check Y-direction
        s = (yS > 0) ? 1 : -1;
        if (abs(yS-L/2) >  L/2) {
            yS = (abs(yS) - L) * s;
        }

        // Check Z-direction
        s = (zS > 0) ? 1 : -1;
        if (abs(zS-L/2) >  L/2) {
            zS = (abs(zS) - L) * s;
        }

        #else
        // Check X-direction
        if (xS < 0)        xS = -xS;
        else if (xS > L)   xS = 2*L - xS;

        // Check Y-direction
        if (yS < 0)        yS = -yS;
        else if (yS > L)   yS = 2*L - yS;

        // Check Z-direction
        if (zS < 0)        zS = -zS;
        else if (zS > L)   zS = 2*L - zS;

        #endif

        // Store the new coordinates
        phages[3*i]     = xS;
        phages[3*i+1]   = yS;
        phages[3*i+2]   = zS;
    }
}

__global__ void CellUpdateKernel(numtype* cells, numtype* cells_new, numtype len, numtype rad, numtype k_rep, numtype k_att, numtype k_int, numtype k_pull, numtype dT, numtype L, int N){

    int i = blockDim.x*blockIdx.x + threadIdx.x;

    if (i < N) {

        // Store the cell coordinates
        numtype xP = cells[8*i];
        numtype yP = cells[8*i+1];
        numtype zP = cells[8*i+2];

        numtype xQ = cells[8*i+3];
        numtype yQ = cells[8*i+4];
        numtype zQ = cells[8*i+5];

        // Define QP vector
        numtype xQP = xQ - xP;
        numtype yQP = yQ - yP;
        numtype zQP = zQ - zP;

        // Compute cell length
        numtype l = gpu_sqrt( xQP * xQP + yQP * yQP + zQP * zQP );

        // Define the gradients of the potential:
        numtype xdV_dP = 0.0;
        numtype ydV_dP = 0.0;
        numtype zdV_dP = 0.0;

        numtype xdV_dQ = 0.0;
        numtype ydV_dQ = 0.0;
        numtype zdV_dQ = 0.0;

        // Compute internal spring forces /////////////////////////////////////////
        // V_i = 0.5 * k_int * (norm(QP)-len)^2
        xdV_dP -=  k_int * (l - len) * (xQP / l);
        ydV_dP -=  k_int * (l - len) * (yQP / l);
        zdV_dP -=  k_int * (l - len) * (zQP / l);
        xdV_dQ +=  k_int * (l - len) * (xQP / l);
        ydV_dQ +=  k_int * (l - len) * (yQP / l);
        zdV_dQ +=  k_int * (l - len) * (zQP / l);

        // Compute repulsion between neighbors ///////////////////////////////////
        #pragma unroll
        for (int n = 0; n < N; n++) {

            // Skip the current cell ( cell I )
            if ( n != i ) {

                // Define end points of the neighbor cell
                numtype xR = cells[8*n];
                numtype yR = cells[8*n+1];
                numtype zR = cells[8*n+2];

                numtype xS = cells[8*n+3];
                numtype yS = cells[8*n+4];
                numtype zS = cells[8*n+5];

                numtype xU, yU, zU, u;
                Shortest_Distance_Vector(&u, &xU, &yU, &zU, xP, yP, zP, xQ, yQ, zQ, xR, yR, zR, xS, yS, zS );

                // Calculate the distance
                numtype d = gpu_sqrt( xU * xU + yU * yU + zU * zU );

                // Check for overlap between the cells
                if ( d < 2 * rad ) {

                    // Compute gradients on the cell
                    xdV_dP += k_rep * (d - 2 * rad) * ( xU / d ) * ( 1 - u);
                    ydV_dP += k_rep * (d - 2 * rad) * ( yU / d ) * ( 1 - u);
                    zdV_dP += k_rep * (d - 2 * rad) * ( zU / d ) * ( 1 - u);
                    xdV_dQ += k_rep * (d - 2 * rad) * ( xU / d ) * u;
                    ydV_dQ += k_rep * (d - 2 * rad) * ( yU / d ) * u;
                    zdV_dQ += k_rep * (d - 2 * rad) * ( zU / d ) * u;

                }
            }
        }

        // Calculate attraction between connected poles  //////////////////////////
        if (k_att > 0) {

            // Get the cell number of the first connected pole
            numtype xC_P, yC_P, zC_P;
            int ID = cells[8*i+6];
            int J = abs(ID) - 1;

            // Check if selected pole has neighbor
            if (J != -1) {

                // Get coordinates of connected point
                if (ID > 0) {
                    xC_P = cells[8*J];
                    yC_P = cells[8*J+1];
                    zC_P = cells[8*J+2];
                } else {
                    xC_P = cells[8*J+3];
                    yC_P = cells[8*J+4];
                    zC_P = cells[8*J+5];
                }
            } else {
                // Set "connected" point to be itself
                xC_P = xP;
                yC_P = yP;
                zC_P = zP;
            }

            // Get the cell number of second connected pole
            numtype xC_Q, yC_Q, zC_Q;

            ID = cells[8*i+7];
            J = abs(ID) - 1;

            // Check if selected pole has neighbor
            if (J != -1) {

                // Get coordinates of connected point
                if (ID > 0) {
                    xC_Q = cells[8*J];
                    yC_Q = cells[8*J+1];
                    zC_Q = cells[8*J+2];
                } else {
                    xC_Q = cells[8*J+3];
                    yC_Q = cells[8*J+4];
                    zC_Q = cells[8*J+5];
                }
            } else {
                // Set "connected" point to be itself
                xC_Q = xQ;
                yC_Q = yQ;
                zC_Q = zQ;
            }

            // Define the distance vectors D
            numtype xD_P = xP - xC_P;
            numtype yD_P = yP - yC_P;
            numtype zD_P = zP - zC_P;

            numtype xD_Q = xQ - xC_Q;
            numtype yD_Q = yQ - yC_Q;
            numtype zD_Q = zQ - zC_Q;

            // Compute the distances to the connected points
            numtype d_P = gpu_sqrt( xD_P * xD_P + yD_P * yD_P + zD_P * zD_P );
            numtype d_Q = gpu_sqrt( xD_Q * xD_Q + yD_Q * yD_Q + zD_Q * zD_Q );

            // Compute the derivatives of the spring potential
            // V = 0.5 * k_att * (norm(P-C_P) - 2*R)^2 + 0.5 * k_att * (norm(Q-C_Q) - 2*R)^2
            if (d_P > 2 * rad) {
                xdV_dP += k_att * (d_P - 2 * rad) * xD_P / d_P;
                ydV_dP += k_att * (d_P - 2 * rad) * yD_P / d_P;
                zdV_dP += k_att * (d_P - 2 * rad) * zD_P / d_P;
            }
            if (d_Q > 2 * rad) {
                xdV_dQ += k_att * (d_Q - 2 * rad) * xD_Q / d_Q;
                ydV_dQ += k_att * (d_Q - 2 * rad) * yD_Q / d_Q;
                zdV_dQ += k_att * (d_Q - 2 * rad) * zD_Q / d_Q;
            }

        }

        // Calculate attraction to the center  ////////////////////////////////////
        if ((k_pull > 0) and (N > 1)) {

            // Define center of mass
            numtype xCOM = (xQ + xP) / 2.0;
            numtype yCOM = (yQ + yP) / 2.0;
            numtype zCOM = (zQ + zP) / 2.0;

            // Compute force pulling to center
            numtype xF = xCOM - L/2;
            numtype yF = yCOM - L/2;
            numtype zF = zCOM - L/2;

            numtype f = gpu_sqrt( xF * xF + yF * yF + zF * zF );

            // Move points towards centre
            if (f > 0.1) {
                xdV_dP += k_pull * xF / f;
                ydV_dP += k_pull * yF / f;
                zdV_dP += k_pull * zF / f;

                xdV_dQ += k_pull * xF / f;
                ydV_dQ += k_pull * yF / f;
                zdV_dQ += k_pull * zF / f;
            }
        }

        // Update positions based on potentials
        xP -= xdV_dP * dT;
        yP -= ydV_dP * dT;
        zP -= zdV_dP * dT;
        xQ -= xdV_dQ * dT;
        yQ -= ydV_dQ * dT;
        zQ -= zdV_dQ * dT;


        // Update the values of the cell
        cells_new[8*i]   = xP;
        cells_new[8*i+1] = yP;
        cells_new[8*i+2] = zP;
        cells_new[8*i+3] = xQ;
        cells_new[8*i+4] = yQ;
        cells_new[8*i+5] = zQ;

    }
}

__global__ void CellOverlaps(numtype* cells, numtype* overlaps, numtype rad, int N){


    int i = blockDim.x*blockIdx.x + threadIdx.x;

    if (i < N - 1) {

        // Store the cell coordinates
        numtype xP = cells[8*i];
        numtype yP = cells[8*i+1];
        numtype zP = cells[8*i+2];

        numtype xQ = cells[8*i+3];
        numtype yQ = cells[8*i+4];
        numtype zQ = cells[8*i+5];

        // Largest overlap detected
        numtype overlap = 0;

        // Compute overlap with neighbors ///////////////////////////////////
        for (int n = i + 1; n < N; n++) {

            // Define end points of the neighbor cell
            numtype xR = cells[8*n];
            numtype yR = cells[8*n+1];
            numtype zR = cells[8*n+2];

            numtype xS = cells[8*n+3];
            numtype yS = cells[8*n+4];
            numtype zS = cells[8*n+5];

            numtype xU, yU, zU, u;
            Shortest_Distance_Vector(&u, &xU, &yU, &zU, xP, yP, zP, xQ, yQ, zQ, xR, yR, zR, xS, yS, zS );

            // Calculate the distance
            numtype d = gpu_sqrt( xU * xU + yU * yU + zU * zU );

            // Check for overlap between the cells
            if ( 2 * rad - d > overlap ) {

                // Store overlap
                overlap =  2 * rad - d;

            }
        }

        // Write overlap to vector
        overlaps[i] = overlap;
    }
}

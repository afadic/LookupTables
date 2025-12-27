/**
 * kraehnert_ammonia_oxidation.c
 * 
 * Kraehnert DAE Model for ammonia oxidation over Pt.
 * C implementation using SUNDIALS IDA 7.x
 * 
 * Compile with:
 * gcc -O3 -march=native -o krahnert_model krahnert.c \
 *   -L/usr/local/sundials/lib -lsundials_ida -lsundials_nvecserial \
 *   -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_core -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_context.h>

#include <time.h>  // For clock_gettime
#include <stdint.h> // For int64_t

#define COMPILER_BARRIER() asm volatile("" ::: "memory")

#define NEQ 6  // Number of equations

// Data structure containing all model parameters
typedef struct {
    double T, P, pNH3, pNO, pO2;
    double R, surface_density, relaxation_factor, T_ref;
    double T_factor, surface_relax_inv;
} KrahnertData;

// Activation energies (kJ/mol)
static const double E[10] = {
    0.0, 60.9, 0.0, 181.0, 99.5, 154.8, 63.5, 139.0, 135.4, 155.2
};

// Pre-exponential factors
static const double k0[10] = {
    6.38e-1, 2.17e0, 2.94e-1, 1.09e-10, 5.91e2,
    1.24e0, 2.63e-1, 6.42e1, 9.34e0, 5.2e0
};

// Forward declarations
static int residual(sunrealtype t, N_Vector y_vec, N_Vector yp_vec, 
                    N_Vector res_vec, void *user_data);

/**
 * Compute rate constants k1-k10 based on Arrhenius expression
 */
static void compute_rate_constants(KrahnertData *data, double *k) {
    double T_factor = data->T_factor;
    double factor = data->surface_relax_inv;
    
    k[0] = k0[0] * exp(-E[0] * T_factor) * factor;
    k[1] = k0[1] * exp(-E[1] * T_factor) * factor;
    k[2] = k0[2] * exp(-E[2] * T_factor) * factor;
    k[3] = k0[3] * exp(-E[3] * T_factor) * factor;
    k[4] = k0[4] * exp(-E[4] * T_factor) * factor;
    k[5] = k0[5] * exp(-E[5] * T_factor) * factor;
    k[6] = k0[6] * exp(-E[6] * T_factor) * factor;
    k[7] = k0[7] * exp(-E[7] * T_factor) * factor;
    k[8] = k0[8] * exp(-E[8] * T_factor) * factor;
    k[9] = k0[9] * exp(-E[9] * T_factor) * factor;
}

/**
 * DAE residual function: res = f(t, y, y')
 */
static int residual(sunrealtype t, N_Vector y_vec, N_Vector yp_vec, 
                    N_Vector res_vec, void *user_data) {
    KrahnertData *data = (KrahnertData*)user_data;
    sunrealtype *y = NV_DATA_S(y_vec);
    sunrealtype *yp = NV_DATA_S(yp_vec);
    sunrealtype *res = NV_DATA_S(res_vec);
    
    double k[10];
    compute_rate_constants(data, k);
    
    double R1 = k[0] * data->pNH3 * y[0];
    double R2 = k[1] * y[1];
    double R3 = k[2] * data->pO2 * y[2] * y[2];
    double R4 = k[3] * y[3] * y[3];
    double R5 = k[4] * y[1] * y[3];
    double R6 = k[5] * y[4];
    double R7 = k[6] * data->pNO * y[2];
    double R8 = k[7] * y[5] * y[5];
    double R9 = k[8] * y[5] * y[3];
    double R10 = k[9] * y[4] * y[5];
    
    res[0] = y[0] + y[1] - 1.0;
    res[1] = yp[1] - (R1 - R2 - R5);
    res[2] = y[2] + y[3] + y[4] + y[5] - 1.0;
    res[3] = yp[3] - (2.0*R3 - 2.0*R4 - 1.5*R5 - R9);
    res[4] = yp[4] - (-R6 + R7 + R9 - R10);
    res[5] = yp[5] - (R5 - 2.0*R8 - R9 - R10);
    
    return 0;
}

/**
 * Compute reaction rates for output
 */
static void compute_reaction_rates(KrahnertData *data, double *y, double *rates) {
    double k[10];
    compute_rate_constants(data, k);
    
    // Scale back for rate calculation
    for (int i = 0; i < 10; i++) k[i] *= data->surface_density * data->relaxation_factor;
    
    double R1 = k[0] * data->pNH3 * y[0];
    double R2 = k[1] * y[1];
    double R3 = k[2] * data->pO2 * y[2] * y[2];
    double R4 = k[3] * y[3] * y[3];
    double R5 = k[4] * y[1] * y[3];
    double R6 = k[5] * y[4];
    double R7 = k[6] * data->pNO * y[2];
    double R8 = k[7] * y[5] * y[5];
    double R9 = k[8] * y[5] * y[3];
    double R10 = k[9] * y[4] * y[5];
    
    rates[0] = (-R1 + R2);
    rates[1] = (R8) ;
    rates[2] = (R6 - R7) ;
    rates[3] = (-R3 + R4) ;
    rates[4] = (1.5 * R5) ;
    rates[5] = (R10) ;
    
    (void)R9; 
}

static int jacobian(sunrealtype tt, sunrealtype cj, N_Vector yy, N_Vector yp, 
                    N_Vector /* not used */, SUNMatrix JJ, void *user_data,
                    N_Vector /* not used */, N_Vector /* not used */) {
    KrahnertData *data = (KrahnertData*)user_data;
    sunrealtype *y = NV_DATA_S(yy);
    
    double k[10];
    compute_rate_constants(data, k);
    
    // Zero the Jacobian matrix first
    SUNMatZero(JJ);
    
    // Row 0: res[0] = y[0] + y[1] - 1.0
    SM_ELEMENT_D(JJ, 0, 0) = 1.0;
    SM_ELEMENT_D(JJ, 0, 1) = 1.0;
    
    // Row 1: res[1] = yp[1] - (R1 - R2 - R5) = yp[1] - R1 + R2 + R5
    SM_ELEMENT_D(JJ, 1, 0) = -k[0] * data->pNH3;
    SM_ELEMENT_D(JJ, 1, 1) = k[1] + k[4] * y[3] + cj;
    SM_ELEMENT_D(JJ, 1, 3) = k[4] * y[1];
    
    // Row 2: res[2] = y[2] + y[3] + y[4] + y[5] - 1.0
    SM_ELEMENT_D(JJ, 2, 2) = 1.0;
    SM_ELEMENT_D(JJ, 2, 3) = 1.0;
    SM_ELEMENT_D(JJ, 2, 4) = 1.0;
    SM_ELEMENT_D(JJ, 2, 5) = 1.0;
    
    // Row 3: res[3] = yp[3] - (2*R3 - 2*R4 - 1.5*R5 - R9) = yp[3] - 2*R3 + 2*R4 + 1.5*R5 + R9
    SM_ELEMENT_D(JJ, 3, 1) = 1.5 * k[4] * y[3];
    SM_ELEMENT_D(JJ, 3, 2) = -4.0 * k[2] * data->pO2 * y[2];  // -2 * (2*k[2]*pO2*y[2])
    SM_ELEMENT_D(JJ, 3, 3) = 6.0 * k[3] * y[3] * y[3] + 1.5 * k[4] * y[1] + k[8] * y[5] + cj;  // 2*(3*k[3]*y[3]^2)
    SM_ELEMENT_D(JJ, 3, 5) = k[8] * y[3];
    
    // Row 4: res[4] = yp[4] - (-R6 + R7 + R9 - R10) = yp[4] + R6 - R7 - R9 + R10
    SM_ELEMENT_D(JJ, 4, 2) = -k[6] * data->pNO;
    SM_ELEMENT_D(JJ, 4, 3) = -k[8] * y[5];
    SM_ELEMENT_D(JJ, 4, 4) = k[5] + k[9] * y[5] + cj;
    SM_ELEMENT_D(JJ, 4, 5) = -k[8] * y[3] + k[9] * y[4];
    
    // Row 5: res[5] = yp[5] - (R5 - 2*R8 - R9 - R10) = yp[5] - R5 + 2*R8 + R9 + R10
    SM_ELEMENT_D(JJ, 5, 1) = -k[4] * y[3];
    SM_ELEMENT_D(JJ, 5, 3) = -k[4] * y[1] + k[8] * y[5];
    SM_ELEMENT_D(JJ, 5, 4) = k[9] * y[5];
    SM_ELEMENT_D(JJ, 5, 5) = 4.0 * k[7] * y[5] + k[8] * y[3] + k[9] * y[4] + cj;  // 2*(2*k[7]*y[5])
    
    return 0;
}

static int solve_krahnert(KrahnertData *data, double *y_final, double t_final) {
    void *mem = NULL;
    N_Vector y = NULL, yp = NULL, constraints = NULL, id = NULL;
    SUNMatrix J = NULL;
    SUNLinearSolver LS = NULL;
    SUNContext sunctx = NULL;
    int flag;
    
    // Create SUNDIALS context
    flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
    if (flag != 0) return -1;
    
    // Initialize vectors
    y = N_VNew_Serial(NEQ, sunctx);
    yp = N_VNew_Serial(NEQ, sunctx);
    if (!y || !yp) {
        N_VDestroy(y); N_VDestroy(yp);
        SUNContext_Free(&sunctx);
        return -1;
    }
    
    sunrealtype *y_data = NV_DATA_S(y);
    sunrealtype *yp_data = NV_DATA_S(yp);
    
    // Initial conditions (algebraic variables)
    y_data[0] = 1.0; y_data[1] = 0.0; y_data[2] = 1.0;
    y_data[3] = 0.0; y_data[4] = 0.0; y_data[5] = 0.0;
    
    // Initialize yp to zero (IDA will correct this)
    for (int i = 0; i < NEQ; i++) yp_data[i] = 0.0;
    
    // Create IDA memory
    mem = IDACreate(sunctx);
    if (mem == NULL) {
        N_VDestroy(y); N_VDestroy(yp);
        SUNContext_Free(&sunctx);
        return -1;
    }
    
    // Initialize IDA
    flag = IDAInit(mem, residual, 0.0, y, yp);
    if (flag != IDA_SUCCESS) {
        IDAFree(&mem); N_VDestroy(y); N_VDestroy(yp);
        SUNContext_Free(&sunctx);
        return -1;
    }
    
    // Set user data and tolerances
    IDASetUserData(mem, data);
    IDASStolerances(mem, 1e-4, 1e-6);
    
    // Set ID vector (1:differential vs 0:algebraic)
    id = N_VNew_Serial(NEQ, sunctx);
    if (id == NULL) {
        IDAFree(&mem); N_VDestroy(y); N_VDestroy(yp);
        SUNContext_Free(&sunctx);
        return -1;
    }
    sunrealtype *id_data = NV_DATA_S(id);
    id_data[0] = 0.0; id_data[1] = 1.0; id_data[2] = 0.0;
    id_data[3] = 1.0; id_data[4] = 1.0; id_data[5] = 1.0;
    IDASetId(mem, id);
    
    // Create matrix and linear solver
    J = SUNDenseMatrix(NEQ, NEQ, sunctx);
    LS = SUNLinSol_Dense(y, J, sunctx);
    if (!J || !LS) {
        SUNMatDestroy(J); SUNLinSolFree(LS);
        IDAFree(&mem); N_VDestroy(y); N_VDestroy(yp); N_VDestroy(id);
        SUNContext_Free(&sunctx);
        return -1;
    }
    
    // Attach linear solver
    IDASetLinearSolver(mem, LS, J);
    IDASetMaxOrd(mem, 2);

    IDASetJacFn(mem, jacobian);
    
    // Set constraints (y_i >= 0)
    constraints = N_VNew_Serial(NEQ, sunctx);
    if (constraints == NULL) {
        SUNLinSolFree(LS); SUNMatDestroy(J);
        IDAFree(&mem); N_VDestroy(y); N_VDestroy(yp); N_VDestroy(id);
        SUNContext_Free(&sunctx);
        return -1;
    }
    sunrealtype *cons_data = NV_DATA_S(constraints);
    //for (int i = 0; i < NEQ; i++) cons_data[i] = 2.0;  // 2.0 = y_i >= 0
    cons_data[1] = 2.0;
    cons_data[3] = 2.0;
    cons_data[4] = 2.0;
    cons_data[5] = 2.0;
    //IDASetConstraints(mem, constraints);
    
    // CRITICAL: Compute consistent initial conditions
    // Use a small time as the first output for IC calculation
    flag = IDACalcIC(mem, IDA_YA_YDP_INIT , 1e-8);
    if (flag != IDA_SUCCESS) {
        fprintf(stderr, "ERROR: IDACalcIC failed with flag %d\n", flag);
        SUNLinSolFree(LS); SUNMatDestroy(J); N_VDestroy(constraints);
        IDAFree(&mem); N_VDestroy(y); N_VDestroy(yp); N_VDestroy(id);
        SUNContext_Free(&sunctx);
        return -1;
    }
    
    // Now solve from t=0 to t_final
    sunrealtype tret;
    for(int w=0; w<1000; w++){
        flag = IDASolve(mem, t_final, &tret, y, yp, IDA_NORMAL);
    }
    
    // Copy results
    for (int i = 0; i < NEQ; i++) y_final[i] = y_data[i];
    
    // Cleanup
    IDAFree(&mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(J);
    N_VDestroy(y); N_VDestroy(yp);
    N_VDestroy(constraints); N_VDestroy(id);
    SUNContext_Free(&sunctx);
    
    return (flag == IDA_SUCCESS) ? 0 : -1;
}


double r_rates(double inputs[5]) {
    // Usage: T(K), <x_NH3>, x_O2, x_NO, relaxation
    double outputs[6];

    KrahnertData data = {
        .T = inputs[0],
        .P = 500, //kPa
        .pNH3 = inputs[1] * 500,
        .pO2 = inputs[2] * 500,
        .pNO = inputs[3] * 500,
        .R = 8.314e-3,
        .surface_density = 2.72e-5,
        .relaxation_factor = inputs[4],
        .T_ref = 385.0 + 273.15
    };

    // Precompute constants
    data.T_factor = (1.0/data.T - 1.0/data.T_ref) / data.R;
    data.surface_relax_inv = 1.0 / (data.surface_density * data.relaxation_factor);    
    
    // Solve

    double y_final[NEQ];
    int status = 0;
    status = solve_krahnert(&data, y_final, 1e1*data.relaxation_factor);
       
    if (status != 0) {
        fprintf(stderr, "Solver failed with status %d\n", status);
        return -1;
    }
       
    // Compute rates
    //double rates[6];
    compute_reaction_rates(&data, y_final, outputs);
    
    //printf("\nReaction rates (mol/m²/s):\n");
    //printf("  NH3: %.6e, N2: %.6e, NO: %.6e, N2O: %.6e\n", 
    //       outputs[0], outputs[1], outputs[2], outputs[5]);

    //printf("output %.6e \n", outputs[0]);

    return outputs[0];
}

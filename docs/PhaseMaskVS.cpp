#include <windows.h>
#include <math.h>
#include <string.h>
#include "usersurf.h"

#pragma warning(disable : 4996)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// Exported User Defined Surface Function
extern "C" {
    int __declspec(dllexport) APIENTRY UserDefinedSurface(USER_DATA* UD, FIXED_DATA* FD);
}

/* a generic Snells law refraction routine */
int Refract(double thisn, double nextn, double* l, double* m, double* n, double ln, double mn, double nn);

/*
BOOL WINAPI DllMain(HANDLE hInst, ULONG ul_reason_for_call, LPVOID lpReserved)
{
    return TRUE;
}
*/

// User Defined Surface Function
int __declspec(dllexport) APIENTRY UserDefinedSurface(USER_DATA* UD, FIXED_DATA* FD) {
    double power, x, y, z, tp, t, loop, r2, r1, sag, dz, mx, my;
    double beta, chi, R_pupil, x_norm, y_norm, exp_term_x, exp_term_y;

    switch (FD->type) {
    case 0: // General Information
        switch (FD->numb) {
        case 0:
            strcpy(UD->string,"ExpPhaseMask");
            break;
        case 1:
            UD->string[0] = '\0'; // Not rotationally symmetric
            break;
        case 2:
            UD->string[0] = '\0'; // Not a gradient index medium
            break;
        }
        break;

    case 1: // Parameter Names
        switch (FD->numb) {
        case 0:
            strcpy(UD->string, "Beta");
            break;
        case 1:
            strcpy(UD->string, "Chi");
            break;
        default:
            UD->string[0] = '\0';
        }
        break;

    case 3: { // Sag Calculation
        // Extract parameters Beta and Chi
        beta = FD->param[0]; // Correct parameter indices
        chi = FD->param[1];
        R_pupil = FD->sdia; // Semi-diameter of the surface


        

        // Coordinates
        x = UD->x;
        y = UD->y;

        // Improved Sag Calculation with a convergence safeguard
        exp_term_x = exp(chi * x * x); // Prevent overflow
        exp_term_y = exp(chi * y * y);

        sag = beta * x * exp_term_x + beta * y * exp_term_y;

        UD->sag1 = sag;
        UD->sag2 = sag; // Same sag for alternate


        break;
    }

    case 4: // Paraxial Ray Trace
        UD->ln = 0.0;
        UD->mn = 0.0;
        UD->nn = -1.0; // Normal points downward

        power = (FD->n2 - FD->n1) * FD->cv; // Optical power

        if ((UD->n) != 0.0)
        {
            (UD->l) = (UD->l) / (UD->n);
            (UD->m) = (UD->m) / (UD->n);

            (UD->l) = (FD->n1 * (UD->l) - (UD->x) * power) / (FD->n2);
            (UD->m) = (FD->n1 * (UD->m) - (UD->y) * power) / (FD->n2);

            /* normalize */
            (UD->n) = sqrt(1 / (1 + (UD->l) * (UD->l) + (UD->m) * (UD->m)));
            /* de-paraxialize */
            (UD->l) = (UD->l) * (UD->n);
            (UD->m) = (UD->m) * (UD->n);
        }

        break;
    case 5: // Real Ray Trace

        // do not allow vertical rays
        if (fabs(UD->n) < 1E-5) return -1;


        x = UD->x;
        y = UD->y;
        z = UD->z;
        tp = 0.0;

        t = 9e9;
        loop = 0;


        while (fabs(t) > 1e-10)
        {
            // what is the radial coordinate? Note this might change as we iterate
            r2 = x * x + y * y;
            r1 = sqrt(r2);

            /////////////////////// Compute the sag //////////////////////////////////////////////////////////
            // Extract parameters Beta and Chi
            beta = FD->param[0]; // Correct parameter indices
            chi = FD->param[1];

            R_pupil = FD->sdia; // Semi-diameter of the surface


           

            // Coordinates
            x = UD->x;
            y = UD->y;

            // Improved Sag Calculation with a convergence safeguard
            exp_term_x = exp(chi * x * x); // Prevent overflow
            exp_term_y = exp(chi * y * y);

            sag = beta * x * exp_term_x + beta * y * exp_term_y;
            ///////////////////////////////////////////////////////////////////////////////////////

            // how far are we away in z?
            dz = sag - z;

            /* now compute how far along the ray this is */
            //t = dz / (UD->n); correct but may be too aggressive and may not converge
            t = dz * fabs(UD->n); // I will use this one

            /* propagate the additional "t" distance */
            x += UD->l * t;
            y += UD->m * t;
            z += UD->n * t;

            /* add in the optical path */
            tp += t;

            /* prevent infinte loop if no convergence */
            loop++;
            if (loop > 1000) return(-1);
        }

        // okay, we should be at the intercept coordinates now
        UD->x = x;
        UD->y = y;
        UD->z = z;

        // don't forget the path!
        UD->path = tp;

        // now do the normals
        if (x == 0.0 && y == 0.0)
        {
            UD->ln = 0;
            UD->mn = 0;
            UD->nn = -1;
        }
        else
        {

            beta = FD->param[0]; // Correct parameter indices
            chi = FD->param[1];
            R_pupil = FD->sdia; // Semi-diameter of the surface

            

            // sag derivatives 

            // Improved Sag Calculation with a convergence safeguard
            exp_term_x = exp(chi * x * x); // Prevent overflow
            exp_term_y = exp(chi * y * y);

            mx = beta * (1 + 2 * chi * x * x) * exp_term_x;

            my = beta * (1 + 2 * chi * y * y) * exp_term_y;

            UD->nn = -sqrt(1 / (1 + (mx * mx) + (my * my)));
            UD->ln = -mx * UD->nn;
            UD->mn = -my * UD->nn;
        }

        if (Refract(FD->n1, FD->n2, &UD->l, &UD->m, &UD->n, UD->ln, UD->mn, UD->nn)) return(-FD->surf);
        break;

    case 6: // GRIN Index Data
        /* ZEMAX wants the index, dn/dx, dn/dy, and dn/dz at the given x, y, z. */
        /* This is only required for gradient index surfaces, so return dummy values */
        UD->index = FD->n2;
        UD->dndx = 0.0;
        UD->dndy = 0.0;
        UD->dndz = 0.0;
        break;

    case 7: // Safe Data Initialization
        /* ZEMAX wants the "safe" data. */
        /* this is used by ZEMAX to set the initial values for all parameters and extra data */
        /* when the user first changes to this surface type. */
        /* this is the only time the DLL should modify the data in the FIXED_DATA FD structure */


        FD->param[0] = 0.1; // Default Beta
        FD->param[1] = 0.01;  // Default Chi
        break;
    }
    return 0;
}
int Refract(double thisn, double nextn, double* l, double* m, double* n, double ln, double mn, double nn)
{
    double nr, cosi, cosi2, rad, cosr, gamma;
    if (thisn != nextn)
    {
        nr = thisn / nextn;
        cosi = fabs((*l) * ln + (*m) * mn + (*n) * nn);
        cosi2 = cosi * cosi;
        if (cosi2 > 1) cosi2 = 1;
        rad = 1 - ((1 - cosi2) * (nr * nr));
        if (rad < 0) return(-1);
        cosr = sqrt(rad);
        gamma = nr * cosi - cosr;
        (*l) = (nr * (*l)) + (gamma * ln);
        (*m) = (nr * (*m)) + (gamma * mn);
        (*n) = (nr * (*n)) + (gamma * nn);
    }
    return 0;
}

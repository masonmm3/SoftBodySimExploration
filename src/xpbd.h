#include "raylib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ---------------------------------------------------------
//  Tunables
// ---------------------------------------------------------
#define CIRCLE_VERTS    24
#define MAX_PARTICLES   (CIRCLE_VERTS + 1)
#define MAX_CONSTRAINTS 256

#define SCREEN_W    900
#define SCREEN_H    700
#define DT          0.016f      // master timestep  (subdivided by SUBSTEPS)
#define SUBSTEPS    4           // XPBD substeps — smaller dt per step
#define GRAVITY     980.0f
#define DAMPING     0.990f      // per-substep velocity damping

// ---------------------------------------------------------
//  Compliance presets
//  α = 1/k in SI units; we work in pixels so scale accordingly.
//  Low  α → stiff spring   (α → 0 is perfectly rigid)
//  High α → soft / jelly
// ---------------------------------------------------------
#define ALPHA_MIN   1e-7f       // nearly rigid
#define ALPHA_MAX   5e-2f       // very soft

// ---------------------------------------------------------
//  Data structures
// ---------------------------------------------------------
typedef struct {
    float x,  y;        // confirmed position (xⁿ)
    float px, py;       // predicted position
    float vx, vy;       // velocity
    float w;            // inverse mass (0 = pinned)
} Particle;

typedef struct {
    int   a, b;
    float rest_len;
    float lambda;       // Lagrange multiplier — reset each substep
} Constraint;

// ---------------------------------------------------------
//  Math helpers
// ---------------------------------------------------------
inline float vlen(float dx, float dy) { return sqrtf(dx*dx + dy*dy); }

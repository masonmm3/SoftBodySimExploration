#include "raylib.h"
#include <math.h>
#include <stdio.h>

#include <stdlib.h>
#include <string.h>

// ---------------------------------------------------------
//  Simulation parameters
// ---------------------------------------------------------
#define MAX_PARTICLES   64
#define MAX_CONSTRAINTS 256
#define CIRCLE_VERTS    24          // perimeter particles
#define INNER_VERTS     1           // 1 center particle

#define TOTAL_PARTICLES (CIRCLE_VERTS + INNER_VERTS)

#define SCREEN_W  900
#define SCREEN_H  700
#define DT        0.016f            // ~60 Hz timestep
#define GRAVITY   980.0f            // px / s^2
#define DAMPING   0.98f

// ---------------------------------------------------------
//  Data structures
// ---------------------------------------------------------
typedef struct {
    float x, y;       // current position
    float px, py;     // predicted position
    float vx, vy;     // velocity
    float w;          // inverse mass  (0 = pinned)
} Particle;

typedef enum {
    CONSTRAINT_DISTANCE,
    CONSTRAINT_SHAPE   // shape-matching style — keep as distance for now
} ConstraintType;

typedef struct {
    int    a, b;
    float  rest_len;
    float  stiffness;  // [0..1]
} Constraint;

// ---------------------------------------------------------
//  Helper math
// ---------------------------------------------------------
float dist2(float ax, float ay, float bx, float by) {
    float dx = bx - ax, dy = by - ay;
    return sqrtf(dx*dx + dy*dy);
}

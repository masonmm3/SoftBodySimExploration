/*
 * XPBD 2D Softbody Circle Simulation
 * Based on Müller et al. "XPBD: Position-Based Simulation of Compliant
 * Constrained Dynamics" (2016), and Servin et al. (2006).
 *
 * Key XPBD departure from PBD:
 *   - Each constraint carries a Lagrange multiplier λ (init 0 each substep)
 *   - Compliance α = 1/k  (physical spring stiffness, NOT a [0..1] scalar)
 *   - Scaled compliance α̃ = α / dt²   (eq. 5 in notes: U = ½ Cᵀ α⁻¹ C)
 *   - Solver update per constraint:
 *       Δλ = -(C + α̃·λ) / (∇C·M⁻¹·∇Cᵀ + α̃)      (from eq.12 in notes)
 *       Δx = w · ∇C · Δλ
 *   - Stiffness is NOW TIMESTEP-INDEPENDENT — the defining win over PBD
 *
 * Controls:
 *   Mouse drag       - grab and pull vertices
 *   UP / DOWN        - increase / decrease compliance (softer ↑, stiffer ↓)
 *   LEFT / RIGHT     - increase / decrease solver iterations
 *   R                - reset simulation
 *   G                - toggle gravity
 *   SPACE            - pause / unpause
 *
 *  my compile commands gcc xpbd.c -o xpbd.exe -I C:\raylib\raylib\src -L C:\raylib\raylib\src -lraylib -lopengl32 -lgdi32 -lwinmm -D__MSVCRT_VERSION__=0x0700
 */

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
//  Globals
// ---------------------------------------------------------
static Particle   particles[MAX_PARTICLES];
static Constraint constraints[MAX_CONSTRAINTS];
static int        num_particles   = 0;
static int        num_constraints = 0;

// α (compliance): user-controlled, maps to physical spring constant k = 1/α
// We expose log-scale so the slider feels linear perceptually
static float global_alpha   = 1e-4f;   // start moderately stiff
static int   solver_iters   = 8;
static bool  gravity_on     = true;
static bool  paused         = false;

// Mouse drag
static int   drag_idx = -1;
static float drag_ox, drag_oy;

// Bounds
static float floor_y, wall_l, wall_r;

// ---------------------------------------------------------
//  Math helpers
// ---------------------------------------------------------
static inline float vlen(float dx, float dy) { return sqrtf(dx*dx + dy*dy); }

// ---------------------------------------------------------
//  Build the soft circle
// ---------------------------------------------------------
static void build_circle(float cx, float cy, float radius) {
    num_particles   = 0;
    num_constraints = 0;

    // Particle 0: center
    particles[0] = (Particle){ cx, cy, cx, cy, 0, 0, 1.0f };
    num_particles = 1;

    // Particles 1..CIRCLE_VERTS: perimeter
    for (int i = 0; i < CIRCLE_VERTS; i++) {
        float angle = (float)i / CIRCLE_VERTS * 2.0f * 3.14159265f;
        float px = cx + cosf(angle) * radius;
        float py = cy + sinf(angle) * radius;
        particles[num_particles++] = (Particle){ px, py, px, py, 0, 0, 1.0f };
    }

    // Constraint helper (lambda always starts 0; α applied at solve time)
#define ADD_DIST(ia, ib) do { \
    int _a = (ia), _b = (ib); \
    float _l = vlen(particles[_a].x - particles[_b].x, \
                    particles[_a].y - particles[_b].y); \
    constraints[num_constraints++] = (Constraint){ _a, _b, _l, 0.0f }; \
} while(0)

    // 1. Perimeter edges
    for (int i = 0; i < CIRCLE_VERTS; i++)
        ADD_DIST(1 + i, 1 + (i + 1) % CIRCLE_VERTS);

    // 2. Spokes: center ↔ rim
    for (int i = 0; i < CIRCLE_VERTS; i++)
        ADD_DIST(0, 1 + i);

    // 3. Diameter cross-braces (volumetric incompressibility)
    for (int i = 0; i < CIRCLE_VERTS / 2; i++)
        ADD_DIST(1 + i, 1 + (i + CIRCLE_VERTS / 2) % CIRCLE_VERTS);

#undef ADD_DIST
}

// ---------------------------------------------------------
//  XPBD constraint projection  (eq. 12 from notes)
//
//  For a single distance constraint C(p) = |pₐ - p_b| - L₀:
//    ∇C w.r.t. pₐ =  n̂  (unit vec a→b)
//    ∇C w.r.t. p_b = -n̂
//
//  The XPBD update (one Gauss-Seidel pass):
//    α̃ = α / (dt_sub²)
//    Δλ = -(C + α̃·λ) / (wₐ + w_b + α̃)
//    λ  += Δλ
//    pₐ += wₐ · n̂ · Δλ
//    p_b -= w_b · n̂ · Δλ
// ---------------------------------------------------------
static void xpbd_project(float dt_sub) {
    // α̃ = α / dt²  — scaled compliance (from eq. 5 / 12 in your notes)
    float alpha_tilde = global_alpha / (dt_sub * dt_sub);

    for (int iter = 0; iter < solver_iters; iter++) {
        for (int ci = 0; ci < num_constraints; ci++) {
            Constraint *c = &constraints[ci];
            Particle   *a = &particles[c->a];
            Particle   *b = &particles[c->b];

            float dx = b->px - a->px;
            float dy = b->py - a->py;
            float d  = vlen(dx, dy);
            if (d < 1e-9f) continue;

            float nx = dx / d;
            float ny = dy / d;

            // C(p) = d - rest_len
            float C = d - c->rest_len;

            float w_sum = a->w + b->w;
            if (w_sum < 1e-9f) continue;

            // Δλ = -(C + α̃·λ) / (∑w + α̃)
            float delta_lambda = -(C + alpha_tilde * c->lambda) / (w_sum + alpha_tilde);
            c->lambda += delta_lambda;

            // Δx = w · n̂ · Δλ
            a->px -= a->w * nx * delta_lambda;
            a->py -= a->w * ny * delta_lambda;
            b->px += b->w * nx * delta_lambda;
            b->py += b->w * ny * delta_lambda;
        }
    }
}

// ---------------------------------------------------------
//  Collision (AABB deprojection, C(p) = (p - q_c)·n_c ≥ 0)
// ---------------------------------------------------------
static void resolve_bounds(float dt_sub) {
    for (int i = 0; i < num_particles; i++) {
        Particle *p = &particles[i];
        if (p->w == 0.0f) continue;

        // Floor
        if (p->py > floor_y) {
            p->py = floor_y;
            // friction-style velocity clamp
            p->vx *= 0.85f;
            p->vy  = 0.0f;
        }
        if (p->py < 22.0f)  { p->py = 22.0f;  p->vy = fabsf(p->vy) * 0.3f; }
        if (p->px < wall_l) { p->px = wall_l;  p->vx = fabsf(p->vx) * 0.4f; }
        if (p->px > wall_r) { p->px = wall_r;  p->vx = -fabsf(p->vx) * 0.4f; }
    }
    (void)dt_sub;
}

// ---------------------------------------------------------
//  One XPBD step (with internal substepping)
//  The substep loop is the proper XPBD integration:
//    for each substep:
//      1. predict positions
//      2. reset λ for all constraints  ← critical XPBD detail
//      3. project constraints
//      4. update velocities from Δx/dt_sub
// ---------------------------------------------------------
static void sim_step(float dt) {
    float dt_sub = dt / (float)SUBSTEPS;

    for (int sub = 0; sub < SUBSTEPS; sub++) {

        // --- Predict ---
        for (int i = 0; i < num_particles; i++) {
            Particle *p = &particles[i];
            if (p->w == 0.0f) continue;

            if (gravity_on) p->vy += GRAVITY * dt_sub;
            p->vx *= DAMPING;
            p->vy *= DAMPING;

            p->px = p->x + dt_sub * p->vx;
            p->py = p->y + dt_sub * p->vy;
        }

        // Mouse drag: move predicted pos directly
        if (drag_idx >= 0) {
            Vector2 mouse = GetMousePosition();
            particles[drag_idx].px = mouse.x + drag_ox;
            particles[drag_idx].py = mouse.y + drag_oy;
        }

        // --- Reset λ for this substep (XPBD requirement) ---
        for (int ci = 0; ci < num_constraints; ci++)
            constraints[ci].lambda = 0.0f;

        // --- Solve ---
        xpbd_project(dt_sub);

        // --- Bounds ---
        resolve_bounds(dt_sub);

        // --- Velocity update: v = (x_new - x_old) / dt_sub ---
        for (int i = 0; i < num_particles; i++) {
            Particle *p = &particles[i];
            if (p->w == 0.0f) continue;
            p->vx = (p->px - p->x) / dt_sub;
            p->vy = (p->py - p->y) / dt_sub;
            p->x  = p->px;
            p->y  = p->py;
        }
    }
}

// ---------------------------------------------------------
//  Drawing
// ---------------------------------------------------------
static void draw_sim(void) {
    // Internal structure lines
    for (int ci = 0; ci < num_constraints; ci++) {
        Constraint *c = &constraints[ci];
        Particle   *a = &particles[c->a];
        Particle   *b = &particles[c->b];
        bool spoke = (c->a == 0 || c->b == 0);
        Color col  = spoke
            ? (Color){ 255, 120, 60,  70 }
            : (Color){ 80,  200, 160, 90 };
        DrawLineEx((Vector2){a->x, a->y}, (Vector2){b->x, b->y}, 1.0f, col);
    }

    // Filled polygon (fan from center)
    for (int i = 0; i < CIRCLE_VERTS; i++) {
        int a = 1 + i;
        int b = 1 + (i + 1) % CIRCLE_VERTS;
        DrawTriangle(
            (Vector2){ particles[0].x, particles[0].y },
            (Vector2){ particles[b].x, particles[b].y },
            (Vector2){ particles[a].x, particles[a].y },
            (Color){ 255, 120, 40, 45 }
        );
    }

    // Perimeter
    for (int i = 0; i < CIRCLE_VERTS; i++) {
        int a = 1 + i;
        int b = 1 + (i + 1) % CIRCLE_VERTS;
        DrawLineEx(
            (Vector2){ particles[a].x, particles[a].y },
            (Vector2){ particles[b].x, particles[b].y },
            2.2f, (Color){ 255, 160, 60, 230 }
        );
    }

    // Particles
    for (int i = 0; i < num_particles; i++) {
        Color c = (i == drag_idx) ? YELLOW
                : (i == 0)        ? (Color){255,80,80,255}
                :                   WHITE;
        DrawCircleV((Vector2){ particles[i].x, particles[i].y },
                    i == 0 ? 4.5f : 3.5f, c);
    }

    // Bounds
    Color wall_col = (Color){ 180, 180, 180, 180 };
    DrawLineEx((Vector2){ wall_l, floor_y }, (Vector2){ wall_r, floor_y }, 2.f, wall_col);
    DrawLineEx((Vector2){ wall_l, 20 },      (Vector2){ wall_l, floor_y }, 2.f, wall_col);
    DrawLineEx((Vector2){ wall_r, 20 },      (Vector2){ wall_r, floor_y }, 2.f, wall_col);
}

// ---------------------------------------------------------
//  HUD
// ---------------------------------------------------------
static void draw_hud(void) {
    DrawRectangle(10, 10, 340, 200, (Color){ 10, 10, 16, 210 });
    DrawRectangleLines(10, 10, 340, 200, (Color){ 255, 160, 60, 120 });

    DrawText("XPBD SOFTBODY  //  2D CIRCLE", 20, 18, 14, (Color){255,160,60,255});

    // Alpha bar (log-scale display)
    int bx = 20, by = 52, bw = 310, bh = 14;

    // Normalize alpha to [0,1] on log scale for the bar
    float log_min   = log10f(ALPHA_MIN);
    float log_max   = log10f(ALPHA_MAX);
    float log_cur   = log10f(global_alpha);
    float bar_frac  = (log_cur - log_min) / (log_max - log_min);  // 0=stiff,1=soft

    DrawText("COMPLIANCE α  [UP=stiffer / DOWN=softer]", bx, by - 14, 10, LIGHTGRAY);
    DrawRectangle(bx, by, bw, bh, (Color){30,20,10,255});
    // Color: blue (stiff) → orange (soft)
    Color bar_col = ColorFromHSV(30.0f - bar_frac * 180.0f, 0.85f, 0.95f);
    DrawRectangle(bx, by, (int)(bw * bar_frac), bh, bar_col);
    DrawRectangleLines(bx, by, bw, bh, GRAY);

    char buf[160];
    // Show α and derived k = 1/α
    snprintf(buf, sizeof(buf), "α = %.2e   k = 1/α = %.1f", global_alpha, 1.0f / global_alpha);
    DrawText(buf, bx, by + bh + 4, 11, WHITE);

    // Stiffness label
    const char *label;
    if (bar_frac < 0.2f)      label = "RIGID";
    else if (bar_frac < 0.4f) label = "STIFF";
    else if (bar_frac < 0.6f) label = "MODERATE";
    else if (bar_frac < 0.8f) label = "SOFT";
    else                       label = "JELLY";
    DrawText(label, bx + bw - MeasureText(label,12) - 4, by + bh + 4, 12,
             bar_col);

    // Solver iters
    int si_y = by + bh + 34;
    DrawText("SOLVER ITERS  [LEFT / RIGHT]", bx, si_y - 14, 10, LIGHTGRAY);
    DrawRectangle(bx, si_y, bw, bh, (Color){20,20,40,255});
    DrawRectangle(bx, si_y, (int)(bw * (solver_iters / 30.0f)), bh,
                  (Color){120, 200, 255, 200});
    DrawRectangleLines(bx, si_y, bw, bh, GRAY);
    snprintf(buf, sizeof(buf), "iters = %d   substeps = %d", solver_iters, SUBSTEPS);
    DrawText(buf, bx, si_y + bh + 4, 11, WHITE);

    // XPBD note: stiffness independent of iterations
    DrawText("* stiffness is timestep-independent (XPBD)", bx, si_y + bh + 20, 10,
             (Color){180,180,100,200});

    // Status
    int st_y = si_y + bh + 38;
    snprintf(buf, sizeof(buf), "gravity %s  |  %s  |  R=reset  G=gravity  SPACE=pause",
             gravity_on ? "ON" : "OFF", paused ? "PAUSED" : "RUNNING");
    DrawText(buf, bx, st_y, 10, (Color){140,140,140,255});
}

// ---------------------------------------------------------
//  Main
// ---------------------------------------------------------
int main(void) {
    InitWindow(SCREEN_W, SCREEN_H, "XPBD Softbody — 2D Circle");
    SetTargetFPS(60);

    floor_y = SCREEN_H - 60.0f;
    wall_l  = 40.0f;
    wall_r  = SCREEN_W - 40.0f;

    build_circle(SCREEN_W / 2.0f, SCREEN_H / 2.0f - 60.0f, 100.0f);

    while (!WindowShouldClose()) {

        // --- Input ---
        if (IsKeyPressed(KEY_R))
            build_circle(SCREEN_W / 2.0f, SCREEN_H / 2.0f - 60.0f, 100.0f);
        if (IsKeyPressed(KEY_G))     gravity_on = !gravity_on;
        if (IsKeyPressed(KEY_SPACE)) paused     = !paused;

        // Compliance: UP = more stiff (lower α), DOWN = softer (higher α)
        // Multiply/divide by factor each frame for smooth log-scale feel
        if (IsKeyDown(KEY_UP)) {
            global_alpha /= 1.04f;
            if (global_alpha < ALPHA_MIN) global_alpha = ALPHA_MIN;
        }
        if (IsKeyDown(KEY_DOWN)) {
            global_alpha *= 1.04f;
            if (global_alpha > ALPHA_MAX) global_alpha = ALPHA_MAX;
        }

        // Solver iterations
        if (IsKeyPressed(KEY_RIGHT)) { solver_iters++; if (solver_iters > 30) solver_iters = 30; }
        if (IsKeyPressed(KEY_LEFT))  { solver_iters--; if (solver_iters <  1) solver_iters  = 1; }

        // Mouse drag
        Vector2 mouse = GetMousePosition();
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            float best = 500.0f;
            drag_idx = -1;
            for (int i = 0; i < num_particles; i++) {
                float d = vlen(mouse.x - particles[i].x, mouse.y - particles[i].y);
                if (d < best) { best = d; drag_idx = i; }
            }
            if (drag_idx >= 0) {
                drag_ox = particles[drag_idx].x - mouse.x;
                drag_oy = particles[drag_idx].y - mouse.y;
            }
        }
        if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) drag_idx = -1;

        // --- Sim ---
        if (!paused) sim_step(DT);

        // --- Draw ---
        BeginDrawing();
            ClearBackground((Color){ 10, 8, 14, 255 });
            draw_sim();
            draw_hud();
        EndDrawing();
    }

    CloseWindow();
    return 0;
}

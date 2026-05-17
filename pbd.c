/*
 * PBD 2D Softbody Circle Simulation
 * Based on Müller et al. "Position Based Dynamics" (2007)
 *
 * Controls:
 *   Mouse drag     - grab and pull vertices
 *   UP/DOWN arrow  - increase/decrease stiffness
 *   LEFT/RIGHT     - increase/decrease solver iterations
 *   R              - reset simulation
 *   G              - toggle gravity
 *   SPACE          - pause/unpause
 *
 * my compile commands. gcc pbd.c -o pbd.exe -I C:\raylib\raylib\src -L C:\raylib\raylib\src -lraylib -lopengl32 -lgdi32 -lwinmm -D__MSVCRT_VERSION__=0x0700
 */

#include "pbd.h"

// ---------------------------------------------------------
//  Globals
// ---------------------------------------------------------
static Particle    particles[MAX_PARTICLES];
static Constraint  constraints[MAX_CONSTRAINTS];
static int         num_particles   = 0;
static int         num_constraints = 0;

static float       global_stiffness = 0.8f;
static int         solver_iters     = 10;
static bool        gravity_on       = true;
static bool        paused           = false;

// Mouse drag
static int   drag_idx  = -1;
static float drag_ox, drag_oy;    // offset from particle to mouse at grab

// Floor / wall bounds
static float floor_y;
static float wall_l, wall_r;



// ---------------------------------------------------------
//  Build the soft circle
// ---------------------------------------------------------
static void build_circle(float cx, float cy, float radius) {
    num_particles   = 0;
    num_constraints = 0;

    // Center particle
    particles[0] = (Particle){ cx, cy, cx, cy, 0, 0, 1.0f };
    num_particles = 1;

    // Perimeter particles
    for (int i = 0; i < CIRCLE_VERTS; i++) {
        float angle = (float)i / CIRCLE_VERTS * 2.0f * 3.14159265f;
        float px = cx + cosf(angle) * radius;
        float py = cy + sinf(angle) * radius;
        particles[num_particles++] = (Particle){ px, py, px, py, 0, 0, 1.0f };
    }

    // --- Constraints ---

    // 1. Perimeter edge constraints (gives the circle its outline)
    for (int i = 0; i < CIRCLE_VERTS; i++) {
        int a = 1 + i;
        int b = 1 + (i + 1) % CIRCLE_VERTS;
        constraints[num_constraints++] = (Constraint){
            a, b, dist2(particles[a].x, particles[a].y,
                        particles[b].x, particles[b].y),
            global_stiffness
        };
    }

    // 2. Spoke constraints: center ↔ perimeter (volumetric stiffness)
    for (int i = 0; i < CIRCLE_VERTS; i++) {
        int a = 0;          // center
        int b = 1 + i;
        constraints[num_constraints++] = (Constraint){
            a, b, dist2(particles[a].x, particles[a].y,
                        particles[b].x, particles[b].y),
            global_stiffness
        };
    }

    // 3. Cross bracing: every other perimeter ↔ opposite (pressure-like)
    for (int i = 0; i < CIRCLE_VERTS / 2; i++) {
        int a = 1 + i;
        int b = 1 + (i + CIRCLE_VERTS / 2) % CIRCLE_VERTS;
        constraints[num_constraints++] = (Constraint){
            a, b, dist2(particles[a].x, particles[a].y,
                        particles[b].x, particles[b].y),
            global_stiffness
        };
    }
}

// ---------------------------------------------------------
//  PBD project constraints  (Müller §3)
//  Δp = -( C(p) / |∇C|² ) * ∇C
//  with stiffness k' = 1 - (1-k)^(1/n_s)
// ---------------------------------------------------------
static void project_constraints(void) {
    for (int iter = 0; iter < solver_iters; iter++) {
        // Compute per-iteration stiffness correction k'
        // k' = 1 - (1-k)^(1/n_s)
        float kprime = 1.0f - powf(1.0f - global_stiffness, 1.0f / (float)solver_iters);

        for (int ci = 0; ci < num_constraints; ci++) {
            Constraint *c = &constraints[ci];
            Particle   *a = &particles[c->a];
            Particle   *b = &particles[c->b];

            float dx = b->px - a->px;
            float dy = b->py - a->py;
            float d   = sqrtf(dx*dx + dy*dy);
            if (d < 1e-6f) continue;

            // C(p) = |p_a - p_b| - rest_len
            float C = d - c->rest_len;

            // ∇C is the unit vector along the edge
            float nx = dx / d;
            float ny = dy / d;

            // Sum of inverse masses
            float w_sum = a->w + b->w;
            if (w_sum < 1e-6f) continue;

            // Δp = -(C / |∇C|²) * k' — |∇C|² = 1 for distance constraint
            float delta = kprime * C / w_sum;

            a->px += a->w * delta * nx;
            a->py += a->w * delta * ny;
            b->px -= b->w * delta * nx;
            b->py -= b->w * delta * ny;
        }
    }
}

// ---------------------------------------------------------
//  Collision response: axis-aligned bounding box
// ---------------------------------------------------------
static void resolve_bounds(void) {
    for (int i = 0; i < num_particles; i++) {
        Particle *p = &particles[i];
        if (p->w == 0.0f) continue;

        if (p->py > floor_y) { p->py = floor_y; p->vy *= -0.3f; }
        if (p->py < 20.0f)   { p->py = 20.0f;  p->vy *= -0.3f; }
        if (p->px < wall_l)  { p->px = wall_l;  p->vx *= -0.3f; }
        if (p->px > wall_r)  { p->px = wall_r;  p->vx *= -0.3f; }
    }
}

// ---------------------------------------------------------
//  One PBD simulation step
// ---------------------------------------------------------
static void sim_step(float dt) {
    // 5. Velocity + force integration
    for (int i = 0; i < num_particles; i++) {
        Particle *p = &particles[i];
        if (p->w == 0.0f) continue;

        if (gravity_on) p->vy += GRAVITY * dt;
        p->vx *= DAMPING;
        p->vy *= DAMPING;

        // 7. Estimate predicted position
        p->px = p->x + dt * p->vx;
        p->py = p->y + dt * p->vy;
    }

    // Mouse drag: teleport grabbed particle to mouse
    if (drag_idx >= 0) {
        Vector2 mouse = GetMousePosition();
        particles[drag_idx].px = mouse.x + drag_ox;
        particles[drag_idx].py = mouse.y + drag_oy;
    }

    // 8+9. Constraint projection (solver iterations inside)
    project_constraints();

    // Bounds
    resolve_bounds();

    // 12. Derive velocity from position delta
    for (int i = 0; i < num_particles; i++) {
        Particle *p = &particles[i];
        if (p->w == 0.0f) continue;
        p->vx = (p->px - p->x) / dt;
        p->vy = (p->py - p->y) / dt;
        p->x  = p->px;
        p->y  = p->py;
    }
}

// ---------------------------------------------------------
//  Drawing
// ---------------------------------------------------------
static void draw_sim(void) {
    // Draw constraint lines (internal structure)
    for (int ci = 0; ci < num_constraints; ci++) {
        Constraint *c = &constraints[ci];
        Particle   *a = &particles[c->a];
        Particle   *b = &particles[c->b];
        // Color by constraint type: spokes = dark, perimeter = brighter
        Color col = (c->a == 0 || c->b == 0)
            ? (Color){ 60, 120, 200, 80 }
            : (Color){ 80, 220, 180, 120 };
        DrawLineEx((Vector2){a->x, a->y}, (Vector2){b->x, b->y}, 1.0f, col);
    }

    // Fill circle with a fan polygon using perimeter particles
    for (int i = 0; i < CIRCLE_VERTS; i++) {
        int    a   = 1 + i;
        int    b   = 1 + (i + 1) % CIRCLE_VERTS;
        // Semi-transparent fill triangle: center → a → b
        DrawTriangle(
            (Vector2){ particles[0].x, particles[0].y },
            (Vector2){ particles[b].x, particles[b].y },
            (Vector2){ particles[a].x, particles[a].y },
            (Color){ 40, 160, 255, 55 }
        );
    }

    // Perimeter outline
    for (int i = 0; i < CIRCLE_VERTS; i++) {
        int a = 1 + i;
        int b = 1 + (i + 1) % CIRCLE_VERTS;
        DrawLineEx(
            (Vector2){ particles[a].x, particles[a].y },
            (Vector2){ particles[b].x, particles[b].y },
            2.2f,
            (Color){ 80, 210, 255, 220 }
        );
    }

    // Particles
    for (int i = 0; i < num_particles; i++) {
        Color c = (i == drag_idx) ? YELLOW : (i == 0 ? RED : WHITE);
        float r = (i == 0) ? 4.0f : 3.5f;
        DrawCircleV((Vector2){ particles[i].x, particles[i].y }, r, c);
    }

    // Floor line
    DrawLineEx((Vector2){ wall_l, floor_y }, (Vector2){ wall_r, floor_y },
               2.0f, (Color){ 180, 180, 180, 200 });
    DrawLineEx((Vector2){ wall_l, 20 }, (Vector2){ wall_l, floor_y },
               2.0f, (Color){ 180, 180, 180, 200 });
    DrawLineEx((Vector2){ wall_r, 20 }, (Vector2){ wall_r, floor_y },
               2.0f, (Color){ 180, 180, 180, 200 });
}

// ---------------------------------------------------------
//  HUD
// ---------------------------------------------------------
static void draw_hud(void) {
    // Panel background
    DrawRectangle(10, 10, 310, 170, (Color){ 10, 12, 20, 200 });
    DrawRectangleLines(10, 10, 310, 170, (Color){ 80, 210, 255, 100 });

    // Stiffness bar
    int bar_x = 20, bar_y = 50, bar_w = 280, bar_h = 14;
    DrawText("PBD SOFTBODY  //  2D CIRCLE", 20, 18, 14, (Color){80,210,255,255});
    DrawText("STIFFNESS  [UP / DOWN]", bar_x, bar_y - 14, 11, LIGHTGRAY);
    DrawRectangle(bar_x, bar_y, bar_w, bar_h, (Color){30,30,50,255});
    DrawRectangle(bar_x, bar_y, (int)(bar_w * global_stiffness), bar_h,
                  ColorFromHSV(200 + global_stiffness * 60, 0.85f, 0.9f));
    DrawRectangleLines(bar_x, bar_y, bar_w, bar_h, GRAY);

    char buf[128];
    snprintf(buf, sizeof(buf), "k = %.3f", global_stiffness);
    DrawText(buf, bar_x, bar_y + bar_h + 4, 12, WHITE);

    // Solver iterations bar
    int si_y = bar_y + bar_h + 28;
    DrawText("SOLVER ITERS  [LEFT / RIGHT]", bar_x, si_y - 14, 11, LIGHTGRAY);
    DrawRectangle(bar_x, si_y, bar_w, bar_h, (Color){30,30,50,255});
    DrawRectangle(bar_x, si_y, (int)(bar_w * (solver_iters / 50.0f)), bar_h,
                  (Color){255,180,50,200});
    DrawRectangleLines(bar_x, si_y, bar_w, bar_h, GRAY);
    snprintf(buf, sizeof(buf), "iterations = %d", solver_iters);
    DrawText(buf, bar_x, si_y + bar_h + 4, 12, WHITE);

    // Status
    int st_y = si_y + bar_h + 26;
    snprintf(buf, sizeof(buf), "gravity %s  |  %s  |  R=reset  G=gravity",
             gravity_on ? "ON" : "OFF",
             paused    ? "PAUSED" : "RUNNING");
    DrawText(buf, bar_x, st_y, 11, (Color){160,160,160,255});

    DrawText("drag vertices with mouse", bar_x, st_y + 16, 11, (Color){120,120,120,255});
}

// ---------------------------------------------------------
//  Main
// ---------------------------------------------------------
int main(void) {
    InitWindow(SCREEN_W, SCREEN_H, "PBD Softbody — 2D Circle");
    SetTargetFPS(60);

    floor_y = SCREEN_H - 60.0f;
    wall_l  = 40.0f;
    wall_r  = SCREEN_W - 40.0f;

    build_circle(SCREEN_W / 2.0f, SCREEN_H / 2.0f - 50.0f, 100.0f);

    while (!WindowShouldClose()) {
        // ---- Input ----
        if (IsKeyPressed(KEY_R)) {
            build_circle(SCREEN_W / 2.0f, SCREEN_H / 2.0f - 50.0f, 100.0f);
        }
        if (IsKeyPressed(KEY_G)) gravity_on = !gravity_on;
        if (IsKeyPressed(KEY_SPACE)) paused = !paused;

        // Stiffness
        if (IsKeyDown(KEY_UP)) {
            global_stiffness += 0.005f;
            if (global_stiffness > 1.0f) global_stiffness = 1.0f;
        }
        if (IsKeyDown(KEY_DOWN)) {
            global_stiffness -= 0.005f;
            if (global_stiffness < 0.01f) global_stiffness = 0.01f;
        }

        // Solver iterations
        if (IsKeyPressed(KEY_RIGHT)) { solver_iters++; if (solver_iters > 50) solver_iters = 50; }
        if (IsKeyPressed(KEY_LEFT))  { solver_iters--; if (solver_iters < 1)  solver_iters  = 1; }

        // Mouse drag — find closest particle on press
        Vector2 mouse = GetMousePosition();
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            float best = 400.0f;
            drag_idx = -1;
            for (int i = 0; i < num_particles; i++) {
                float d = dist2(mouse.x, mouse.y, particles[i].x, particles[i].y);
                if (d < best) { best = d; drag_idx = i; }
            }
            if (drag_idx >= 0) {
                drag_ox = particles[drag_idx].x - mouse.x;
                drag_oy = particles[drag_idx].y - mouse.y;
            }
        }
        if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) drag_idx = -1;

        // ---- Simulation ----
        if (!paused) sim_step(DT);

        // ---- Draw ----
        BeginDrawing();
            ClearBackground((Color){ 8, 10, 18, 255 });
            draw_sim();
            draw_hud();
        EndDrawing();
    }

    CloseWindow();
    return 0;
}

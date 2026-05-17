# why PBD
---
Integration produces error porportional to the timestep. Solutions such as substepping **CAN** solve the problem however they are computationally expensive and not a viable solution for real time.

PBD gets around this by controlling the position hence the name, Position Based Dynamics. This allows for precise control of positional limits **AND** collision behavior making it ideal for soft body and high load dynamics.

> [!INFO]- Previous work
>  The paper mentions `Verlet Integration` as a precursor
# Basic Equation
---

> [!info]+ Paramters
>  `N` is vertices count
>  `M` is constraints 
>  `Constraint` is j ∈ [`1` … `M`]
>  Cardinality $N_{j}$
>  Function $C_j$ : $\mathbb{R}^{3nj}$ → $\mathbb{R}$
>  Indices {$i_1$ … $i_{N_j}$}, $i_k$ ∈ [0….1]

Constraint j with type equality is satisfied if $C_j$($x_{i_1}$ , . . . , $x_{i_{n_j}}$ ) = 0. If its type is inequality then it is satisfied if $C_j$($x_{i_1}$ , . . . , $x_{i_{n_j}}$ ) ≥ 0. The stiffness parameter $k_j$ defines the strength of the constraint in a range from zero to one.

based on this data and a timestep $\Delta$t the object is simulated as follows

The following is psuedo code based on the paper
```c
// 1. Initialization phase (This part is often done before the main loop)
for (int i = 0; i < num_vertices; i++) {
    xi[i] = x0[i];
    v[i] = v0[i];
    w[i] = 1.0 / mi[i];
}

// 4. Main simulation loop
while (1) {
    // 5. Velocity Update (Integration step)
    for (int i = 0; i < num_vertices; i++) {
        v[i] = v[i] + dt * twifext(xi[i]);
    }

    // 6. Damp Velocities
    dampVelocities(v1, ..., vN); // Assuming v1..vN are arrays/vectors

    // 7. Position Update
    for (int i = 0; i < num_vertices; i++) {
        pi[i] = xi[i] + dt * v[i];
    }

    // 8. Collision Constraint Generation
    for (int i = 0; i < num_vertices; i++) {
        generateCollisionConstraints(xi[i], pi[i]);
    }

    // 9. Solver Iterations
    for (int solver_iteration = 0; solver_iteration < times; solver_iteration++) {
        //10. Project Constraints
        projectConstraints(C1, C2, ..., CM+Mcoll, p1, p2, ..., pN);
    }

    // 12. Final Position Update for the iteration
    for (int i = 0; i < num_vertices; i++) {
        v[i] = (pi[i] - xi[i]) / dt;
        xi[i] = pi[i];
    }
    
    // 16. Final Velocity Update (If this is meant to be a final step after the main loop)
	VelocityUpdate(v1, v2, ..., vN);
}
```

The core of the PBD solver is 
```C
//7. Position Update
for (int i = 0; i < num_vertices; i++) {
    pi[i] = xi[i] + dt * v[i];
}

    // 8. Collision Constraint Generation
for (int i = 0; i < num_vertices; i++) {
    generateCollisionConstraints(xi[i], pi[i]);
}

    // 9. Solver Iterations
for (int solver_iteration = 0; solver_iteration < times; solver_iteration++) {
    // 10. Project Constraints
    projectConstraints(C1, C2, ..., CM+Mcoll, p1, p2, ..., pN);
}
```

---
```C
for (int i = 0; i < num_vertices; i++) {
	Pi[i] = xi[i] + dt * v[i];
}
```
The above estimates the new $P_i$ using explicit Euler integration.

---
```C
for (int solver_iteration = 0; solver_iteration < times; solver_iteration++) {
    // 10. Project Constraints
    projectConstraints(C1, C2, ..., CM+Mcoll, p1, p2, ..., pN);
}
```
The above is the iterative solver. Its job is to manipulate the Euler integration to satisfy the required constraints of the system. It achieves this by completing `Guessien Seidal` style solver repeatedly. It moves the predictions and velocities accordingly. This is pulled from `Verlet integration`

---

>[!quote]+ Muller et al
>“The scheme is unconditionally stable. This is because the integration steps (13) and (14) do not extrapolate blindly into the future as traditional explicit schemes do but move the vertices to a physically valid configuration pi computed by the constraint solver. The only possible source for instabilities is the solver itself which uses the Newton-Raphson method to solve for valid positions”
>
>“ increasing the number of iterations shifts the bottleneck from collision detection to the solver“

> [!important]- Gaussian Seidel 
> [[Gaussian Seidel]] Becomes a very important throughline.  The general Idea is that for matrices Ax = b you perform A * b iteratively starting a a 0 matrix. Each step you write to the 0 matrix making it non 0. You then perform the step by step multiplication for term 1-3 repeatedly until the $\Delta$ of the term is below a threshold for all three.

Gaussian Seidel is only for linear equations *however* the core idea of iteratively guessing repeatedly transfers to the solver used here.

> [!info]- $\nabla$ 
> This operator is used to indicate a [[Gradients]] calculation
> - **Mathematical Representation:** Expressed as a vector of partial derivatives:
>    
   > $$\nabla f = \left[ \frac{\partial f}{\partial x_1}, \frac{\partial f}{\partial x_2}, \dots, \frac{\partial f}{\partial x_n} \right]^T$$
 >   
>- **Geometric Behavior:** * Points in the direction of steepest ascent.
 >   
>    - Its magnitude signifies the maximum rate of change.
> 

# Solver Theory
---
$\nabla{_P}C$ gives us the direction of maximal positional change *relative* to the Constraint term. Using this we can get our positional change equation to include the Constraint. P($\Delta{P}$ + P) $\approx$ C(P) + $\nabla{C(p)}$ * $\Delta{P}$ = 0

This requires a $\lambda$ such that $\Delta{P}$ = $\lambda{\nabla{_p}}C_p$ holds true.

> [!important]- Figure 2
> ![[IMG_0302.jpeg]]

We are then able to resolve for $\Delta{P}$ after substituting $\Delta{P}$ = $\lambda{\nabla{_p}C(p)}$ into P($\Delta{P}$ + P) $\approx$ C(P) + $\nabla{C(p)}$ * $\Delta{P}$ = 0 solving for $\lambda$ and substituting back into $\Delta{P}$ = $\lambda{\nabla{_p}C(p)}$ gives us $\Delta{P}$ = -$\frac{C(p)}{|\nabla{_p}C(p)|^2}{\nabla{p}}C(p)$ 

Finishing up with how the `k` constant of stiffness is affected it is applied to the $\Delta{P}$ if C($P_1$…$P_n$) < 0. To avoid solver buildup multiplying by k’ not k generates a lower error and is a linear relationship because k’= 1 - $(1-k)^{\frac{1}{n_s}}$ minimizing the error propagation over time

# Collisions
---
For continuous style of collisions we simply use raycast entry de-projection relying on ray detection to shift to static calculations if continuous fails by computing the entry point $q_c$ and the surface normal $n_c$ at this position. Creating an inequality constraint with constraint function C(p) = (p− $q_c$) * $N_c$ and stiffness k= 1 is added to the list of constraints.

# XPBD  and the short faults of PBD
---
while PBD created a stable simulation environment that could handle soft body without sub-step iterations it was still reliant on an iterative solver that could be slow in many scenarios. Additionally the solvers stiffnes was timestep **AND** solver dependent.

# The XPBD changes
---

XPBD starts from Newtons equations specifically for the energy potential behavior `(U(x))`
> [!note]- Servin et al. (2006)
> this paper is noted as the way we move from F to $\lambda$ 

* Starting with Mx = -$\nabla{U^t}$(x)
* we then discretize it M($\frac{x^{n+1}-2x^n+x^{n-1}}{\Delta{t^2}}$) = -$\nabla{U^yt}$($x^{n+1}$) (4)
* Establishing U(x): U(x) = $\frac{1}{2}$$C(x)^T{a^{-1}}C(x)$ (5)
* We can use the negative gradient to find $F_{elastic}$ = -$\nabla{_X}U^T$ = $-\nabla{C^T}a^{-1}{C}$ (6)
* $\lambda{_{elastic}}=-a^{-1}{C(x)}$ (7)
* M($x^{n+1}$ - x) = $\nabla$C$(X^{n+1})^T{\lambda}^{n+1}$ = 0 (8)
* $C(x^{n+1})^{T}{\lambda}^{n+1}$ = 0 (9)
* 8 and 9 become g and h respectively meaning we can create the below equation for solving
* $$\begin{bmatrix} \mathbf{K} & -\nabla\mathbf{C}^T(\mathbf{x}_i) \\ \nabla\mathbf{C}(\mathbf{x}_i) & \tilde{\boldsymbol{\alpha}} \end{bmatrix} \begin{bmatrix} \Delta\mathbf{x} \\ \Delta\boldsymbol{\lambda} \end{bmatrix} = -\begin{bmatrix} \mathbf{g}(\mathbf{x}_i, \boldsymbol{\lambda}_i) \\ \mathbf{h}(\mathbf{x}_i, \boldsymbol{\lambda}_i) \end{bmatrix} \qquad (12)$$

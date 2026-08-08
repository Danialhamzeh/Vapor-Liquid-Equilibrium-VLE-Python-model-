# VLE Python Model

A Python-based thermodynamic solver for Vapor-Liquid Equilibrium (VLE) modeling of the quaternary system: $CO_2 - N_2 - Ethanol - Water$.

## Thermodynamic Framework

*   **Vapor Phase:** Peng-Robinson Equation of State (PR-EOS).
*   **Liquid Phase:** 
    *   Activity coefficient model for condensable components ($Ethanol$, $Water$).
    *   Henry's Law for supercritical/non-condensable components ($CO_2$, $N_2$).

## Phase Fugacity Formulation

*   **Vapor Fugacity ($\hat{f}_i^V$):** Calculated directly via PR-EOS.
*   **Liquid Fugacity ($\hat{f}_i^L$):** 
    *   Condensables: $\hat{f}_i^L = x_i \gamma_i f_i^{L,pure}$
    *   The pure liquid fugacity ($f_i^{L,pure}$) is computed within the vapor phase subroutine using the saturation fugacity coefficient and Poynting correction factor, then passed to the liquid phase subroutine.

## Solution Algorithm

1.  **Initialization:** Provide initial guesses for equilibrium ratios ($K_i$).
2.  **Rachford-Rice Solution:** Solve the Rachford-Rice equation to determine vapor fraction ($V$) and mole fractions ($x_i$, $y_i$).
3.  **Vapor Phase Subroutine:**
    *   Input: $y_i$, $P$, $T$
    *   Compute vapor fugacities: $\hat{f}_i^V$.
    *   Compute pure liquid fugacities for condensables: $f_i^{L,pure} = \phi_i^{sat} P_i^{sat} \exp \left( \frac{V_i^L (P - P_i^{sat})}{RT} \right)$.
4.  **Liquid Phase Subroutine:**
    *   Input: $x_i$, $T$, and $f_i^{L,pure}$ (imported from the vapor subroutine).
    *   Compute activity coefficients ($\gamma_i$).
    *   Compute liquid fugacities: $\hat{f}_i^L$.
5.  **Convergence Check:** Evaluate the objective function based on the isofugacity criterion:
    $$ Obj = \sum_{i} \left( 1 - \frac{\hat{f}_i^L}{\hat{f}_i^V} \right)^2 < 10^{-12} $$
6.  **K-Value Update:** If the convergence criterion is not met, update $K_i$ values and iterate from Step 2 until $\hat{f}_i^L = \hat{f}_i^V$ for all components.

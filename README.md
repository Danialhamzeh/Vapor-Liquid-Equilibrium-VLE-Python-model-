# VLE Python Model

A Python-based thermodynamic solver for Vapor-Liquid Equilibrium (VLE) modeling of the quaternary system: $CO_2 - N_2 - Ethanol - Water$. 
The final results have been validated against Aspen HYSYS.

## Thermodynamic Framework

*   **Vapor Phase:** Peng-Robinson Equation of State (PR-EOS).
*   **Liquid Phase:** 
    *   Activity coefficient model for condensable components ($Ethanol$, $Water$).
    *   Henry's Law for supercritical/non-condensable components ($CO_2$, $N_2$).

## Phase Fugacity Formulation

*   **Vapor Fugacity ($\hat{f}_i^V$):** 
    Calculated at system pressure $P$ (e.g., 50 bar):

$$ \hat{f}_i^V = y_i \phi_i P $$

*   **Liquid Fugacity ($\hat{f}_i^L$):** 
    *   **Condensables:** 

$$ \hat{f}_i^L = x_i \gamma_i \phi_i^{sat} P_i^{sat} \exp \left( \frac{V_i^L (P - P_i^{sat})}{RT} \right) $$

    *   **Non-condensables:** 
        Using Henry's Law with Poynting correction ($\bar{V}_i^{\infty}$ is the partial molar volume at infinite dilution):

$$ \hat{f}_i^L = x_i H_i \exp \left( \frac{\bar{V}_i^{\infty} (P - P_{solvent}^{sat})}{RT} \right) $$

## Solution Algorithm

1.  **Initialization (Highly Sensitive):** 
    The problem is highly sensitive to the initial guess for equilibrium ratios ($K_i$). At the operating conditions (50 °C, 50 bar):
    *   $CO_2$ & $N_2$: Low solubility in the liquid phase (mostly vapor). Initial guess must be $K_i > 1$.
    *   $Ethanol$ & $Water$: High solubility in the liquid phase. Initial guess must be $K_i < 1$.
2.  **Rachford-Rice Solution:** Solve the Rachford-Rice equation to determine vapor fraction ($V$) and mole fractions ($x_i$, $y_i$).
3.  **Vapor Phase Subroutine:** Compute vapor fugacities ($\hat{f}_i^V$) and pure liquid fugacities for condensables.
4.  **Liquid Phase Subroutine:** Compute activity coefficients ($\gamma_i$) and liquid fugacities ($\hat{f}_i^L$).
5.  **Convergence Check:** Evaluate the objective function based on the isofugacity criterion:

$$ Obj = \sum_{i} \left( 1 - \frac{\hat{f}_i^L}{\hat{f}_i^V} \right)^2 < 10^{-12} $$

6.  **K-Value Update:** If the convergence criterion is not met, update $K_i$ values and iterate from Step 2 until phase equilibrium is reached.


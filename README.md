# Thermodynamic VLE Model: CO2-N2-EtOH-H2O System

A Python-based thermodynamic solver for Vapor-Liquid Equilibrium (VLE) modeling of the quaternary system: $CO_2 (1) - N_2 (2) - EtOH (3) - H_2O (4)$. The custom model results have been rigorously validated against Aspen HYSYS.

## Operating Conditions & Feed Specifications

*   **Temperature ($T$):** 50 °C
*   **Pressure ($P$):** 50 bar
*   **Feed Mole Fractions ($z_i$):** 
    *   $z_{CO_2} = 0.25$
    *   $z_{N_2} = 0.05$
    *   $z_{EtOH} = 0.40$
    *   $z_{H_2O} = 0.30$

## Thermodynamic Framework

*   **Vapor Phase:** Soave-Redlich-Kwong Equation of State (SRK-EOS).
*   **Liquid Phase:** 
    *   Extended NRTL Activity Coefficient Model for condensable components.
    *   Henry's Law for supercritical/non-condensable components.

## Phase Fugacity Formulation

### Vapor Fugacity ($\hat{f}_i^V$)
Calculated at system pressure $P$:

$$ \hat{f}_i^V = y_i \phi_i P $$

### Liquid Fugacity ($\hat{f}_i^L$)

**Condensables ($EtOH$, $H_2O$):**

$$ \hat{f}_i^L = x_i \gamma_i \phi_i^{sat} P_i^{sat} \exp \left( \frac{V_i^L (P - P_i^{sat})}{RT} \right) $$

**Non-condensables ($CO_2$, $N_2$):**

$$ \hat{f}_i^L = x_i H_i \exp \left( \frac{\bar{V}_i^{\infty} (P - P_{solvent}^{sat})}{RT} \right) $$

**Mixed Solvent Henry's Constant:**
For the mixed solvent system, the effective Henry's constant ($H_i$) for non-condensable gases is calculated using a logarithmic mixing rule:

$$ \ln(H_i) = x_{EtOH} \ln(H_{i, EtOH}) + x_{H_2O} \ln(H_{i, H_2O}) $$

## Solution Algorithm

1.  **Initialization (Highly Sensitive):** 
    The problem is highly sensitive to the initial guess for equilibrium ratios ($K_i$).
    *   $CO_2$ & $N_2$: Low solubility in the liquid phase. Initial guess must be $K_i > 1$.
    *   $EtOH$ & $H_2O$: High solubility in the liquid phase. Initial guess must be $K_i < 1$.
2.  **Rachford-Rice Solution:** Solve the Rachford-Rice equation to determine vapor fraction ($\beta$) and mole fractions ($x_i$, $y_i$).
3.  **Vapor Phase Subroutine:** Compute vapor fugacities ($\hat{f}_i^V$) using SRK-EOS.
4.  **Liquid Phase Subroutine:** Compute activity coefficients ($\gamma_i$) via Extended NRTL and liquid fugacities ($\hat{f}_i^L$).
5.  **Convergence Check:** Evaluate the objective function based on the isofugacity criterion:

$$ Obj = \sum_{i} \left( 1 - \frac{\hat{f}_i^L}{\hat{f}_i^V} \right)^2 < 10^{-12} $$

6.  **K-Value Update:** If the convergence criterion is not met, update $K_i$ values and iterate from Step 2 until phase equilibrium is reached.
7.  **HYSYS Validation:** Once equilibrium is achieved, the Python outputs are directly compared against Aspen HYSYS simulation results. Key evaluated parameters include Vapor Fraction ($\beta$), K-values ($K_i$), Phase Compositions ($x_i, y_i$), Compressibility Factors ($Z_{mix}, Z_{pure,i}$), and Henry's Constants ($H_i$).
8.  **Error Analysis:** The model accuracy is quantified using (ARE):

$$ ARE (\%) = \left| \frac{Value_{Python} - Value_{HYSYS}}{Value_{Python}} \right| \times 100 $$

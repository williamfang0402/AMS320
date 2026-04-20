# AMS 320 — Chapter 14 Study Guide
## Wiener Processes, Itô's Lemma & HW5 Solutions

---

## 📚 Table of Contents

1. [Core Concepts](#core-concepts)
2. [Key Formulas](#key-formulas)
3. [Common Pitfalls](#common-pitfalls)
4. [Homework Solutions](#homework-solutions)
5. [Quick Reference Card](#quick-reference-card)

---

## Core Concepts

### 1. The Markov Property

A stochastic process is **Markov** if only the *current* value matters for predicting the future — the past is irrelevant.

- Stock prices are assumed Markov.
- Consistent with the **weak form of market efficiency**: current price reflects all info in past prices.

### 2. Wiener Process (Brownian Motion)

A variable $z$ follows a **Wiener process** if:
- $\Delta z = \epsilon\sqrt{\Delta t}$, where $\epsilon \sim N(0,1)$
- $\Delta z$ values over non-overlapping intervals are **independent**

**Key consequence:** $\Delta z \sim N(0, \Delta t)$

**Why $\sqrt{\Delta t}$, not $\Delta t$?** Randomness partially cancels — like a drunk taking 100 random steps only wanders ~10 feet. This is why **financial volatility scales with $\sqrt{T}$**, not $T$.

### 3. Generalized Wiener Process

$$dx = a\,dt + b\,dz$$

- $a$ = **drift rate** (mean change per unit time)
- $b^2$ = **variance rate** (variance per unit time)
- ⚠️ $b$ itself is NOT the variance rate — $b^2$ is

**Metaphor:** A sailboat crossing the ocean — $a\,dt$ is the steady current, $b\,dz$ is the random waves.

### 4. Itô Process

$$dx = a(x,t)\,dt + b(x,t)\,dz$$

Drift and volatility can depend on the current state $x$ and time $t$. Still Markov.

### 5. Geometric Brownian Motion (GBM) — Stock Prices

$$dS = \mu S\,dt + \sigma S\,dz$$

- $\mu$ = expected rate of return
- $\sigma$ = volatility
- Percentage changes (not absolute changes) are constant in expectation.

### 6. Discrete vs. Continuous Notation

| Symbol | Meaning |
|--------|---------|
| $\Delta z$ | Finite change over finite interval $\Delta t$ (use for simulation) |
| $dz$ | Infinitesimal change as $\Delta t \to 0$ (use for SDEs) |

Technical subtlety: $(dz)^2 = dt$ (not zero!) — this is why Itô's Lemma has an extra term.

### 7. Correlated Wiener Processes

If $dz_1$ and $dz_2$ have correlation $\rho$:
$$\text{Cov}(\Delta z_1, \Delta z_2) = \rho\,\Delta t$$

To simulate: $\epsilon_1 = u$, $\epsilon_2 = \rho u + \sqrt{1-\rho^2}\,v$ where $u, v \sim N(0,1)$ independent.

### 8. Itô's Lemma — "Chain Rule for Randomness"

Given $dx = a(x,t)\,dt + b(x,t)\,dz$, any function $G(x,t)$ follows:

$$dG = \left(\frac{\partial G}{\partial t} + a\frac{\partial G}{\partial x} + \frac{1}{2}b^2\frac{\partial^2 G}{\partial x^2}\right)dt + b\frac{\partial G}{\partial x}\,dz$$

The extra term $\frac{1}{2}b^2\frac{\partial^2 G}{\partial x^2}$ is the **Itô correction** — it's what ordinary calculus misses.

**Why does it appear?** Because $(\Delta x)^2$ has a piece of order $\Delta t$ (from $b^2\epsilon^2\Delta t$), which survives as $\Delta t \to 0$.

**Jensen's correction intuition:**
- **Convex** $G$ ($G'' > 0$): noise pushes drift **up** 📈
- **Concave** $G$ ($G'' < 0$): noise pushes drift **down** 📉
- **Linear** $G$ ($G'' = 0$): no correction — ordinary chain rule works

---

## Key Formulas

### Distribution Results

| Quantity | Distribution |
|----------|--------------|
| $\Delta z$ over interval $\Delta t$ | $N(0, \Delta t)$ |
| Change $x_T - x_0$ in generalized Wiener | $N(aT, b^2 T)$ |
| Value $x_T$ in generalized Wiener | $N(x_0 + aT, b^2 T)$ |
| $\ln S_T$ under GBM | $N\left(\ln S_0 + (\mu - \tfrac{1}{2}\sigma^2)T,\; \sigma^2 T\right)$ |

### Variance Identity

$$\text{Var}(\sigma\,dz) = \sigma^2\,dt$$

### Combining Processes (Key for 14.4)

For $X = X_1 + X_2$ where each follows a generalized Wiener:

| Case | Drift | Variance Rate | Volatility |
|------|-------|---------------|------------|
| Independent | $\mu_1 + \mu_2$ | $\sigma_1^2 + \sigma_2^2$ | $\sqrt{\sigma_1^2 + \sigma_2^2}$ |
| Correlated (ρ) | $\mu_1 + \mu_2$ | $\sigma_1^2 + \sigma_2^2 + 2\rho\sigma_1\sigma_2$ | $\sqrt{\sigma_1^2 + \sigma_2^2 + 2\rho\sigma_1\sigma_2}$ |

**Sanity checks for correlated case:**
- $\rho = 0$: reduces to independent case ✓
- $\rho = +1$: $\sigma = \sigma_1 + \sigma_2$ (volatilities add linearly — perfect co-movement)
- $\rho = -1$: $\sigma = |\sigma_1 - \sigma_2|$ (partial cancellation — perfect hedge)

### Itô's Lemma Applied to GBM

For $dS = \mu S\,dt + \sigma S\,dz$:

$$dG = \left(\frac{\partial G}{\partial t} + \mu S\frac{\partial G}{\partial S} + \frac{1}{2}\sigma^2 S^2\frac{\partial^2 G}{\partial S^2}\right)dt + \sigma S\frac{\partial G}{\partial S}\,dz$$

### Common Function Results

| Function | Process |
|----------|---------|
| $G = \ln S$ | $dG = (\mu - \tfrac{1}{2}\sigma^2)\,dt + \sigma\,dz$ (generalized Wiener!) |
| $G = S^n$ | $dG = [\mu n + \tfrac{1}{2}\sigma^2 n(n-1)]\,G\,dt + \sigma n\,G\,dz$ (still GBM!) |
| $G = 2S$ | $dG = \mu G\,dt + \sigma G\,dz$ (same GBM params) |
| $G = S^2$ | $dG = (2\mu + \sigma^2)\,G\,dt + 2\sigma G\,dz$ |

---

## Common Pitfalls

⚠️ **The Top 5 Ways to Lose Points**

1. **Variance vs. Standard Deviation**
   - Variance is $b^2 T$, std dev is $b\sqrt{T}$
   - Z-score formula uses **std dev** in denominator, not variance

2. **Change vs. Value**
   - Change formula: $N(aT, b^2 T)$ — no $x_0$!
   - Value formula: $N(x_0 + aT, b^2 T)$ — add $x_0$!

3. **Sign errors in Itô's Lemma**
   - Especially when $G$ involves $1/S$ or $e^{-something}$
   - Always double-check signs in partial derivatives

4. **Exponent rules require same base**
   - $S^2 \cdot S^{n-2} = S^n$ ✓ (same base)
   - $n(n-1) \cdot S^n$ stays as is — you CANNOT add $n$-type exponents to $S$-type exponents!

5. **Confusing $b$ with $b^2$ in SDE**
   - Problem says "variance rate 4.0" → $b^2 = 4$, so $b = 2$ in the SDE
   - Problem says "volatility 0.30" → $b = 0.30$ directly

---

## Homework Solutions

### Problem 14.3 — Cash Position Probability

**Setup:** Generalized Wiener with drift $a = 0.5$/quarter, variance rate $b^2 = 4$/quarter. $T = 4$ quarters (1 year).

**Distribution:** $x_T \sim N(x_0 + 2, 16)$, so std dev = 4.

**Standardize:**
$$P(x_T < 0) = P\left(Z < \frac{-(x_0 + 2)}{4}\right) < 0.05$$

**Solve:** $\frac{-(x_0 + 2)}{4} < -1.645 \implies x_0 > 4.58$

**Answer:** Initial cash must exceed **$4.58 million**.

---

### Problem 14.4 — Sums of Wiener Processes

Given $dX_1 = \mu_1\,dt + \sigma_1\,dz_1$ and $dX_2 = \mu_2\,dt + \sigma_2\,dz_2$.

**(a) Uncorrelated:**
$$d(X_1 + X_2) = (\mu_1 + \mu_2)\,dt + \sqrt{\sigma_1^2 + \sigma_2^2}\,dz$$

**(b) Correlation ρ:**
$$d(X_1 + X_2) = (\mu_1 + \mu_2)\,dt + \sqrt{\sigma_1^2 + \sigma_2^2 + 2\rho\sigma_1\sigma_2}\,dz$$

**Derivation (b):** Match variance of $\sigma_1\Delta z_1 + \sigma_2\Delta z_2$ to $\sigma\Delta z$:
$$\text{Var}(\sigma_1\Delta z_1 + \sigma_2\Delta z_2) = (\sigma_1^2 + \sigma_2^2 + 2\rho\sigma_1\sigma_2)\Delta t$$
Setting equal to $\sigma^2\Delta t$ and solving gives $\sigma$.

**Big Picture:** This is the foundation of Modern Portfolio Theory — diversification benefits depend on correlation.

---

### Problem 14.5 — Piecewise Wiener Process

**Setup:** $dS = \mu\,dt + \sigma\,dz$ with $S_0 = 5$.
- Years 0–3: $\mu = 2$, $\sigma = 3$
- Years 3–6: $\mu = 3$, $\sigma = 4$

**Strategy:** Decompose $S_6 = S_0 + (S_3 - S_0) + (S_6 - S_3)$. Each change is independent.

**Individual changes:**
- $S_3 - S_0 \sim N(aT, b^2 T) = N(2 \cdot 3,\; 9 \cdot 3) = N(6, 27)$
- $S_6 - S_3 \sim N(3 \cdot 3,\; 16 \cdot 3) = N(9, 48)$

**Combine (add means, add variances since independent, add $S_0$ as constant):**
$$S_6 \sim N(5 + 6 + 9,\; 27 + 48) = \boxed{N(20, 75)}$$

---

### Problem 14.10 — Process for $S^n$

**Partials of $y = S^n$:**
$$\frac{\partial y}{\partial t} = 0, \quad \frac{\partial y}{\partial S} = nS^{n-1}, \quad \frac{\partial^2 y}{\partial S^2} = n(n-1)S^{n-2}$$

**Apply Itô with $a = \mu S$, $b = \sigma S$:**
$$dy = \left(\mu S \cdot nS^{n-1} + \tfrac{1}{2}\sigma^2 S^2 \cdot n(n-1)S^{n-2}\right)dt + \sigma S \cdot nS^{n-1}\,dz$$

**Simplify** (combining exponents: $S^2 \cdot S^{n-2} = S^n$):
$$\boxed{dy = \left[\mu n + \tfrac{1}{2}\sigma^2 n(n-1)\right]y\,dt + \sigma n\,y\,dz}$$

**Conclusion:** $S^n$ follows GBM with drift $\mu n + \tfrac{1}{2}\sigma^2 n(n-1)$ and volatility $\sigma n$.

---

### Problem 14.16 — Itô's Lemma on Various Functions

**(a) $y = 2S$:**
- Partials: $0, 2, 0$
- Itô correction vanishes (linear function!)
- **Answer:** $dy = \mu y\,dt + \sigma y\,dz$ (same GBM as $S$)

**(b) $y = S^2$:**
- Partials: $0, 2S, 2$
- **Answer:** $dy = (2\mu + \sigma^2)\,y\,dt + 2\sigma y\,dz$

**(c) $y = e^S$:**
- Partials: $0, e^S, e^S$
- **Answer:** $dy = \left(\mu S + \tfrac{1}{2}\sigma^2 S^2\right)y\,dt + \sigma S\,y\,dz$
- Note: Coefficients depend on $S$ — this is NOT standard GBM

**(d) $y = e^{r(T-t)}/S = e^{r(T-t)} \cdot S^{-1}$:**
- Partials (rewritten in terms of $y$):
  - $\frac{\partial y}{\partial t} = -ry$ (chain rule: $\frac{d}{dt}(T-t) = -1$)
  - $\frac{\partial y}{\partial S} = -y/S$
  - $\frac{\partial^2 y}{\partial S^2} = 2y/S^2$
- **Answer:** $dy = (-r - \mu + \sigma^2)\,y\,dt - \sigma y\,dz$

---

### Problem 14.17 — Stock Price Probability (Capstone)

**Given:** $S_0 = 50$, $\mu = 0.12$, $\sigma = 0.30$, $T = 2$. Find $P(S_T > 80)$.

**Step 1 — Apply Itô's Lemma to $G = \ln S$:**
$$d(\ln S) = (\mu - \tfrac{1}{2}\sigma^2)\,dt + \sigma\,dz$$

**Step 2 — Distribution of $\ln S_T$:**
$$\ln S_T \sim N\left(\ln 50 + (0.12 - 0.045)(2),\; (0.09)(2)\right) = N(4.062,\; 0.18)$$

**Step 3 — Convert to log question:**
$$P(S_T > 80) = P(\ln S_T > \ln 80) = P(\ln S_T > 4.382)$$

**Step 4 — Standardize:**
$$Z = \frac{4.382 - 4.062}{\sqrt{0.18}} = \frac{0.320}{0.4243} \approx 0.7543$$

**Step 5 — Look up probability:**
$$P(Z > 0.7543) = 1 - P(Z < 0.7543) = 1 - 0.7747 = \boxed{0.2253}$$

**Answer:** About **22.53%** probability that $S_T > 80$.

**Sanity check:** Expected value after 2 years is $50 \cdot e^{0.24} \approx 63.56$. Reaching 80 exceeds the expected value by a reasonable amount given 30% annual volatility — roughly 1-in-5 odds feels right.

---

## Quick Reference Card

### The 4-Step Recipe for Itô's Lemma

1. **Compute** three partial derivatives: $\frac{\partial G}{\partial t}, \frac{\partial G}{\partial x}, \frac{\partial^2 G}{\partial x^2}$
2. **Identify** $a$ and $b$ from $dx = a\,dt + b\,dz$
3. **Plug** into Itô's formula
4. **Simplify** and (if asked) express coefficients in terms of $G$

### The 5-Step Recipe for Wiener Probability Problems

1. Identify $a$, $b^2$, $T$ (check units!)
2. Write distribution: $x_T \sim N(x_0 + aT,\; b^2 T)$
3. Set up the probability statement (convert to log if GBM!)
4. Standardize: $Z = \frac{x_T - \mu}{\sigma}$ (note: $\sigma$ = std dev, not variance)
5. Look up in standard normal table; use $P(Z > z) = 1 - P(Z < z)$ if needed

### Critical Standard Normal Values

| Tail Probability | Z-value |
|:---:|:---:|
| $P(Z < z) = 0.05$ | $z = -1.645$ |
| $P(Z < z) = 0.95$ | $z = +1.645$ |
| $P(Z < z) = 0.025$ | $z = -1.960$ |
| $P(Z < z) = 0.01$ | $z = -2.326$ |

### Mental Checklist Before Submitting

- [ ] Did I use variance $(b^2 T)$ or std dev $(b\sqrt{T})$ where appropriate?
- [ ] Did I add $x_0$ if asked for value (not change)?
- [ ] Did I include the Itô correction for non-linear $G$?
- [ ] Did I double-check all signs?
- [ ] Did I simplify coefficients in terms of the new variable if asked?
- [ ] Does my answer pass a sanity check?

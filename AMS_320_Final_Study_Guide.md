# AMS 320 Final Exam Study Guide

**Course:** AMS 320 — Quantitative Finance, Stony Brook University
**Textbook:** Hull, *Options, Futures, and Other Derivatives*, 9th Edition

---

## Table of Contents

1. [Exam Structure & Strategy](#exam-structure--strategy)
2. [Forwards, Futures & Cost of Carry](#1-forwards-futures--cost-of-carry-problems-12)
3. [Bonds: Zero Rates, Bootstrapping, Par Yields, Forward Rates](#2-bonds-zero-rates-bootstrapping-par-yields-forward-rates-problems-35)
4. [Option Strategies & No-Arbitrage Inequalities](#3-option-strategies--no-arbitrage-inequalities-problems-67)
5. [Binomial Trees](#4-binomial-trees-problems-89)
6. [Volatility Estimation](#5-volatility-estimation-problem-10)
7. [Wiener Processes & Itô's Lemma](#6-wiener-processes--itôs-lemma-problems-11-12)
8. [Black-Scholes-Merton Pricing](#7-black-scholes-merton-pricing-problems-13-14)
9. [Hull 15.12: Derivative Paying S^n](#8-hull-1512-derivative-paying-s_tn-final-problem-11)
10. [Capped Options](#9-capped-options-final-problem-10)
11. [Master Formula Sheet](#master-formula-sheet)
12. [Common Traps & Final Tips](#common-traps--final-tips)

---

## Exam Structure & Strategy

- **14 problems total, 160 points** (graded out of 120 — anything above is extra credit)
- Problems 1–10: 10 points each
- Problems 11–14: 15 points each
- Problems 1–7: similar to practice final
- Problems 8–9: Itô's Lemma (similar to Midterm 2)
- Problem 10: Capped option payoff replication
- Problem 11: Hull 15.12 ($S_T^n$ derivative)
- Problems 12–14: intentionally more challenging (within scope)

**Strategy:** You only need 120/160 for full credit. Nail Problems 1–11 (95 pts) and pick up partial credit on 12–14.

---

## 1. Forwards, Futures & Cost of Carry (Problems 1–2)

### Cost-of-Carry Formula

For an asset with continuous dividend yield $q$:

$$F_0 = S_0 \cdot e^{(r-q)T}$$

**Variants by underlying:**
- Non-dividend-paying stock: $q = 0 \Rightarrow F_0 = S_0 e^{rT}$
- Stock index: $q$ = continuous dividend yield
- Currency: $q$ = foreign risk-free rate
- Commodity with storage cost $u$: $F_0 = S_0 e^{(r+u)T}$

**Intuition:** Higher $r$ → higher $F_0$ (financing costs more). Higher $q$ → lower $F_0$ (you're paid to hold).

### Forward Price vs. Forward Value

These are **two different concepts** — never confuse them.

| Concept | Symbol | Formula |
|---------|--------|---------|
| Forward **price** (delivery price) | $F$ | $F = S e^{(r-q)T}$ |
| Forward contract **value** | $f$ | $f = (F - K) e^{-rT}$ |

where $K$ is the **originally agreed** delivery price.

**At initiation:** $f = 0$ by construction (this defines $F_0$).

**Equivalent for non-dividend stock:** $f = S - K e^{-rT}$

### Worked Example (Practice Problem 1)

Given: $S_0 = 100$, $F_0 = 105$, $r = 0.1864$, $T = 1/3$. Find $q$.

$$q = r - \frac{1}{T}\ln\left(\frac{F_0}{S_0}\right) = 0.1864 - 3\ln(1.05) = 0.0400$$

**Answer:** $q = 4.00\%$

### Worked Example (Practice Problem 2)

(a) Stock at \$100, $r = 8\%$, $T = 1$ year, no dividends.
- $F_0 = 100 e^{0.08} = \$108.33$
- $f_0 = \$0$

(b) 9 months later, $S = 110$, 3 months left.
- New forward price: $F = 110 e^{0.08 \times 0.25} = \$112.22$
- Value of original contract: $f = (112.22 - 108.33) e^{-0.08 \times 0.25} = \$3.82$

**Common trap:** Forgetting to discount the value back to today.

---

## 2. Bonds: Zero Rates, Bootstrapping, Par Yields, Forward Rates (Problems 3–5)

### Zero Rates via Bootstrapping

**Step 1 — Use zero-coupon bonds (T-bills) first:**

$$\text{Price} = \text{Face} \cdot e^{-RT} \implies R = -\frac{1}{T}\ln\left(\frac{\text{Price}}{\text{Face}}\right)$$

**Step 2 — Then coupon bonds, in order of maturity:**

For an $n$-period bond with coupon $C$ and price $P$:

$$P = C e^{-R_1 t_1} + C e^{-R_2 t_2} + \dots + (C + \text{Face}) e^{-R_n t_n}$$

All $R_i$ for $i < n$ are known. Solve for $R_n$:

$$R_n = -\frac{1}{t_n}\ln\left(\frac{P - \sum_{i<n} C e^{-R_i t_i}}{C + \text{Face}}\right)$$

**Critical reminders:**
- Last cash flow is **coupon + face value**, not just the coupon
- Solve in order: shortest maturity → longest
- Each new bond gives exactly **one** new unknown — that's the magic

### Par Yield

**Definition:** The coupon rate $c$ that makes price = face value.

$$c = \frac{(\text{Face} - \text{Face} \cdot d_n) \cdot m}{\sum_{i=1}^n d_i}$$

where:
- $d_i = e^{-R_i t_i}$ is the discount factor at time $t_i$
- $m$ = coupons per year (semiannual = 2)
- $d_n$ = discount factor at final maturity

**Setup from scratch (if you forget the formula):**
$$\text{Face} = \frac{c}{m}\sum_i d_i + \text{Face} \cdot d_n$$

then isolate $c$.

**Sanity check:** Par yield sits near the longest-maturity zero rate (the principal at maturity is the dominant cash flow).

### Forward Rates

The forward rate $R_F$ between times $T_1$ and $T_2$ ($T_2 > T_1$) is:

$$R_F = \frac{R_2 T_2 - R_1 T_1}{T_2 - T_1}$$

**Intuition:** It's a "marginal" rate — what rate is the curve implying between $T_1$ and $T_2$?

**Rule of thumb:** When zero rates are linear in maturity, forward rates rise *twice as fast*.

### Compounding Conversions

**Continuous → $m$-times-per-year:**
$$R_m = m\left(e^{R_c/m} - 1\right)$$

**$m$-times-per-year → Continuous:**
$$R_c = m\ln\left(1 + \frac{R_m}{m}\right)$$

---

## 3. Option Strategies & No-Arbitrage Inequalities (Problems 6–7)

### Strangle (Problem 6)

**Setup:** Long call (high strike $K_2$) + long put (low strike $K_1$). Bets on volatility.

**Three regions:**

| Range of $S_T$ | Call Payoff | Put Payoff | Profit |
|---------------|-------------|------------|--------|
| $S_T \leq K_1$ | $0$ | $K_1 - S_T$ | $K_1 - S_T - \text{premium}$ |
| $K_1 < S_T \leq K_2$ | $0$ | $0$ | $-\text{premium}$ (max loss) |
| $S_T > K_2$ | $S_T - K_2$ | $0$ | $S_T - K_2 - \text{premium}$ |

**Breakevens:**
- Lower: $S_T = K_1 - \text{premium}$
- Upper: $S_T = K_2 + \text{premium}$

**Workflow for any spread:** Identify each leg's payoff → sum payoffs → subtract net premium paid → sketch profit, find breakevens.

### No-Arbitrage Inequalities (Problem 7)

**The universal proof technique** for any "show that $X \leq Y$" on option prices:

1. **Construct a portfolio** whose initial cost equals the inequality
2. **Show its payoff at expiration is $\geq 0$ in every state**
3. **Conclude initial cost must be $\geq 0$** (no-arbitrage)

### Convexity of Call Prices

For $K_3 > K_2 > K_1$ with $K_3 - K_2 = K_2 - K_1$:

$$c_2 \leq \frac{1}{2}(c_1 + c_3)$$

**Proof:** Long 1 call at $K_1$ + Long 1 call at $K_3$ − 2 calls at $K_2$ (a **butterfly spread**). This portfolio has nonnegative payoff in every state, so its initial cost $c_1 + c_3 - 2c_2 \geq 0$.

**Same result for puts:** $p_2 \leq 0.5(p_1 + p_3)$ (use a put butterfly).

---

## 4. Binomial Trees (Problems 8–9)

### Cox-Ross-Rubinstein Parameters

If volatility $\sigma$ is given:
$$u = e^{\sigma\sqrt{\Delta t}}, \quad d = 1/u, \quad \tilde{p} = \frac{e^{(r-q)\Delta t} - d}{u - d}$$

If $u, d$ are given directly (like Problem 9), use them as-is.

### Backward Induction

At each non-terminal node:

$$f = e^{-r\Delta t}\left[\tilde{p} f_{\text{up}} + (1-\tilde{p}) f_{\text{down}}\right]$$

### European Option Workflow (Problem 8)

1. Build stock price tree forward from $S_0$
2. Compute terminal payoffs at leaves: $\max(S_T - K, 0)$ for call, $\max(K - S_T, 0)$ for put
3. Backward-induct using the formula above
4. Final answer is the option value at the root

### American Option Twist (Problem 9)

At every non-terminal node, take the max of:
- **Continuation value:** $e^{-r\Delta t}[\tilde{p} f_{\text{up}} + (1-\tilde{p}) f_{\text{down}}]$
- **Immediate exercise value:** $\max(K - S, 0)$ for put, $\max(S - K, 0)$ for call

$$f_{\text{node}} = \max(\text{continuation}, \text{exercise})$$

**Critical theoretical facts:**
- **American calls on non-dividend stocks:** Never optimal to exercise early. Treat as European.
- **American puts:** Often optimal to exercise early when deep in the money.

### Delta (Hedge Ratio)

At any node:
$$\Delta = \frac{f_{\text{up}} - f_{\text{down}}}{S_{\text{up}} - S_{\text{down}}}$$

**Hedge sign rules:**
- Long calls → **short** $\Delta$ shares
- Short calls → **buy** $\Delta$ shares
- Long puts → delta is negative → **buy** stock
- Short puts → **short** stock

**Common trap:** At each backward step, use the previously computed **option values**, not the original terminal payoffs.

---

## 5. Volatility Estimation (Problem 10)

### Recipe

1. Compute log returns: $u_i = \ln(S_i / S_{i-1})$
2. Compute sample standard deviation $s$ with divisor $n-1$:
$$s = \sqrt{\frac{1}{n-1}\sum_{i=1}^n (u_i - \bar{u})^2}$$
3. Annualize: $\hat{\sigma} = s / \sqrt{\tau}$ where $\tau$ is the period length in years

### Annualization Factors (Square-Root-of-Time Rule)

| Data frequency | Periods/year | Multiplier |
|---|---|---|
| Daily (trading days) | 252 | $\sqrt{252} \approx 15.87$ |
| Weekly | 52 | $\sqrt{52} \approx 7.21$ |
| Monthly | 12 | $\sqrt{12} \approx 3.46$ |
| Quarterly | 4 | $\sqrt{4} = 2$ |

### Why Log Returns?

Because $d(\ln S) = (\mu - \sigma^2/2)dt + \sigma\,dz$, log returns have *constant* variance per unit time. Simple returns break down for large moves.

### Why $\sqrt{T}$?

**Variance scales linearly with time, std dev scales with $\sqrt{T}$.**

---

## 6. Wiener Processes & Itô's Lemma (Problems 11, 12)

### Wiener Process Properties

$z_t$ is a Wiener process if:

1. $\Delta z = \epsilon \sqrt{\Delta t}$ where $\epsilon \sim N(0,1)$
2. Increments over non-overlapping intervals are independent

**Consequences:**
- $\Delta z \sim N(0, \Delta t)$
- $E[z_T] = 0$, $\text{Var}(z_T) = T$
- $z_T \sim N(0, T)$

### The Stochastic Calculus Multiplication Table

| Product | Result |
|---------|--------|
| $(dt)^2$ | $0$ |
| $dt \cdot dz$ | $0$ |
| $(dz)^2$ | $dt$ |

**Key insight:** $(dz)^2 = dt$ is the rule that creates the Itô correction. It survives because $\Delta z \sim \sqrt{\Delta t}$, so $(\Delta z)^2 \sim \Delta t$ — the same order as $dt$ itself.

### Itô's Lemma — General Form

If $dx = a(x,t)\,dt + b(x,t)\,dz$ and $G = G(x,t)$:

$$dG = \left(\frac{\partial G}{\partial t} + a\frac{\partial G}{\partial x} + \frac{1}{2}b^2\frac{\partial^2 G}{\partial x^2}\right)dt + b\frac{\partial G}{\partial x}\,dz$$

### Itô's Lemma — GBM Specialization

If $dS = \mu S\,dt + \sigma S\,dz$ and $G = G(S,t)$:

$$dG = \left(\frac{\partial G}{\partial t} + \mu S\frac{\partial G}{\partial S} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 G}{\partial S^2}\right)dt + \sigma S\frac{\partial G}{\partial S}\,dz$$

### Workflow

1. Identify $G$ and underlying process (stock = GBM)
2. Compute three partials: $\frac{\partial G}{\partial t}$, $\frac{\partial G}{\partial S}$, $\frac{\partial^2 G}{\partial S^2}$
3. Plug into template
4. Simplify and express in terms of $G$ (or $y$) if asked

### Important Applications

| $G(S,t)$ | $dG$ |
|----------|------|
| $\ln S$ | $\left(\mu - \frac{\sigma^2}{2}\right)dt + \sigma\,dz$ |
| $S^n$ | $n\left(\mu + \frac{(n-1)\sigma^2}{2}\right)G\,dt + n\sigma G\,dz$ |
| $S^2$ | $(2\mu + \sigma^2)G\,dt + 2\sigma G\,dz$ |
| $S^3$ | $(3\mu + 3\sigma^2)G\,dt + 3\sigma G\,dz$ |
| $S^4$ | $(4\mu + 6\sigma^2)G\,dt + 4\sigma G\,dz$ |
| $1/S^2 = S^{-2}$ | $(-2\mu + 3\sigma^2)G\,dt - 2\sigma G\,dz$ |
| $Se^{r(T-t)}$ | $(\mu - r)G\,dt + \sigma G\,dz$ |
| $\cos(t)$ | $-\sin(t)\,dt$ (no $dz$ term — deterministic) |

### Geometric Brownian Motion — Key Result

$$S_T = S_0 \exp\left[\left(\mu - \frac{\sigma^2}{2}\right)T + \sigma\sqrt{T}\,\epsilon\right], \quad \epsilon \sim N(0,1)$$

So $\ln(S_T/S_0) \sim N\left((\mu - \sigma^2/2)T, \sigma^2 T\right)$.

The $-\sigma^2/2$ term is the **volatility drag** — the gap between arithmetic and geometric mean returns.

### Problem 12 Type: Direct Stochastic Integration

When drift and diffusion don't depend on $S$, integrate $dS$ directly.

**Example:** $dS_t = \ln(t)\,dt + dz$, $S_0 = 0$.

$$S_T = \int_0^T \ln(t)\,dt + z_T = (T\ln T - T) + z_T$$

- $E(S_T) = T\ln T - T$ (deterministic part)
- $\text{Var}(S_T) = T$ (variance of $z_T$)

**Standard limit:** $\lim_{t \to 0^+} t\ln t = 0$.

### Itô Isometry (for variance with non-constant volatility)

$$\text{Var}\left(\int_0^T b(t)\,dz_t\right) = \int_0^T b(t)^2\,dt$$

---

## 7. Black-Scholes-Merton Pricing (Problems 13, 14)

### The BSM PDE (Universal)

Applies to **any** derivative $f(S,t)$ on the same underlying stock:

$$\frac{\partial f}{\partial t} + rS\frac{\partial f}{\partial S} + \frac{1}{2}\sigma^2 S^2\frac{\partial^2 f}{\partial S^2} = rf$$

What distinguishes different derivatives is the **boundary condition** at $t = T$ (the payoff).

### BSM Pricing Formulas (Provided on Exam)

For a European call/put on a non-dividend stock:

$$c = S_0 N(d_1) - K e^{-rT} N(d_2)$$
$$p = K e^{-rT} N(-d_2) - S_0 N(-d_1)$$

where:
$$d_1 = \frac{\ln(S_0/K) + (r + \sigma^2/2)T}{\sigma\sqrt{T}}, \quad d_2 = d_1 - \sigma\sqrt{T}$$

### Dividend-Yield Version

For a stock with continuous dividend yield $q$:
$$c = S_0 e^{-qT} N(d_1) - K e^{-rT} N(d_2)$$
$$d_1 = \frac{\ln(S_0/K) + (r - q + \sigma^2/2)T}{\sigma\sqrt{T}}, \quad d_2 = d_1 - \sigma\sqrt{T}$$

### Problem 13 — Direct BSM Application

> Calculate the price of a 9-month European call on a non-dividend stock with $K = 100$, $S_0 = 105$, $r = 0.10$, $\sigma = 0.30$.

**Direct formula method:**

Given: $S_0 = 105$, $K = 100$, $r = 0.10$, $\sigma = 0.30$, $T = 0.75$.

Step 1: Compute $d_1$ and $d_2$:

$$d_1 = \frac{\ln(105/100) + (0.10 + 0.09/2)(0.75)}{0.30\sqrt{0.75}}$$

- $\ln(105/100) = \ln(1.05) = 0.04879$
- $(0.10 + 0.045)(0.75) = 0.145 \times 0.75 = 0.10875$
- Numerator: $0.04879 + 0.10875 = 0.15754$
- Denominator: $0.30 \times 0.8660 = 0.2598$
- $d_1 = 0.15754 / 0.2598 = 0.6065$
- $d_2 = 0.6065 - 0.2598 = 0.3467$

Step 2: Plug into BSM formula:

$$c = 105 \cdot N(0.6065) - 100 e^{-0.10 \times 0.75} N(0.3467)$$
$$c = 105 \cdot N(0.6065) - 92.77 \cdot N(0.3467)$$

Leave as $N(d_1)$, $N(d_2)$ unless the exam asks for the numeric value.

### Problem 13 — Alternative: No-Arbitrage / Risk-Neutral Method

The risk-neutral pricing formula:

$$c = e^{-rT} \cdot \tilde{E}[\max(S_T - K, 0)]$$

Under the risk-neutral measure, $S_T = S_0 \exp\left[(r - \sigma^2/2)T + \sigma\sqrt{T}\,Z\right]$ where $Z \sim N(0,1)$.

So:
$$c = e^{-rT} \int_{-\infty}^{\infty} \max(S_T(z) - K, 0) \phi(z)\,dz$$

where $\phi(z)$ is the standard normal density.

The integrand is nonzero only when $S_T > K$, i.e., when:
$$Z > \frac{\ln(K/S_0) - (r - \sigma^2/2)T}{\sigma\sqrt{T}} = -d_2$$

Evaluating the integral splits into two pieces:

**Piece 1:** $e^{-rT} \cdot \tilde{E}[S_T \cdot \mathbf{1}_{S_T > K}] = S_0 N(d_1)$
**Piece 2:** $e^{-rT} \cdot K \cdot \tilde{P}(S_T > K) = K e^{-rT} N(d_2)$

Subtracting:
$$c = S_0 N(d_1) - K e^{-rT} N(d_2)$$

This is the no-arbitrage derivation of the BSM formula. It shows *why* the formula has the structure it does.

### Problem 14 — BSM with Non-Standard Payoff

Payoff: $\max(S_T^3 - K, 0)$.

**Part 1: The BSM PDE.** The PDE is the **same** as any derivative — only the boundary condition changes:

$$\frac{\partial f}{\partial t} + rS\frac{\partial f}{\partial S} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 f}{\partial S^2} = rf$$

with boundary condition $f(S,T) = \max(S^3 - K, 0)$.

**Part 2: Price.** Use the substitution trick.

1. Let $Y_t = S_t^3$. Then payoff = $\max(Y_T - K, 0)$ — a standard call on $Y$.

2. Apply Itô's Lemma to find the process for $Y$:
$$dY = 3(\mu + \sigma^2)Y\,dt + 3\sigma Y\,dz$$

3. Under the risk-neutral measure ($\mu \to r$):
$$dY = 3(r + \sigma^2)Y\,dt + 3\sigma Y\,dz$$

4. Treat $Y$ as an asset with:
   - Spot: $Y_0 = S_0^3$
   - Volatility: $\sigma_Y = 3\sigma$
   - Dividend yield: $q_Y = r - 3(r + \sigma^2) = -2r - 3\sigma^2$ (negative!)

5. Apply dividend-yield BSM:
$$f(0) = S_0^3 e^{(2r + 3\sigma^2)T} N(d_1) - K e^{-rT} N(d_2)$$

where:
$$d_1 = \frac{\ln(S_0^3/K) + (3r + \frac{15\sigma^2}{2})T}{3\sigma\sqrt{T}}, \quad d_2 = d_1 - 3\sigma\sqrt{T}$$

**Why negative dividend yield?** It's a bookkeeping device meaning "$Y$ grows faster than $r$ under risk-neutral measure."

### General Result for $G = S^n$

If $S$ is GBM, then under the risk-neutral measure:

$$\tilde{E}[S_T^n] = S_0^n \exp\left[\left(nr + \frac{n(n-1)\sigma^2}{2}\right)T\right]$$

This is the formula that drives Hull 15.12.

---

## 8. Hull 15.12: Derivative Paying $S_T^n$ (Final Problem 11)

> A derivative pays $S_T^n$ at time $T$. Show its price at time $t$ has the form $h(t,T) S^n$ where:
> $$h(t,T) = \exp\left\{\left[\frac{1}{2}\sigma^2 n(n-1) + r(n-1)\right](T-t)\right\}$$

### Setup

Assume the price has the form $f(S,t) = h(t,T) S^n$.

**Compute partial derivatives:**
- $\frac{\partial f}{\partial t} = h'(t) S^n$ where $h'(t) = \frac{dh}{dt}$
- $\frac{\partial f}{\partial S} = h(t) \cdot n S^{n-1}$
- $\frac{\partial^2 f}{\partial S^2} = h(t) \cdot n(n-1) S^{n-2}$

### Substitute into BSM PDE

$$\frac{\partial f}{\partial t} + rS\frac{\partial f}{\partial S} + \frac{1}{2}\sigma^2 S^2\frac{\partial^2 f}{\partial S^2} = rf$$

$$h'(t) S^n + rS \cdot h(t) n S^{n-1} + \frac{1}{2}\sigma^2 S^2 \cdot h(t) n(n-1) S^{n-2} = r h(t) S^n$$

Factor out $S^n$ (it's nonzero):

$$h'(t) + rn \cdot h(t) + \frac{1}{2}\sigma^2 n(n-1) h(t) = r h(t)$$

Rearrange:
$$h'(t) = h(t)\left[r - rn - \frac{1}{2}\sigma^2 n(n-1)\right] = -h(t)\left[\frac{1}{2}\sigma^2 n(n-1) + r(n-1)\right]$$

### Boundary Condition

At maturity $t = T$: $f(S,T) = S^n$, so $h(T,T) S^n = S^n$, giving:
$$h(T,T) = 1$$

### Solve the ODE

Let $\alpha = \frac{1}{2}\sigma^2 n(n-1) + r(n-1)$. The ODE is:

$$\frac{dh}{dt} = -\alpha h \implies h(t) = C e^{-\alpha t}$$

Applying $h(T) = 1$:
$$1 = C e^{-\alpha T} \implies C = e^{\alpha T}$$

So:
$$h(t,T) = e^{\alpha T} e^{-\alpha t} = e^{\alpha(T-t)}$$

$$\boxed{h(t,T) = \exp\left\{\left[\frac{1}{2}\sigma^2 n(n-1) + r(n-1)\right](T-t)\right\}}$$ ✓

### Workflow Summary

1. Assume $f = h(t)S^n$ and compute partials
2. Substitute into BSM PDE
3. Factor out $S^n$ → ODE for $h(t)$
4. Apply boundary condition $h(T) = 1$
5. Solve the ODE

---

## 9. Capped Options (Final Problem 10)

### Definition

A **capped option** has payoff:

$$\text{Payoff} = \min\left(\max(S_T - K, 0), L\right)$$

The payoff is a call, but **capped at maximum value $L$**.

### Decomposition into Bull Spread

A capped option is equivalent to a **bull spread** with strikes $K$ and $K + L$:

$$\min\left(\max(S_T - K, 0), L\right) = \max(S_T - K, 0) - \max(S_T - (K+L), 0)$$

**Verification (three regions):**

| Region | Long call at $K$ | Short call at $K+L$ | Net payoff |
|--------|------------------|---------------------|------------|
| $S_T \leq K$ | $0$ | $0$ | $0$ |
| $K < S_T \leq K+L$ | $S_T - K$ | $0$ | $S_T - K$ |
| $S_T > K+L$ | $S_T - K$ | $-(S_T - K - L)$ | $L$ |

### Pricing

Since the payoff decomposes into two calls:

$$\text{Capped call price} = c(K) - c(K+L)$$

where $c(K)$ is the BSM price of a standard European call with strike $K$.

### Generalized Payoff Replication (e.g., Hull Problem 8)

Given payoff:
$$Z = \begin{cases} 0, & S_T \leq K_1 \\ S_T - K_1, & K_1 < S_T \leq K_2 \\ K_2 - K_1, & S_T > K_2 \end{cases}$$

**Replicate with European calls:**

$$Z = \max(S_T - K_1, 0) - \max(S_T - K_2, 0)$$

**Time-0 value:** $c(K_1) - c(K_2)$, where $c(K)$ is the call price at strike $K$.

**Nonnegative because:** $K_2 > K_1 \Rightarrow c(K_1) \geq c(K_2)$ (call prices decrease in strike).

---

## Master Formula Sheet

### Forwards & Futures
- $F_0 = S_0 e^{(r-q)T}$ — cost of carry
- $f = (F - K)e^{-rT}$ — forward contract value
- $f_0 = 0$ at initiation (defines $F_0$)

### Bonds
- Zero rate: $R = -\frac{1}{T}\ln(P/\text{Face})$
- Par yield: $c = \frac{(\text{Face} - \text{Face} \cdot d_n) \cdot m}{\sum d_i}$
- Forward rate: $R_F = \frac{R_2 T_2 - R_1 T_1}{T_2 - T_1}$
- Compounding: $R_m = m(e^{R_c/m} - 1)$

### Binomial Trees
- $u = e^{\sigma\sqrt{\Delta t}}$, $d = 1/u$
- $\tilde{p} = \frac{e^{(r-q)\Delta t} - d}{u - d}$
- $f = e^{-r\Delta t}[\tilde{p} f_u + (1-\tilde{p}) f_d]$
- American: take max with exercise value
- $\Delta = (f_u - f_d)/(S_u - S_d)$

### Volatility Estimation
- $u_i = \ln(S_i/S_{i-1})$
- $s = $ sample SD with $n-1$ divisor
- $\hat{\sigma} = s/\sqrt{\tau}$
- Monthly: multiply by $\sqrt{12}$

### Stochastic Calculus
- $(dz)^2 = dt$
- $\Delta z = \epsilon\sqrt{\Delta t}$, $\epsilon \sim N(0,1)$
- $E[z_T] = 0$, $\text{Var}(z_T) = T$
- Itô isometry: $\text{Var}\left(\int_0^T b(t)\,dz\right) = \int_0^T b(t)^2\,dt$

### Itô's Lemma (GBM)
$$dG = \left(\frac{\partial G}{\partial t} + \mu S\frac{\partial G}{\partial S} + \frac{1}{2}\sigma^2 S^2\frac{\partial^2 G}{\partial S^2}\right)dt + \sigma S\frac{\partial G}{\partial S}\,dz$$

### Black-Scholes-Merton
- PDE: $\frac{\partial f}{\partial t} + rS\frac{\partial f}{\partial S} + \frac{1}{2}\sigma^2 S^2\frac{\partial^2 f}{\partial S^2} = rf$
- Call: $c = S_0 N(d_1) - K e^{-rT} N(d_2)$
- Put: $p = K e^{-rT} N(-d_2) - S_0 N(-d_1)$
- $d_1 = \frac{\ln(S_0/K) + (r + \sigma^2/2)T}{\sigma\sqrt{T}}$
- $d_2 = d_1 - \sigma\sqrt{T}$
- With dividend yield: replace $S_0$ with $S_0 e^{-qT}$ in $c$, replace $r$ with $r-q$ in $d_1$ numerator

### Risk-Neutral Pricing
$$f(0) = e^{-rT} \cdot \tilde{E}[\text{Payoff}_T]$$

### For $G = S^n$ (Under Risk-Neutral Measure)
- $\tilde{E}[S_T^n] = S_0^n \exp\left[\left(nr + \frac{n(n-1)\sigma^2}{2}\right)T\right]$
- $\sigma_Y = n\sigma$
- $q_Y = (1-n)r - \frac{n(n-1)\sigma^2}{2}$

---

## Common Traps & Final Tips

### Top Errors to Avoid

1. **Forgetting to discount** future values back to today
2. **Using $n$ instead of $n-1$** in sample variance
3. **Mixing up forward price and forward value** ($F$ vs. $f$)
4. **Using terminal payoffs instead of node values** in backward induction
5. **Forgetting principal in final bond cash flow** ($C + \text{Face}$, not just $C$)
6. **Sign errors in $q_Y$** when converting to dividend-yield BSM
7. **Wrong annualization factor** (use $\sqrt{12}$ for monthly, not 12)
8. **Skipping the early-exercise check** for American puts
9. **Using ordinary chain rule** instead of Itô's Lemma (missing the second-order term)
10. **Confusing risk-neutral $\mu = r$ with real-world $\mu$** in derivations

### Strategic Tips

- **Show your work clearly.** Partial credit is your friend on Problems 11–14.
- **Sanity-check answers.** Bond zero rates should be roughly comparable. Volatility above ~100% should make you pause.
- **Leave $N(d_1)$, $N(d_2)$ unevaluated** unless the exam asks for numerics.
- **Box your final answers.** Grader-friendly.
- **Time budget:** 14 problems in (likely) 2 hours = 8-9 minutes per problem on average. Spend less on early problems, more on 11-14.
- **If you blank on Itô's Lemma:** re-derive it from the Taylor expansion — it's faster than panicking.

### Day-of-Exam Checklist

- ✅ Calculator (allowed)
- ✅ Pencils, erasers
- ❌ No textbook, no notes, no hand-held computer
- BSM formulas will be provided on the exam booklet

### The Big Picture

Almost every problem reduces to one of these workflows:

1. **Cost-of-carry** → forwards, futures
2. **Discounting + bootstrapping** → bonds, rates
3. **Payoff decomposition** → spreads, capped options, replication
4. **Backward induction** → binomial trees
5. **Itô's Lemma machinery** → finding the process of any function $G(S,t)$
6. **BSM PDE/formula** → pricing any derivative

Master these six workflows and you've mastered the exam.

---

*Good luck — you've earned this preparation.*

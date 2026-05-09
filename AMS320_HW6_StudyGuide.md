# AMS320 HW6 — Comprehensive Study Guide

**Topics:** Wiener Processes, Itô's Lemma, and the Black-Scholes-Merton Model
**Reference:** Hull, *Options, Futures, and Other Derivatives* (9th ed.), Chapters 14–15

---

## Table of Contents

1. [Foundational Concepts](#1-foundational-concepts)
2. [Wiener Process & Itô Process](#2-wiener-process--itô-process)
3. [Stock Price Process: GBM](#3-stock-price-process-geometric-brownian-motion-gbm)
4. [Itô's Lemma](#4-itôs-lemma)
5. [Continuously Compounded Returns](#5-continuously-compounded-rate-of-return)
6. [Volatility Drag](#6-volatility-drag-variance-drain)
7. [Black-Scholes-Merton Formulas](#7-black-scholes-merton-formulas)
8. [The BSM PDE](#8-the-bsm-partial-differential-equation)
9. [Greeks](#9-the-greeks-partial-derivatives-of-option-prices)
10. [Put-Call Parity](#10-put-call-parity)
11. [American vs European Options](#11-american-vs-european-options)
12. [Worked Examples Recap](#12-worked-examples-recap)
13. [Common Pitfalls Checklist](#13-common-pitfalls-checklist)
14. [Key Identities Cheat Sheet](#14-key-identities-cheat-sheet)

---

## 1. Foundational Concepts

### Volatility vs. Variance

🔑 **Volatility = Standard Deviation, NOT variance**

| Term | Symbol | Units |
|------|--------|-------|
| Volatility (= standard deviation) | σ | Same units as the variable (%, $, etc.) |
| Variance | σ² | Squared units |

If a problem says "volatility is 30% per annum":
- σ = 0.30 (standard deviation per year)
- σ² = 0.09 (variance per year)

**Why volatility is defined as std dev, not variance:** Standard deviation has the same units as returns, making it directly interpretable.

---

### The Square-Root-of-Time Rule

🔑 **Volatility scales with the square root of time.**

If annual volatility is σ, then over a time interval Δt (expressed as a fraction of a year):

$$\text{Std. dev. over } \Delta t = \sigma \sqrt{\Delta t}$$

| Period | Δt | Volatility |
|--------|-----|-----------|
| Daily | 1/252 | σ/√252 |
| Weekly | 1/52 | σ/√52 |
| Monthly | 1/12 | σ/√12 |

**Convention:** Hull uses **252 trading days** per year (not 365 calendar days), since stock prices don't change on weekends/holidays.

This rule comes from the Wiener process property that **Var(Δz) = Δt**.

---

## 2. Wiener Process & Itô Process

### Basic Wiener Process

A variable z follows a Wiener process if:

$$\Delta z = \epsilon \sqrt{\Delta t}, \quad \epsilon \sim N(0,1)$$

- Drift rate = 0
- Variance rate = 1
- $\Delta z \sim N(0, \Delta t)$
- The Markov property holds (future depends only on the present)

### Generalized Wiener Process

$$dx = a\, dt + b\, dz$$

where a and b are constants.
- Drift rate = a
- Variance rate = b²
- Over time T: change has mean **aT** and variance **b²T**

### Itô Process

$$dx = a(x,t)\, dt + b(x,t)\, dz$$

Drift and variance rates are functions of x and t. Still Markov, since drift and diffusion at time t depend only on the current value $x$, not on the history.

---

## 3. Stock Price Process: Geometric Brownian Motion (GBM)

$$\frac{dS}{S} = \mu\, dt + \sigma\, dz \quad\Longleftrightarrow\quad dS = \mu S\, dt + \sigma S\, dz$$

Discrete version:

$$\Delta S = \mu S \Delta t + \sigma S \epsilon \sqrt{\Delta t}$$

This implies:

$$\frac{\Delta S}{S} \sim N(\mu \Delta t, \sigma^2 \Delta t)$$

Where:
- μ = expected rate of return
- σ = volatility

**In a risk-neutral world**, μ = r (the risk-free rate). This is the world used for option pricing.

---

## 4. Itô's Lemma

If $dx = a(x,t)\, dt + b(x,t)\, dz$ and G = G(x, t), then:

$$dG = \left(\frac{\partial G}{\partial t} + a\frac{\partial G}{\partial x} + \frac{1}{2}b^2\frac{\partial^2 G}{\partial x^2}\right)dt + b\frac{\partial G}{\partial x}\, dz$$

🔑 The extra $\frac{1}{2}b^2 \frac{\partial^2 G}{\partial x^2}$ term is what distinguishes stochastic from ordinary calculus. It comes from the fact that $(dz)^2 = dt$ (heuristically), so second-order terms in $dz$ contribute to the drift.

### Key Application: Log of Stock Price

Let G = ln(S) where dS = μS dt + σS dz. Then:

$$d(\ln S) = \left(\mu - \frac{\sigma^2}{2}\right)dt + \sigma\, dz$$

Therefore over time T:

$$\ln S_T - \ln S_0 \sim N\left(\left(\mu - \frac{\sigma^2}{2}\right)T,\; \sigma^2 T\right)$$

🔑 The $-\sigma^2/2$ term is the **Itô correction** — it's not just μ!

### Itô's Lemma for GBM (Quick Reference)

If $dS = \mu S\, dt + \sigma S\, dz$ and $G = G(S, t)$:

$$dG = \left(\frac{\partial G}{\partial t} + \mu S\frac{\partial G}{\partial S} + \frac{1}{2}\sigma^2 S^2\frac{\partial^2 G}{\partial S^2}\right)dt + \sigma S \frac{\partial G}{\partial S}\, dz$$

This is the form most often used for derivatives pricing.

---

## 5. Continuously Compounded Rate of Return

Define rate x such that $S_T = S_0 e^{xT}$. Then:

$$x = \frac{1}{T}\ln\left(\frac{S_T}{S_0}\right)$$

**Distribution of x:**

$$x \sim N\left(\mu - \frac{\sigma^2}{2},\; \frac{\sigma^2}{T}\right)$$

**Key observations:**
- Mean does NOT depend on T
- Variance shrinks as T grows (longer-term rates are more predictable)

🔑 **Two different "expected returns":**

| Quantity | Value | Meaning |
|----------|-------|---------|
| E[ΔS/S] per unit time | μ | Expected *arithmetic* return |
| E[Δ(ln S)] per unit time | μ − σ²/2 | Expected *log* (continuously compounded) return |

### Lognormal Moments

From the distribution of $\log S_T$:

$$E(S_T) = S_0 e^{\mu T}, \qquad \text{Var}(S_T) = S_0^2 e^{2\mu T}\left(e^{\sigma^2 T} - 1\right)$$

---

## 6. Volatility Drag (Variance Drain)

The gap σ²/2 between arithmetic and log returns is called **volatility drag**.

**Concrete example:** Stock goes up 50%, then down 50%:
- $100 → $150 → $75 (lost 25% over 2 periods!)

**In log space:**
- ln(1.5) ≈ +0.405
- ln(0.5) ≈ −0.693
- Sum ≈ −0.288 → $100 × e^(-0.288) ≈ $75 ✓

The log-return of a loss is larger in magnitude than the log-return of an equivalent-sized gain. Higher volatility → bigger drag on compounded returns.

### Why are stock prices lognormal?

```
Stock prices can't go negative
            ↓
Model percentage changes (GBM) → S is lognormal
            ↓
        Apply Itô's lemma
            ↓
d(ln S) has drift μ − σ²/2 (not μ)
            ↓
    Volatility drag exists
```

Volatility drag is a **consequence** of lognormality (via Itô's lemma & Jensen's inequality), not the cause.

---

## 7. Black-Scholes-Merton Formulas

### European Call & Put on Non-Dividend-Paying Stock

$$c = S_0 N(d_1) - Ke^{-rT}N(d_2)$$

$$p = Ke^{-rT}N(-d_2) - S_0 N(-d_1)$$

where:

$$d_1 = \frac{\ln(S_0/K) + (r + \sigma^2/2)T}{\sigma\sqrt{T}}$$

$$d_2 = d_1 - \sigma\sqrt{T} = \frac{\ln(S_0/K) + (r - \sigma^2/2)T}{\sigma\sqrt{T}}$$

🔑 **Memorization tip:** The "S side" gets d₁, the "K side" gets d₂.

### Symmetry Property of Normal CDF

$$N(-x) = 1 - N(x)$$

⚠️ **Calculator tip:** Use **normCdf** (cumulative), NOT **normPdf** (density)!
- normPdf(x) = height of bell curve at x
- normCdf(-∞, x) = N(x) = area under curve up to x ← what we want

**TI-Nspire syntax:** `normCdf(-1E99, x, 0, 1)` for $N(x)$.

### Boundary Behavior

| Limit | Call price $c$ | Put price $p$ |
|-------|---------------|---------------|
| $S_0 \to \infty$ | $S_0 - Ke^{-rT}$ | 0 |
| $S_0 \to 0$ | 0 | $Ke^{-rT} - S_0$ |
| $\sigma \to 0$ | $\max(S_0 - Ke^{-rT}, 0)$ | $\max(Ke^{-rT} - S_0, 0)$ |
| $T \to \infty$ | $S_0$ | 0 |
| $t \to T$ | $\max(S - K, 0)$ | $\max(K - S, 0)$ |

---

## 8. The BSM Partial Differential Equation

Any derivative $f(S, t)$ on a non-dividend-paying stock following GBM must satisfy:

$$\boxed{\frac{\partial f}{\partial t} + rS\frac{\partial f}{\partial S} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 f}{\partial S^2} = rf}$$

### Key ideas behind the derivation

1. **Same source of uncertainty:** Stock and option both driven by the same Wiener process $dz$.
2. **Risk-elimination via hedging:** Construct portfolio $\Pi = -f + \frac{\partial f}{\partial S} \cdot S$ (short option, long $\Delta = \partial f/\partial S$ shares).
3. **Riskless portfolio earns risk-free rate:** $d\Pi = r\Pi\, dt$.
4. **Equating gives the PDE.**

### Boundary conditions for different derivatives

| Derivative | Boundary condition (at maturity $T$) |
|------------|--------------------------------------|
| European call | $f = \max(S - K, 0)$ |
| European put | $f = \max(K - S, 0)$ |
| Forward contract | $f = S - K$ → solution: $f = S - Ke^{-r(T-t)}$ |
| Power payoff $S^n$ | $f = S^n$ → solution: $f = h(t,T)\cdot S^n$ |

### Power-payoff derivative (Problem 15.12)

For a derivative paying $S_T^n$ at time $T$, with ansatz $f(S,t) = h(t,T)\cdot S^n$:

**ODE for $h$:**
$$h'(t,T) + \left[\tfrac{1}{2}\sigma^2 n(n-1) + r(n-1)\right] h(t,T) = 0$$

**Boundary condition:** $h(T, T) = 1$

**Solution:**
$$h(t,T) = \exp\left\{\left[\tfrac{1}{2}\sigma^2 n(n-1) + r(n-1)\right](T-t)\right\}$$

---

## 9. The Greeks (Partial Derivatives of Option Prices)

For a European call on a non-dividend-paying stock:

| Greek | Symbol | Formula | Meaning |
|-------|--------|---------|---------|
| Delta | $\Delta$ | $\dfrac{\partial c}{\partial S} = N(d_1)$ | Sensitivity to stock price |
| Theta | $\Theta$ | $\dfrac{\partial c}{\partial t} = -rKe^{-r(T-t)}N(d_2) - \dfrac{SN'(d_1)\sigma}{2\sqrt{T-t}}$ | Time decay |
| Gamma | $\Gamma$ | $\dfrac{\partial^2 c}{\partial S^2} = \dfrac{N'(d_1)}{S\sigma\sqrt{T-t}}$ | Convexity / sensitivity of delta |

### Useful Identity (from Problem 15.17b)

$$\boxed{S N'(d_1) = K e^{-r(T-t)} N'(d_2)}$$

This identity lets many calculations cancel beautifully — used heavily in proving that $c$ satisfies the BSM PDE.

### Useful Partial Derivatives

$$\frac{\partial d_1}{\partial S} = \frac{\partial d_2}{\partial S} = \frac{1}{S\sigma\sqrt{T-t}}$$

$$\frac{\partial d_1}{\partial t} = \frac{d_1}{2(T-t)}, \qquad \frac{\partial d_2}{\partial t} = \frac{d_1}{2(T-t)} + \frac{\sigma}{2\sqrt{T-t}}$$

**Note:** Since $d_2 - d_1$ doesn't depend on $S$, the partials w.r.t. $S$ are equal. But $d_2 - d_1 = -\sigma\sqrt{T-t}$ does depend on $t$, so the partials w.r.t. $t$ differ by $\sigma/(2\sqrt{T-t})$.

### Standard Normal PDF/CDF Relationship

$$N(x) = \int_{-\infty}^x \frac{1}{\sqrt{2\pi}} e^{-u^2/2}\, du$$

$$N'(x) = \frac{1}{\sqrt{2\pi}} e^{-x^2/2}$$

By the Fundamental Theorem of Calculus, $N'$ is the standard normal PDF.

---

## 10. Put-Call Parity

For European options on a non-dividend-paying stock:

$$\boxed{c + Ke^{-rT} = p + S_0}$$

**Intuition:** Both portfolios pay $\max(S_T, K)$ at maturity:
- LHS: Buy call, deposit $Ke^{-rT}$ in cash
  - If $S_T > K$: exercise call, get stock worth $S_T$. Cash grows to $K$. Net: keep stock + use cash to do nothing. Value = $S_T$. Wait — actually you exercise and pay $K$ from the cash. So you end with the stock, worth $S_T$. ✓
  - If $S_T < K$: don't exercise. Cash grows to $K$. Value = $K$. ✓
- RHS: Buy put, buy stock
  - If $S_T > K$: don't exercise put. Hold stock. Value = $S_T$. ✓
  - If $S_T < K$: exercise put, sell stock for $K$. Value = $K$. ✓

Both portfolios pay $\max(S_T, K)$, so by no-arbitrage they cost the same today.

---

## 11. American vs European Options

🔑 **Key result for non-dividend-paying stocks:**

$$C_{\text{American call}} = c_{\text{European call}}$$

**Why?** Two reasons:
1. **Time value of strike:** Paying $K$ at maturity (PV = $Ke^{-rT}$) is cheaper than paying $K$ now.
2. **Insurance value:** Holding the option preserves the right to walk away if the stock crashes.

Together, these give the lower bound:
$$C_{\text{American}} \geq c_{\text{European}} \geq S - Ke^{-r(T-t)} > S - K$$

So you'd always prefer to **sell** the option ($\geq S - Ke^{-r(T-t)}$) rather than **exercise** it ($S - K$). Early exercise is never optimal, so the American premium adds nothing.

⚠️ **This argument FAILS for:**
- American **puts** (early exercise can be optimal even without dividends)
- American **calls on dividend-paying stocks** (just-before-ex-dividend exercise can be optimal)

---

## 12. Worked Examples Recap

### Problem 15.2 — Daily Std Dev of Returns

Given σ = 30% per annum. Find std dev of % price change in one trading day.

$$\sigma\sqrt{\Delta t} = 0.30 \times \sqrt{1/252} \approx 0.0189 = 1.89\%$$

---

### Problem 15.4 — European Put Pricing

S₀ = K = $50, r = 10%, σ = 30%, T = 0.25 (3 months)

| Step | Value |
|------|-------|
| ln(S₀/K) | 0 (at-the-money!) |
| d₁ | 0.2417 |
| d₂ | 0.0917 |
| N(−d₁) | 0.4045 |
| N(−d₂) | 0.4634 |
| Ke^(−rT) | 48.7654 |
| **p** | **≈ $2.38** |

---

### Problem 15.7 — Rate of Return Distribution

S₀ = $40, μ = 15%, σ = 25%, T = 2 years.

$$x \sim N\left(\mu - \frac{\sigma^2}{2},\; \frac{\sigma^2}{T}\right) = N\left(0.15 - \frac{0.0625}{2},\; \frac{0.0625}{2}\right) = N(0.11875,\; 0.03125)$$

- Mean ≈ **11.875%** per year
- Std dev ≈ **17.68%** per year

---

### Problem 15.12 — Power-Payoff Derivative

Given derivative paying $S_T^n$, ansatz $f(S,t) = h(t,T)\cdot S^n$.

**(a) Substitute into BSM PDE:**

Compute partials:
$$\frac{\partial f}{\partial t} = h'(t,T)\cdot S^n, \quad \frac{\partial f}{\partial S} = nh(t,T)\cdot S^{n-1}, \quad \frac{\partial^2 f}{\partial S^2} = n(n-1)h(t,T)\cdot S^{n-2}$$

After substitution and dividing by $S^n$:
$$h'(t,T) + \left[\tfrac{1}{2}\sigma^2 n(n-1) + r(n-1)\right] h(t,T) = 0$$

**(b) Boundary condition:** $h(T,T) = 1$ (since $f(S,T) = S^n$ must equal $h(T,T)\cdot S^n$).

**(c) Solution by separation of variables:**

Let $k = \tfrac{1}{2}\sigma^2 n(n-1) + r(n-1)$. The ODE is $h'(t) = -k h(t)$.

Separate: $\frac{dh}{h} = -k\, dt$. Integrate from $t$ to $T$:

$$\ln h(T) - \ln h(t) = -k(T - t)$$

Using $h(T) = 1$:
$$h(t,T) = \exp\{k(T-t)\} = \exp\left\{\left[\tfrac{1}{2}\sigma^2 n(n-1) + r(n-1)\right](T-t)\right\}$$

---

### Problem 15.13 — European Call Pricing

S₀ = $52, K = $50, r = 12%, σ = 30%, T = 0.25 (3 months)

| Step | Value |
|------|-------|
| ln(S₀/K) = ln(1.04) | 0.03922 |
| (r + σ²/2)T = (0.165)(0.25) | 0.04125 |
| Numerator of d₁ | 0.08047 |
| σ√T | 0.15 |
| **d₁** | **0.5365** |
| **d₂** | **0.3865** |
| N(d₁) | ≈ 0.7042 |
| N(d₂) | ≈ 0.6504 |
| **c** | **≈ $5.06** |

**Sanity check:** Intrinsic value = $2; time value = $3.06. ✓

---

### Problem 15.17 — Verifying BSM Formula Satisfies the PDE

**(a)** $N'(x) = \frac{1}{\sqrt{2\pi}} e^{-x^2/2}$ (standard normal PDF, by FTC)

**(b)** Show $SN'(d_1) = Ke^{-r(T-t)} N'(d_2)$:

Form the ratio $\frac{N'(d_1)}{N'(d_2)} = e^{(d_2^2 - d_1^2)/2}$.

Use difference of squares: $d_2^2 - d_1^2 = (d_2 - d_1)(d_2 + d_1)$, where:
- $d_2 - d_1 = -\sigma\sqrt{T-t}$
- $d_1 + d_2 = \frac{2\ln(S/K) + 2r(T-t)}{\sigma\sqrt{T-t}}$

After multiplying: $d_2^2 - d_1^2 = -2\ln(S/K) - 2r(T-t)$.

So $\frac{N'(d_1)}{N'(d_2)} = e^{-\ln(S/K) - r(T-t)} = \frac{K}{S} e^{-r(T-t)}$. ✓

**(c)** $\frac{\partial d_1}{\partial S} = \frac{\partial d_2}{\partial S} = \frac{1}{S\sigma\sqrt{T-t}}$

**(d)** Computing $\frac{\partial c}{\partial t}$ requires:
- $\frac{\partial d_1}{\partial t} = \frac{d_1}{2(T-t)}$
- $\frac{\partial d_2}{\partial t} = \frac{d_1}{2(T-t)} + \frac{\sigma}{2\sqrt{T-t}}$

After product rule and using identity (b) to combine terms:
$$\frac{\partial c}{\partial t} = -rKe^{-r(T-t)} N(d_2) - \frac{SN'(d_1)\sigma}{2\sqrt{T-t}}$$

**(e)** $\frac{\partial c}{\partial S} = N(d_1)$ — using parts (b) and (c), the messy terms cancel.

**(f)** Plugging $\partial c/\partial t$, $\partial c/\partial S$, and $\partial^2 c/\partial S^2 = \frac{N'(d_1)}{S\sigma\sqrt{T-t}}$ into the BSM PDE:

The two "messy" $\frac{SN'(d_1)\sigma}{2\sqrt{T-t}}$ terms cancel (one from $\partial c/\partial t$, one from $\frac{1}{2}\sigma^2 S^2 \partial^2 c/\partial S^2$).

What remains: $-rKe^{-r(T-t)}N(d_2) + rSN(d_1) = r[SN(d_1) - Ke^{-r(T-t)}N(d_2)] = rc$ ✓

**(g)** As $t \to T$, three cases:

| Case | $d_1, d_2$ | $N(d_1), N(d_2)$ | $c \to$ | $\max(S-K, 0)$ |
|------|-----------|------------------|---------|----------------|
| $S > K$ | $+\infty$ | 1 | $S - K$ | $S - K$ ✓ |
| $S < K$ | $-\infty$ | 0 | 0 | 0 ✓ |
| $S = K$ | 0 | 1/2 | 0 | 0 ✓ |

---

### Problem 15.30 — Multi-Part Option Pricing

S₀ = $30, K = $29, r = 5%, σ = 25%, T = 1/3 year (4 months)

| Step | Value |
|------|-------|
| ln(S₀/K) = ln(30/29) | 0.0339 |
| (r + σ²/2)T | 0.02708 |
| σ√T | 0.1443 |
| **d₁** | **0.4226** |
| **d₂** | **0.2783** |

**(a) European call:** $c \approx \$2.52$

**(b) American call:** $C = c \approx \$2.52$ (no dividends → early exercise never optimal)

**(c) European put:** $p \approx \$1.05$

**(d) Put-call parity check:**
- LHS: $c + Ke^{-rT} = 2.52 + 29 e^{-0.05/3} = 2.52 + 28.52 = \$31.04$
- RHS: $p + S_0 = 1.05 + 30 = \$31.05$
- Match (within rounding) ✓

---

## 13. Common Pitfalls Checklist

- [ ] Don't confuse volatility (σ) with variance (σ²)
- [ ] Use 252 trading days, not 365
- [ ] Volatility scales with √Δt (not Δt)
- [ ] Use **normCdf**, not normPdf, on your calculator
- [ ] In BSM put formula: K-side pairs with d₂, S-side pairs with d₁
- [ ] Expected log return is μ − σ²/2, NOT μ
- [ ] Option prices are always ≥ 0 (sanity check!)
- [ ] N(−x) = 1 − N(x) (symmetry of standard normal)
- [ ] When applying chain rule on $N(d_i)$, don't forget the $\partial d_i / \partial \cdot$ factor!
- [ ] $\sqrt{T-t}$ — the square root covers the *whole* expression, not just $T$
- [ ] Don't try to factor out $h(t,T)$ from $h'(t,T)$ — they're different quantities
- [ ] American calls = European calls only on **non-dividend-paying** stocks
- [ ] Put-call parity: it's $Ke^{-rT}$ on the call side, plain $S_0$ on the put side

---

## 14. Key Identities Cheat Sheet

### Stochastic Calculus
- **Itô product:** $(dz)^2 = dt$, $(dt)^2 = 0$, $dz \cdot dt = 0$
- **GBM solution:** $S_T = S_0 \exp\left[(\mu - \sigma^2/2)T + \sigma(z_T - z_0)\right]$
- **Lognormal moments:** $E(S_T) = S_0 e^{\mu T}$, $\text{Var}(S_T) = S_0^2 e^{2\mu T}(e^{\sigma^2 T} - 1)$

### BSM Formulas
- **Call:** $c = S_0 N(d_1) - Ke^{-rT}N(d_2)$
- **Put:** $p = Ke^{-rT}N(-d_2) - S_0 N(-d_1)$
- **Put-call parity:** $c + Ke^{-rT} = p + S_0$

### Useful Greeks Identities
- $SN'(d_1) = Ke^{-r(T-t)}N'(d_2)$ ← workhorse identity
- $\partial d_1 / \partial S = \partial d_2 / \partial S = 1/(S\sigma\sqrt{T-t})$
- $\partial c/\partial S = N(d_1)$ (delta of call)
- $\partial p/\partial S = N(d_1) - 1 = -N(-d_1)$ (delta of put)
- $\partial^2 c/\partial S^2 = \partial^2 p/\partial S^2 = N'(d_1)/(S\sigma\sqrt{T-t})$ (gamma)

### Standard Normal
- $N'(x) = \frac{1}{\sqrt{2\pi}}e^{-x^2/2}$
- $N(-x) = 1 - N(x)$
- $N(0) = 1/2$
- $N(\infty) = 1$, $N(-\infty) = 0$

---

*End of Study Guide*

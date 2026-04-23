# Binomial Trees — AMS 320 Reference Guide
*Based on Hull, "Options, Futures, and Other Derivatives," 9th Edition*

---

## 1. Core Concepts

### What is a Binomial Tree?
A binomial tree models stock price movements over discrete time steps. At each step, the stock price can either:
- Move **up** by a factor of **u > 1**
- Move **down** by a factor of **d < 1**

### Key Parameters
| Symbol | Meaning |
|--------|---------|
| S₀ | Current stock price |
| u | Up factor (e.g. u = 1.1 means 10% up) |
| d | Down factor (e.g. d = 0.9 means 10% down) |
| r | Risk-free interest rate (per annum) |
| Δt | Length of each time step (in years) |
| T | Total time to expiration |
| K | Strike price |
| σ | Volatility of the stock |

### Option Payoffs at Expiration
- **European Call:** f = max(S_T − K, 0)
- **European Put:** f = max(K − S_T, 0)

---

## 2. Matching Volatility: u and d (Cox-Ross-Rubinstein)

Rather than choosing u and d arbitrarily, we calibrate them to real market volatility σ:

$$u = e^{\sigma\sqrt{\Delta t}}, \quad d = e^{-\sigma\sqrt{\Delta t}} = \frac{1}{u}$$

**Key insight:** 
- **Low σ** → u and d are close to 1 → small stock movements (stable stock)
- **High σ** → u and d are far from 1 → large stock movements (volatile stock)

**Useful check:** u × d = 1 always holds.

---

## 3. Risk-Neutral Valuation

### Why Risk-Neutral?
Pricing with real-world probabilities requires knowing the correct **risk-adjusted discount rate** for the option — which depends on investor risk preferences and is very difficult to determine. Risk-neutral valuation sidesteps this by using a fictitious world where:
- The expected return on every asset is the **risk-free rate r**
- We discount all payoffs at the **risk-free rate r**

### Risk-Neutral Probability
$$p = \frac{e^{r\Delta t} - d}{u - d}$$

The probability of a down move is simply **(1 − p)**.

**Verification:** In the risk-neutral world, the expected stock price grows at the risk-free rate:
$$p \cdot S_0u + (1-p) \cdot S_0d = S_0 e^{r\Delta t}$$

### Recursive Pricing Formula (at each node)
$$f = e^{-r\Delta t}\left[p \cdot f_u + (1-p) \cdot f_d\right]$$

Work **backwards** from expiration to today.

---

## 4. No-Arbitrage Approach

### Core Idea
If two portfolios have **identical future payoffs**, they must have the **same price today**. Otherwise, a riskless profit (arbitrage) would exist.

### Setting Up the Riskless Portfolio
Construct a portfolio: **long Δ shares + short 1 option**

Set the up and down portfolio values equal:
$$\Delta \cdot S_0u - f_u = \Delta \cdot S_0d - f_d$$

Solving for Δ (the **hedge ratio**):
$$\Delta = \frac{f_u - f_d}{S_0(u - d)}$$

### Pricing the Option
Since the riskless portfolio must earn the risk-free rate:
$$(\Delta \cdot S_0 - f) \cdot e^{rT} = \Delta \cdot S_0u - f_u$$

Solving for f gives exactly the same result as risk-neutral valuation:
$$f = e^{-rT}\left[p \cdot f_u + (1-p) \cdot f_d\right]$$

**Both methods always give the same answer.**

---

## 5. Step-by-Step Procedure

### European Option (any number of steps)

**Step 1:** Calculate Δt = T / (number of steps)

**Step 2:** Calculate u, d (if not given): $u = e^{\sigma\sqrt{\Delta t}}$, $d = 1/u$

**Step 3:** Calculate risk-neutral probability:
$$p = \frac{e^{r\Delta t} - d}{u - d}$$

**Step 4:** Build the stock price tree **forward**:
$$S_{u^i d^j} = S_0 \cdot u^i \cdot d^j$$

**Step 5:** Calculate option payoffs at expiration (rightmost nodes)

**Step 6:** Work **backwards** through the tree:
$$f = e^{-r\Delta t}\left[p \cdot f_u + (1-p) \cdot f_d\right]$$

---

### American Option (additional step)

At **every intermediate node**, compare:
1. **Discounted value:** $e^{-r\Delta t}[p \cdot f_u + (1-p) \cdot f_d]$
2. **Early exercise value:** max(K − S, 0) for put, or max(S − K, 0) for call

Take the **maximum** of the two:
$$f = \max\left(e^{-r\Delta t}[p \cdot f_u + (1-p) \cdot f_d],\ \text{early exercise value}\right)$$

---

## 6. Delta (Δ) — The Hedge Ratio

$$\Delta = \frac{f_u - f_d}{S_u - S_d}$$

- Δ represents the number of shares needed to hedge one option
- Δ **changes at every node** — so the riskless portfolio must be continuously rebalanced
- This is called **dynamic hedging** or **delta hedging**
- In practice, rebalancing incurs transaction costs not captured by the basic model

---

## 7. Three-Step Tree — Stock Price Nodes

For a three-step tree, there are **4 final stock prices** and **10 total nodes**:

### Final Stock Prices (t = T)
| Node | Price |
|------|-------|
| S_uuu | S₀ · u³ |
| S_uud | S₀ · u² · d |
| S_udd | S₀ · u · d² |
| S_ddd | S₀ · d³ |

### General Formula for n-Step Trees
$$S_{i,j} = S_0 \cdot u^j \cdot d^{i-j}$$
where i = time step, j = number of up moves

---

## 8. Worked Example: Three-Step European Call

**Parameters:** S₀ = $100, u = 1.1, d = 0.9, r = 5%, Δt = 1/3, K = $100

**Step 1: Risk-neutral probability**
$$p = \frac{e^{0.05/3} - 0.9}{1.1 - 0.9} = 0.5840$$

**Step 2: Final stock prices**
- S_uuu = 133.1, S_uud = 108.9, S_udd = 89.1, S_ddd = 72.9

**Step 3: Option payoffs**
- f_uuu = 33.1, f_uud = 8.9, f_udd = 0, f_ddd = 0

**Step 4: Work backwards** (using $e^{-r\Delta t} = 0.9835$)

*At t = 2Δt:*
- f_uu = 0.9835 × [0.5840 × 33.1 + 0.4160 × 8.9] = **22.65**
- f_ud = 0.9835 × [0.5840 × 8.9 + 0.4160 × 0] = **5.11**
- f_dd = 0

*At t = Δt:*
- f_u = 0.9835 × [0.5840 × 22.65 + 0.4160 × 5.11] = **15.10**
- f_d = 0.9835 × [0.5840 × 5.11 + 0.4160 × 0] = **2.935**

*At t = 0:*
$$f_0 = 0.9835 \times [0.5840 \times 15.10 + 0.4160 \times 2.935] = \mathbf{\$9.874}$$

---

## 9. Real World vs. Risk-Neutral World

| | Real World | Risk-Neutral World |
|--|------------|-------------------|
| Expected return | μ (hard to estimate) | r (risk-free rate) |
| Discount rate | Risk-adjusted (hard to determine) | r (risk-free rate) |
| Probability of up move | p* | p |
| Volatility σ | Same | Same |

**Girsanov's Theorem:** When moving between the real world and risk-neutral world, the **volatility stays the same** — only the expected return changes.

---

## 10. Key Reminders

- Always work **forward** to build the stock price tree, then **backward** to price the option
- For **American options**, check early exercise at **every intermediate node**
- **Δt = T / n** where n is the number of steps — don't confuse Δt with T
- The no-arbitrage and risk-neutral methods **always give the same price**
- u × d = 1 when using Cox-Ross-Rubinstein parameters (useful sanity check)

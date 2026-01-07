---
layout: post
title: The 2D strain rate tensor and its invariants
date: 2026-01-06 18:10:00-0000
description: Commonly overlooked conceptual inaccuracies in strain-rate studies 
tags: strain tensors linear-algebra 
categories: methods 
related_posts: false
---

**"Errare humanum est, perseverare autem diabolicum."**
<br> 
*"Making mistakes is human; persisting in them is diabolical."*

---

Before diving into the formal definitions, I want to note a few conceptual and
notational inaccuracies that I have run into over time when working with strain-rate
tensors. Some common in the literature, and some of my own making:

1. Implicitly **associating $\dot\varepsilon_1$ with extension and
   $\dot\varepsilon_2$ with shortening**, overlooking that these indices reflect
   eigenvalue ordering rather than deformation style, which depends on the sign of the
   principal strain rates.

2. Referring to the **Frobenius norm** of the strain-rate tensor as the “**second invariant**,”
   rather than as a distinct (and perfectly valid) scalar invariant.

The notes below are my best attempt to write these points down carefully, both to clarify
them for myself and to avoid having to rediscover them again later.

---

Given a set of geodetically derived horizontal surface velocities $\mathbf{v}$, the horizontal strain-rate tensor is defined as

$$
\dot{\varepsilon}_{ij}
=
\tfrac12\!\left(
\frac{\partial v_i}{\partial x_j}
+
\frac{\partial v_j}{\partial x_i}
\right)
=
\begin{pmatrix}
\frac{\partial v_1}{\partial x_1}
&
\tfrac12\!\left(\frac{\partial v_1}{\partial x_2}+\frac{\partial v_2}{\partial x_1}\right)
\\
\tfrac12\!\left(\frac{\partial v_2}{\partial x_1}+\frac{\partial v_1}{\partial x_2}\right)
&
\frac{\partial v_2}{\partial x_2}
\end{pmatrix},
\label{eq:strain_tensor}
$$

where $i,j$ denote spatial coordinates $x$ and velocity components $v$.  
Here we consider **horizontal surface strain rates**, so $i,j = 1,2$.

---

## Principal strain rates

The eigenvalues of $\dot{\varepsilon}$ correspond to the principal strain rates $\dot\varepsilon_1$ and $\dot\varepsilon_2$,

$$
\dot{\varepsilon}_{1,2}
=
\frac{\dot{\varepsilon}_{11}+\dot{\varepsilon}_{22}}{2}
\pm
\sqrt{
\left(\frac{\dot{\varepsilon}_{11}-\dot{\varepsilon}_{22}}{2}\right)^2
+
\dot{\varepsilon}_{12}^2
},
\label{eq:principal_strain}
$$

and the associated eigenvectors define the principal strain-rate orientations.

---

## Tensor invariants

In continuum mechanics, coordinate-system-independent quantities (tensor invariants) are often used to quantify geodetically derived horizontal strain rates. Under this definition, **invariants are not unique**.
The first and second invariants, $I_1$ and $I_2$, are distinguished because they appear as coefficients of the characteristic polynomial of the strain-rate tensor:

$$
\lambda^2 - I_1^{2D}(\dot\varepsilon)\,\lambda + I_2^{2D}(\dot\varepsilon) = 0,
$$

where $I_1^{2D}$ and $I_2^{2D}$ are the first and second invariants of the 2D strain rate tensor, respectively. 

---

## First invariant (dilatation rate)

The first invariant is the trace of the strain-rate tensor and represents the areal strain rate (dilatation rate),

$$
I_1^{2D}(\dot\varepsilon)
=
\dot\Delta
=
\nabla\cdot\mathbf{v}
=
\operatorname{tr}(\dot\varepsilon)
=
\dot\varepsilon_{11}+\dot\varepsilon_{22}
=
\dot\varepsilon_1+\dot\varepsilon_2 .
\label{eq:I1}
$$

A related quantity is the mean horizontal strain rate,

$$
\dot\varepsilon_m^{2D}
=
\tfrac12 I_1^{2D}
=
\tfrac12(\dot\varepsilon_1+\dot\varepsilon_2).
\label{eq:mean_strain}
$$

---

## Second invariant

The formal definition of the second invariant is

$$
I_2(\dot\varepsilon)
=
\tfrac12\!\left[
\operatorname{tr}(\dot\varepsilon)^2
-
\operatorname{tr}(\dot\varepsilon^2)
\right].
\label{eq:I2_def}
$$

For a 2-D tensor, this reduces to the determinant,

$$
I_2^{2D}(\dot\varepsilon)
=
\det(\dot\varepsilon)
=
\dot\varepsilon_{11}\dot\varepsilon_{22}
-
\dot\varepsilon_{12}^2,
$$

or equivalently,

$$
I_2^{2D}(\dot\varepsilon)
=
\dot\varepsilon_1\,\dot\varepsilon_2 .
$$

Note that in three dimensions the determinant is not the second invariant, but the third invariant of the tensor.

---

## Frobenius norm (commonly misidentified as $I_2$)

In the literature, the following quantity is often (incorrectly) referred to as the second invariant:

$$
\begin{aligned}
\|\dot\varepsilon\|_F^{2D} = \sqrt{\sum_{i,j}\dot\varepsilon_{ij}^2} 
    = \sqrt{\dot\varepsilon : \dot\varepsilon} 
    = \sqrt{\operatorname{tr}(\dot\varepsilon^2)} 
    = \sqrt{\dot\varepsilon_{11}^2 + \dot\varepsilon_{22}^2 + 2\dot\varepsilon_{12}^2}
    = \sqrt{\dot\varepsilon_{1}^2 + \dot\varepsilon_{2}^2} 
\end{aligned}
$$

This is the **Frobenius (Euclidean) norm** of the 2D strain-rate tensor, not the formal second invariant $I_2$.  
The Frobenius norm is always non-negative, whereas $I_2$ can be negative when
$\dot\varepsilon_{11}\dot\varepsilon_{22}-\dot\varepsilon_{12}^2<0$.  

The Frobenius norm is itself an invariant under orthogonal coordinate transformations.  
For a rotation $\mathbf{Q}$ satisfying $\mathbf{Q}^T\mathbf{Q}=\mathbf{I}$,

$$
\dot\varepsilon^{r} = \mathbf{Q}^T \dot\varepsilon \mathbf{Q}
\quad\Rightarrow\quad
\operatorname{tr}\!\left((\dot\varepsilon^{r})^2\right)
=
\operatorname{tr}\!\left(\dot\varepsilon^2\right).
$$

Thus, although the Frobenius norm is **not** the formal second invariant $I_2$, it is a valid scalar invariant of the strain-rate tensor.

---

## Infinitely many invariants

More generally, **any symmetric scalar function of the eigenvalues**
$f(\dot\varepsilon_1,\dot\varepsilon_2)$
**is an invariant**.  
Examples include

- $\dot\varepsilon_1+\dot\varepsilon_2$ (trace, $I_1$),
- $\dot\varepsilon_1\dot\varepsilon_2$ (determinant, $I_2$),
- $\sqrt{\dot\varepsilon_1^2+\dot\varepsilon_2^2}$ (Frobenius norm),
- $\frac{1}{2}\left(\dot\varepsilon_1-\dot\varepsilon_2\right)$ (maximum shear),
- $\frac{1}{2}\left(\dot\varepsilon_1+\dot\varepsilon_2\right)$ (mean horizontal strain rate),
- higher-order combinations (e.g., $\dot\varepsilon_1^4+\dot\varepsilon_2^4$).

This implies that there exists an **infinite family of invariants**, each capturing a different aspect of the strain-rate tensor.

---

## Maximum shear strain rate

Another commonly used measure is the maximum shear strain rate, which relates to the second invariant of the deviatoric strain-rate tensor $J_2$,

$$
\begin{aligned}
|\dot\gamma_{\max}^{2D}| &= \sqrt{J_2}, \\
\dot\gamma_{\max}^{2D}
&=
\tfrac12
\sqrt{(\dot\varepsilon_{11}-\dot\varepsilon_{22})^2 + 4\dot\varepsilon_{12}^2}
=
\tfrac12(\dot\varepsilon_1-\dot\varepsilon_2).
\end{aligned}
\label{eq:gamma_max}
$$

The principal strain rates $\dot\varepsilon_{1,2}$ can be expressed in terms of the dilatation rate and max. shear strain rate:

$$
\dot\varepsilon_{1,2} = \tfrac12\left(\dot\Delta \pm {\color{red}{2}}\dot\gamma_{\max}\right),
$$

**Note (erratum).** The factor of ${\color{red}{2}}$ was inadvertently omitted in the corresponding equation reported in the appendix of [this paper](https://doi.org/10.1029/2025JB031738) on strain rates across the Alpine–Himalayan belt. This was an oversight on my part. 

---

## Relationship between invariants and norms

The Frobenius norm relates to $I_1$ and $J_2$ as

$$
\|\dot\varepsilon\|_F^{2D}
=
\sqrt{2J_2+\tfrac12 I_1^2}.
\label{eq:frob_relation}
$$

For a purely deviatoric tensor $\dot\varepsilon'$, the Frobenius norm is just a scaled version of the max. shear strain rate:

$$
\|\dot\varepsilon'\|_F^{2D}
=
\sqrt{2J_2}
=
\sqrt{2}\,\dot\gamma_{\max}^{2D}.
$$

---

## A common confusion: principal strain rates vs. maximum shortening and extension

A recurring conceptual mistake in strain-rate studies is the implicit conflation of the *principal strain rates* $(\dot\varepsilon_1,\dot\varepsilon_2)$ with the notions of *maximum extension* and *maximum shortening*.

By definition, $\dot\varepsilon_1$ and $\dot\varepsilon_2$ are simply the eigenvalues of the strain-rate tensor, conventionally ordered such that
$\dot\varepsilon_1 \ge \dot\varepsilon_2$.
This ordering is mathematical, not physical: it does not, by itself, specify whether a principal strain rate direction is extensional or compressional.

In particular:

- $\dot\varepsilon_1$ represents the **maximum principal strain rate**, not necessarily the maximum *extensional* strain rate.
- $\dot\varepsilon_2$ represents the **minimum principal strain rate**, not necessarily the maximum *shortening* strain rate.

**Whether a given principal strain rate direction corresponds to extension or shortening depends on the sign of its associated eigenvalue, not on its index.**

Specifically:
- If $\dot\varepsilon_1 > 0$ and $\dot\varepsilon_2 < 0$, the principal direction associated with $\dot\varepsilon_1$ corresponds to maximum extension, and that associated with $\dot\varepsilon_2$ corresponds to maximum shortening.
- If $\dot\varepsilon_1 > 0$ and $\dot\varepsilon_2 > 0$, the strain rate tensor is purely extensional. Both principal strain rates correspond to extension, and no shortening direction exists.
- If $\dot\varepsilon_1 < 0$ and $\dot\varepsilon_2 < 0$, the strain rate tensor is purely compressional. Both principal strain rates correspond to shortening; $\dot\varepsilon_2$ represents the maximum shortening direction, while $\dot\varepsilon_1$ represents the least compressive direction, and no extensional direction exists. 

**Note (erratum)** In the appendix of [this paper](https://doi.org/10.1029/2025JB031738), I incorrectly stated that when $\dot\varepsilon_{\max}<0$, the azimuth $Az_{\dot\varepsilon_{\max}}$ represents the direction of maximum shortening. In a purely compressional regime, $\dot\varepsilon_{\max}=\dot\varepsilon_1$ corresponds to the least compressive direction, whereas maximum shortening is associated with $\dot\varepsilon_2$. 

In all cases, associating $\dot\varepsilon_1$ with extension and $\dot\varepsilon_2$ with shortening without explicitly evaluating their signs can lead to incorrect physical interpretations.

A robust practice is therefore to:
1. compute the principal strain rates,
2. examine their signs,
3. and only then assign the labels *maximum extension* and *maximum shortening* accordingly.

**Note (erratum)**. A confusion led me to implicitly associate $\dot\varepsilon_1$ with maximum extension and $\dot\varepsilon_2$ with maximum shortening in [this paper](https://doi.org/10.1029/2025JB031738). Fortunately, these issues were purely notational rather than methodological. In the analysis itself, I explicitly evaluated the signs of the eigenvalues to determine whether the principal directions represented extension or shortening. Therefore, all results presented are correct despite the unfortunate notational confusion throughout the paper. 

---

## Principal strain-rate azimuths

The azimuth of the eigenvector associated with the **maximum principal strain rate**
(i.e., the largest eigenvalue) is given by

$$
Az_{\dot\varepsilon_{\max}}
=
90^\circ
-
\tfrac12
\arctan\!\left(
\frac{2\dot\varepsilon_{12}}
{\dot\varepsilon_{11}-\dot\varepsilon_{22}}
\right).
\label{eq:azimuth}
$$

This azimuth describes the orientation of a **principal axis** (an eigenvector of the
strain-rate tensor) and is therefore a *geometric* property of the tensor itself.
By construction, it does not encode whether deformation along this axis corresponds
to extension or shortening. Consequently, $Az_{\dot\varepsilon_{\max}}$ should be interpreted strictly
as the orientation of the principal axis associated with the largest eigenvalue.
Assigning it to either extension or shortening requires an explicit evaluation of
the signs of the principal strain rates.

<br>

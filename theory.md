---
author:
- Charles Turner
bibliography:
- library.bib
date: 2024-11-08
title: Detecting an ocean transient response to anthropogenic
  atmospheric forcings
---

# Trend Emergence {#A:trendEmergence}

At present, there are two major techniques used in the statistical
evaluation of trend emergence:

-   An autocorrelation and observational uncertainty based approach,
    formulated by [@Tiao1990] and [@Weatherhead1998].

-   An ensemble standard deviation based approach. Examples include
    [@McKinley2016], [@Deser2012].

Both approaches are useful for their respective use cases, but neither
are suitable work which seeks to blend observations and modelling: we do
not have ensembles of observations, nor do models output measurement
uncertainties. In addition, neither is of any use whatsoever if we wish
to estimate trend emergence from a single model run.

To get around this, I have developed a new approach, which aims to be
consistent with both.

## The Ensemble Method

Consider a Gaussian distributed variable, $X$, with standard deviation
$\sigma$ and mean $\mu$. If we have $N$ observations of $X$ with mean
$\bar{x}$ and standard deviation $s$, we can say with 95% confidence
that $\mu \neq 0$ if: 

$$\bar{x} \ge \pm \frac{2s}{\sqrt{N-1}}$$

 This is simply the formula for the standard
error of the mean, using Central Limit Theorem. This can be rearranged
as follows:

$$\frac{\bar{x}}{s} = \mathrm{SNR} \ge \pm 2 \frac{1}{\sqrt{N - 1}}$$

where the ratio of mean to standard deviation where the ratio of mean to standard deviation
is the signal-to-noise ration, SNR.

This SNR statistic forms the core of the large ensemble approach. The
ensemble standard deviation, $s$, is assumed to be perfectly
representative of the variable standard deviation, $\sigma$, and
likewise for the ensemble mean. $\mu$ is therefore significantly
different from 0 at a $2\sigma$ level if SNR exceeds
$\pm 2 / \sqrt{N - 1}$, and the trend has therefore emerged from
statistical noise.

It is not possible to compute Equation
[2](#eqn-AstdErrSNR) for $N=1$, ie. a single run. As such, a
bootstrap technique, which preserves the power spectrum of the variable
timeseries, and thus its autocorrelation, is used.

## Bootstrap procedure

The standard bootstrapping technique is to take a sample of observations
of a variable $X$, labelled, $\{x_i\}$, and repeatedly subsample this
series of observations to obtain a new mean. By subsampling enough times
(nominally 10,000), the standard deviation of the mean can be obtained,
assuming that the initial sample is representative of the true
distribution of the variable. In the context of determining trend
emergence the process is more involved, but follows the same principle.

It is performed as follows, for a single continuous variable $y(t)$,
though of course observations and models are both discrete. Extension to
the discrete case is straightforward.

$y(t)$ is assumed to be composed of two components: one which increases
linearly with time, and a noise component which may be decomposed as a
Fourier series. Hence,

$$y(t) = \gamma t+ f(t) = \gamma t + \int_{-\infty}^{\infty} \frac{d \omega}{2 \pi} f (\omega) e ^{ - i \omega t} ,
  \label{eqn:AFourierDecomp1}$$ 

where $\gamma$ is the gradient of the linearly increasing component, and $f(t)$ the noise component. Each frequency component, $f(\omega)$, is given by

$$f(\omega) = \int_{-\infty}^{\infty} dt f(t) e ^{i \omega t}
  \label{eqn:AfOmegaDef}$$ 

Now, each frequency component has its phase
shuffled by a random amount $\phi(\omega)$, where $\phi(\omega)$ is
uniformly distributed in the range $[-\pi,\pi]$, for a new noise
component in frequency space, $f'(\omega)$:

$$f'(\omega) = \int_{-\infty}^{\infty} dt f(t) e ^{i(\omega t + \phi (\omega))}
  \label{eqn:AfPrimeOmegaDef}$$ 

Now, inverse transforming $f'(\omega)$ to obtain a phase suffled noise component, $f'(t)$, we see

$$f'(t) = \mathcal{F}^{-1}\big(f'(\omega)\big) = \int_{-\infty}^{\infty} \frac{dt d\omega'}{2 \pi} f(t) e ^{
  i ( \omega t + \phi (\omega) - \omega' t)} ,
  \label{eqn:AInvFourierTransform}$$ 

which, noting the identity

$$\int_{-\infty}^{\infty} dk e ^ {ik ( x - x')} = 2\pi \delta (x - x'),
  \label{eqn:ADiracDeltaIdendity}$$ 

  gives

$$f'(t) = \int_{-\infty}^{\infty} \frac{d \omega}{2 \pi} f (\omega) e ^{ - i \omega t} e ^{i \phi (\omega)},
  \label{eqn:AfPrimeTime}$$ 

  and so

$$y'(t) = \gamma t+ f'(t) = \gamma t + \int_{-\infty}^{\infty} \frac{d \omega}{2 \pi} f (\omega) e ^{ - i \omega t} e ^{i \phi (\omega)}
  \label{eqn:AyPrimeTime}$$ 

We now have a new variable of time $y'(t)$, which has the same linearly increasing component. 
The power spectrum of $f'(t)$ and $f(t)$ are identical, as the power spectrum is phase
independent.

This procedure forms the basis of the bootstrap. By performing this
procedure a number of times on the original timeseries, it is possible
to obtain a large number of new timeseries, with the same linear trend
but different noise, $y_i'(t)$.

Now consider sampling $y_i'(t)$ over between times $t_0,t_1$, before
estimating the mean gradient of $y_i'(t)$ over the interval. In the
limit $t_1 \rightarrow -\infty$ and $t_2 \rightarrow \infty$, the mean
gradient will always be $\gamma$, as all periodic variations will
average out. However, as the time interval $\Delta t = t_1 - t_0$
decreases, Fourier components with a period longer than $\Delta t$ will
cause the mean gradient to diverge from $\gamma$.

In the case that $\gamma$ is not known, we can estimate it by simply
linearly detrending $y$ using a least squares regression<sup>1</sup>. In this
procedure, the signal we wish to test for statisitical significance is
then $\gamma\Delta t$: that is, we wish to test whether our estimate of
the beginning to end difference can be explained solely by noise.

To do so, we fit linear trends to each member of the noise timeseries
ensemble: $f_i'(t)$. Unlike $f(t)$, these linear trends will not be by
definition zero, but will be Gaussian distributed about zero. Then, we
take the width of this distribution, $\sigma$, and compare it to the
value of $\gamma \Delta t$. This gives us statistical significance at
95% confidence if $$\gamma \Delta t \ge 2 \sigma $$ Note that we now have no factor of
$\sqrt{N - 1}$ as Central Limit Theorem is not used to estimate the
error of the mean. Thus, once the initial timeseries has been
bootstrapped to obtain a large ensemble, the procedure is essentially
the same as the large ensemble method. This section seeks to describe
the algorithm qualitatively, whilst giving mathematical detail. The original
implementation in MATLAB can be found at
<https://github.com/charles-turner-1/FFT_bootstrap>.

## Consistency with Autocorrelation

The autocorrelation based approach to determining trend emergence is not
described in detail here. However, it seeks to account for the reduction
of new information in successive measurements due to autocorrelation.
Consider the correlation between two points of the same function:

$$\langle A(x_2)A(x_1)\rangle = C_{AA}(x_2,x_1) = C_{AA}(x_1,x_2)
  \label{eqn:AAutoCorr}$$ 

If we now insert the Fourier representations for $A(x_1),A(x_2)$, we obtain

$$C_{AA}(x_1,x_2) =  \Big\langle \sum_{k_2} A_{k_2} e^{ik_2 x_2 } \sum_{k_1} A_{k_1} e^{ik_1 x_1 }\Big\rangle 
  \label{eqn:AAutoCorrFourier1}$$ 

  Combining the sums, we obtain

$$C_{AA}(x_1,x_2) =  \Big\langle \sum_{k_2} \sum_{k_1} A_{k_2} A_{k_1} e^{ik_2 x_2 }  e^{ik_1 x_1 }\Big\rangle ,
  \label{eqn:AAutoCorrFourier2}$$ 

   which we may rewrite using $r = x_2 - x_1$ as

$$C_{AA}(x_1,x_2) =  \Big\langle \sum_{k_2} \sum_{k_1} A_{k_2} A_{k_1} e^{ik_2r} e^{i(k_1+k_2)x_1}\Big\rangle
  \label{eqn:AAutoCorrFourier3}$$ 

In order to recover an autocorrelation, the correlation must only depend on the distance
between the two points, not the points themselves. Therefore, it must be
the case that if $k_1 +k_2 \neq 0$, $\langle A_{k_1}A_{k_2}\rangle = 0$,
and therefore:

$$\langle A_{k_1}A_{k_2} \rangle = \langle A_{k_1}A_{k_{-1}}\rangle \delta_{k_1 + k_2}=\langle A_{k_1}A^*_{k_{1}}\rangle \delta_{k_1 + k_2},
  \label{eqn:ADiscreteDiracDelta}$$ 

  which may also be observed as a corollary of Equation
[](eqn-diracdelta). Therefore, Fourier coefficients at
differing wavenumbers must be uncorrelated: the random phase shuffling
of the bootstrap respects this requirement.

Now, rewriting Equation
[\[eqn:AAutoCorrFourier3\]](#eqn:AAutoCorrFourier3){reference-type="ref"
reference="eqn:AAutoCorrFourier3"} observing the constraints of Equation
[\[eqn:ADiscreteDiracDelta\]](#eqn:ADiscreteDiracDelta){reference-type="ref"
reference="eqn:ADiscreteDiracDelta"}, we obtain 

$$\begin{split}
  C_{AA}(r) =&  \sum_{k_1} \sum_{k_{-1}} \langle A_{k_1}A_{k_{-1}} \rangle e ^{-i k_1 r} \\
  =& \sum_{k}\sum_{k}\langle A_k A_k^*\rangle e ^{ikr} = \sum_{k}\sum_{k}\langle |A_k|^2  \rangle e ^{ikr} ,
\end{split}
  \label{eqn:APowerSpectrum} $$

having now moved the angled brackets
inside the sum and removed the unecessary numbered subscripts. From
this, we can see that the 1<sup>st</sup> order autocorrelation function is the
Fourier transform of the power spectrum.

Therefore, as the bootstrap method preserves the power spectrum, it also
preserves the autocorrelation, and so this method will give results
comparable to the autocorrelation based estimates of the SNR. However,
they will not be identical, as the bootstrap method does not consider
observational uncertainties as the autocorrelation method does.

___

1: It is therefore the case that the expected gradient of $f(t)$ is 0.

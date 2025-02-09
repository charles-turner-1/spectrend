"""
A toolbox for determining the presence of trends using spectral analysis
"""

import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft

__all__ = ["fft", "ifft"]


class TimeSeries:
    """
    Class for time series data.
    """

    def __init__(self, data: np.ndarray):
        """
        Create a spectrend.TimeSeries object from a one-dimensional numpy array.

        Parameters
        ----------
        data : np.ndarray
            One-dimensional numpy array containing the time series data.
        """
        self._input_ts = np.squeeze(data)
        self._detrended_ts = self._detrend()

        self._freq_domain: np.ndarray | None = None

    @property
    def input_ts(self) -> np.ndarray:
        """
        View the input time series data. No setter is provided for this data -
        time series classes can hold only one time series at a time.
        """
        return self._input_ts

    @property
    def detrended_ts(self) -> np.ndarray:
        """
        View the detrended time series data. No setter is provided for this data -
        detrended time series data is generated from the input time series data.
        """
        try:
            getattr(self, "_detrended_ts")
        except AttributeError:
            self._detrended_ts = self._detrend()

        return self._detrended_ts

    @property
    def freq_domain(self) -> np.ndarray:
        """
        Get the frequency domain representation of the time series data.
        """
        try:
            getattr(self, "_freq_domain")
        except AttributeError:
            self._freq_domain = fft(self.detrended_ts)

        return self._freq_domain

    @property
    def power_spectrum(self) -> np.ndarray:
        """
        Get the power spectrum of the time series data.
        """
        try:
            getattr(self, "_power_spectrum")
        except AttributeError:
            self._power_spectrum = np.abs(self.freq_domain) ** 2

        return self._power_spectrum

    def _detrend(self) -> np.ndarray:
        """
        Linearly detrend a time series, making it stationary & so amenable to spectral
            analysis. Additionally, cache the trend for later use.
        """
        detrended = signal.detrend(self._input_ts)
        start_diff = self._input_ts[0] - detrended[0]

        self.linear_slope = np.diff(self.input_ts - detrended).mean()

        return detrended + start_diff

    def plot_ts(
        self,
        kwargs_input: dict | None = None,
        kwargs_detrended: dict | None = None,
        **kwargs,
    ) -> tuple[plt.Figure, plt.Axes]:
        """
        Plot the time series data.

        Parameters
        ----------
        **kwargs : dict
            Keyword arguments to be passed to the plot function.
        **kwargs_input: dict
            Keyword arguments to be passed to the plot function for the raw time series.
        **kwargs_detrended : dict
            Keyword arguments to be passed to the plot function for the detrended time series.

        """
        kwargs_input, kwargs_detrended = kwargs_input or {}, kwargs_detrended or {}

        kwargs_input = {**kwargs, **kwargs_input}
        kwargs_detrended = {**kwargs, **kwargs_detrended}

        fig, ax = plt.subplots()
        ax.plot(self._input_ts, label="Input Time Series", **kwargs_input)
        ax.plot(self.detrended_ts, label="Detrended Time Series", **kwargs_detrended)

        return fig, ax

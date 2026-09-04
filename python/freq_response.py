import numpy as np
from scipy import signal


def _resolve_deskew(deskew, fs):
    """Turn the ``deskew`` argument into a delay in seconds to *advance* H by.

    Accepts a number (samples), ``'zoh'``/``True`` (half a sample, the
    zero-order-hold group delay of a sample-and-held stimulus), or a falsey /
    ``'off'`` value for no correction. There is deliberately no 'auto' mode:
    near a resonance the sampling transport delay and the circuit's own group
    delay are the same measurement, so the artifact can only be removed from
    its known analytic value (T/2 for a once-per-timestep sample and hold)."""
    if deskew is None or deskew is False:
        return 0.0
    if deskew is True:
        return 0.5 / fs
    if isinstance(deskew, str):
        key = deskew.strip().lower()
        if key in ('', 'off', 'none', 'false', '0'):
            return 0.0
        if key in ('zoh', 'half', 'true'):
            return 0.5 / fs
        return float(key) / fs
    return float(deskew) / fs


def freq_response(x, y, fs=1.0, nperseg=1024, window='flattop',
                      scaling='spectrum', detrend='constant', deskew=0.0):
    """
    Estimate frequency response H(f) = Y(f)/X(f) using cross power spectral density

    Parameters:
    -----------
    x : array_like
        Input signal (reference / stimulus)
    y : array_like
        Output signal (response)
    fs : float, optional
        Sampling frequency (Hz). Default = 1.0
    nperseg : int, optional
        Length of each segment. Default = 1024
    window : str or tuple or array_like, optional
        Desired window type. Default = 'hann'
    scaling : {'spectrum', 'density'}, optional
        Selects between power spectrum ('spectrum') and power spectral density ('density')
    detrend : {'constant', 'linear', False}, optional
        Detrending applied to each segment
    deskew : number or {'zoh', 'off'}, optional
        Remove a bulk transport delay from the estimate (a phase that droops
        linearly with frequency). A fixed-step ODE testbench that samples and
        holds its stimulus once per timestep feeds the Welch estimator a
        half-sample (T/2) zero-order-hold group delay that the true circuit
        does not have; ``'zoh'`` (or ``True``) removes exactly that, a number
        removes that many samples. Default 0.0 (no correction).

    Returns:
    --------
    f : ndarray
        Array of frequency bins
    H : ndarray (complex)
        Frequency response estimate H(f) = Sxy(f) / Sxx(f)
    coh : ndarray
        Magnitude-squared coherence (0...1)
    """
    # scipy's spectral estimators index inputs with tuples internally
    # (e.g. x[..., i0:i1]), which pandas Series/DataFrame don't support.
    x = np.asarray(x)
    y = np.asarray(y)

    # Compute cross power spectral density and auto power spectral density
    f, Pxx = signal.csd(x,x, fs=fs, window=window, nperseg=nperseg,
                         scaling=scaling, detrend=detrend, axis=-1)
    
    f, Pxy = signal.csd(x, y, fs=fs, window=window, nperseg=nperseg,
                       scaling=scaling, detrend=detrend, axis=-1)
    
    # Frequency response H1 estimator: H(f) = Sxy(f) / Sxx(f)
    H = Pxy / Pxx
    
    # Coherence (useful for quality assessment)
    f, Pyy = signal.welch(y, fs=fs, window=window, nperseg=nperseg,
                         scaling=scaling, detrend=detrend)
    coh = np.abs(Pxy)**2 / (Pxx * Pyy)

    tau = _resolve_deskew(deskew, fs)
    if tau:
        H = H * np.exp(2j * np.pi * f * tau)
        print(f"    freq_response: de-skewed {tau * fs:+.3f} sample "
              f"({tau * 1e9:+.0f} ns) bulk delay out of the phase")

    return f, H, coh


# ─────────────────────────────────────────────────────────────────────
# Quick usage example:
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    # Create example signals
    fs = 1000.0
    t = np.arange(0, 10000, 1/fs)
    x = np.random.randn(len(t))
    # System: 2nd order low-pass with resonance at ~80 Hz
    y = signal.lfilter(*signal.butter(8, 150/(fs/2), btype='low'), x) + np.random.randn(len(t))
    
    # Calculate frequency response
    f, H, _ = freq_response(x, y, fs=fs, nperseg=250)

    # Plot
    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, sharex=True, figsize=(12, 8))

    ax_mag.semilogx(f, 20*np.log10(np.abs(H)))
    ax_mag.grid(True, alpha=0.3)
    ax_mag.set_title("Magnitude |H(f)|")
    ax_mag.set_ylabel("Gain")
    # ax_mag.set_xlim(0, 300)

    ax_phase.semilogx(f, np.angle(H, deg=True))
    ax_phase.grid(True, alpha=0.3)
    ax_phase.set_title("Phase")
    ax_phase.set_xlabel("Frequency [Hz]")
    ax_phase.set_ylabel("Phase [°]")

    plt.tight_layout()
    plt.show()

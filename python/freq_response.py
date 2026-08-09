import numpy as np
from scipy import signal


def freq_response(x, y, fs=1.0, nperseg=1024, window='flattop', 
                      scaling='spectrum', detrend='constant'):
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

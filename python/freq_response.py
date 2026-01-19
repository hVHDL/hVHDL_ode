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
    y = signal.lfilter(*signal.butter(8, 150/(fs/2), btype='low'), x) + 0.2*np.random.randn(len(t))
    
    # Calculate frequency response
    f, H, coh = freq_response(x, y, fs=fs, nperseg=1000)
    
    # Plot
    plt.figure(figsize=(12, 8))
    
    plt.subplot(3,1,1)
    plt.semilogy(f, np.abs(H))
    plt.grid(True, alpha=0.3)
    plt.title("Magnitude |H(f)|")
    plt.ylabel("Gain")
    plt.xlim(0, 300)
    
    plt.subplot(3,1,2)
    plt.plot(f, np.angle(H, deg=True))
    plt.grid(True, alpha=0.3)
    plt.title("Phase")
    plt.ylabel("Phase [°]")
    plt.xlim(0, 300)
    
    plt.subplot(3,1,3)
    plt.semilogy(f, coh)
    plt.grid(True, alpha=0.3)
    plt.title("Coherence")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Coherence")
    plt.xlim(0, 300)
    plt.ylim(0, 1.1)
    
    plt.tight_layout()
    plt.show()

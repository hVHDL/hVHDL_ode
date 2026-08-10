import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def read_plot_config(filename):
    # VHDL testbenches may write lines like "#CONFIG T_title=Currents (A)"
    # to configure titles/labels. Recognized keys: title, xlabel,
    # <prefix>_title, <prefix>_ylabel (e.g. T_title, B_ylabel),
    # label_<column_name> to rename a signal's legend entry
    # (e.g. label_T_i0=Inductor current (A)), xlim/T_ylim/B_ylim for the
    # time-domain plot, freq_pair_<name>, freq_fs, freq_nperseg (or
    # freq_num_windows to split the data into N segments instead of
    # specifying segment length directly), freq_xlim, phase_ylim,
    # freq_title, freq_unwrap_phase=true for frequency response plots
    # (see parse_freq_pairs/plot_freq_response), and combined_layout=true
    # to put the time-domain and frequency response plots in one 2x2 figure.
    config = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#CONFIG"):
                continue
            entry = line[len("#CONFIG"):].strip()
            if "=" not in entry:
                continue
            key, value = entry.split("=", 1)
            config[key.strip()] = value.strip()
    return config

def parse_freq_pairs(config):
    # "freq_pair_<name>=<x_column>,<y_column>" declares a signal pair to
    # compute the frequency response H(f) = Y(f)/X(f) for.
    pairs = []
    for key, value in config.items():
        if not key.startswith("freq_pair_"):
            continue
        name = key[len("freq_pair_"):]
        parts = [p.strip() for p in value.split(",")]
        if len(parts) == 2:
            pairs.append((name, parts[0], parts[1]))
        else:
            print(f"Warning: freq_pair_{name} must be 'x_column,y_column', got '{value}'. Skipping.")
    return pairs

def is_true(value):
    return str(value).strip().lower() in ("1", "true", "yes")

def parse_limit_pair(config, key):
    # "<key>=low,high", e.g. "xlim=0,1e-3" or "T_ylim=-5,5"
    if key not in config:
        return None
    try:
        lo, hi = (float(v) for v in config[key].split(","))
        return lo, hi
    except ValueError:
        print(f"Warning: {key} must be 'low,high', got '{config[key]}'. Ignoring.")
        return None

def plot_time_domain(loaded, config, signal_labels, ax_top, ax_bottom):
    for filename, df in loaded:
        # Split columns into top and bottom based on the prefix
        top_columns = [col for col in df.columns if col.startswith("T_")]
        bottom_columns = [col for col in df.columns if col.startswith("B_")]

        # Plot top columns
        if top_columns:
            df[top_columns].rename(columns=signal_labels).plot(ax=ax_top, title=config.get("T_title", "Top Data Plot"), legend=True)
            ax_top.set_ylabel(config.get("T_ylabel", "Top Data Values"))
            ax_top.set_xlabel("")
            ax_top.grid(True)

        # Plot bottom columns
        if bottom_columns:
            df[bottom_columns].rename(columns=signal_labels).plot(ax=ax_bottom, title=config.get("B_title", "Bottom Data Plot"), legend=True)
            ax_bottom.set_ylabel(config.get("B_ylabel", "Bottom Data Values"))
            ax_bottom.grid(True)

    ax_bottom.set_xlabel(config.get("xlabel", "Time"))

    xlim = parse_limit_pair(config, "xlim")
    if xlim:
        ax_bottom.set_xlim(*xlim)

    t_ylim = parse_limit_pair(config, "T_ylim")
    if t_ylim:
        ax_top.set_ylim(*t_ylim)

    b_ylim = parse_limit_pair(config, "B_ylim")
    if b_ylim:
        ax_bottom.set_ylim(*b_ylim)

def plot_freq_response(loaded, config, freq_pairs, signal_labels, ax_mag, ax_phase):
    from freq_response import freq_response

    fs = float(config["freq_fs"])
    unwrap_phase = is_true(config.get("freq_unwrap_phase", "false"))

    for filename, df in loaded:
        # freq_num_windows=N splits this file's data into N segments;
        # takes priority over freq_nperseg (default 1024) when set.
        if "freq_num_windows" in config:
            num_windows = max(1, int(float(config["freq_num_windows"])))
            nperseg = max(1, len(df) // num_windows)
        else:
            nperseg = int(config.get("freq_nperseg", 1024))

        for name, x_col, y_col in freq_pairs:
            if x_col not in df.columns or y_col not in df.columns:
                continue
            try:
                f, h, _ = freq_response(df[x_col], df[y_col], fs=fs, nperseg=nperseg)
            except Exception as e:
                print(f"Could not compute frequency response for '{name}' in '{filename}': {e}. Skipping.")
                continue

            label = signal_labels.get(name, name)
            if len(loaded) > 1:
                label = f"{label} ({filename})"
            phase = np.angle(h)
            if unwrap_phase:
                phase = np.unwrap(phase)
            ax_mag.semilogx(f, 20 * np.log10(np.abs(h)), label=label)
            ax_phase.semilogx(f, np.degrees(phase), label=label)

    freq_xlim = parse_limit_pair(config, "freq_xlim")
    if freq_xlim:
        ax_mag.set_xlim(*freq_xlim)

    phase_ylim = parse_limit_pair(config, "phase_ylim")
    if phase_ylim:
        ax_phase.set_ylim(*phase_ylim)

    ax_mag.set_title(config.get("freq_title", "Frequency Response"))
    ax_mag.set_ylabel("Magnitude [dB]")
    ax_mag.grid(True)
    ax_mag.legend()

    ax_phase.set_ylabel("Phase [deg]")
    ax_phase.set_xlabel("Frequency [Hz]")
    ax_phase.grid(True)
    ax_phase.legend()

def plot_data(filenames):
    loaded = []
    config = {}
    suptitle = None

    # Load each file's data and #CONFIG lines
    for filename in filenames:
        try:
            file_config = read_plot_config(filename)
            config.update(file_config)
            if "title" in file_config:
                suptitle = file_config["title"]

            # Load the data file into a DataFrame, ignoring #CONFIG lines
            df = pd.read_csv(filename, sep=r'\s+', comment='#')

            # Check if 'time' column is present
            if 'time' not in df.columns:
                print(f"Error: 'time' column is missing in file '{filename}'. Skipping this file.")
                continue

            # Set 'time' as the index for the plot
            df.set_index('time', inplace=True)
            loaded.append((filename, df))

        except FileNotFoundError:
            print(f"Error: File '{filename}' not found. Skipping this file.")
        except pd.errors.EmptyDataError:
            print(f"Error: The file '{filename}' is empty. Skipping this file.")
        except Exception as e:
            print(f"An error occurred with file '{filename}': {e}. Skipping this file.")

    # Signal display names, e.g. "label_T_i0=Inductor current (A)"
    signal_labels = {
        key[len("label_"):]: value
        for key, value in config.items()
        if key.startswith("label_")
    }

    freq_pairs = parse_freq_pairs(config)
    has_freq = bool(freq_pairs)
    if freq_pairs and "freq_fs" not in config:
        print("Warning: freq_pair_* config found but freq_fs is not set. Skipping frequency response plot.")
        has_freq = False

    combined = has_freq and is_true(config.get("combined_layout", "false"))

    if combined:
        # 2x2 figure: time-domain on the left column, frequency response on the right
        fig, ((ax_top, ax_mag), (ax_bottom, ax_phase)) = plt.subplots(2, 2, sharex='col', figsize=(12, 6))
    else:
        fig, (ax_top, ax_bottom) = plt.subplots(2, 1, sharex=True, figsize=(7, 5))

    plot_time_domain(loaded, config, signal_labels, ax_top, ax_bottom)

    if suptitle:
        fig.suptitle(suptitle)

    if has_freq and combined:
        plot_freq_response(loaded, config, freq_pairs, signal_labels, ax_mag, ax_phase)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
    else:
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        if has_freq:
            fig2, (ax_mag, ax_phase) = plt.subplots(2, 1, sharex=True, figsize=(7, 5))
            plot_freq_response(loaded, config, freq_pairs, signal_labels, ax_mag, ax_phase)
            plt.tight_layout()

    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_data.py <filename1> <filename2> ...")
    else:
        filenames = sys.argv[1:]  # List of all provided filenames
        plot_data(filenames)

import os
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
    # freq_title, freq_unwrap_phase=true, freq_save_<name>=<path> to save
    # that pair's computed response to a file (see save_freq_response,
    # is_saved_freq_response) for overlaying onto a later plot, for
    # frequency response plots (see parse_freq_pairs/plot_freq_response),
    # and combined_layout=true to put the time-domain and frequency
    # response plots in one 2x2 figure.
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

def is_saved_freq_response(filename):
    # Detects a saved frequency response, as opposed to a simulation .dat
    # file: either the marker written by save_freq_response, or a bare
    # "frequency,magnitude_db,phase_deg" CSV (e.g. hand-made or exported
    # from elsewhere), skipping any leading #-comment lines.
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line == "#FREQRESPONSE":
                return True
            if line.startswith("#"):
                continue
            return line.replace(" ", "").lower() == "frequency,magnitude_db,phase_deg"
    return False

def save_freq_response(path, label, f, mag_db, phase_deg):
    with open(path, "w") as out:
        out.write("#FREQRESPONSE\n")
        out.write(f"#CONFIG label={label}\n")
        out.write("frequency,magnitude_db,phase_deg\n")
        for freq, mag, ph in zip(f, mag_db, phase_deg):
            out.write(f"{freq},{mag},{ph}\n")

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

def plot_freq_response(loaded, config, freq_pairs, signal_labels, saved_responses, ax_mag, ax_phase):
    from freq_response import freq_response

    unwrap_phase = is_true(config.get("freq_unwrap_phase", "false"))

    if freq_pairs:
        for filename, df in loaded:
            if "freq_fs" in config:
                fs = float(config["freq_fs"])
            else:
                # No freq_fs configured: infer the sample rate from the
                # time column itself, so a bode plot can be requested
                # (e.g. via --freq-pair) with no #CONFIG lines at all.
                fs = 1.0 / np.median(np.diff(df.index.to_numpy()))

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
                mag_db = 20 * np.log10(np.abs(h))
                phase_deg = np.degrees(phase)
                ax_mag.semilogx(f, mag_db, label=label)
                ax_phase.semilogx(f, phase_deg, label=label)

                save_path = config.get(f"freq_save_{name}")
                if save_path:
                    save_freq_response(save_path, label, f, mag_db, phase_deg)
                    print(f"Saved frequency response for '{name}' to '{save_path}'")

    # Previously saved responses (see save_freq_response), overlaid for comparison
    for label, f, mag_db, phase_deg in saved_responses:
        ax_mag.semilogx(f, mag_db, label=label, linestyle="--")
        ax_phase.semilogx(f, phase_deg, label=label, linestyle="--")

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

def plot_data(filenames, config_overrides=None, freq_only=False):
    loaded = []
    saved_responses = []
    config = {}
    suptitle = None

    # Load each file's data and #CONFIG lines
    for filename in filenames:
        try:
            if is_saved_freq_response(filename):
                file_config = read_plot_config(filename)
                fr_df = pd.read_csv(filename, comment='#')
                label = file_config.get("label", os.path.basename(filename))
                saved_responses.append((label, fr_df["frequency"].to_numpy(), fr_df["magnitude_db"].to_numpy(), fr_df["phase_deg"].to_numpy()))
                continue

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

    # Command-line overrides (e.g. --save-freq) take precedence over the
    # testbench's own #CONFIG lines.
    if config_overrides:
        config.update(config_overrides)

    # Signal display names, e.g. "label_T_i0=Inductor current (A)"
    signal_labels = {
        key[len("label_"):]: value
        for key, value in config.items()
        if key.startswith("label_")
    }

    freq_pairs = parse_freq_pairs(config)
    has_freq = bool(freq_pairs) or bool(saved_responses)

    # Nothing time-domain to show (e.g. only saved frequency-response
    # files were given): skip straight to the frequency-only layout.
    if not loaded and has_freq:
        freq_only = True

    if freq_only:
        if not has_freq:
            print("Warning: --freq-only was given but no frequency response data (freq_pair_* config or saved response files) was found.")
            return
        fig, (ax_mag, ax_phase) = plt.subplots(2, 1, sharex=True, figsize=(7, 5))
        plot_freq_response(loaded, config, freq_pairs, signal_labels, saved_responses, ax_mag, ax_phase)
        if suptitle:
            fig.suptitle(suptitle)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.show()
        return

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
        plot_freq_response(loaded, config, freq_pairs, signal_labels, saved_responses, ax_mag, ax_phase)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
    else:
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        if has_freq:
            fig2, (ax_mag, ax_phase) = plt.subplots(2, 1, sharex=True, figsize=(7, 5))
            plot_freq_response(loaded, config, freq_pairs, signal_labels, saved_responses, ax_mag, ax_phase)
            plt.tight_layout()

    plt.show()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plot hVHDL_ode simulation output files.")
    parser.add_argument("filenames", nargs="+", help="Simulation .dat files and/or saved frequency-response files")
    parser.add_argument("--save-freq", action="append", default=[], metavar="NAME=PATH",
                         help="Save the freq_pair_NAME frequency response to PATH, so it can be "
                              "overlaid on a later plot (repeatable). Same effect as the testbench's "
                              "freq_save_NAME config, but set from the command line instead.")
    parser.add_argument("--freq-only", action="store_true",
                         help="Skip the time-domain plot and show only the frequency response.")
    parser.add_argument("--freq-pair", action="append", default=[], metavar="NAME=X_COL,Y_COL",
                         help="Add a frequency response pair from X_COL to Y_COL (repeatable). Same "
                              "effect as the testbench's freq_pair_NAME config, but set from the "
                              "command line instead, so it works even without a #CONFIG line.")
    parser.add_argument("--freq-fs", metavar="HZ",
                         help="Sampling frequency for the frequency response. Optional: if omitted, "
                              "it's inferred from the file's time column.")
    args = parser.parse_args()

    config_overrides = {}
    for item in args.save_freq:
        if "=" not in item:
            print(f"Warning: --save-freq expects NAME=PATH, got '{item}'. Ignoring.")
            continue
        name, path = item.split("=", 1)
        config_overrides[f"freq_save_{name.strip()}"] = path.strip()

    for item in args.freq_pair:
        if "=" not in item:
            print(f"Warning: --freq-pair expects NAME=X_COL,Y_COL, got '{item}'. Ignoring.")
            continue
        name, cols = item.split("=", 1)
        config_overrides[f"freq_pair_{name.strip()}"] = cols.strip()

    if args.freq_fs:
        config_overrides["freq_fs"] = args.freq_fs

    plot_data(args.filenames, config_overrides, freq_only=args.freq_only)

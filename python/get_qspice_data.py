from PyQSPICE import clsQSPICE as pqs
import pandas as pd
import os
import sys
import ast  # For safely parsing string representation of a list

def get_qspice_data(qsh_file, variable_names, output_csv_path):
    """
    Simulates a QSPICE circuit and saves the results to a CSV file.

    Parameters:
        qsh_file (str): Path to the QSPICE `.qsh` file.
        variable_names (list of str): A list of strings specifying variables to extract (e.g., ["V(dc)", "I(L1)", "V(sampled_current)"]).
        output_csv_path (str): The path to save the output CSV file.
    """
    # Get the current working directory
    current_dir = os.getcwd()

    # Get the directory containing the .qsh file
    qsh_directory = os.path.dirname(os.path.realpath(qsh_file))

    # Change directory to the directory containing the .qsh file
    os.chdir(qsh_directory)

    # Initialize and prepare QSPICE simulation
    run = pqs(os.path.basename(qsh_file))  # Load the provided `.qsh` file
    run.qsch2cir()
    run.cir2qraw()
    run.setNline(4999)  # Adjust Nline if required for simulation length

    # Run QSPICE simulation and load specified variables
    df = run.LoadQRAW(variable_names)

    os.chdir(current_dir)


    # Save the DataFrame to a CSV file with space as the delimiter
    df.to_csv(output_csv_path, sep=' ', index=False)

    print(f"Simulation results saved to: {output_csv_path}")


if __name__ == "__main__":
    # Parse command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python get_qspice_data.py <qsh_file> <variable_names> <output_csv_path>")
        print("Example: python get_qspice_data.py ./boost_ref.qsh \"['V(dc)', 'I(L1)', 'V(sampled_current)']\" simulation_results.csv")
        sys.exit(1)

    qsh_file = sys.argv[1]  # First argument: Path to the .qsh file
    variable_names_str = sys.argv[2]  # Second argument: List of variables as a string
    output_csv_path = sys.argv[3]  # Third argument: Path to save the CSV file

    # Safely parse the variable names string into a list
    try:
        variable_names = ast.literal_eval(variable_names_str)
        if not isinstance(variable_names, list):
            raise ValueError
    except (ValueError, SyntaxError):
        print("Error: <variable_names> must be a valid Python list (e.g., \"['V(dc)', 'I(L1)', 'V(sampled_current)']\").")
        sys.exit(1)

    # Run the simulation and save the results
    get_qspice_data(qsh_file, variable_names, output_csv_path)


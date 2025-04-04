import matplotlib.pyplot as plt
import re
from collections import defaultdict

TESTS = {
    #'5893 0.05 3 10 10': ['0.002 0.035', '2'],
    #'8555 0.05 3 10 10': ['0.016 0.049', '1'],
    #'12 100 5 10000 10000': ['76.732 61.943', '2209'],
    #'-11 3500 20 500000 10': ['1984.878 1625.992', '35'],
    #'1 5000 100 1000000 4': ['3936.506 131.472', '4'],
    #'1 5000 100 1000000 100': ['3899.787 156.291', '163'],
    '1 5000 20 1000000 10': ['3918.912 143.364', '19'],
    '1 1000 3 10000 10000': ['287.788 261.446', '31'],
    '3 5000 50 1000000 300': ['3819.032 25.659', '469']
    #'3 5000 50 1000000 500': ['3738.436 58.743', '804'],
    #'-1 1000 30 100000 1000': ['575.878 370.663', '1203'],
    #'32 100 6 10000 100000': ['29.530 30.082', '4448']
}

NTASKS=[1, 2, 4, 8, 16, 32]
CPUS_PER_TASK=[1, 2, 4, 6, 8]

def parse_speedup_summary(file_path):
    """Parse the speedup_summary.txt file and extract execution times."""
    data = defaultdict(lambda: defaultdict(list))  # {test_params: {config: mpi_time}}
    current_config = None

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("Configuration:"):
                # Extract ntasks and ncpu from the configuration line
                match = re.search(r"(\d+) tasks, (\d+) CPUs per task", line)
                if match:
                    ntasks = int(match.group(1))
                    ncpu = int(match.group(2))
                    current_config = (ntasks, ncpu)
            elif line.startswith("Test:"):
                # Extract test parameters
                test_params = line.split(":")[1].strip()
            elif line.startswith("MPI time:"):
                # Extract MPI execution time
                mpi_time = float(line.split(":")[1].strip()[:-1])  # Remove 's' and convert to float
                if current_config and test_params:
                    data[test_params][current_config].append(mpi_time)

    return data

def plot_execution_time(data):
    """Plot execution time trends for each test."""
    for test_params, configs in data.items():
        plt.figure(figsize=(10, 6))
        configs_sorted = sorted(configs.items(), key=lambda x: (x[0][0], x[0][1]))  # Sort by ntasks, ncpu
        x_labels = [f"{nt} tasks, {ncpu} CPUs" for (nt, ncpu), _ in configs_sorted]
        y_values = [times[0] for _, times in configs_sorted]  # Extract the first MPI time for each config

        plt.plot(x_labels, y_values, marker='o', label=f"Test: {test_params}")
        plt.xlabel("Configuration (ntasks, ncpu)")
        plt.ylabel("Execution Time (s)")
        plt.title(f"Execution Time Trend for Test: {test_params}")
        plt.xticks(rotation=45, ha="right")
        plt.grid(True)
        plt.tight_layout()
        plt.legend()
        plt.savefig(f"results/plots/execution_time_{test_params.replace(' ', '_')}.png")
        plt.show()

def plot_execution_time_combined(data, selected_tests):
    """Plot execution time trends for selected tests in a single graph."""
    plt.figure(figsize=(12, 8))
    
    for test_params, configs in data.items():
        if test_params not in selected_tests:
            continue  # Skip tests that are not selected
        
        configs_sorted = sorted(configs.items(), key=lambda x: (x[0][0], x[0][1]))  # Sort by ntasks, ncpu
        x_labels = [f"{nt} tasks, {ncpu} CPUs" for (nt, ncpu), _ in configs_sorted]
        y_values = [times[0] for _, times in configs_sorted]  # Extract the first MPI time for each config

        plt.plot(x_labels, y_values, marker='o', label=f"Test: {test_params}")

    plt.xlabel("Configuration (ntasks, ncpu)")
    plt.ylabel("Execution Time (s)")
    plt.title("Execution Time Trends for Selected Tests")
    plt.xticks(rotation=45, ha="right")
    plt.grid(True)
    plt.tight_layout()
    plt.legend()
    plt.savefig("results/plots/execution_time_selected_tests.png")
    plt.show()

def plot_fixed_tasks(data, test_params, fixed_ntasks):
    """Plot execution time trend for a specific test with a fixed number of tasks while varying CPUs per task."""
    plt.figure(figsize=(10, 6))
    
    # Ensure the fixed_ntasks is valid
    if fixed_ntasks not in NTASKS:
        print(f"Invalid number of tasks: {fixed_ntasks}. Valid options are: {NTASKS}")
        return

    # Check if the test exists in the data
    if test_params not in data:
        print(f"Test '{test_params}' not found in the data.")
        return

    # Filter configurations for the specific test and fixed number of tasks
    configs = {
        (ntasks, ncpu): times
        for (ntasks, ncpu), times in data[test_params].items()
        if ntasks == fixed_ntasks
    }
    
    if not configs:
        print(f"No data found for test '{test_params}' with {fixed_ntasks} tasks.")
        return

    # Sort configurations by number of CPUs per task
    configs_sorted = sorted(configs.items(), key=lambda x: x[0][1])  # Sort by ncpu
    x_labels = []
    y_values = []

    # Iterate over CPUS_PER_TASK and check if data exists
    for ncpu in CPUS_PER_TASK:
        if (fixed_ntasks, ncpu) in configs:
            x_labels.append(f"{ncpu} CPUs")
            y_values.append(configs[(fixed_ntasks, ncpu)][0])  # Extract the first MPI time
        else:
            print(f"Skipping missing configuration: {fixed_ntasks} tasks, {ncpu} CPUs")

    if not x_labels:
        print(f"No valid configurations found for test '{test_params}' with {fixed_ntasks} tasks.")
        return

    # Plot the trend
    plt.plot(x_labels, y_values, marker='o', label=f"Test: {test_params}, {fixed_ntasks} tasks")
    plt.xlabel("CPUs per Task")
    plt.ylabel("Execution Time (s)")
    plt.title(f"Execution Time Trend for Test: {test_params} with {fixed_ntasks} Tasks")
    plt.xticks(rotation=45, ha="right")
    plt.grid(True)
    plt.tight_layout()
    plt.legend()
    plt.savefig(f"results/plots/execution_time_{test_params.replace(' ', '_')}_ntasks_{fixed_ntasks}.png")
    plt.show()

def main():
    # Parse the speedup summary file
    file_path = "results/speedup_summary.txt"
    data = parse_speedup_summary(file_path)

    # Dynamically select tests based on the TESTS dictionary
    selected_tests = {key for key in TESTS.keys() if not key.startswith('#')}  # Exclude commented-out tests

    # Create a combined plot for execution time trends
    plot_execution_time_combined(data, selected_tests)

    # Plot the trend for a specific test with a fixed number of tasks
    specific_test = '3 5000 50 1000000 500'  # Change this to the test you want
    fixed_ntasks = 2  # Changed to a valid number of tasks (can be 1, 2, 4, 8, 16, or 32)
    plot_fixed_tasks(data, specific_test, fixed_ntasks)

if __name__ == "__main__":
    main()
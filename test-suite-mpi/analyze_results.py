import os
import re

# Checks MPI correctness and compute the speedup respect to the serial one

# TODO put this TEST list in and external common file for all the scripts
TESTS = {
    '5893 0.05 3 10 10': ['0.002 0.035', '2'],
    '8555 0.05 3 10 10': ['0.016 0.049', '1'],
    '12 100 5 10000 10000': ['76.732 61.943', '2209'],
    '-11 3500 20 500000 10': ['1984.878 1625.992', '35'],
    '1 5000 100 1000000 4': ['3936.506 131.472', '4'],
    '1 5000 100 1000000 100': ['3899.787 156.291', '163'],
    '1 5000 20 1000000 10': ['3918.912 143.364', '19'],
    '1 1000 3 10000 10000': ['287.788 261.446', '31'],
    '3 5000 50 1000000 300': ['3819.032 25.659', '469'],
    '3 5000 50 1000000 500': ['3738.436 58.743', '804'],
    '-1 1000 30 100000 1000': ['575.878 370.663', '1203'],
    '32 100 6 10000 100000': ['29.530 30.082', '4448']
}

def get_serial_outputs():
    SERIAL_OUTPUT = {}
    with open('results/serial_results.txt') as serial_file:
        for line in serial_file:
            if line.startswith('$ ./parsim'):
                test_params = ' '.join(line.split()[2:])
                SERIAL_OUTPUT[test_params] = ['', '', 0.0]  # Initialize with proper types
                coords = serial_file.readline().strip()
                collisions = serial_file.readline().strip()
                exec_time = serial_file.readline().strip()
                SERIAL_OUTPUT[test_params][0] = coords 
                SERIAL_OUTPUT[test_params][1] = collisions
                # Remove 's' and convert to float
                SERIAL_OUTPUT[test_params][2] = float(exec_time.rstrip('s'))
    return SERIAL_OUTPUT

# TODO run this correctness checker here for the serial too, and remove it from run_serial.sh, makes more sense
def check_correctness(mpi_dic):
    all_correct = True
    for config, tests in mpi_dic.items():
        ntasks, ncpu = config
        for test_params, (coords, collisions, _) in tests.items():
            if test_params in TESTS:
                expected_coords, expected_collisions = TESTS[test_params]
                if coords != expected_coords or collisions != expected_collisions:
                    all_correct = False
                    print(f"Warning: Mismatch for test {test_params}")
                    print(f"Expected: {expected_coords}, {expected_collisions}")
                    print(f"Got: {coords}, {collisions}")
    return all_correct

def compute_speedup(serial_dic, mpi_dic):
    speedups = {}
    with open('results/speedup_summary.txt', 'w') as f:
        f.write("Speedup Summary\n")
        
        for config, tests in mpi_dic.items():
            ntasks, ncpu = config
            f.write(f"\nConfiguration: {ntasks} tasks, {ncpu} CPUs per task\n")
            speedups[config] = {}
            
            for test_params, (_, _, mpi_time) in tests.items():
                if test_params in serial_dic:
                    serial_time = serial_dic[test_params][2]
                    if (serial_time > 0):
                        if (mpi_time > 0):
                            speedup = serial_time / mpi_time
                        else:
                            speedup = 9999999 # mpi_time 0s
                        speedups[config][test_params] = speedup
                        
                        f.write(f"\nTest: {test_params}\n")
                        f.write(f"Serial time: {serial_time:.1f}s\n")
                        f.write(f"MPI time: {mpi_time:.1f}s\n")
                        f.write(f"Speedup: {speedup:.2f}x\n")
    
    return speedups

def get_mpi_results_file():
    mpi_outputs_dir = 'results/mpi-outputs/'
    mpi_results_file = 'results/mpi_results.txt'
    pattern = re.compile(r'output_nt(\d+)_ncpu(\d+)_(\d+)\.txt')
    is_created = False

    # Collect and sort all matching files
    matching_files = []
    for filename in os.listdir(mpi_outputs_dir):
        match = pattern.match(filename)
        if match:
            nt = int(match.group(1))
            ncpu = int(match.group(2))
            file_id = int(match.group(3))
            matching_files.append((nt, ncpu, file_id, filename))
    
    # Sort by ntasks, ncpu, and file_id
    matching_files.sort()

    with open(mpi_results_file, 'w') as results_file:
        for nt, ncpu, file_id, filename in matching_files:
            test_id = file_id // 3  # normalize the id number
            params = list(TESTS.keys())[test_id]
            results_file.write(f"\n$ cat output_nt{nt}_ncpu{ncpu}_id{test_id}\n")
            results_file.write(f"params {params}\n")
            is_n_collisions = False
            full_path = os.path.join(mpi_outputs_dir, filename)
            with open(full_path) as job_output:
                for line in job_output:
                    line = line.strip()
                    if line.endswith('s'):
                        exec_time = float(line[:-1]) 
                    else:
                        if (is_n_collisions): 
                            collisions = line
                            results_file.write(f"{collisions}\n")
                            is_n_collisions = False
                        else: 
                            coords = line
                            results_file.write(f"{coords}\n")
                            is_n_collisions = True # always follows coords
                results_file.write(f"{exec_time}s\n")                    

        is_created = True

    return is_created

def get_mpi_outputs(): 
    MPI_OUTPUT = {}   
    with open('results/mpi_results.txt') as mpi_file:
        current_config = None
        current_params = None
        for line in mpi_file:
            line = line.strip()
            if not line:  # Skip empty lines
                continue
                
            if line.startswith('$ cat output_nt'):
                parts = line.split('_')
                ntasks = int(parts[1][2:])
                ncpu = int(parts[2][4:])
                current_config = (ntasks, ncpu)
                if current_config not in MPI_OUTPUT:
                    MPI_OUTPUT[current_config] = {}

            elif line.startswith('params'):
                current_params = ' '.join(line.split()[1:])
                coords = mpi_file.readline().strip()
                collisions = mpi_file.readline().strip()
                time_line = mpi_file.readline().strip()
                
                if time_line and time_line.endswith('s'):
                    # Remove 's' and convert to float
                    exec_time = float(time_line.rstrip('s'))
                    MPI_OUTPUT[current_config][current_params] = (coords, collisions, exec_time)

    return MPI_OUTPUT

def main():
    print("\nCollecting the Serial outputs..")
    SERIAL_OUTPUT = get_serial_outputs() 

    print("\nCollecting all the MPI outputs..")
    if get_mpi_results_file():
        MPI_OUTPUT = get_mpi_outputs()
    
    print("\nChecking correctness of MPI results...")
    if check_correctness(MPI_OUTPUT):
        print("All MPI tests passed!")
    else:
        print("Some tests failed - check warnings above")
    
    print("\nComputing speedups...")
    speedups = compute_speedup(SERIAL_OUTPUT, MPI_OUTPUT)
    print("Speedup summary written to results/speedup_summary.txt")

if __name__ == "__main__":
    main()
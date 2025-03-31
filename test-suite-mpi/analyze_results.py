
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
                SERIAL_OUTPUT[test_params] = ['', '', '']
                coords = serial_file.readline().strip()
                collisions = serial_file.readline().strip()
                exec_time = serial_file.readline().strip()
                SERIAL_OUTPUT[test_params][0] = coords 
                SERIAL_OUTPUT[test_params][1] = collisions
                SERIAL_OUTPUT[test_params][2] = float(exec_time[:-1])  
    return SERIAL_OUTPUT

def check_correctness(mpi_dic):
    all_correct = True
    for config, tests in mpi_dic.items():
        ntasks, ncpu = config
        print(f"\nChecking config: {ntasks} tasks, {ncpu} CPUs per task")
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
                    speedup = serial_time / mpi_time
                    speedups[config][test_params] = speedup
                    
                    f.write(f"\nTest: {test_params}\n")
                    f.write(f"Serial time: {serial_time:.1f}s\n")
                    f.write(f"MPI time: {mpi_time:.1f}s\n")
                    f.write(f"Speedup: {speedup:.2f}x\n")
    
    return speedups

def get_mpi_outputs(): 
    MPI_OUTPUT = {}   
    with open('results/mpi_outputs/byhand.txt') as mpi_file:
        current_config = None
        for line in mpi_file:
            line = line.strip()
            if line.startswith('cat output_nt'):
                # Extract ntasks and ncpu from filename
                parts = line.split('_')
                ntasks = int(parts[1][2:])
                ncpu = int(parts[2][4:].split('.')[0])
                current_config = (ntasks, ncpu)
                MPI_OUTPUT[current_config] = {}

            elif line.startswith('$ ./parsim-mpi'):
                test_params = ' '.join(line.split()[2:])
                coords = mpi_file.readline().strip()
                collisions = mpi_file.readline().strip()
                exec_time = float(mpi_file.readline().strip()[:-1])
                
                MPI_OUTPUT[current_config][test_params] = (coords, collisions, exec_time)

    return MPI_OUTPUT

def main():
    SERIAL_OUTPUT = get_serial_outputs() 
    MPI_OUTPUT = get_mpi_outputs()
    
    # Check correctness
    print("\nChecking correctness of MPI results...")
    if check_correctness(MPI_OUTPUT):
        print("All tests passed!")
    else:
        print("Some tests failed - check warnings above")
    
    # Compute and save speedups
    print("\nComputing speedups...")
    speedups = compute_speedup(SERIAL_OUTPUT, MPI_OUTPUT)
    print("Speedup summary written to results/speedup_summary.txt")

if __name__ == "__main__":
    main()
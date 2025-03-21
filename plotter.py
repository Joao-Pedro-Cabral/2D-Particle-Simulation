import matplotlib.pyplot as plt

def parse_speedup_data(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    speedup_data = {}
    current_test = None

    for line in lines:
        line = line.strip()
        if line.startswith("Test ["):
            current_test = line.split("[")[1].split("]")[0]
            speedup_data[current_test] = {}
        elif line.startswith("1 threads:") or line.startswith("2 threads:") or line.startswith("4 threads:") or line.startswith("8 threads:"):
            threads, speedup = line.split(":")
            threads = int(threads.split()[0])
            speedup = float(speedup.strip()) if speedup.strip() else None
            speedup_data[current_test][threads] = speedup

    return speedup_data

def plot_speedup(speedup_data):
    plt.figure(figsize=(12, 8))

    for test, data in speedup_data.items():
        threads = sorted(data.keys())
        speedups = [data[t] for t in threads]
        plt.plot(threads, speedups, marker='o', label=test)

    plt.xlabel('Number of Threads', fontsize=18)
    plt.ylabel('Speedup', fontsize=18)
    plt.title('Speedup vs Number of Threads', fontsize=20)
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    filename = 'output.txt'
    speedup_data = parse_speedup_data(filename)
    plot_speedup(speedup_data)
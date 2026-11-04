import h5py 
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from scipy.fft import rfft, rfftfreq
import os 
import concurrent.futures

def main(i, j):
    print(f"Processing file {i}.{j}")
    

    yf_file = f"output/{i}.{j}/fft.npy"
    xf_file = f"output/{i}.{j}/freq.npy"


    if os.path.exists(f"output/{i}.{j}/fft.npy"):

        xf = np.load(xf_file)
        yf = np.load(yf_file)


    else:
        h5_file = f"output/{i}.{j}/data.h5"
        
        with h5py.File(h5_file, 'r') as f:
            mx = f['m'][:-1, :, :, 0, 0]
            t = f['t'][:-1]
            
        mx = np.mean(mx, axis=(1, 2))
        
        
        xf = rfftfreq(len(mx), t[1] - t[0])
        yf = rfft(mx)
        
        
        np.save(f"output/{i}.{j}/fft.npy", yf)
        np.save(f"output/{i}.{j}/freq.npy", xf)


    max_freq = xf[np.argmax(np.abs(yf))]

    return max_freq
    

if __name__ == '__main__':
    j_list = []
    max_freq_list = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        futures = []
        for i in range(1, 20):
            for j in range(0, 10):
                j_val = float(f"{i}.{j}")
                futures.append((j_val, executor.submit(main, i, j)))

        for j_val, future in futures:
            j_list.append(j_val)
            max_freq_list.append(future.result())

    # plt.plot(j_list, max_freq_list, 'o')
    # plt.savefig("j_vs_freq.png")

    j_list = np.array(j_list)
    max_freq_list = np.array(max_freq_list)

    j_list_sorted = np.sort(j_list)
    max_freq_list_sorted = max_freq_list[np.argsort(j_list)]

    plt.plot(j_list_sorted, max_freq_list_sorted/1e9, 'o')
    plt.xlabel("J (10^11 A/m^2)")
    plt.ylabel("Frequency (GHz)")
    plt.savefig("j_vs_freq.png", dpi=300)
    plt.close()

    df = pd.DataFrame({'J': j_list_sorted, 'Frequency': max_freq_list_sorted/1e9})
    df.to_csv("j_vs_freq.csv", index=False)

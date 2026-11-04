import h5py 
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from scipy.fft import rfft, rfftfreq
import os 
import concurrent.futures

def main(i):
    print(f"Processing file {i}")
    

    yf_file = f"output/{i}/fft.npy"
    xf_file = f"output/{i}/freq.npy"


    if os.path.exists(f"output/{i}/fft.npy"):

        xf = np.load(xf_file)
        yf = np.load(yf_file)


    else:
        h5_file = f"output/{i}/data.h5"
        
        with h5py.File(h5_file, 'r') as f:
            mx = f['m'][:-1, :, :, 0, 0]
            t = f['t'][:-1]
            
        mx = np.mean(mx, axis=(1, 2))
        
        
        xf = rfftfreq(len(mx), t[1] - t[0])
        yf = rfft(mx)
        
        
        np.save(f"output/{i}/fft.npy", yf)
        np.save(f"output/{i}/freq.npy", xf)


    max_freq = xf[np.argmax(np.abs(yf))]

    return max_freq
    

if __name__ == '__main__':
    diameter_list = []
    max_freq_list = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        futures = []
        for i in range(100, 701, 10):
            futures.append((i, executor.submit(main, i)))

        for i, future in futures:
            diameter_list.append(i)
            max_freq_list.append(future.result())



    diameter_list = np.array(diameter_list)
    max_freq_list = np.array(max_freq_list)

    diameter_list_sorted = np.sort(diameter_list)
    max_freq_list_sorted = max_freq_list[np.argsort(diameter_list)]

    plt.plot(diameter_list_sorted, max_freq_list_sorted/1e9, 'o')
    plt.xlabel("Diameter (nm)")
    plt.ylabel("Frequency (GHz)")
    plt.savefig("diameter_vs_freq.png", dpi=300)
    plt.close()

    df = pd.DataFrame({'Diameter': diameter_list_sorted, 'Frequency': max_freq_list_sorted/1e9})
    df.to_csv("diameter_vs_freq.csv", index=False)

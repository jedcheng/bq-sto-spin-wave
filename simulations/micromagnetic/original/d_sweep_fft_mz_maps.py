import numpy as np 
from matplotlib import pyplot as plt

from scipy.fft import rfftn, rfftfreq
from scipy.signal import find_peaks

from scipy.interpolate import interp1d
import h5py
import argparse 

parser = argparse.ArgumentParser(description="FFT of M_z for CoFe J sweep")
parser.add_argument("--cpu", type=int, default=8, help="Number of CPUs to use")
args = parser.parse_args()


def main(i):
    print(f"Processing {i}")
    with h5py.File(f"output/{i}/data.h5", 'r') as hf:
        m = hf['m'][:-1, :, :, 0, 2]
        # t = hf['t'][:-1]
    #m = np.load(f"{diameter}/m_save_1.0e11.npy")

    vaccum = np.where(m[0, :, :] == 0)


    yf_file = f"output/{i}/fft.npy"
    xf_file = f"output/{i}/freq.npy"

    xf = np.load(xf_file)
    yf = np.load(yf_file)

    max_freq_id = np.argmax(np.abs(yf))
    max_freq = xf[max_freq_id]

    m_shape = m.shape


    # print(m_shape)



    yf_full = rfftn(m, axes=(0))

    # print(yf_full.shape, yf_full.dtype)

    avg = np.mean(np.abs(yf_full), axis=(1, 2))

    # print(avg.shape)


    peaks = find_peaks(avg, height=10)

    plt.plot(xf/1e9, avg)
    plt.plot(xf[peaks[0]]/1e9, avg[peaks[0]], "x")
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    plt.ylabel("Amplitude (a.u.)")
    plt.xlabel("Frequency (GHz)")
    plt.title(f"FFT of M_z of {i}")

    plt.savefig(f"output/{i}/fft_avg.png", dpi=300)
    plt.close()

    np.save(f"output/{i}/fft_avg.npy", avg)

    heat_map = yf_full[max_freq_id]
    heat_map[vaccum] = np.nan
    
    np.save(f"output/{i}/heat_map.npy", heat_map)

    phase_map = np.angle(heat_map)
    heat_map = np.abs(heat_map)




    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    im0 = axes[0].imshow(heat_map, cmap='jet', origin='lower')
    axes[0].set_title(f"Heatmap of M_z of {i}")
    im1 = axes[1].imshow(phase_map, cmap=plt.cm.hsv, origin='lower')
    axes[1].set_title(f"Phase map of M_z of {i}")
    plt.colorbar(im0, ax=axes[0])
    plt.colorbar(im1, ax=axes[1])
    plt.savefig(f"output/{i}/heatmap.png", dpi=300)
    plt.close()
    
    return

# import multiprocess  as mp
from joblib import Parallel, delayed

if __name__ == "__main__":
    tasks = [i for i in range(100, 701, 10)]
    Parallel(n_jobs=args.cpu)(delayed(main)(i) for i in tasks)

    # main(320)

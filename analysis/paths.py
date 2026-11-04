"""Where the scripts look for data (one setting: the environment variable DATA_DIR).

DATA_DIR       simulation outputs, default <repo>/data/sim.  Same layout as the public Nextcloud share
               (https://genkai-nextcloud.hpc.kyushu-u.ac.jp/s/iHs5kSBoLckLR6a), so a full or partial copy made by
               data/download.sh (or a local mount/copy of the share) can be used directly:
                   J_sweep/bq/<J:.1f>/            current sweep with BQ, D = 320 nm, J in 1e11 A/m^2
                   J_sweep/no_bq/<J:.1f>/         current sweep without BQ, D = 320 nm
                   D_sweep/bq/<D>/                diameter sweep with BQ, J = 6.0, D in nm
                   D_sweep/no_bq_5.5/<D>/         diameter sweep without BQ, J = 5.5
               each run directory holds freq.npy, fft.npy, fft_avg.npy, heat_map.npy, final.npy, initial.npy and
               data.h5 (see data/README.md).
DERIVED        derived tables shipped in git: <repo>/data/derived (wave_table.csv, power_balance.csv,
               power_balance_J.csv, j_vs_freq_bq.csv, j_vs_freq_no_bq.csv); override with BQSTO_DERIVED if needed.
"""
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.expanduser(os.environ.get("DATA_DIR", os.path.join(REPO, "data", "sim")))
DERIVED = os.environ.get("BQSTO_DERIVED", os.path.join(REPO, "data", "derived"))
D_SWEEP_J = {"bq": 6.0, "no_bq": 5.5}          # current density of the diameter sweeps (1e11 A/m^2)
D_SWEEP_DIR = {"bq": "bq", "no_bq": "no_bq_5.5"}  # folder names of the diameter sweeps in DATA_DIR/D_sweep


def run_dir(tag, J=None, D=None):
    """Run directory of the current sweep (J given, D = 320 nm) or of the diameter sweep (D given); tag = bq | no_bq."""
    if D is None:
        return os.path.join(DATA, "J_sweep", tag, f"{float(J):.1f}")
    return os.path.join(DATA, "D_sweep", D_SWEEP_DIR[tag], str(int(D)))

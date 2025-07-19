import pandas as pd
import numpy as np

def load_hz_data(include_farooq=True, deduplicate=True, return_raw=False):
    # BAO data
    bao = pd.DataFrame({
        "z": [0.106, 0.15, 0.32, 0.57, 0.61, 0.73, 2.34],
        "Hz": [69.0, 67.0, 79.2, 96.8, 97.3, 97.9, 222.0],
        "sigma": [19.6, 12.0, 5.6, 3.4, 2.1, 2.7, 7.0],
        "source": "BAO"
    })

    # Cosmic chronometers
    cc = pd.DataFrame({
        "z": [0.07, 0.12, 0.17, 0.179, 0.199, 0.27, 0.4, 0.48, 0.88, 1.3, 1.75, 2.0],
        "Hz": [69, 68.6, 83, 75, 75, 77, 95, 97, 90, 168, 202, 222],
        "sigma": [19.6, 26.2, 8, 4.9, 5, 14, 17, 60, 40, 17, 40, 41],
        "source": "CC"
    })

    # Farooq & Ratra
    farooq = pd.DataFrame({
        "z": [0.07, 0.09, 0.12, 0.17, 0.179, 0.199, 0.2, 0.27, 0.28, 0.352, 0.4,
              0.44, 0.48, 0.57, 0.593, 0.6, 0.68, 0.73, 0.781, 0.875, 0.88, 0.9,
              1.037, 1.3, 1.363, 1.43, 1.53, 1.75, 1.965, 2.34, 2.36],
        "Hz": [69, 69, 68.6, 83, 75, 75, 72.9, 77, 88.8, 83, 95,
               82.6, 97, 96.8, 104, 87.9, 92, 97.3, 105, 125, 90, 117,
               154, 168, 160, 177, 140, 202, 186.5, 222, 226],
        "sigma": [19.6, 12, 26.2, 8, 4.9, 5, 29.6, 14, 36.6, 14, 17,
                  7.8, 60, 3.4, 13, 6.1, 8, 7, 12, 17, 40, 23,
                  20, 17, 33.6, 18, 14, 40, 50.4, 7, 8],
        "source": "Farooq"
    })

    # Merge with priority: Farooq > CC > BAO
    dfs = []
    if include_farooq:
        dfs.append(farooq)
    dfs.extend([cc, bao])

    hz_all = pd.concat(dfs, ignore_index=True)

    if deduplicate:
        hz_all = hz_all.drop_duplicates(subset="z", keep="first")

    hz_all = hz_all.sort_values("z").reset_index(drop=True)

    if return_raw:
        return hz_all

    return hz_all["z"].values, hz_all["Hz"].values, hz_all["sigma"].values

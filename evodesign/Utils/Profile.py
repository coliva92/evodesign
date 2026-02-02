from .AminoAcids import AMINO_ACIDS, MAP_AMINO_ACID_TO_INT
import numpy as np
import numpy.typing as npt


def load_profile(filename: str) -> npt.NDArray[np.float64]:
    profile = None
    i = 0
    for line in open(filename, "rt", encoding="utf-8"):
        line = line.strip()
        if i == 0:
            sequence_length = int(float(line))
            profile = np.zeros((sequence_length, len(AMINO_ACIDS)), dtype=np.float64)
            i += 1
            continue
        columns = line.split(" ")
        res_idx = int(float(columns[0])) - 1
        for col in columns[1:]:
            tmp = col.split(":")
            aa_idx = MAP_AMINO_ACID_TO_INT[tmp[0]]
            prob = float(tmp[1])
            profile[res_idx][aa_idx] = prob
    for res_idx, total in enumerate(np.sum(profile, axis=1).tolist()):
        if total > 1.0:
            raise ValueError
        if total == 0.0:
            profile[res_idx, :] = 1.0 / len(AMINO_ACIDS)
            continue
        p = 1.0 - total
        mask = profile[res_idx, :] == 0.0
        profile[res_idx, mask] = p / np.sum(mask)
    assert profile is not None
    return profile


def save_profile(profile: npt.NDArray[np.float64], filename: str):
    with open(filename, "wt", encoding="utf-8") as txt:
        txt.write(f"{profile.shape[0]}\n")
        for i in range(profile.shape[0]):
            txt.write(f"{i}")
            for j in range(profile.shape[1]):
                txt.write(f" {AMINO_ACIDS[j]}:{profile[i][j]}")
            txt.write("\n")

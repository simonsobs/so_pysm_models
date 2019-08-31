import healpy as hp
import sys
from pathlib import Path

if __name__ == "__main__":

    folder = Path(sys.argv[1])

    out_folder = folder / "equatorial"

    rot = hp.Rotator(coord=["G","C"])

    for f in folder.glob("*.fits"):
        print(f)
        f = str(f)
        m = hp.read_map(f)
        m_equatorial = rot.rotate_map_alms(m)
        out_file = f.replace(str(folder), str(out_folder))
        print(out_file)
        hp.write_map(out_file, m_equatorial)

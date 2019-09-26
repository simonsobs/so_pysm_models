import healpy as hp
import sys
from pathlib import Path

if __name__ == "__main__":

    folder = Path(sys.argv[1])

    out_folder = folder / "512"


    for f in folder.glob("[cib,tsz,ksz,kap]*.fits"):
        print(f)
        f = str(f)
        m = hp.read_map(f)
        m_out = hp.ud_grade(m, 512)
        out_file = f.replace(str(folder), str(out_folder))
        print(out_file)
        hp.write_map(out_file, m_out)

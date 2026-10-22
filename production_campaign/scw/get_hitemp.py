import os, radis
DLP = os.path.join(os.environ["DATA"], "radisdb")
radis.config["DEFAULT_DOWNLOAD_PATH"] = DLP
from radis.io.hitemp import fetch_hitemp
df = fetch_hitemp("H2O",
                  local_databases=os.path.join(DLP, "hitemp"),
                  load_wavenum_min=10, load_wavenum_max=11000,
                  parallel=False)          # serial parse -> ~1/3 peak memory, avoids the login-node kill
print("HITEMP H2O lines in 10-11000 cm-1:", len(df))

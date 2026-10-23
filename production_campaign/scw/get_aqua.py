"""get_aqua.py -- download the AQUA rho-T table without git-lfs (via the LFS batch API)."""
import json, urllib.request, zipfile, os
OID="e731d43bff1905b6821d80e295688b0cedac95247c8df6485e13429281505224"; SIZE=19311061
req=urllib.request.Request("https://github.com/mnijh/AQUA.git/info/lfs/objects/batch",
    data=json.dumps({"operation":"download","transfers":["basic"],
                     "objects":[{"oid":OID,"size":SIZE}]}).encode(),
    headers={"Accept":"application/vnd.git-lfs+json","Content-Type":"application/vnd.git-lfs+json"})
href=json.load(urllib.request.urlopen(req,timeout=60))["objects"][0]["actions"]["download"]["href"]
print("downloading AQUA rhoT table ...")
urllib.request.urlretrieve(href,"aqua_eos_rhot_v1_0.zip")
zipfile.ZipFile("aqua_eos_rhot_v1_0.zip").extractall(".")
print("extracted:", [f for f in os.listdir(".") if f.startswith("aqua_eos_rhot")])

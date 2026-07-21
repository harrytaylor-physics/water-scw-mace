import json, os, pathlib, pytest

ROOT = pathlib.Path(__file__).parent.parent
POT = ROOT / "models" / "water_scw_v4.model"
DIP = ROOT / "models" / "water_dipole_v2.model"
GOLDEN = pathlib.Path(__file__).parent / "golden"
DEVICE = os.environ.get("WATER_SCW_DEVICE", "cpu")

@pytest.fixture(scope="session")
def calc():
    if not POT.exists():
        pytest.skip(f"model absent: {POT} (see models/README.md)")
    from mace.calculators import MACECalculator
    return MACECalculator(model_paths=str(POT), device=DEVICE,
                          default_dtype="float64")

@pytest.fixture(scope="session")
def dipole_calc():
    if not DIP.exists():
        pytest.skip(f"dipole model absent: {DIP}")
    from mace.calculators import MACECalculator
    try:
        return MACECalculator(model_paths=str(DIP), device=DEVICE,
                              default_dtype="float64", model_type="DipoleMACE")
    except Exception as e:
        pytest.skip(f"dipole model failed to load as DipoleMACE: {e}")

@pytest.fixture()
def water_dimer():
    from ase import Atoms
    return Atoms("OHHOHH",
        positions=[[0,0,0],[0.76,0.59,0],[-0.76,0.59,0],
                   [0,0,2.9],[0.76,0.59,2.9],[-0.76,0.59,2.9]],
        cell=[12,12,12], pbc=True)

@pytest.fixture()
def single_water_small_box():
    """One molecule in a 5 A box: L/2 < r_max, so the model MUST handle
    multiple periodic images of the same neighbour (as in the dense
    production cells, where L ~ 9.5 A < 2*r_max)."""
    from ase import Atoms
    return Atoms("OHH",
        positions=[[2.5,2.5,2.5],[3.26,3.09,2.5],[1.74,3.09,2.5]],
        cell=[5,5,5], pbc=True)

def golden_check(name, values, tol):
    """Self-pinning golden values: first run writes tests/golden/<name>.json
    (commit it); later runs compare. Returns (ref, is_first_run)."""
    GOLDEN.mkdir(exist_ok=True)
    f = GOLDEN / f"{name}.json"
    if not f.exists():
        f.write_text(json.dumps(values, indent=2))
        return values, True
    ref = json.loads(f.read_text())
    return ref, False

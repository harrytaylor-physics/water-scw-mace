import os, pathlib, pytest

MODEL = pathlib.Path(__file__).parent.parent / "models" / "water_scw_v4.model"
DEVICE = os.environ.get("WATER_SCW_DEVICE", "cpu")

@pytest.fixture(scope="session")
def calc():
    if not MODEL.exists():
        pytest.skip(f"model file absent: {MODEL} (copy via models/README.md)")
    from mace.calculators import MACECalculator
    return MACECalculator(model_paths=str(MODEL), device=DEVICE,
                          default_dtype="float64")

@pytest.fixture()
def water_dimer():
    from ase import Atoms
    return Atoms("OHHOHH",
        positions=[[0,0,0],[0.76,0.59,0],[-0.76,0.59,0],
                   [0,0,2.9],[0.76,0.59,2.9],[-0.76,0.59,2.9]],
        cell=[12,12,12], pbc=True)

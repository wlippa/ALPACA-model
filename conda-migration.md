# Conda Migration

This note describes the changes needed to make ALPACA installable with conda or mamba while keeping Gurobi support available as an optional second step.

## Target packaging model

- Publish the base `alpaca` package without a hard dependency on `gurobipy`.
- Make the Pyomo backend the default solver for a fresh install.
- Keep Gurobi as an optional add-on for users with a license.
- Support this user flow:

```bash
mamba create -n alpaca -c conda-forge alpaca
mamba install -n alpaca -c gurobi gurobi>=13
conda activate alpaca
grbgetkey YOUR-LICENSE-KEY
```

Users who install only the base package should be able to run ALPACA with the Pyomo backend and an open-source solver such as SCIP, CBC, or GLPK.

## Required code changes

### 1. Remove Gurobi as a hard runtime dependency

In `pyproject.toml`:

- Remove `gurobipy` from `[project].dependencies`.
- Keep the core scientific Python dependencies in the base package.
- Optionally define an extra such as `gurobi`, but this is not required for conda packaging.

Reason:

- The public conda package should be installable without requiring a proprietary solver package.
- Users without Gurobi access should still get a working installation.

### 2. Make Pyomo the default solver

Update the default solver from `gurobi` to `pyomo` in the runtime defaults.

Places to review:

- `alpaca/ALPACA_model_class.py`
- `alpaca/solvers/factory.py`
- Any CLI argument defaults or config builders that currently assume `gurobi`
- Any README examples that imply Gurobi is the default backend

Reason:

- After removing the hard Gurobi dependency, the default behavior must still work in a clean conda environment.

### 3. Keep Gurobi import paths optional

The code already lazy-imports `gurobipy` in `alpaca/ALPACA_model_class.py`, which is the right direction. Preserve that behavior and ensure:

- Importing the package does not fail when Gurobi is absent.
- Selecting `solver=gurobi` without Gurobi installed raises a clear error message.
- Tests cover the no-Gurobi base installation path.

## Required packaging work

### 4. Add a conda recipe

Create a `recipe/` directory with a `meta.yaml` file.

The recipe should:

- Build the package from the source tree using `pip install . --no-deps`.
- Use the Python build backend requirements from `pyproject.toml`.
- Package ALPACA as `noarch: python` if no platform-specific code is required.
- Include a small smoke test such as:

```bash
python -c "import alpaca"
alpaca --help
```

### 5. Translate dependencies to conda package names

Check each Python dependency against conda-forge naming. In particular:

- `kaleido` may need to become `python-kaleido` in the conda recipe.
- Solver executables used by Pyomo should be listed explicitly if they are meant to work out of the box.

For the base package, prefer open-source solver support only.

### 6. Review versioning for source builds

The project currently uses `hatch-vcs` for dynamic versioning. Confirm that conda builds from source archives get a valid version.

If needed:

- add a fallback version strategy, or
- ensure builds are performed from tagged sources with version metadata available.

## Documentation updates

### 7. Rewrite installation instructions

Update `README.md` to document two installation paths.

Base install:

```bash
mamba create -n alpaca -c conda-forge alpaca
```

Optional Gurobi install:

```bash
mamba install -n alpaca -c gurobi gurobi>=13
```

The README should explicitly state:

- Pyomo is the default backend.
- Gurobi is optional and requires a valid license.
- Users who want Gurobi must install it separately into the same environment.

### 8. Document solver expectations clearly

Document which solver is used in each case:

- Base install: Pyomo backend with an open-source solver.
- Gurobi-enabled install: either direct Gurobi backend or Pyomo configured to use Gurobi, depending on the intended supported path.

If ALPACA supports both, pick one documented recommendation to avoid ambiguity.

## Test and release checklist

Before publishing:

1. Build the wheel with `hatch build`.
2. Build the conda package locally with `conda-build` or `mambabuild`.
3. Create a fresh environment and install only the base package.
4. Verify import, CLI startup, and one small example using the Pyomo backend.
5. Add Gurobi to the same environment and verify the Gurobi backend works.
6. Confirm that the error message is clear when a user requests Gurobi without having it installed.

## Recommended migration order

1. Remove `gurobipy` from the base package dependencies.
2. Change default solver selection to Pyomo.
3. Update README installation and solver documentation.
4. Add the conda recipe.
5. Test both base and Gurobi-enabled environments.
6. Submit the package to the target conda channel.

## Scope decision recorded here

The agreed migration plan is:

- base conda package uses Pyomo by default
- Gurobi is installed separately as a second step for licensed users
- no requirement to hide Gurobi support, only to decouple it from the default installation
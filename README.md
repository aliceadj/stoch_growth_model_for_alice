# Stochastic Growth Model

## Setup

To set up the project environment and install dependencies:

**For local development (using root project files):**
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

**For server usage:**
```julia
using Pkg
Pkg.activate("stoch_growth_model_server")
Pkg.instantiate()
```

## Usage

### Using on a different machine

### Running Simulations

The main simulation script is located at `./scripts/simulate/simulate.jl`. This script runs parameter sweeps across different `ns` (nutrient sensitivity) and `abx` (antibiotic concentration) values using distributed computing.

**To run simulations:**

On the server, open a new tmux session:

```bash
tmux
```

Navigate to the project directory:

```bash
cd stoch_growth_model
```

Run the simulation script:
```bash
julia scripts/simulate/simulate.jl
```

**What the script does:**
- Runs simulations across parameter combinations of `ns_vec` and `abx_vec`
- Uses parallel workers for distributed computing
- Saves results as Arrow files organized by parameter values
- Generates a detailed log file tracking simulation performance
- Creates output directory structure: `simulation_data/[timestamp]/ns[value]/sim_abx[value].arrow`

**Key parameters to modify:**
- `ns_vec`: Array of nutrient sensitivity values
- `abx_vec`: Array of antibiotic concentration values  
- `addprocs(n)`: Number of parallel workers (adjust based on your system)
- Time span: Currently set to `(0.0, 100000)` in the DiscreteProblem but simulation will end once 500 cell cycles have occurred

**Output:**
- **Data files**: `simulation_data/[timestamp]/ns[value]/sim_abx[value].arrow`
- **Log file**: `simulation_data/[timestamp]/simulation.log` with timing information
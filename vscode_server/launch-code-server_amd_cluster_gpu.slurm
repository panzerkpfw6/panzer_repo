
#!/bin/bash
#SBATCH --job-name=cs_mi300a
#SBATCH --partition=MI300A_SPX
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=1
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j-slurm.out
#SBATCH --error=%x-%j-slurm.err

set -euo pipefail

# Environment
module load rocm

export PATH="$HOME/.local/bin:$PATH"
export CODE_SERVER_EXTENSIONS="$HOME/code-server/extensions"
export XDG_CONFIG_HOME="$HOME/.config"

mkdir -p "$CODE_SERVER_EXTENSIONS"

PROJECT_DIR="$HOME/panzer_repo"
cd "$PROJECT_DIR"

# Verify tools
command -v code-server || {
    echo "ERROR: code-server is not installed"
    exit 1
}

echo "Compute node: $(hostname)"
echo "GPU allocation: $SLURM_JOB_GPUS"

amd-smi list || true
hipcc --version || true

# Select an available port
COMPUTE_NODE=$(hostname -s)
CODE_SERVER_PORT=$(python3 -c \
'import socket; s=socket.socket(); s.bind(("",0)); print(s.getsockname()[1]); s.close()')

echo ""
echo "============================================"
echo "AMD AAC6 VS Code Server"
echo "============================================"
echo "Node: $COMPUTE_NODE"
echo "Port: $CODE_SERVER_PORT"
echo ""
echo "On your MacBook, run:"
echo ""
echo "ssh -N -L ${CODE_SERVER_PORT}:${COMPUTE_NODE}:${CODE_SERVER_PORT} paplotnitskii@aac6.amd.com"
echo ""
echo "Then open:"
echo "http://localhost:${CODE_SERVER_PORT}"
echo "============================================"

# Use password authentication because the server
# listens on the compute-node network interface.
# Configure your password in ~/.config/code-server/config.yaml
exec code-server \
    --auth password \
    --bind-addr "${COMPUTE_NODE}:${CODE_SERVER_PORT}" \
    --extensions-dir "$CODE_SERVER_EXTENSIONS" \
    "$PROJECT_DIR"

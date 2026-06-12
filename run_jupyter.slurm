#!/bin/bash -l

#SBATCH --job-name=jupyter
#SBATCH --output=jupyter-%j.log
#SBATCH --time=12:00:00
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mail-type=end
#SBATCH --mail-user=abarne65@jhu.edu
##SBATCH --account=mfalk1

# Load environment
source /data/apps/go.sh
source /home/abarne65/scratch16-mfalk1/annie/jpytr/myenv_ovito311/bin/activate

# Choose a port
PORT=8888

# Get the hostname of the compute node
HOSTNAME=$(hostname)

echo -e "\n"
echo "===================================================="
echo "Jupyter Notebook running on node: $HOSTNAME"
echo "To connect, open a terminal on your local machine and run:"
echo -e "\n"
echo "ssh -N -L ${PORT}:${HOSTNAME}:${PORT} your_username@your_cluster_login_node"
echo -e "\nThen open in your browser: http://localhost:${PORT}"
echo "===================================================="
echo -e "\n"

# Launch Jupyter
jupyter-notebook --no-browser --port=$PORT --ip=0.0.0.0

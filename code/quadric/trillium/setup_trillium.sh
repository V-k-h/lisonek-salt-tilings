#!/bin/bash
# Run once on a Trillium login node: builds the virtualenv used by m2_dim4.sbatch.
# numpy comes from the Alliance wheelhouse; python-sat (imported by cross_tiling.py) from there if it is
# listed, otherwise from PyPI (login nodes have internet, compute nodes do not).
set -e
module load python/3.13
virtualenv --no-download $HOME/.virtualenvs/ks
source $HOME/.virtualenvs/ks/bin/activate
pip install --no-index --upgrade pip
pip install --no-index numpy
pip install --no-index python-sat || pip install python-sat
python -c "import numpy, pysat; print('ok: numpy', numpy.__version__, 'pysat', pysat.__version__)"

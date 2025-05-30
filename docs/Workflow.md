# Workflow Overview
Here we discuss the workflow in more detail, describing the scripts, their roles and some important details.


### Testing
Note the awkward code below is needed at the top of testing scripts. It allows testing scripts to find files in all dirs of the project.

'''import sys
from pathlib import Path

# Add project root to sys.path
sys.path.append(str(Path(__file__).resolve().parents[1]))'''




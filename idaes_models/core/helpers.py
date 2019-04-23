import os

# FIXME (from Qi): is this necessary anymore? Should it be moved to util?
def cmd_exists(cmd):
    """
    Returns true if command exists in the users PATH env variable
    """
    return any(
        os.access(os.path.join(path, cmd), os.X_OK)
        for path in os.environ["PATH"].split(os.pathsep)
    )

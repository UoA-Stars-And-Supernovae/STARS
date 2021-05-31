from datetime import datetime
from pathlib import Path

def exp_to_int(exp):
    exp = exp.lower().split("e")
    return float(exp[0])*10**int(exp[1])

with open('modout', 'r') as f:
    file = f.readlines()

    n_models = len(file) / 399

    # First, we create a directory to store the output models
    T = datetime.now().isoformat().split(".")[0]
    Path(f"model_outputs/{T}").mkdir(parents=True, exist_ok=True)

    for i in range(int(n_models)):
        """
        Each model in modout is 399 lines long. The first line is the
        model information:
            1:  stellar mass
            2:  current timestep size
            3:  model age
            4:  period of the binary
            5:  total mass of the binary
            6:  artificial energy generation time
            7:  number of mesh points in the model
            8:  desired number of models to be computed
            9:  starting model number
            10: number indicating whether this is star 1 or star 2 of binary
            11: pressure in H-burning shell
            12: pressure in He-burning shell
        """
        model = file[i*399:i*399+399]

        model_info = []
        info = model[0].strip().split(" ")

        # Fixed-length formatting is a pain
        for inf in info:
            if inf.strip() != '':
                model_info.append(inf)

        # Convert 1.435E+00 to 1.435 for example
        mass = exp_to_int(model_info[0])

        # output to file, tag everything by stellar mass
        fname = f"model_outputs/{T}/modout_{mass}.mod"

        with open(fname, 'w+') as f2:
            f2.writelines(model)
